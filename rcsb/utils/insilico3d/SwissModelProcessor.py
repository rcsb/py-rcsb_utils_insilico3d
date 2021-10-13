##
# File:    SwissModelProcessor.py
# Author:  Dennis Piehl
# Date:    11-Oct-2021
#
# Update:
#
# To Do:
# - pylint: disable=fixme
# - Add mkdssp calculation
# - Add check that converted files are consistent with mmCIF dictionaries
##

"""
Processors for converting SwissModel PDB models into mmCIF format.

"""

__docformat__ = "google en"
__author__ = "Dennis Piehl"
__email__ = "dennis.piehl@rcsb.org"
__license__ = "Apache 2.0"

import logging
import os.path
import time
import collections

from rcsb.utils.insilico3d.SwissModelPdbToCifConverter import SwissModelPdbToCifConverter
from rcsb.utils.insilico3d import __version__
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.multiproc.MultiProcUtil import MultiProcUtil

logger = logging.getLogger(__name__)

SwissModelEntry = collections.namedtuple("SwissModelEntry", ["name", "model", "uniProtId"])


class SwissModelWorker(object):
    """A skeleton class that implements the interface expected by the multiprocessing
    for converting SwissModel PDB files to mmCIF files.
    """

    def __init__(self, **kwargs):
        _ = kwargs
        self.__workPath = kwargs.get("workPath", None)
        self.__fU = FileUtil(workPath=self.__workPath)

    def convert(self, dataList, procName, optionsD, workingDir):
        """Enumerate SwissModel PDB model files for conversion to mmCIF.

        Args:
            dataList (list): list of model files to convert, each stored as a dict
            procName (str): worker process name
            optionsD (dict): dictionary of additional options that worker can access
            workingDir (str): path to working directory

        Returns:
            successList (list): list of input data items that were successfully processed
            retList (list): list of processed results
            diagList (list): list of unique diagnostics
        """
        _ = workingDir
        _ = optionsD
        successList = []
        failList = []
        retList = []
        diagList = []

        try:
            converterObj = SwissModelPdbToCifConverter(cachePath=workingDir)
            for modelE in dataList:
                convertedCifFileZ = None
                # calculatedDssp = False
                modelNameRoot, pdbFile, uniProtId = modelE.name, modelE.model, modelE.uniProtId
                modelDir = os.path.dirname(self.__fU.getFilePath(pdbFile))
                if pdbFile.endswith(".pdb"):
                    cF = os.path.join(modelDir, modelNameRoot + ".cif")
                    # Then attempt the conversion with alignment
                    convertedCifFileZ = self.__convertPdbToCif(procName=procName, pdbFile=pdbFile, mmCifOutFile=cF,
                                                               uniProtId=uniProtId, swissModelConverter=converterObj, optionsD=optionsD)
                    if convertedCifFileZ:
                        successList.append(modelE)
                    print('\n', convertedCifFileZ)
                    self.__fU.remove(cF)
                retList.append((modelE, convertedCifFileZ, uniProtId))

            failList = sorted(set(dataList) - set(successList))
            if failList:
                logger.info("%s returns %d definitions with failures: %r", procName, len(failList), failList)

            logger.debug("%s converted %d/%d models, failures %d", procName, len(retList), len(dataList), len(failList))
        except Exception as e:
            logger.exception("Failing %s for %d data items %s", procName, len(dataList), str(e))

        return successList, retList, diagList

    def __convertPdbToCif(self, procName, pdbFile, mmCifOutFile, uniProtId, swissModelConverter, optionsD):
        """Internal method to convert a SwissModel PDB file to mmCIF format.

        Args:
            procName (str): worker process name
            pdbFile (str): path to PDB model file to convert
            mmCifOutFile (str): path to which to write mmCIF model file
            uniProtId (str): UniProt ID corresponding to the given model file
            swissModelConverter (obj): SwissModelPdbToCifConverter class object
            optionsD (dict): additional options/parameters to use in conversion

        Returns:
            str: path to converted (and gzipped) mmCIF file if successful; otherwise None
        """

        try:
            species = optionsD.get("species")
            ok = swissModelConverter.convertPdbToCif(pdbFileIn=pdbFile, cifFileOut=mmCifOutFile, organism=species, uniProtId=uniProtId)
            assert ok
            mmCifOutFileZ = mmCifOutFile + ".gz"
            ok = self.__fU.compress(inpPath=mmCifOutFile, outPath=mmCifOutFileZ)
            if ok:
                return mmCifOutFileZ
        except Exception as e:
            logger.exception("%s failing on PDB file %s with %s", procName, pdbFile, str(e))
            # logger.debug("%s failing on PDB file %s, with %s", procName, pdbFile, str(e))


class SwissModelProcessor(object):
    """Generators and accessors for SwissModel model files."""

    def __init__(self, cachePath=None, useCache=False, speciesD=None, **kwargs):
        """Initialize SwissModelProcessor object.

        Args:
            cachePath (str): path to species-specific cache file containing list of processed model files.
            useCache (bool): whether to use the existing data cache or re-run conversion process
            speciesD (dict): dictionary containing the following necessary key:value pairs for conversion:
                             "speciesModelDir": path to species data directory
                             "lastModified": last modified date of the downloaded species archive tarball
                             "speciesName": name of the species as it is stored in the SwissModelProvider cache
                             "speciesPdbModelFileList": list of the PDB model files to convert
        """

        try:
            self.__version = __version__
            self.__numProc = kwargs.get("numProc", 2)
            self.__chunkSize = kwargs.get("chunkSize", 10)

            self.__speciesModelDir = speciesD.get("speciesModelDir")
            self.__speciesModDate = speciesD.get("lastModified", None)
            self.__speciesName = speciesD.get("speciesName")
            self.__speciesPdbModelFileList = speciesD.get("speciesPdbModelFileList", [])

            self.__cachePath = cachePath if cachePath else self.__getModelCachePath()

            self.__mU = MarshalUtil(workPath=self.__speciesModelDir)
            self.__fU = FileUtil(workPath=self.__speciesModelDir)

            self.__modelD = self.__reload(fmt="pickle", useCache=useCache)

        except Exception as e:
            logger.exception("Failing with %s", str(e))

    def testCache(self, minCount=1):
        try:
            if minCount == 0:
                return True
            if self.__modelD and minCount and "modelsCif" in self.__modelD and len(self.__modelD["modelsCif"]) >= minCount:
                logger.info("mmCIF models converted for (%d) models, created %r, version %r", len(self.__modelD["modelsCif"]), self.__modelD["created"], self.__modelD["version"])
                return True
        except Exception:
            pass
        return False

    def getModelsCif(self):
        """Return a list of model root names for which models that have been successfully converted to mmCIF.

        Returns:
            (list): [modelName1, modelName2, ... ]
        """
        try:
            return list(self.__modelD["modelsCif"].keys())
        except Exception:
            pass
        return []

    def generate(self, updateOnly=False, fmt="pickle", indent=0):
        """Generate converted mmCIF models from SwissModel PDB files.

        Args:
            updateOnly (bool): only convert new or previously-failed models.  Defaults to False.
            fmt (str, optional): export file format. Defaults to "pickle".
            indent (int, optional): json format indent. Defaults to 0.

        Returns:
            bool: True for success or False otherwise
        """
        ok = False
        try:
            tS = time.strftime("%Y %m %d %H:%M:%S", time.localtime())
            mD, failD = self.__convertSwissModelPdb(numProc=self.__numProc, chunkSize=self.__chunkSize, updateOnly=updateOnly)
            self.__modelD = {
                "version": self.__version,
                "created": tS,
                "species": self.__speciesName,
                "archiveModDate": self.__speciesModDate,
                "speciesModelDir": self.__speciesModelDir,
                "modelsCif": mD,
                "modelsFailed": failD,
            }
            kwargs = {"indent": indent} if fmt == "json" else {"pickleProtocol": 4}
            modelCachePath = self.__getModelCachePath(fmt=fmt)
            ok = self.__mU.doExport(modelCachePath, self.__modelD, fmt=fmt, **kwargs)
            logger.info("Wrote %r status %r", modelCachePath, ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return ok

    def reload(self, fmt="pickle", useCache=True):
        self.__modelD = self.__reload(fmt=fmt, useCache=useCache)
        return self.__modelD is not None

    def __reload(self, fmt="pickle", useCache=True):
        """Reload from the current cache directory."""
        try:
            modelCachePath = self.__getModelCachePath(fmt=fmt)
            tS = time.strftime("%Y %m %d %H:%M:%S", time.localtime())
            modelD = {"version": self.__version, "created": tS, "modelsCif": {}}
            logger.debug("useCache %r modelCachePath %r", useCache, modelCachePath)
            #
            if useCache and self.__mU.exists(modelCachePath):
                modelD = self.__mU.doImport(modelCachePath, fmt=fmt)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        #
        return modelD

    def getCachePath(self):
        return self.__cachePath

    def __getModelCachePath(self, fmt="pickle"):
        ext = "pic" if fmt == "pickle" else "json"
        speciesNameNoSpace = self.__speciesName.replace(" ", "_")
        pth = os.path.join(self.__speciesModelDir, speciesNameNoSpace + "-model-data." + ext)
        return pth

    def __convertSwissModelPdb(self, numProc=2, chunkSize=10, updateOnly=False):
        """Prepare multiprocessor queue and workers for converting all SwissModel PDB model files to mmCIF.

        Args:
            numProc (int, optional): number of processes to use. Defaults to 2.
            chunkSize (int, optional): incremental chunk size used for distribute work processes. Defaults to 10.
            updateOnly (bool): only convert new or previously-failed models.  Defaults to False.

        Returns:
            mD (dict): dictionary of successfully-converted mmCIF models, in the following structure:
                       {modelNameRoot: {"model": pathToGzippedCifFile, "uniProtId": uniProtId}, ...}
            failD (dict): dictionary of SwissModel models for which conversion failed, in the following structure:
                          {modelNameRoot: {"model": pathToCompressedPDBModelFile, "uniProtId": uniProtId}, ...}
        """
        mD = {}
        failD = {}
        exD = {}
        #
        # updateOnly - will reuse any existing data loaded when this is instantiated
        #              otherwise the cache context is cleared before the calculation.
        if updateOnly:
            exD = {k: True for k in self.getModelsCif()}
            logger.info("Reusing (%d) already converted mmCIF models", len(exD))
            mD = self.__modelD["modelsCif"] if "modelsCif" in self.__modelD else {}
        #
        pdbModelFileList = self.__speciesPdbModelFileList
        modelList = []
        for pF in pdbModelFileList:
            modelNameRoot = self.__fU.getFileName(pF).split(".pdb")[0]
            if modelNameRoot not in mD:
                uniProtId = "".join(pF.split("/")[-5:-2])
                modelE = SwissModelEntry(name=modelNameRoot, model=pF, uniProtId=uniProtId)
                modelList.append(modelE)

        logger.info("Starting with %d models, numProc %d, updateOnly (%r)", len(modelList), self.__numProc, updateOnly)
        #
        rWorker = SwissModelWorker(workPath=self.__speciesModelDir)
        mpu = MultiProcUtil(verbose=True)
        optD = {"species": self.__speciesName, "speciesModDate": self.__speciesModDate}
        mpu.setOptions(optD)
        mpu.set(workerObj=rWorker, workerMethod="convert")
        mpu.setWorkingDir(workingDir=self.__speciesModelDir)
        ok, failList, resultList, _ = mpu.runMulti(dataList=modelList, numProc=numProc, numResults=1, chunkSize=chunkSize)
        if failList:
            logger.info("mmCIF conversion failures (%d): %r", len(failList), failList)
        #
        for (model, convertedCifFileZ, uniProtId) in resultList[0]:
            if convertedCifFileZ:
                mD[model.name] = {"model": convertedCifFileZ, "uniProtId": uniProtId}
            else:
                failD[model.name] = {"model": model.model, "uniProtId": uniProtId}
        #
        logger.info("Completed with multi-proc status %r, failures %r, total models with data (%d)", ok, len(failList), len(mD))
        return mD, failD

    def convertJsonPickle(self, fmt1="json", fmt2="pickle"):
        # Keep method name as "convertJsonPickle" instead of "convert", to avoid complications from using the other "convert" for MPU method above
        modelCachePath = self.__getModelCachePath(fmt=fmt1)
        self.__modelD = self.__mU.doImport(modelCachePath, fmt=fmt1)
        #
        modelCachePath = self.__getModelCachePath(fmt=fmt2)
        ok = self.__mU.doExport(modelCachePath, self.__modelD, fmt=fmt2, pickleProtocol=4)
        return ok
