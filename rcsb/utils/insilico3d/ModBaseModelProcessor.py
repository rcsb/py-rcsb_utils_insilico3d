##
# File:    ModBaseModelProcessor.py
# Author:  Dennis Piehl
# Date:    27-Sep-2021
#
# Updates:
#
# To Do:
# - pylint: disable=fixme
# - Add mkdssp calculation
# - Add check that converted files are consistent with mmCIF dictionaries
##

"""
Processors for converting ModBase PDB models into mmCIF format.

"""

__docformat__ = "google en"
__author__ = "Dennis Piehl"
__email__ = "dennis.piehl@rcsb.org"
__license__ = "Apache 2.0"

import logging
import os.path
import time
import collections

import modelcif

import rcsb.utils.modbase_utils.modbase_pdb_to_cif as modbase
from rcsb.utils.insilico3d import __version__
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.multiproc.MultiProcUtil import MultiProcUtil

logger = logging.getLogger(__name__)

ModBaseEntry = collections.namedtuple("ModBaseEntry", ["name", "model", "alignment"])


class ModBaseModelWorker(object):
    """A skeleton class that implements the interface expected by the multiprocessing
    for converting ModBase PDB files to mmCIF files.
    """

    def __init__(self, **kwargs):
        _ = kwargs
        self.__workPath = kwargs.get("workPath", None)
        self.__fU = FileUtil(workPath=self.__workPath)

    def convert(self, dataList, procName, optionsD, workingDir):
        """Enumerate ModBase PDB model files for conversion to mmCIF.

        Args:
            dataList (list): list of model (and alignment) files to convert, each stored as a dict
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
            for modelE in dataList:
                convertedCifFileZ = None
                alignmentFileUsed = None
                modelNameRoot, pdbFileZ, alignFileZ = modelE.name, modelE.model, modelE.alignment
                # First unzip pdb and alignmet files
                pF = self.__fU.uncompress(pdbFileZ, outputDir=workingDir)
                aF = self.__fU.uncompress(alignFileZ, outputDir=workingDir)
                if pF.endswith(".pdb") and aF.endswith(".xml"):
                    cF = os.path.join(workingDir, modelNameRoot + ".cif")
                    # Then attempt the conversion with alignment
                    convertedCifFileZ = self.__convertPdbToCif(procName=procName, pdbFile=pF, alignmentFile=aF, mmCifOutFile=cF, optionsD=optionsD)
                    if convertedCifFileZ:
                        alignmentFileUsed = alignFileZ
                        successList.append(modelE)
                        self.__fU.remove(cF)  # Remove the unzipped cif file
                    else:  # Still append it to success list, since it was at least capable of being converted (this only occurs if it doesn't meet the quality score criteria)
                        successList.append(modelE)
                    # Remove the unzipped pdb and alignment files
                    self.__fU.remove(pF)
                    self.__fU.remove(aF)
                retList.append((modelE, convertedCifFileZ, alignmentFileUsed))

            failList = sorted(set(dataList) - set(successList))
            if failList:
                logger.info("%s returns %d definitions with failures: %r", procName, len(failList), failList)

            logger.debug("%s converted %d/%d models, failures %d", procName, len(retList), len(dataList), len(failList))
        except Exception as e:
            logger.exception("Failing %s for %d data items %s", procName, len(dataList), str(e))

        return successList, retList, diagList

    def __convertPdbToCif(self, procName, pdbFile, alignmentFile, mmCifOutFile, optionsD):
        """Internal method to convert a ModBase PDB file to mmCIF format.

        Args:
            procName (str): worker process name
            pdbFile (str): path to PDB model file to convert
            alignmentFile (str): path to ModBase alignment file for the corresponding model file
            mmCifOutFile (str): path to which to write mmCIF model file
            optionsD (dict): additional options/parameters to use in conversion

        Returns:
            str: path to converted (and gzipped) mmCIF file if successful; otherwise None
        """

        # Running in command line:
        # python ./modbase_pdb_to_cif.py   XP_039771084.1_1.pdb   XP_039771084.1_1.cif   -a XP_039771084.1_1.ali.xml   -r /Volumes/ftp.wwpdb.org/pub/pdb/data/structures/divided/mmCIF

        try:
            mmCifRepo = modbase.Repository(optionsD.get("pdbxRepoPath"))
            mmCifOutFileZ = None
            with open(pdbFile, "r", encoding="utf-8") as fh:
                sF = modbase.read_pdb(fh, mmCifRepo)
            systemWithAlign = sF.get_system(alignmentFile)
            # Create quality metric dictionary to only write models with MPQS > 1.1 and ZDOPE < -1.0
            qMD = {}
            for qM in systemWithAlign.model_groups[0][0].qa_metrics:
                qMD.update({qM.name: qM.value})
            #
            if (float(qMD['MPQS']) > 1.1 and float(qMD['zDOPE']) < -1.0):
                with open(mmCifOutFile, "w", encoding="utf-8") as fh:
                    modelcif.dumper.write(fh, [systemWithAlign], format="mmCIF")
                mmCifOutFileZ = mmCifOutFile + ".gz"
                self.__fU.compress(inpPath=mmCifOutFile, outPath=mmCifOutFileZ)
            else:
                logger.debug("Skipping PDB file %s - Does not meet quality score criteria", pdbFile)
            #
            return mmCifOutFileZ
            #
        except Exception as e:
            logger.exception("%s failing on PDB file %s and alignment file %s, with %s", procName, pdbFile, alignmentFile, str(e))
            # logger.debug("%s failing on PDB file %s and alignment file %s, with %s", procName, pdbFile, alignmentFile, str(e))


class ModBaseModelProcessor(object):
    """Generators and accessors for ModBase model files."""

    def __init__(self, cachePath=None, useCache=False, speciesModelDir=None, speciesName=None, speciesPdbModelFileList=None, pdbxRepoPath=None, **kwargs):
        """Initialize ModBaseModelProcessor object.

        Args:
            cachePath (str): path to species-specific cache file containing list of processed model files.
            useCache (bool): whether to use the existing data cache or re-run conversion process
            speciesModelDir (str): path to species data directory
            speciesName (str): name of the species as it is stored in the ModBaseModelProvider cache
            speciesPdbModelFileList (list): list of the PDB model files to convert
            pdbxRepoPath (str): path to repository containing divided mmCIF structure files (in current wwpdb archive)
                                (e.g., /Volumes/ftp.wwpdb.org/pub/pdb/data/structures/divided/mmCIF)
        """

        try:
            self.__version = __version__
            self.__numProc = kwargs.get("numProc", 2)
            self.__chunkSize = kwargs.get("chunkSize", 10)

            self.__speciesModelDir = speciesModelDir
            self.__speciesName = speciesName
            self.__speciesPdbModelFileList = speciesPdbModelFileList if speciesPdbModelFileList else []
            self.__pdbxRepoPath = pdbxRepoPath

            self.__cachePath = cachePath if cachePath else self.__speciesModelDir
            self.__cacheFormat = kwargs.get("cacheFormat", "pickle")
            self.__workPath = kwargs.get("workPath", self.__speciesModelDir)

            self.__mU = MarshalUtil(workPath=self.__workPath)
            self.__fU = FileUtil(workPath=self.__workPath)

            self.__modelD = self.__reload(useCache=useCache)

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

    def generate(self, updateOnly=False, indent=4):
        """Generate converted mmCIF models from ModBase PDB and alignment files.

        Args:
            updateOnly (bool): only convert new or previously-failed models.  Defaults to False.
            indent (int, optional): json format indent. Defaults to 0.

        Returns:
            bool: True for success or False otherwise
        """
        ok = False
        try:
            tS = time.strftime("%Y %m %d %H:%M:%S", time.localtime())
            mD, failD = self.__convertModBasePdb(numProc=self.__numProc, chunkSize=self.__chunkSize, updateOnly=updateOnly)
            self.__modelD = {
                "version": self.__version,
                "created": tS,
                "species": self.__speciesName,
                "speciesModelDir": self.__speciesModelDir,
                "modelsCif": mD,
                "modelsFailed": failD,      # Will contain failed models as well as models that didn't meet the minimum quality score requirments
            }
            kwargs = {"indent": indent} if self.__cacheFormat == "json" else {"pickleProtocol": 4}
            modelCachePath = self.__getModelCachePath()
            ok = self.__mU.doExport(modelCachePath, self.__modelD, fmt=self.__cacheFormat, **kwargs)
            logger.info("Wrote %r status %r", modelCachePath, ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return ok

    def getModelD(self):
        return self.__modelD

    def reload(self, useCache=True):
        self.__modelD = self.__reload(useCache=useCache)
        return self.__modelD is not None

    def __reload(self, useCache=True):
        """Reload from the current cache directory."""
        try:
            modelCachePath = self.__getModelCachePath()
            tS = time.strftime("%Y %m %d %H:%M:%S", time.localtime())
            modelD = {"version": self.__version, "created": tS, "modelsCif": {}}
            logger.debug("useCache %r modelCachePath %r", useCache, modelCachePath)
            #
            if useCache and self.__mU.exists(modelCachePath):
                modelD = self.__mU.doImport(modelCachePath, fmt=self.__cacheFormat)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        #
        return modelD

    def getCachePath(self):
        return self.__cachePath

    def __getModelCachePath(self):
        ext = "pic" if self.__cacheFormat == "pickle" else "json"
        speciesNameNoSpace = self.__speciesName.replace(" ", "_")
        pth = os.path.join(self.__cachePath, speciesNameNoSpace + "-model-data." + ext)
        return pth

    def __convertModBasePdb(self, numProc=2, chunkSize=10, updateOnly=False):
        """Prepare multiprocessor queue and workers for converting all ModBase PDB model files to mmCIF, using alignment files.

        Args:
            numProc (int, optional): number of processes to use. Defaults to 2.
            chunkSize (int, optional): incremental chunk size used for distribute work processes. Defaults to 10.
            updateOnly (bool): only convert new or previously-failed models.  Defaults to False.

        Returns:
            mD (dict): dictionary of successfully-converted mmCIF models, in the following structure:
                       {modelNameRoot: {"model": pathToGzippedCifFile, "alignment": pathToCompressedAlignmentFile}, ...}
            failD (dict): dictionary of ModBase models for which conversion failed, in the following structure:
                          {modelNameRoot: {"model": pathToCompressedPDBModelFile, "alignment": pathToCompressedAlignmentFile}, ...}
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
        for pFZ in pdbModelFileList:
            modelNameRoot = self.__fU.getFileName(pFZ).split(".pdb.xz")[0]
            if modelNameRoot not in mD:
                aFZ = os.path.join(self.__speciesModelDir, "alignment", modelNameRoot + ".ali.xml.xz")
                modelE = ModBaseEntry(name=modelNameRoot, model=pFZ, alignment=aFZ)
                modelList.append(modelE)

        logger.info("Starting with %d models, numProc %d, updateOnly (%r)", len(modelList), self.__numProc, updateOnly)
        #
        rWorker = ModBaseModelWorker(workPath=self.__workPath)
        mpu = MultiProcUtil(verbose=True)
        optD = {"pdbxRepoPath": self.__pdbxRepoPath}
        mpu.setOptions(optD)
        mpu.set(workerObj=rWorker, workerMethod="convert")
        mpu.setWorkingDir(workingDir=self.__workPath)
        ok, failList, resultList, _ = mpu.runMulti(dataList=modelList, numProc=numProc, numResults=1, chunkSize=chunkSize)
        if failList:
            logger.info("mmCIF conversion failures (%d): %r", len(failList), failList)
        #
        for (model, convertedCifFileZ, alignmentFileUsed) in resultList[0]:
            if convertedCifFileZ:
                mD[model.name] = {"model": convertedCifFileZ, "alignment": alignmentFileUsed}
            else:
                failD[model.name] = {"model": model.model, "alignment": model.alignment}
        #
        logger.info("Completed with multi-proc status %r, failures %r, total models with data (%d)", ok, len(failList), len(mD))
        return mD, failD

    # def convertJsonPickle(self, fmt1="json", fmt2="pickle"):
    #     # Keep method name as "convertJsonPickle" instead of "convert", to avoid complications from using the other "convert" for MPU method above
    #     modelCachePath = self.__getModelCachePath(fmt=fmt1)
    #     self.__modelD = self.__mU.doImport(modelCachePath, fmt=fmt1)
    #     #
    #     modelCachePath = self.__getModelCachePath(fmt=fmt2)
    #     ok = self.__mU.doExport(modelCachePath, self.__modelD, fmt=fmt2, pickleProtocol=4)
    #     return ok
