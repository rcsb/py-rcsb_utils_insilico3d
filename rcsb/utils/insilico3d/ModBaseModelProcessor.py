##
# File:    ModBaseModelProcessor.py
# Author:  Dennis Piehl
# Date:    16-Sep-2021
#
# Update:
#
#
##

"""
Processors for converting ModBase PDB models into mmCIF format.

"""
# pylint: disable=fixme

import datetime
import logging
import os.path
import time
import collections

# from rcsb.utils.insilico3d import __version__
# import rcsb.utils.modbase_utils.modbase_pdb_to_cif as modbase
import modbase_pdb_to_cif as modbase
# from ModBaseModelProvider import ModBaseModelProvider

# from rcsb.utils.insilico3d.ModBaseModelProvider import ModBaseModelProvider
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.StashableBase import StashableBase
from rcsb.utils.multiproc.MultiProcUtil import MultiProcUtil

logger = logging.getLogger(__name__)

ModBaseEntry = collections.namedtuple("ModBaseEntry", ["name", "model", "alignment"])

class ModBaseModelWorker(object):
    """A skeleton class that implements the interface expected by the multiprocessing
    for converting ModBase PDB files to mmCIF files.
    """

    def __init__(self, **kwargs):
        _ = kwargs
        self.__fU = FileUtil()

    def convert(self, dataList, procName, optionsD, workingDir):
        """Enumerate ModBase PDB model files for conversion to mmCIF.

        Args:
            dataList (list): list of model (and alignment) files to convert, each stored as a dict
            procName (str): worker process name (?)
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
                modelNameRoot, pdbFileZ, alignFileZ = modelE.name, modelE.model, modelE.alignment
                # fnRoot = self.__fU.getFileName(pdbFileZ).split(".pdb")[0]
                # print('\n',modelNameRoot,'\n')
                # First unzip pdb and alignmet files
                pF = self.__fU.uncompress(pdbFileZ, outputDir=workingDir)
                aF = self.__fU.uncompress(alignFileZ, outputDir=workingDir)
                if pF.endswith(".pdb") and aF.endswith(".xml"):
                    # fnRoot = self.__fU.getFileName(pF).split(".pdb")[0]
                    # print('\n',fnRoot,'\n')
                    cF = os.path.join(workingDir, modelNameRoot + ".cif")

                    # Then perform conversion
                    # ok = self.__convertPdbToCif(procName=procName, pdbFile=pF, alignmentFile=aF, mmCifOutFile=cF)
                    try:
                        with open(pF) as fh:
                            sF = modbase.read_pdb(fh)
                        with open(cF, "w") as fh:
                            sF.write_mmcif(fh, aF)
                        # retList.append((modelE, True))
                        # Next zip the cif file
                        cifFileZ = cF + ".gz"
                        ok = self.__fU.compress(inpPath=cF, outPath=cifFileZ)
                        if ok:
                            successList.append(modelE)
                            convertedCifFileZ = cifFileZ
                    except Exception as e:
                        logger.exception("Failing for data items %s, with %s", modelE, str(e))

                    # Last remove the unzipped pdb, alignment and cif files
                    ok = self.__fU.remove(pF)
                    ok = self.__fU.remove(aF)
                    ok = self.__fU.remove(cF)

                retList.append((modelE, convertedCifFileZ))

            failList = sorted(set(dataList) - set(successList))
            if failList:
                logger.info("%s returns %d definitions with failures: %r", procName, len(failList), failList)

            logger.debug("%s converted %d/%d models, failures %d", procName, len(retList), len(dataList), len(failList))
        except Exception as e:
            logger.exception("Failing %s for %d data items %s", procName, len(dataList), str(e))

        return successList, retList, diagList

    def __convertPdbToCif(self, procName, pdbFile, alignmentFile, mmCifOutFile):
        """Internal method to convert a ModBase PDB file to mmCIF format, using the alignment information.

        Args:
            procName (str): worker process name (?)
            pdbFile (str): path to PDB model file to convert
            alignmentFile (str): path to ModBase alignment file for the corresponding model file
            mmCifOutFile (str): path to which to write mmCIF model file

        Returns:
            bool: True for success or False otherwise
        """
        try:
            with open(pdbFile) as fh:
                sF = modbase.read_pdb(fh)
            with open(mmCifOutFile, "w") as fh:
                sF.write_mmcif(fh, alignmentFile)
            return True
        except Exception as e:
            logger.exception("%s failing on file %s with %s", procName, pdbFile, str(e))
            return False


class ModBaseModelProcessor(object):
    """Generators and accessors for ModBase model files."""

    def __init__(self, cachePath=None, useCache=False, speciesModelDir=None, speciesPdbModelFileList=None, **kwargs):
        try:
            # self.__version = __version__
            self.__version = "0.12"
            self.__numProc = kwargs.get("numProc", 2)
            self.__chunkSize = kwargs.get("chunkSize", 10)

            self.__speciesModelDir = speciesModelDir
            self.__speciesName = self.__speciesModelDir.split("/")[-1]

            self.__cachePath = cachePath
            # print("modBaseModelProcessor__cachePath", self.__cachePath)

            # print("modBaseModelProcessor__speciesModelDir", self.__speciesModelDir)
            self.__speciesPdbModelFileList = speciesPdbModelFileList if speciesPdbModelFileList else []

            # This was commented out in previous exmaple usage
            # useCache = kwargs.get("useCache", True)

            #  - Configuration for stash services -
            #    Local target directory name to be stashed.  (subdir of dirPath)
            self.__mU = MarshalUtil(workPath=self.__speciesModelDir)
            self.__fU = FileUtil(workPath=self.__speciesModelDir)
            
            ok = self.reload(fmt="pickle", useCache=useCache)

        except Exception as e:
            # May occur if 'speciesModelDir' is not specified in instantiation
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

    def getEntries(self):
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
        """Generate converted mmCIF models from ModBase PDB and alignment files.

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
            mD, failD = self.__convertModBasePdb(numProc=self.__numProc, chunkSize=self.__chunkSize, updateOnly=updateOnly)
            self.__modelD = {
                "version": self.__version,
                "created": tS,
                "species": self.__speciesName,
                "speciesModelDir": self.__speciesModelDir,
                "modelsCif": mD,
                "modelFailed": failD,
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
                # if fmt != "pickle":
                #     for _, cD in modelD["modelsCif"].items():
                #         cD["nearestNeighbors"] = [LigandTargetInstance(*neighbor) for neighbor in nD["nearestNeighbors"]]
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        #
        return modelD

    def __getModelCachePath(self, fmt="pickle"):
        ext = "pic" if fmt == "pickle" else "json"
        pth = os.path.join(self.__speciesModelDir, self.__speciesName + "-model-data." + ext)
        return pth

    # def getEntries(self):
    #     return []

    def __convertModBasePdb(self, numProc=2, chunkSize=10, updateOnly=False):
        """Prepare multiprocessor queue and workers for converting all ModBase PDB model files to mmCIF, using alignment files.

        Args:
            numProc (int, optional): number of processes to use. Defaults to 2.
            chunkSize (int, optional): incremental chunk size used for distribute work processes. Defaults to 10.
            updateOnly (bool): only convert new or previously-failed models.  Defaults to False.

        Returns:
            (dict): dictionary of successfully-converted mmCIF models, in the following structure:
                    {modelNameRoot: pathToGzippedCifFile, ...}
            (dict): dictionary of ModBase models for which conversion failed, in the following structure:
                    {modelNameRoot: {"model": pathToCompressedPDBModelFile, "alignment": pathToCompressedAlignmentFile}, ...}
        """
        # updateDate = datetime.datetime.now().isoformat()
        mD = {}
        failD = {}
        exD = {}
        # fU = FileUtil()
        #
        # updateOnly - will reuse any existing data loaded when this is instantiated
        #              otherwise the cache context is cleared before the calculation.
        if updateOnly:
            exD = {k: True for k in self.getEntries()}
            logger.info("Reusing (%d) already converted mmCIF models", len(exD))
            mD = self.__modelD["modelsCif"] if "modelsCif" in self.__modelD else {}
        #
        pdbModelFileList = self.__speciesPdbModelFileList
        modelList = []
        for pFZ in pdbModelFileList:
            modelNameRoot = self.__fU.getFileName(pFZ).split(".pdb.xz")[0]
            if modelNameRoot not in mD:
                aFZ = os.path.join(self.__speciesModelDir,  "alignment", modelNameRoot + ".ali.xml.xz")
                modelE = ModBaseEntry(name=modelNameRoot, model=pFZ, alignment=aFZ)
                modelList.append(modelE)

        logger.info("Starting with %d models, numProc %d, updateOnly (%r)", len(modelList), self.__numProc, updateOnly)
        #
        rWorker = ModBaseModelWorker()
        mpu = MultiProcUtil(verbose=True)
        # optD = {"updateDate": updateDate}  # Not needed, unless storing cache...
        # mpu.setOptions(optD)  # Not needed, unless storing cache...
        mpu.set(workerObj=rWorker, workerMethod="convert")
        mpu.setWorkingDir(workingDir=self.__speciesModelDir)
        ok, failList, resultList, _ = mpu.runMulti(dataList=modelList, numProc=numProc, numResults=1, chunkSize=chunkSize)
        # The resultList from above is a list of lists (# of sublists determined by numResults, and the return statement in convert method above)
        if failList:
            logger.info("mmCIF conversion failures (%d): %r", len(failList), failList)
        #
        for (model, convertedCifFileZ) in resultList[0]:
            if convertedCifFileZ:
                mD[model.name] = convertedCifFileZ
            else:
                failD[model.name] = {"model": model.model, "alignment": model.alignment}
            # mD[(entry.model, entry.alignment)] = eResult
        #
        logger.info("Completed with multi-proc status %r, failures %r, total models with data (%d)", ok, len(failList), len(mD))
        # print("\n\nSUCCESSES")
        # for k,v in mD.items():
        #     print(k,'\t', v)
        # print("\n\nFAILURES")
        # for k,v in failD.items():
        #     print(k,'\n\t', v)
        # for i in failList:
        #     print(i)
        return mD, failD

