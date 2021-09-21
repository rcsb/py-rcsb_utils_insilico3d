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


class ModBaseModelWorker(object):
    """A skeleton class that implements the interface expected by the multiprocessing
    for converting ModBase PDB files to mmCIF files.
    """

    def __init__(self, **kwargs):
        # self.__mP = ModBaseModelProvider()  # Is this the proper place/way to initialize this class? Or initialize in ModBaseModelProcessor class below only?
        # self.__dirPath = self.__mP.getCacheDirPath()
        # self.__cachePath = kwargs.get("cachePath", "./CACHE-insilico3d-models")
        # self.__dirPath = os.path.join(self.__cachePath, "ModBase")
        self.__speciesModelDir = kwargs.get("speciesModelDir", "./CACHE-insilico3d-models/ModBase")
        # self.__workingPath = kwargs.get("speciesDir", self.__dirPath)
        _ = kwargs
        # print("modBaseModelWorker__dirPath", self.__dirPath)
        print("modBaseModelWorker__speciesModelDir", self.__speciesModelDir)

        self.__mU = MarshalUtil(workPath=self.__speciesModelDir)
        # self.__mU = MarshalUtil()  # Should workPath be defined up here, or down below for each species directory?
        self.__fU = FileUtil(workPath=self.__speciesModelDir)

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

        # Ask John about failList in NeighborInteractionProvider, and if it's even needed up here or only down below
        try:
            modelDataList = [d["model"] for d in dataList]
            for modelD in dataList:
                pdbFileZ, alignFileZ = modelD["model"], modelD["alignment"]
                print("pdbFileZ", pdbFileZ)
                print("workingDir", workingDir)
                print("dirPath", self.__speciesModelDir)
                # exit()
                # First unzip pdb and alignmet files
                pF = self.__fU.uncompress(pdbFileZ, outputDir=self.__speciesModelDir)
                # pF = os.path.join()
                aF = self.__fU.uncompress(alignFileZ, outputDir=self.__speciesModelDir)
                if pF.endswith(".pdb") and aF.endswith(".xml"):
                    fnRoot = self.__fU.getFileName(pF).split(".pdb")[0]
                    cF = os.path.join(self.__speciesModelDir, fnRoot + ".cif")
                    # Then perform conversion
                    print("pF", pF)
                    print("aF", aF)
                    print("cF", cF)
                    # exit()
                    ok = self.__convertPdbToCif(procName=procName, pdbFile=pF, alignmentFile=aF, mmCifOutFile=cF)
                    retList.append((modelD, ok))
                    # Next zip the cif file
                    cifFileZ = cF + ".gz"
                    ok = self.__fU.compress(inpPath=cF, outPath=cifFileZ)
                    # Last remove the unzipped pdb, alignment and cif files
                    ok = self.__fU.remove(pF)
                    ok = self.__fU.remove(aF)
                    ok = self.__fU.remove(cF)
                else:
                    retList.append((modelD, False))

            successList = sorted(set(modelDataList) - set(failList))  # this is subtracting the failList from the list of pdb model paths only (not the model and alignment pair dicts)
            if failList:
                logger.info("%s returns %d definitions with failures: %r", procName, len(failList), failList)

            logger.debug("%s converted %d/%d entries failures %d", procName, len(retList), len(dataList), len(failList))
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
                # print(dir(sF))
                # exit()
            with open(mmCifOutFile, "w") as fh:
                sF.write_mmcif(fh, alignmentFile)
            return True

        except Exception as e:
            logger.exception("%s failing on file %s with %s", procName, pdbFile, str(e))
            return False


# class ModBaseModelProcessor(StashableBase):
    # What is this stashable base? Necessary for this application?
class ModBaseModelProcessor(object):
    """Generators and accessors for ModBase model files."""

    def __init__(self, cachePath=None, useCache=False, speciesModelDir=None, speciesPdbModelFileList=None, **kwargs):
        # Set default useCache=True ...?

        # modelGroup: name of the species directory of model files to process

        # TODO: Determine proper way and place in code to initialize ModBaseModelProvider() class
        # TODO: Figure out declaration of cachePath and dirPath
        # TODO: Also determine necessity of using super(class) below
        # TODO: Determine appropriate default setting to use for "useCache"

        try:
            # self.__mP = ModBaseModelProvider(useCache=True)  # Is this the proper place/way to initialize this class?

            # self.__version = __version__
            # self.__fileLimit = kwargs.get("fileLimit", None)
            self.__numProc = kwargs.get("numProc", 2)
            self.__chunkSize = kwargs.get("chunkSize", 10)

            self.__speciesModelDir = speciesModelDir
            self.__speciesName = self.__speciesModelDir.split("/")[-1]

            self.__cachePath = cachePath
            print("modBaseModelProcessor__cachePath", self.__cachePath)
            # if cachePath:
            #     self.__cachePath = cachePath
            # else:
            #     self.__cachePath = self.__mP.getCacheDirPath()

            # self.__dirPath = self.__speciesModelDir
            print("modBaseModelProcessor__speciesModelDir", self.__speciesModelDir)
            self.__speciesPdbModelFileList = speciesPdbModelFileList if speciesPdbModelFileList else []

            # super(ModBaseModelProcessor, self).__init__(self.__cachePath, [self.__speciesName])

            # This was commented out in previous exmaple usage...if deemed OK to set default to True in method(parameter) declaration above, then keep this here too
            # useCache = kwargs.get("useCache", True)

            #  - Configuration for stash services -
            #    Local target directory name to be stashed.  (subdir of dirPath)

            self.__mU = MarshalUtil(workPath=self.__speciesModelDir)

            # self.__modelCifD = self.__reload(fmt="pickle", useCache=useCache)

        except Exception as e:
            # May occur if 'speciesModelDir' is not specified in instantiation
            logger.exception("Failing with %s", str(e))

    # def testCache(self, minCount=1):
    #     try:
    #         if minCount == 0:
    #             return True
    #         if self.__modelCifD and minCount and "entries" in self.__modelCifD and len(self.__modelCifD["entries"]) >= minCount:
    #             logger.info("Model meatadata for (%d) entries created %r version %r", len(self.__modelCifD["entries"]), self.__modelCifD["created"], self.__modelCifD["version"])
    #             return True
    #     except Exception:
    #         pass
    #     return False

    def generate(self, updateOnly=False):
        """Generate and export metadata for all ModBase mmCIF model files.

        Args:
            updateOnly (bool): only extract and update metadata for new entries.  Defaults to False.
            fmt (str, optional): export file format. Defaults to "pickle".
            indent (int, optional): json format indent. Defaults to 0.

        Returns:
            bool: True for success or False otherwise
        """
        ok = False
        try:
            # tS = time.strftime("%Y %m %d %H:%M:%S", time.localtime())
            tD = self.__convertModBasePdb(numProc=self.__numProc, chunkSize=self.__chunkSize, updateOnly=updateOnly)
            # self.__modelCifD = {"version": self.__version, "created": tS, "entries": tD, "species": self.__speciesName}
            # kwargs = {"indent": indent} if fmt == "json" else {"pickleProtocol": 4}
            # targetFilePath = self.__getTargetFilePath(fmt=fmt)
            # ok = self.__mU.doExport(targetFilePath, self.__modelCifD, fmt=fmt, **kwargs)
            # logger.info("Wrote %r status %r", targetFilePath, ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return ok

    # def reload(self, fmt="pickle"):
    #     self.__modelCifD = self.__reload(fmt=fmt, useCache=True)
    #     return self.__modelCifD is not None

    # def __reload(self, fmt="pickle", useCache=True):
    #     """Reload from the current cache directory."""
    #     try:
    #         targetFilePath = self.__getTargetFilePath(fmt=fmt)
    #         tS = time.strftime("%Y %m %d %H:%M:%S", time.localtime())
    #         modelCifD = {"version": self.__version, "created": tS, "entries": {}}
    #         logger.debug("useCache %r targetFilePath %r", useCache, targetFilePath)
    #         #
    #         if useCache and self.__mU.exists(targetFilePath):
    #             modelCifD = self.__mU.doImport(targetFilePath, fmt=fmt)
    #             # if fmt != "pickle":
    #             #     for _, nD in modelCifD["entries"].items():
    #             #         nD["nearestNeighbors"] = [LigandTargetInstance(*neighbor) for neighbor in nD["nearestNeighbors"]]
    #     except Exception as e:
    #         logger.exception("Failing with %s", str(e))
    #     #
    #     return modelCifD

    # def __getTargetFilePath(self, fmt="pickle"):
    #     ext = "pic" if fmt == "pickle" else "json"
    #     pth = os.path.join(self.__dirPath, self.__speciesName + "-model-data." + ext)
    #     return pth

    # def getEntries(self):
    #     return []

    def __convertModBasePdb(self, numProc=2, chunkSize=10, updateOnly=False):
        """Prepare multiprocessor queue and workers for extracting metadata for all ModBase mmCIF model files.
        Called by generate() method -> __convertModBasePdb() -> Specifies use of workerMethod="convert"

        Args:
            numProc (int, optional): number of processes to use. Defaults to 2.
            chunkSize (int, optional): incremental chunk size used for distribute work processes. Defaults to 10.
            updateOnly (bool): only extract and update metadata for new entries.  Defaults to False.

        Returns:
            (dict): dictionary of extracted metadata for the provided model file, in the following structure:
                    {dataContainerEntryName: [{categoryName: [{attributeName: ...}], ...}], ...}
        """
        updateDate = datetime.datetime.now().isoformat()
        rD = {}
        # exD = {}
        fU = FileUtil()
        #
        # updateOnly - will reuse any existing data loaded when this is instantiated
        #              otherwise the cache context is cleared before the calculation.
        # if updateOnly:
        #     existingDataList = self.__mP.getModelFileList(inputPathList=[self.__speciesModelDir])
        #     logger.info("Reusing (%d) entries", len(existingDataList))
        #     rD = 
        #
        pdbModelFileList = self.__speciesPdbModelFileList
        modelList = []
        for pFZ in pdbModelFileList:
            pFZNameRoot = fU.getFileName(pFZ).split(".pdb.xz")[0]
            aFZ = os.path.join(self.__speciesModelDir,  "alignment", pFZNameRoot + ".ali.xml.xz")
            modelD = {"model": pFZ, "alignment": aFZ}
            modelList.append(modelD)

        logger.info("Starting with %d entries numProc %d updateOnly (%r)", len(modelList), self.__numProc, updateOnly)
        #
        rWorker = ModBaseModelWorker(speciesModelDir=self.__speciesModelDir)
        mpu = MultiProcUtil(verbose=True)
        optD = {"updateDate": updateDate}
        mpu.setOptions(optD)
        mpu.set(workerObj=rWorker, workerMethod="convert")
        ok, failList, resultList, _ = mpu.runMulti(dataList=modelList, numProc=numProc, numResults=1, chunkSize=chunkSize)
        if failList:
            logger.info("mmCIF conversion failures (%d): %r", len(failList), failList)
        #
        print(resultList[0])
        # exit()
        for (entry, eResult) in resultList[0]:
            entryModel = entry["model"]
            rD[entryModel] = eResult
        #
        logger.info("Completed with multi-proc status %r failures %r total entries with data (%d)", ok, len(failList), len(rD))
        return rD

