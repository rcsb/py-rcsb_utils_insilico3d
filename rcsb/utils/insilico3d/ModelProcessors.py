##
# File:    ModelProcessors.py
# Author:  Dennis Piehl
# Date:    10-Mar-2022
#
# Update:
#   10-Mar-2022  dwp Start class
#
# To Do:
# - pylint: disable=fixme
# - Add mkdssp calculation
# - Add check that converted files are consistent with mmCIF dictionaries
##

"""
Common worker methods for processing and handling computed model files.

"""

__docformat__ = "google en"
__author__ = "Dennis Piehl"
__email__ = "dennis.piehl@rcsb.org"
__license__ = "Apache 2.0"

import logging
import os.path
import time
import collections

import rcsb.utils.modbase_utils.modbase_pdb_to_cif as modbase
from rcsb.utils.insilico3d import __version__
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.multiproc.MultiProcUtil import MultiProcUtil

logger = logging.getLogger(__name__)


class ModelWorker(object):
    """A skeleton class that implements the interface expected by the multiprocessing
    for converting ModBase PDB files to mmCIF files.
    """

    def __init__(self, **kwargs):
        _ = kwargs
        self.__workPath = kwargs.get("workPath", None)
        self.__fU = FileUtil(workPath=self.__workPath)
        self.__mU = MarshalUtil(workPath=self.__workPath)

    def reorganize(self, dataList, procName, optionsD, workingDir):
        """Enumerate and reorganize/rename model files.

        Args:
            dataList (list): list of model files to work on
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
            modelSourcePrefix = optionsD.get("modelSourcePrefix")  # e.g., "af" or "ma"
            destModelDir = optionsD.get("destModelDir")  # self.__dividedDataPath, which is:  os.path.join(self.__cachePath, "computed-models")

            for modelPath in dataList:
                modelFileOut = None
                modelFileNameIn = self.__fU.getFileName(modelPath)
                containerList = self.__mU.doImport(modelPath, fmt="mmcif")
                if len(containerList) > 1:
                    # Expecting all computed models to have one container per file. When this becomes no longer the case, update this to handle it accordingly.
                    logger.error("Skipping - model file %s has more than one container (%d)", modelPath, len(containerList))
                    continue
                modelEntryId = containerList[0].getObj("entry").getValue("id", 0)
                internalModelId = modelSourcePrefix+"_"+modelEntryId
                internalModelName = internalModelId+".cif.gz"
                modelName = self.__fU.getFileName(modelPath)
                uniProtID = modelName.split(".cif.gz")[0].split("-")[1]
                first2 = uniProtID[0:2]
                mid2 = uniProtID[2:4]
                last2 = uniProtID[4:6]
                destDir = os.path.join(destModelDir, first2, mid2, last2)
                if not self.__fU.exists(destDir):
                    self.__fU.mkdir(destDir)
                modelFileOut = os.path.join(destDir, internalModelName)
                #
                ok = self.__fU.put(modelPath, modelFileOut)  # Copies files (only use for testing)
                # ok = self.__fU.replace(modelPath, modelFileOut)  # Moves files (use for production)
                #
                if ok:
                    successList.append(modelPath)
                #
                retList.append((modelFileNameIn, internalModelId, modelFileOut))
            failList = sorted(set(dataList) - set(successList))
            if failList:
                logger.info("%s returns %d definitions with failures: %r", procName, len(failList), failList)

            logger.debug("%s converted %d/%d models, failures %d", procName, len(retList), len(dataList), len(failList))
        except Exception as e:
            logger.exception("Failing %s for %d data items %s", procName, len(dataList), str(e))

        return successList, retList, diagList


class ModelReorganizer(object):
    """Generators and accessors for model files."""

    def __init__(self, cachePath=None, useCache=False, speciesD=None, **kwargs):
        """Initialize ModelReorganizer object.

        Args:
            cachePath (str): path to species-specific cache file containing list of processed model files;
                             should be provided!; else defaults to some temporary cache file.
            useCache (bool): whether to use the existing data cache or re-run conversion process
            speciesD (dict): dictionary containing the following necessary key:value pairs for processing:
                             "speciesModelDir": path to species data directory, will be used as working directory
                             "lastModified": last modified date of the downloaded species archive tarball
                             "speciesName": name of the species as it is stored in the ModBaseModelProvider cache
                             "speciesModelFileList": list of the model files to process
                             "modelSourcePrefix": internal identifier prefix to use for provided model source type (e.g., "af", "ma")
                             "destModelDir": base destination directory into which to reorganize model files
        """

        try:
            # self.__version = __version__
            self.__numProc = kwargs.get("numProc", 4)
            self.__chunkSize = kwargs.get("chunkSize", 20)

            self.__speciesModelDir = speciesD.get("speciesModelDir")
            self.__workPath = self.__speciesModelDir
            self.__speciesModDate = speciesD.get("lastModified", None)
            self.__speciesName = speciesD.get("speciesName")
            self.__speciesModelFileList = speciesD.get("speciesModelFileList", [])
            self.__modelSourcePrefix = speciesD.get("modelSourcePrefix")
            self.__destModelDir = speciesD.get("destModelDir")

            self.__cachePath = cachePath if cachePath else self.__getModelCachePath()

            self.__mU = MarshalUtil(workPath=self.__workPath)
            self.__fU = FileUtil(workPath=self.__workPath)

        except Exception as e:
            logger.exception("Failing with %s", str(e))

    def getCachePath(self):
        return self.__cachePath

    # def __getModelCachePath(self, fmt="pickle"):
    def __getModelCachePath(self, fmt="json"):
        ext = "pic" if fmt == "pickle" else "json"
        # speciesNameNoSpace = self.__speciesName.replace(" ", "_")
        # pth = os.path.join(self.__speciesModelDir, speciesNameNoSpace + "-model-files-cache." + ext)
        pth = os.path.join(self.__speciesModelDir, "species-model-files-cache." + ext)
        return pth

    # def reorganize(self, fmt="pickle", indent=0):
    def reorganize(self, fmt="json", indent=4):
        """Move model files from organism-wide model listing to hashed directory structure constructed
        from the 6-character UniProt ID (e.g., "P52078" will be moved to "./P5/20/78"). Also rename files
        to follow internal identifier naming convention

        Args:
            fmt (str, optional): export file format. Defaults to "pickle".
            indent (int, optional): json format indent. Defaults to 0.

        Returns:
            bool: True for success or False otherwise
        """
        ok = False
        try:
            # tS = time.strftime("%Y %m %d %H:%M:%S", time.localtime())
            mD, failD = self.__reorganizeModels(numProc=self.__numProc, chunkSize=self.__chunkSize)
            if len(failD) > 0:
                logger.error("Failed to move %d model files.", len(failD))
            # self.__modelD = {
            #     # "version": self.__version,
            #     # "created": tS,
            #     # "species": self.__speciesName,
            #     # "archiveModDate": self.__speciesModDate,
            #     "speciesModelDir": self.__speciesModelDir,
            #     "modelsCif": mD,
            #     "modelsFailed": failD,
            # }
            kwargs = {"indent": indent} if fmt == "json" else {"pickleProtocol": 4}
            ok = self.__mU.doExport(self.__cachePath, mD, fmt=fmt, **kwargs)
            logger.info("Wrote %r status %r", self.__cachePath, ok)
            # modelCachePath = self.__getModelCachePath(fmt=fmt)
            # modelCachePath = self.__cachePath
            # ok = self.__mU.doExport(modelCachePath, self.__modelD, fmt=fmt, **kwargs)
            # ok = self.__mU.doExport(modelCachePath[:-5]+"_failed.json", failD, fmt=fmt, **kwargs)
            # logger.info("Wrote %r status %r", modelCachePath, ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return ok

    def __reorganizeModels(self, numProc=2, chunkSize=10):
        """Prepare multiprocessor queue and workers for reorganizing and renaming all downloaded model files into
        structured directory tree using internal identifiers.

        Args:
            numProc (int, optional): number of processes to use. Defaults to 2.
            chunkSize (int, optional): incremental chunk size used for distribute work processes. Defaults to 10.

        Returns:
            mD (dict): dictionary of successfully processed models, in the following structure:
                       {modelNameRoot: {"model": newPathToModelFile}, ...}
            failD (dict): dictionary of models for which processing failed, in the following structure:
                          {modelNameRoot: {"model": originalPathToModelFile}, ...}
        """
        mD = {}
        failD = {}
        # exD = {}
        #
        # modelFileList = self.__speciesModelFileList
        modelFileList = self.__speciesModelFileList[0:100]
        # print("modelFileList: ", modelFileList)
        #
        logger.info("Starting with %d models, numProc %d", len(modelFileList), self.__numProc)
        #
        rWorker = ModelWorker(workPath=self.__workPath)
        mpu = MultiProcUtil(verbose=True)
        optD = {"species": self.__speciesName, "speciesModDate": self.__speciesModDate, "modelSourcePrefix": self.__modelSourcePrefix, "destModelDir": self.__destModelDir}
        mpu.setOptions(optD)
        mpu.set(workerObj=rWorker, workerMethod="reorganize")
        mpu.setWorkingDir(workingDir=self.__workPath)
        ok, failList, resultList, _ = mpu.runMulti(dataList=modelFileList, numProc=numProc, numResults=1, chunkSize=chunkSize)
        if failList:
            logger.info("model file failures (%d): %r", len(failList), failList)
        #
        for (modelFileNameIn, internalModelId, modelFileOut) in resultList[0]:
            if modelFileOut:
                # internalModelId = self.__fU.getFileName(modelFileOut).split(".cif.gz")[0]   # Provide this as output from MPU method above, don't extract it here
                mD[internalModelId] = {"sourceModelName": modelFileNameIn, "modelPath": modelFileOut}
            else:
                failD[internalModelId] = {"sourceModelName": modelFileNameIn}
        #
        logger.info("Completed with multi-proc status %r, failures %r, total models with data (%d)", ok, len(failList), len(mD))
        return mD, failD

    def convertJsonToPickle(self, fmt1="json", fmt2="pickle"):
        modelCachePath = self.__getModelCachePath(fmt=fmt1)
        self.__modelD = self.__mU.doImport(modelCachePath, fmt=fmt1)
        #
        modelCachePath = self.__getModelCachePath(fmt=fmt2)
        ok = self.__mU.doExport(modelCachePath, self.__modelD, fmt=fmt2, pickleProtocol=4)
        return ok
