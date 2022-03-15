##
# File:    ModelProcessors.py
# Author:  Dennis Piehl
# Date:    10-Mar-2022
#
# Updates:
#
# To Do:
# - pylint: disable=fixme
# - Add mkdssp calculation
# - Add check that model files are consistent with mmCIF dictionaries
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

# import rcsb.utils.modbase_utils.modbase_pdb_to_cif as modbase
# from rcsb.utils.insilico3d import __version__
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.multiproc.MultiProcUtil import MultiProcUtil

logger = logging.getLogger(__name__)


class ModelWorker(object):
    """A skeleton class that implements the interface expected by the multiprocessing
    for working on model files.
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
            destModelDir = optionsD.get("destModelDir")  # base path for all computed models (i.e., "computed-models")
            modelSourceUrlMapD = optionsD.get("modelSourceUrlMapD")  # mapping between source model filenames and accession URLs
            keepSource = optionsD.get("keepSource", False)  # whether to copy files over (instead of moving them)

            for modelPath in dataList:
                modelFileOut = None
                modelFileNameIn = self.__fU.getFileName(modelPath)
                containerList = self.__mU.doImport(modelPath, fmt="mmcif")
                if len(containerList) > 1:
                    # Expecting all computed models to have one container per file. When this becomes no longer the case, update this to handle it accordingly.
                    logger.error("Skipping - model file %s has more than one container (%d)", modelPath, len(containerList))
                    continue
                modelEntryId = containerList[0].getObj("entry").getValue("id", 0)
                internalModelId = modelSourcePrefix + "_" + modelEntryId
                internalModelName = internalModelId + ".cif.gz"
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
                sourceModelUrl = modelSourceUrlMapD[modelFileNameIn]
                #
                if keepSource:
                    ok = self.__fU.put(modelPath, modelFileOut)  # Copies files (only use for testing)
                else:
                    ok = self.__fU.replace(modelPath, modelFileOut)  # Moves files (use for production)
                #
                if ok:
                    successList.append(modelPath)
                #
                retList.append((modelFileNameIn, internalModelId, modelFileOut, sourceModelUrl))
            failList = sorted(set(dataList) - set(successList))
            if failList:
                logger.info("%s returns %d definitions with failures: %r", procName, len(failList), failList)

            logger.debug("%s processed %d/%d models, failures %d", procName, len(retList), len(dataList), len(failList))
        except Exception as e:
            logger.exception("Failing %s for %d data items %s", procName, len(dataList), str(e))

        return successList, retList, diagList


class ModelReorganizer(object):
    """Generators and accessors for model files."""

    def __init__(self, cachePath=None, modelD=None, keepSource=False, **kwargs):
        """Initialize ModelReorganizer object.

        Args:
            cachePath (str): path to species-specific cache file containing list of processed model files;
                             should be provided!; else defaults to some temporary cache file.
            useCache (bool): whether to use the existing data cache or re-run conversion process
            modelD (dict): dictionary containing the following necessary key:value pairs for processing:
                             "modelDir": path to species data directory, will be used as working directory
                             "modelFileList": list of the model files to process
                             "modelSourcePrefix": internal identifier prefix to use for provided model source type (e.g., "af", "ma")
                             "destModelDir": base destination directory into which to reorganize model files
                             "modelSourceUrlMapD": dictionary mapping between source model filenames and accession URLs
            keepSource (bool): whether to copy model files to new directory (instead of moving them). Defaults to False.
        """

        try:
            self.__cachePath = cachePath if cachePath else self.__getModelCachePath()
            self.__numProc = kwargs.get("numProc", 2)
            self.__chunkSize = kwargs.get("chunkSize", 20)
            self.__modelDir = modelD.get("modelDir")
            self.__workPath = self.__modelDir
            self.__modelFileList = modelD.get("modelFileList", [])
            self.__modelSourcePrefix = modelD.get("modelSourcePrefix")
            self.__destModelDir = modelD.get("destModelDir")
            self.__modelSourceUrlMapD = modelD.get("modelSourceUrlMapD")
            self.__keepSource = keepSource

            self.__mU = MarshalUtil(workPath=self.__workPath)

        except Exception as e:
            logger.exception("Failing with %s", str(e))

    def getCachePath(self):
        return self.__cachePath

    # def __getModelCachePath(self, fmt="pickle"):
    def __getModelCachePath(self, fmt="json"):
        ext = "pic" if fmt == "pickle" else "json"
        pth = os.path.join(self.__modelDir, "species-model-files-cache." + ext)
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
            mD, failD = self.__reorganizeModels(numProc=self.__numProc, chunkSize=self.__chunkSize)
            if len(failD) > 0:
                logger.error("Failed to move %d model files.", len(failD))
            kwargs = {"indent": indent} if fmt == "json" else {"pickleProtocol": 4}
            ok = self.__mU.doExport(self.__cachePath, mD, fmt=fmt, **kwargs)
            logger.info("Wrote %r status %r", self.__cachePath, ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return ok

    def __reorganizeModels(self, numProc=2, chunkSize=20):
        """Prepare multiprocessor queue and workers for reorganizing and renaming all downloaded model files into
        structured directory tree using internal identifiers.

        Args:
            numProc (int, optional): number of processes to use. Defaults to 2.
            chunkSize (int, optional): incremental chunk size used for distribute work processes. Defaults to 20.

        Returns:
            mD (dict): dictionary of successfully processed models, in the following structure:
                       {internalModelId: {"sourceModelName": modelFileNameIn, "modelPath": modelFileOut, "sourceModelUrl": sourceModelUrl}, ...}
            failD (dict): dictionary of models for which processing failed, in the following structure:
                          {internalModelId: {"sourceModelName": modelFileNameIn}, ...}
        """
        mD = {}
        failD = {}
        #
        logger.info("Starting with %d models, numProc %d", len(self.__modelFileList), self.__numProc)
        #
        rWorker = ModelWorker(workPath=self.__workPath)
        mpu = MultiProcUtil(verbose=True)
        optD = {
            "modelSourcePrefix": self.__modelSourcePrefix,
            "destModelDir": self.__destModelDir,
            "modelSourceUrlMapD": self.__modelSourceUrlMapD,
            "keepSource": self.__keepSource
        }
        mpu.setOptions(optD)
        mpu.set(workerObj=rWorker, workerMethod="reorganize")
        mpu.setWorkingDir(workingDir=self.__workPath)
        ok, failList, resultList, _ = mpu.runMulti(dataList=self.__modelFileList, numProc=numProc, numResults=1, chunkSize=chunkSize)
        if failList:
            logger.info("model file failures (%d): %r", len(failList), failList)
        #
        for (modelFileNameIn, internalModelId, modelFileOut, sourceModelUrl) in resultList[0]:
            if modelFileOut:
                mD[internalModelId] = {"modelPath": modelFileOut, "sourceModelName": modelFileNameIn, "sourceModelUrl": sourceModelUrl}
            else:
                failD[internalModelId] = {"sourceModelName": modelFileNameIn}
        #
        logger.info("Completed with multi-proc status %r, failures %r, total models with data (%d)", ok, len(failList), len(mD))
        return mD, failD
