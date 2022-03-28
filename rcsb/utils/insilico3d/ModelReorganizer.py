##
# File:    ModelReorganizer.py
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
import copy

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
        """Enumerate and reorganize and rename model files.

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
            destBaseDir = optionsD.get("destBaseDir")  # base path for all computed models (i.e., "computed-models"); Or will be root path at HTTP endpoint
            keepSource = optionsD.get("keepSource", False)  # whether to copy files over (instead of moving them)
            #
            for modelFileIn in dataList:
                modelD = {}
                success = False
                modelFileOut = None
                modelFileNameIn = self.__fU.getFileName(modelFileIn)
                #
                containerList = self.__mU.doImport(modelFileIn, fmt="mmcif")
                if len(containerList) > 1:
                    # Expecting all computed models to have one container per file. When this becomes no longer the case, update this to handle it accordingly.
                    logger.error("Skipping - model file %s has more than one container (%d)", modelFileNameIn, len(containerList))
                    continue
                modelEntryId = containerList[0].getObj("entry").getValue("id", 0)
                # Remove all punctuation and make ALL CAPS
                modelEntryId = "".join(char for char in modelEntryId if char.isalnum()).upper()
                internalModelId = modelSourcePrefix + "_" + modelEntryId
                internalModelName = internalModelId + ".cif.gz"
                # Use last six to last two characters for second-level hashed directory
                firstDir, secondDir = modelEntryId[-6:-4], modelEntryId[-4:-2]
                destModelDir = os.path.join(destBaseDir, modelSourcePrefix, firstDir, secondDir)
                if not self.__fU.exists(destModelDir):
                    self.__fU.mkdir(destModelDir)
                modelFileOut = os.path.join(destModelDir, internalModelName)
                #
                sourceModelUrl = self.__getSourceUrl(modelSourcePrefix, modelFileNameIn, modelEntryId)
                #
                modelD["modelId"] = internalModelId
                modelD["modelPath"] = modelFileOut
                modelD["sourceModelFileName"] = modelFileNameIn
                modelD["sourceModelUrl"] = sourceModelUrl
                #
                try:
                    if keepSource:
                        self.__fU.put(modelFileIn, modelFileOut)  # Copies files (only use for testing)
                    else:
                        self.__fU.replace(modelFileIn, modelFileOut)  # Moves files (use for production)
                    success = True
                    successList.append(modelFileIn)
                except Exception as e:
                    logger.debug("Failing to reorganize %s --> %s, with %s", modelFileIn, modelFileOut, str(e))
                #
                retList.append((modelFileIn, modelD, success))
                #
            failList = sorted(set(dataList) - set(successList))
            if failList:
                logger.info("%s returns %d definitions with failures: %r", procName, len(failList), failList)
            #
            logger.debug("%s processed %d/%d models, failures %d", procName, len(retList), len(dataList), len(failList))
        except Exception as e:
            logger.exception("Failing %s for %d data items %s", procName, len(dataList), str(e))

        return successList, retList, diagList

    def __getSourceUrl(self, modelSourcePrefix, sourceModelFileName, modelEntryId):
        """Construct model accession URL for each model source.

        Args:
            modelFileList (list): List of model file paths.

        Returns:
            (list): list of dictionary mapping betwen source model filenames and accession URLs.
                    Note that file-specific downloads aren't gzipped, unlike model files in species tarball.
                    E.g., "https://alphafold.ebi.ac.uk/files/AF-Q9RQP8-F1-model_v2.cif"
        """
        sourceModelUrl = None
        try:
            if modelSourcePrefix == "AF":
                modelFileNameInUrl = sourceModelFileName.split(".gz")[0]
                sourceModelUrl = os.path.join("https://alphafold.ebi.ac.uk/files/", modelFileNameInUrl)
            elif modelSourcePrefix == "MB":
                modbaseInternalId = modelEntryId.split("model_")[-1]
                # sourceModelUrl = "https://salilab.org/modbase/searchbyid?modelID=" + modbaseInternalId + "&displaymode=moddetail"  # Directs to entry page
                sourceModelUrl = "https://salilab.org/modbase/retrieve/modbase/?modelID=" + modbaseInternalId + "&format=mmcif"  # Directs to mmCIF file displayed in web browser
                # E.g.: https://salilab.org/modbase/retrieve/modbase/?modelID=ecac68b60ee6877ccde36af05cdeac58&format=mmcif
            elif modelSourcePrefix == "MA":
                sourceModelUrl = "https://www.modelarchive.org/api/projects/" + modelEntryId + "?type=basic__model_file_name"
            else:
                logger.error("URL mapping process not ready yet for %s", modelSourcePrefix)
        except Exception as e:
            logger.exception("Failing with %s", str(e))

        return sourceModelUrl


class ModelReorganizer(object):
    """Generators and accessors for model files."""

    def __init__(self, cachePath=None, useCache=True, **kwargs):
        """Initialize ModelReorganizer object.

        Args:
            cachePath (str): directory path for storing cache file containing a dictionary of all reorganized models.
            useCache (bool): whether to use the existing data cache or re-run entire model reorganization process.
            cacheFile (str, optional): filename for cache file; default "computed-models-cache.json".
            cacheFilePath (str, optional): full filepath and name for cache file (will override "cachePath" and "cacheFile" if provided).
            numProc (int, optional): number of processes to use; default 2.
            chunkSize (int, optional): incremental chunk size used for distribute work processes; default 20.
            workPath (str, optional): directory path for workers to operate in; default is cachePath.
            keepSource (bool, optional): whether to copy model files to new directory instead of moving them; default False.
        """

        try:
            self.__cachePath = cachePath if cachePath else "."
            self.__cacheFormat = kwargs.get("cacheFormat", "json")
            # self.__cacheFormat = kwargs.get("cacheFormat", "pickle")
            self.__cacheFile = kwargs.get("cacheFile", "computed-models-cache." + self.__cacheFormat)
            self.__cacheFilePath = kwargs.get("cacheFilePath", os.path.join(self.__cachePath, self.__cacheFile))
            self.__numProc = kwargs.get("numProc", 2)
            self.__chunkSize = kwargs.get("chunkSize", 20)
            self.__workPath = kwargs.get("workPath", self.__cachePath)
            self.__keepSource = kwargs.get("keepSource", False)

            self.__mU = MarshalUtil(workPath=self.__workPath)
            self.__fU = FileUtil(workPath=self.__workPath)

            self.__mD = self.__reload(cacheFilePath=self.__cacheFilePath, useCache=useCache)

        except Exception as e:
            logger.exception("Failing with %s", str(e))

    def testCache(self, minCount=20):
        try:
            if minCount == 0:
                return True
            if self.__mD and len(self.__mD) >= minCount:
                logger.info("Reorganized models in cache (%d)", len(self.__mD))
                return True
        except Exception:
            pass
        return False

    def reload(self, cacheFilePath, useCache=True):
        self.__mD = self.__reload(cacheFilePath, useCache=useCache)
        return len(self.__mD) > 0  # Returns True or False, if reload was successful or not

    def __reload(self, cacheFilePath, useCache=True):
        """Reload from the current cache file."""
        try:
            mD = {}
            logger.info("useCache %r cacheFilePath %r", useCache, cacheFilePath)
            if useCache and self.__mU.exists(cacheFilePath):
                if cacheFilePath != self.__cacheFilePath:
                    self.__cacheFilePath = cacheFilePath
                mD = self.__mU.doImport(cacheFilePath, fmt=self.__cacheFormat)
                logger.info("Reorganized models (%d) in cacheFilePath %r", len(mD), cacheFilePath)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        #
        return mD

    def getCachePath(self):
        return self.__cachePath

    def getCacheFilePath(self):
        return self.__cacheFilePath

    def reorganize(self, inputModelList, modelSource, destBaseDir, useCache=True):
        """Move model files from organism-wide model listing to hashed directory structure and rename files
        to follow internal identifier naming convention.

        Args:
            inputModelList (list): List of input model filepaths to reorganize.
            modelSource (str): Source of model files ("AlphaFold", "ModBase", "ModelArchive", or "SwissModel")
            destBaseDir (str): Base destination directory into which to reorganize model files (e.g., "computed-models")

        Returns:
            bool: True for success or False otherwise
        """
        ok = False
        try:
            mD, failD = self.__reorganizeModels(
                inputModelList=inputModelList,
                modelSource=modelSource,
                destBaseDir=destBaseDir,
                numProc=self.__numProc,
                chunkSize=self.__chunkSize
            )
            if len(failD) > 0:
                logger.error("Failed to process %d model files.", len(failD))
            kwargs = {"indent": 4} if self.__cacheFormat == "json" else {"pickleProtocol": 4}
            if useCache:
                self.__mD = self.__reload(cacheFilePath=self.__cacheFilePath, useCache=useCache)
                for modelId, modelD in mD.items():
                    # if modelId not in self.__mD:
                    #     self.__mD[modelId] = modelD
                    self.__mD.update({modelId: modelD})
            else:
                self.__mD = copy.deepcopy(mD)
            ok = self.__mU.doExport(self.__cacheFilePath, self.__mD, fmt=self.__cacheFormat, **kwargs)
            logger.info("Wrote %r status %r", self.__cacheFilePath, ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return ok

    def __reorganizeModels(self, inputModelList, modelSource, destBaseDir, numProc=2, chunkSize=20):
        """Prepare multiprocessor queue and workers for generating input:output map for model files, to use in reorganizing
        them into a structured directory tree and naming them with internal identifiers.

        Args:
            inputModelList (list): List of input model filepaths to reorganize.
            modelSource (str): Source of model files ("AlphaFold", "ModBase", "ModelArchive", or "SwissModel")
            destBaseDir (str): Base destination directory into which to reorganize model files (e.g., "computed-models")
            numProc (int, optional): number of processes to use; default 2.
            chunkSize (int, optional): incremental chunk size used for distribute work processes; default 20.

        Returns:
            mD (dict): dictionary of successfully processed models, in the following structure:
                        { internalModelId_1 : {"sourceModelFileName": ..., "modelPath": ..., "sourceModelUrl": ...},
                          internalModelId_2 : {...}, ...}
                E.g:    { "AF_AFO25670F1" : {
                            "sourceModelFileName": "AF-O25670-F1-model_v2.cif.gz",
                            "modelPath": "./AF/56/70/AF_AFO25670F1.cif.gz",
                            "sourceModelUrl": "https://alphafold.ebi.ac.uk/files/AF-O25670-F1-model_v2.cif"},
                          "AF_AFO24939F1" : {...},
                        ...}
            failD (dict): dictionary of models for which processing failed and the associated output file path and source details
                          that were attempted, in the following structure:
                          {inputModelFilePath:  {"sourceModelFileName": ..., "modelPath": ..., "sourceModelUrl": ...}, ...}
        """
        mD = {}
        failD = {}
        #
        logger.info("Starting with %d models, numProc %d", len(inputModelList), numProc)
        #
        # Create the base destination directory if it doesn't exist
        if not self.__fU.exists(destBaseDir):
            logger.info("Creating base destination directory for model file reorganization, %s", destBaseDir)
            self.__fU.mkdir(destBaseDir)
        #
        modelSourcePrefixD = {"AlphaFold": "AF", "ModBase": "MB", "ModelArchive": "MA", "SwissModel": "SM"}
        modelSourcePrefix = modelSourcePrefixD[modelSource]
        #
        rWorker = ModelWorker(workPath=self.__workPath)
        mpu = MultiProcUtil(verbose=True)
        optD = {"modelSourcePrefix": modelSourcePrefix, "destBaseDir": destBaseDir, "keepSource": self.__keepSource}
        mpu.setOptions(optD)
        mpu.set(workerObj=rWorker, workerMethod="reorganize")
        mpu.setWorkingDir(workingDir=self.__workPath)
        ok, failList, resultList, _ = mpu.runMulti(dataList=inputModelList, numProc=numProc, numResults=1, chunkSize=chunkSize)
        if failList:
            logger.info("model file failures (%d): %r", len(failList), failList)
        #
        for (modelFileIn, modelD, success) in resultList[0]:
            modelIdD = copy.deepcopy(modelD)
            if success:
                modelId = modelIdD.pop("modelId")
                mD[modelId] = modelIdD
            else:
                failD[modelFileIn] = modelIdD
        #
        logger.info("Completed with multi-proc status %r, failures %r, total models with data (%d)", ok, len(failList), len(mD))
        return mD, failD
