##
# File:    ModelArchiveModelProvider.py
# Author:  Dennis Piehl
# Date:    18-Mar-2022
#
# Updates:
#
##
"""
Accessors for ModelArchive 3D In Silico Models (mmCIF).

"""

__docformat__ = "google en"
__author__ = "Dennis Piehl"
__email__ = "dennis.piehl@rcsb.org"
__license__ = "Apache 2.0"

import datetime
import logging
import os.path
import time
from pathlib import Path
import glob

from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.insilico3d.ModelReorganizer import ModelReorganizer

logger = logging.getLogger(__name__)


class ModelArchiveModelProvider:
    """Accessors for ModelArchive 3D in silico models (mmCIF)."""

    def __init__(self, cachePath, useCache=False, cfgOb=None, configName=None, **kwargs):
        # Use the same root cachePath for all types of insilico3D model sources, but with unique dirPath names (sub-directory)
        self.__cachePath = cachePath
        self.__cfgOb = cfgOb
        self.__configName = configName
        self.__dirPath = os.path.join(self.__cachePath, "ModelArchive")
        self.__dataSetCacheFile = os.path.join(self.__dirPath, "ma-model-data.json")
        self.__computedModelsDataPath = self.__cfgOb.getPath("PDBX_COMP_MODEL_SANDBOX_PATH", sectionName=self.__configName, default=os.path.join(self.__cachePath, "computed-models"))

        self.__mU = MarshalUtil(workPath=self.__dirPath)
        self.__fU = FileUtil(workPath=self.__dirPath)

        self.__oD, self.__createdDate = self.__reload(useCache=useCache, **kwargs)

    def testCache(self, minCount=0):  # Increase minCount once we are consistently downloading more than one data set
        if self.__oD and len(self.__oD) > minCount:
            return True
        else:
            return False

    def getCacheDirPath(self):
        return self.__dirPath

    def __reload(self, **kwargs):
        """Reload cached list of ModelArchive model dataset files and check server for updated data sets,
        or re-download latest versions of all model data sets from server.

        Returns:
            oD (dict): dictionary of cached/downloaded model data sets with general metadata about each data set (local data path, size, # files, ...)
            createdDate (str): timestamp in isoformat of when the data cache was last created
        """
        try:
            oD = None
            createdDate = None
            ok = False
            startTime = time.time()
            startDateTime = datetime.datetime.now().isoformat()
            useCache = kwargs.get("useCache", True)
            baseUrl = kwargs.get("baseUrl", "https://www.modelarchive.org/api/projects")
            serverDataSetPathD = kwargs.get("serverDataSetPathD", {
                "ma-bak-cepc": {
                    "urlEnd": "ma-bak-cepc?type=materials_procedures__accompanying_data_file_name",
                    "fileName": "ma-bak-cepc.zip"
                },
            })

            self.__fU.mkdir(self.__dirPath)

            logger.info("useCache %r self.__dataSetCacheFile %r", useCache, self.__dataSetCacheFile)
            if useCache and self.__mU.exists(self.__dataSetCacheFile):
                logger.info("Loading data cache, %s.", self.__dataSetCacheFile)
                cacheD = self.__mU.doImport(self.__dataSetCacheFile, fmt="json")
                createdDate = cacheD["created"]
                oD = cacheD["data"]

                logger.info("Checking consistency of cached data with data available on server")
                for dataSet, pathD in serverDataSetPathD.items():
                    try:
                        cacheArchiveDir = oD[dataSet]["dataDirectory"]
                        cacheArchiveFileSize = oD[dataSet]["archiveFileSizeBytes"]
                        cacheArchiveFileDownloadDate = oD[dataSet]["lastDownloaded"]
                        cacheArchiveFileDownloadAge = (datetime.datetime.now() - datetime.datetime.fromisoformat(cacheArchiveFileDownloadDate)).days
                        if not os.path.exists(cacheArchiveDir):
                            logger.warning("Missing archive data for dataSet %s from server: %s", dataSet, pathD)
                        # If 120 days old, check consistency of local cache with data on ModelArchive
                        # (requires redownloading the entire archive file, since can't get requests.header info from the ModelArchive download URL)
                        if cacheArchiveFileDownloadAge > 120:
                            logger.info(
                                "Cached archive data for dataset %s last downloaded on %s (%d days ago): %s. Redownloading to check consistency with cached data.",
                                dataSet, cacheArchiveFileDownloadDate, cacheArchiveFileDownloadAge, pathD
                            )
                            dataSetFileName = pathD["fileName"]
                            dataSetFilePath = os.path.join(baseUrl, pathD["urlEnd"])
                            dataSetDataDumpDir = os.path.join(self.__dirPath, dataSet.replace(" ", "_"))
                            self.__fU.mkdir(dataSetDataDumpDir)
                            dataSetFileDumpPath = os.path.join(dataSetDataDumpDir, dataSetFileName)
                            logger.info("Fetching file %s from server to local path %s", dataSetFilePath, dataSetFileDumpPath)
                            ok = self.__fU.get(dataSetFilePath, dataSetFileDumpPath)
                            logger.info("Completed fetch (%r) at %s (%.4f seconds)", ok, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime)
                            dataSetFileSize = int(self.__fU.size(dataSetFileDumpPath))
                            if cacheArchiveFileSize != dataSetFileSize:
                                logger.warning("Cached archive data for dataset %s not up-to-date with data archive on server: %s. You should Rebuild the cache!", dataSet, pathD)
                            logger.info("Deleting compressed dataset download, %s", dataSetFileDumpPath)
                            self.__fU.remove(dataSetFileDumpPath)
                    except Exception as e:
                        logger.info("Failing on checking of cache data for dataSet %s: %s", dataSet, pathD)
                        logger.exception("Failing with %s", str(e))
            else:
                logger.info("Refetching all files from server.")
                cacheD = {}
                cacheD.update({"created": startDateTime, "data": {}})
                for dataSet, pathD in serverDataSetPathD.items():
                    try:
                        dataSetFileName = pathD["fileName"]
                        dataSetFilePath = os.path.join(baseUrl, pathD["urlEnd"])
                        dataSetDataDumpDir = os.path.join(self.__dirPath, dataSet.replace(" ", "_"))
                        self.__fU.mkdir(dataSetDataDumpDir)
                        dataSetFileDumpPath = os.path.join(dataSetDataDumpDir, dataSetFileName)
                        #
                        logger.info("Fetching file %s from server to local path %s", dataSetFilePath, dataSetFileDumpPath)
                        ok = self.__fU.get(dataSetFilePath, dataSetFileDumpPath)
                        logger.info("Completed fetch (%r) at %s (%.4f seconds)", ok, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime)
                        dataSetFileSize = int(self.__fU.size(dataSetFileDumpPath))
                        ok = self.__fU.unbundleZipfile(dataSetFileDumpPath, dirPath=dataSetDataDumpDir)
                        logger.info("Completed unbundle (%r) at %s (%.4f seconds)", ok, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime)
                        #
                        sD = {
                            "dataSetName": dataSet,
                            "archiveFileName": dataSetFileName,
                            "archiveFileSizeBytes": dataSetFileSize,
                            "lastDownloaded": startDateTime,
                            "dataDirectory": dataSetDataDumpDir,
                        }
                        #
                        logger.info("Clearing non-model files from extracted zip bundle...")
                        for nonModelFile in Path(dataSetDataDumpDir).glob("*.a3m"):
                            nonModelFile.unlink()
                        for nonModelFile in Path(dataSetDataDumpDir).glob("*_local_pairwise_qa.cif"):
                            nonModelFile.unlink()
                        #
                        if ok:
                            cacheD["data"].update({dataSet: sD})
                            self.__fU.remove(dataSetFileDumpPath)
                    #
                    except Exception as e:
                        logger.info("Failing on fetching and expansion of file for dataSet %s: %s", dataSet, pathD)
                        logger.exception("Failing with %s", str(e))

                createdDate = cacheD["created"]
                oD = cacheD["data"]
                ok = self.__mU.doExport(self.__dataSetCacheFile, cacheD, fmt="json", indent=3)
                logger.info("Export ModelArchive dataSet model data (%d) status %r", len(oD), ok)

        except Exception as e:
            logger.exception("Failing with %s", str(e))

        return oD, createdDate

    def getArchiveDirList(self):
        archiveDirList = [self.__oD[k]["dataDirectory"] for k in self.__oD]

        return archiveDirList

    def getModelFileList(self, inputPathList=None):
        """Return a list of filepaths for all mmCIF models under the provided set of directories.

        Args:
            inputPathList (list): Optional (but recommended) list of paths to directories containing model mmCIF files to return.
                                  If not provided, method will retrieve all cached dataset paths, so all model files from all datasets will be returned in a single list;
                                  thus, it's recommended to provide a list of specific dataset directory to break the returned model list down into more manageable parts.

        Returns:
            (list): list of model mmCIF file paths (only matches compressed ".cif.gz" files, to ensure only one file per model)
        """

        if not inputPathList:
            inputPathList = self.getArchiveDirList()

        modelFileList = []

        for modelDir in inputPathList:
            try:
                modelFiles = glob.glob(os.path.join(modelDir, "*.cif"))
                modelFileList = [os.path.abspath(f) for f in modelFiles]
            except Exception as e:
                logger.exception("Failing with %s", str(e))

        return modelFileList

    def getArchiveDataDownloadDate(self):
        return self.__createdDate

    def getBaseDataPath(self):
        return self.__dirPath

    def getArchiveDataCacheFilePath(self):
        return self.__dataSetCacheFile

    def getModelReorganizer(self, cachePath=None, useCache=True, **kwargs):
        cachePath = cachePath if cachePath else self.__cachePath
        return ModelReorganizer(cachePath=cachePath, useCache=useCache, **kwargs)

    def getComputedModelsDataPath(self):
        return self.__computedModelsDataPath

    def reorganizeModelFiles(self, cachePath=None, useCache=True, inputModelList=None, **kwargs):
        """Reorganize model files from organism-wide model listing to hashed directory structure and rename files
        to follow internal identifier naming convention.

        Args:
            cachePath (str): Path to cache directory.
            inputModelList (list, optional): List of input model filepaths to reorganize; defaults to all models for all model datasets.
            **kwargs (optional):
                numProc (int): number of processes to use; default 2.
                chunkSize (int): incremental chunk size used for distribute work processes; default 20.
                keepSource (bool): whether to copy files to new directory (instead of moving them); default False.
                cacheFilePath (str): full filepath and name for cache file containing a dictionary of all reorganized models.

        Returns:
            (bool): True if successful; False otherwise.
        """

        try:
            ok = False
            #
            mR = self.getModelReorganizer(cachePath=cachePath, useCache=useCache, **kwargs)
            #
            if inputModelList:  # Only reorganize given list of model files
                ok = mR.reorganize(inputModelList=inputModelList, modelSource="ModelArchive", destBaseDir=self.__computedModelsDataPath, useCache=useCache)
                if not ok:
                    logger.error("Reorganization of model files failed for inputModelList starting with item, %s", inputModelList[0])
            #
            else:  # Reorganize ALL model files for ALL available model datasets
                archiveDirList = self.getArchiveDirList()
                for archiveDir in archiveDirList:
                    inputModelList = self.getModelFileList(inputPathList=[archiveDir])
                    ok = mR.reorganize(inputModelList=inputModelList, modelSource="ModelArchive", destBaseDir=self.__computedModelsDataPath, useCache=useCache)
                    if not ok:
                        logger.error("Reorganization of model files failed for dataset archive %s", archiveDir)
                        break
            return ok
        #
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            return False
