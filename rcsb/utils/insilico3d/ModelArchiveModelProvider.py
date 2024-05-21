##
# File:    ModelArchiveModelProvider.py
# Author:  Dennis Piehl
# Date:    18-Mar-2022
#
# Updates:
#   24-Oct-2022  dwp Add functionality to fetch the release date associated for a given ModelArchive data set
#                    (to use for data loading in case model mmCIF file doesn't contain this information already)
#   16-Nov-2022  dwp Add new default functionality to fetch model files individually instead of the full bulk download
#    9-Jan-2023  dwp Fetch data set model IDs directly (don't try to construct them here), and add more ModelArchive data sets
#    7-Dec-2023  dwp Update base URL for individual model file downloads, which are now served as .cif.gz when using aiohttp
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
import json
from pathlib import Path
import glob
import asyncio
import requests
import aiohttp
import aiofiles

from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.insilico3d.ModelReorganizer import ModelReorganizer

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger(__name__)


class ModelArchiveModelProvider:
    """Accessors for ModelArchive 3D in silico models (mmCIF)."""

    def __init__(self, cachePath, baseWorkPath=None, useCache=False, reload=True, **kwargs):
        """Initialize AlphaFoldModelProvider class.

        Args:
            cachePath (str): Path to directory where model files will be reorganized and stored permanently (also contains computed-model-cache.json file)  (i.e., 'computed-models'),
                             and will serve as the parent directory of where model files will be downloaded and where MA-specific cache file will sit (under 'work-dir/ModelArchive').
            useCache (bool, optional): Start from cache of already downloaded and/or reorganized model files. Defaults to False.
                                       When True, checks if last downloaded set of files is up-to-date and downloads any newly available models.
                                       When False (default), redownloads all model files.
            reload (bool, optional): Peform full reload (i.e., download/update) upon instantiation. Defaults to True.
        """
        # Use the same root cachePath for all types of insilico3D model sources, but with unique dirPath names (sub-directory)
        self.__cachePath = cachePath  # Cache path is where all model files will eventually be reorganized and stored in (i.e. "computed-models")
        self.__baseWorkPath = baseWorkPath if baseWorkPath else self.__cachePath
        self.__workPath = os.path.join(self.__baseWorkPath, "work-dir", "ModelArchive")  # Directory where model files will be downloaded (also contains MA-specific cache file)
        self.__dataSetCacheFile = os.path.join(self.__workPath, "model-download-cache.json")
        self.__dataSetHoldingsFileName = "modelarchive-holdings.json.gz"

        self.__modelArchiveSummaryPageBaseApiUrl = "https://www.modelarchive.org/api/projects/"
        self.__modelArchiveBaseDownloadUrl = "https://www.modelarchive.org/doi/10.5452/"
        # Use above for direct gzipped file downloads (e.g., https://www.modelarchive.org/doi/10.5452/ma-bak-cepc-0001.cif.gz)
        self.__modelArchiveBulkDownloadUrlEnd = "?type=materials_procedures__accompanying_data_file_name"  # E.g., "ma-bak-cepc?type=materials_procedures__accompanying_data_file_name"

        self.__mU = MarshalUtil(workPath=self.__workPath)
        self.__fU = FileUtil(workPath=self.__workPath)

        if reload:
            self.__oD, self.__createdDate = self.__reload(useCache=useCache, **kwargs)

    def testCache(self, minCount=0):  # Increase minCount once we are consistently downloading more than one data set
        if self.__oD and len(self.__oD) > minCount:
            return True
        else:
            return False

    def getCacheDirPath(self):
        return self.__workPath

    def reload(self, useCache, **kwargs):
        self.__oD, self.__createdDate = self.__reload(useCache=useCache, **kwargs)

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
            baseUrl = kwargs.get("baseUrl", self.__modelArchiveSummaryPageBaseApiUrl)
            modelArchiveRequestedDatasetD = kwargs.get("modelArchiveRequestedDatasetD", {})
            if not modelArchiveRequestedDatasetD:  # use default
                modelArchiveRequestedDatasetD = {
                    "ma-bak-cepc": {
                        # "bulkFileName": "ma-bak-cepc.zip",
                    },
                    "ma-ornl-sphdiv": {
                        # "bulkFileName": "ma-ornl-sphdiv.zip",
                    },
                    "ma-coffe-slac": {
                        # "bulkFileName": "",
                    },
                    "ma-asfv-asfvg": {
                        # "bulkFileName": "ma-asfv-asfvg.zip",
                    },
                    "ma-t3vr3": {
                        # "bulkFileName": "ma-t3vr3.zip",
                    },
                }
            #
            self.__fU.mkdir(self.__workPath)
            #
            logger.info("useCache %r self.__dataSetCacheFile %r", useCache, self.__dataSetCacheFile)
            if useCache and self.__mU.exists(self.__dataSetCacheFile):
                logger.info("Loading data cache, %s.", self.__dataSetCacheFile)
                cacheD = self.__mU.doImport(self.__dataSetCacheFile, fmt="json")
                createdDate = cacheD["created"]
                oD = cacheD["data"]

                logger.info("Checking consistency of cached data with data available on server")
                for dataSet, pathD in modelArchiveRequestedDatasetD.items():
                    try:
                        cacheArchiveDir = oD[dataSet]["dataDirectory"]
                        cacheArchiveFileDownloadDate = oD[dataSet]["lastDownloaded"]
                        cacheArchiveFileDownloadAge = (datetime.datetime.now() - datetime.datetime.fromisoformat(cacheArchiveFileDownloadDate)).days
                        if not os.path.exists(cacheArchiveDir):
                            logger.warning("Missing archive data for dataSet %s from server: %s", dataSet, pathD)
                        # If 120 days old, log WARNING about age of archive and possibly being out-of-date
                        if cacheArchiveFileDownloadAge > 120:
                            logger.warning(
                                "Cached archive data for dataset %s last downloaded on %s (%d days ago): %s. Recommend redownloading data.",
                                dataSet, cacheArchiveFileDownloadDate, cacheArchiveFileDownloadAge, pathD
                            )
                        else:
                            logger.info(
                                "Cached archive data for dataset %s last downloaded on %s (%d days ago): %s. Skipping redownload of data.",
                                dataSet, cacheArchiveFileDownloadDate, cacheArchiveFileDownloadAge, pathD
                            )
                    except Exception as e:
                        logger.info("Failing on checking of cache data for dataSet %s: %s", dataSet, pathD)
                        logger.exception("Failing with %s", str(e))
            else:
                logger.info("Refetching all files from server.")
                cacheD = {}
                cacheD.update({"created": startDateTime, "data": {}})
                for dataSet, pathD in modelArchiveRequestedDatasetD.items():
                    try:
                        sD = {}
                        numModelsToDownload = pathD.get("numModels", None)  # Used for testing purposes, defaults to total number of models
                        bulkFileName = pathD.get("bulkFileName", None)
                        numModelsDownloaded = 0
                        if bulkFileName:
                            # Download bulk model archive file (contains associated local pairwise quality data and a3m files)
                            sD.update({"downloadMethod": "bulk"})
                            sD.update({"bulkArchiveFileName": bulkFileName})
                            dataSetFilePath = os.path.join(baseUrl, dataSet + self.__modelArchiveBulkDownloadUrlEnd)
                            dataSetDataDumpDir = os.path.join(self.__workPath, dataSet.replace(" ", "_"))
                            self.__fU.mkdir(dataSetDataDumpDir)
                            dataSetFileDumpPath = os.path.join(dataSetDataDumpDir, bulkFileName)
                            #
                            logger.info("Fetching file %s from server to local path %s", dataSetFilePath, dataSetFileDumpPath)
                            ok = self.__fU.get(dataSetFilePath, dataSetFileDumpPath)
                            logger.info("Completed fetch (%r) at %s (%.4f seconds)", ok, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime)
                            ok = self.__fU.unbundleZipfile(dataSetFileDumpPath, dirPath=dataSetDataDumpDir)
                            logger.info("Completed unbundle (%r) at %s (%.4f seconds)", ok, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime)
                            #
                            logger.info("Clearing non-model files from extracted zip bundle...")
                            for nonModelFile in Path(dataSetDataDumpDir).glob("*.a3m"):
                                nonModelFile.unlink()
                            for nonModelFile in Path(dataSetDataDumpDir).glob("*_local_pairwise_qa.cif"):
                                nonModelFile.unlink()
                            numModelsDownloaded = len(list(Path(dataSetDataDumpDir).glob("*.cif*")))
                        else:
                            # Download model files individually
                            sD.update({"downloadMethod": "individual"})
                            dataSetDataDumpDir = os.path.join(self.__workPath, dataSet.replace(" ", "_"))
                            self.__fU.mkdir(dataSetDataDumpDir)
                            logger.info("Fetching files for %s from server to local path %s", dataSet, dataSetDataDumpDir)
                            ok, numModelsDownloaded = asyncio.run(self.downloadIndividualModelFiles(
                                modelSetName=dataSet,
                                destDir=dataSetDataDumpDir,
                                numModels=numModelsToDownload
                            ))
                            logger.info("Completed fetch (%r) at %s (%.4f seconds)", ok, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime)
                        #
                        sD.update({
                            "dataSetName": dataSet,
                            "numModels": numModelsDownloaded,
                            "lastDownloaded": startDateTime,
                            "dataDirectory": dataSetDataDumpDir,
                        })
                        if ok:
                            cacheD["data"].update({dataSet: sD})
                            if bulkFileName:
                                self.__fU.remove(dataSetFileDumpPath)
                    #
                    except Exception as e:
                        logger.info("Failing on fetching of dataSet %s: %s", dataSet, pathD)
                        logger.exception("Failing with %s", str(e))

                createdDate = cacheD["created"]
                oD = cacheD["data"]
                ok = self.__mU.doExport(self.__dataSetCacheFile, cacheD, fmt="json", indent=3)
                logger.info("Export ModelArchive dataSet model data (%d) status %r", len(oD), ok)

        except Exception as e:
            logger.exception("Failing with %s", str(e))

        return oD, createdDate

    def fetchModelIdList(self, modelSetName):
        """Fetech the list of individual models files for a ModelArchive data set.

        Args:
            modelSetName (str): model set name (e.g., "ma-ornl-sphdiv").

        Returns:
            list: list of individual model IDs
        """
        modelSetResp = requests.get(os.path.join(self.__modelArchiveSummaryPageBaseApiUrl, modelSetName), timeout=600)
        modelSetRespMaterials = modelSetResp.json()["materials_procedures"]["materials"]
        startIdx = modelSetRespMaterials.index("linkData=") + len("linkData=")
        endIdx = modelSetRespMaterials.index("];", startIdx) + 1
        modelSetString = modelSetRespMaterials[startIdx:endIdx]
        modelSetL = json.loads(modelSetString)
        modelIdList = []
        for data in modelSetL:
            modelIdList.append(data["id"])
        return modelIdList

    async def downloadIndividualModelFiles(self, modelSetName=None, destDir=None, limit=100, breakTime=3, numModels=None):
        """Download model files individually, in case bulk download not available or want to avoid downloading (currently) unnecessary associated metdata.

        Args:
            modelSetName (str): model set name (e.g., "ma-ornl-sphdiv").
            destDir (str): destination directory to which to download model files.
            numModels (int): number of models to download from model set.
            limit (int, optional): max number of models to download asynchronously at once. Splits the total set into batches of size(limit),
                                   and forces sleep(breakTime) between batches to spare traffic load on ModelArchive server). Defaults to 100.
            breakTime (int, optional): seconds to wait between subsequent batches of asynchronous requests. Defaults to 3.

        Returns:
            (bool): True if successful; False otherwise.
        """

        if not (modelSetName and destDir):
            return False

        # First, fetch list of model set IDs
        modelSetIdFullList = self.fetchModelIdList(modelSetName)
        numModels = numModels if numModels else len(modelSetIdFullList)
        modelSetIdL = modelSetIdFullList[0:numModels]

        modelUrlList = []
        for mId in modelSetIdL:
            mUrl = os.path.join(self.__modelArchiveBaseDownloadUrl, f"{mId}.cif.gz")
            modelUrlList.append(mUrl)
        logger.info("First few items in modelUrlList %r", modelUrlList[0:5])

        sema = asyncio.BoundedSemaphore(20)

        async def fetchFile(url, session):
            try:
                fname = url.split("/")[-1]
                # fname = url.split("/")[-1].split("?")[0] + ".cif.gz"
                async with sema, session.get(url, timeout=10) as resp:
                    assert resp.status == 200
                    data = await resp.read()
                async with aiofiles.open(os.path.join(destDir, fname), "wb") as outfile:
                    await outfile.write(data)
                return True
            #
            except Exception as e:
                logger.exception("Failing to fetch url %s with %s", url, str(e))
                return url

        def modelUrlBatches(fullUrlList, batchSize):
            for i in range(0, len(fullUrlList), batchSize):
                yield fullUrlList[i:i + batchSize]
        #
        resultList, failList = [], []
        maxRetries = 10
        for batchNum, batchUrls in enumerate(modelUrlBatches(modelUrlList, limit)):
            logger.info("Downloading batch %d", batchNum + 1)
            async with aiohttp.ClientSession() as session:
                tasks = []
                for modelUrl in batchUrls:
                    tasks.append(fetchFile(modelUrl, session))
                resL = await asyncio.gather(*tasks)
                failL = [i for i in resL if i is not True]
                resultList += resL
                time.sleep(breakTime)
            # Re-run any failed model file downloads
            if len(failL) > 0:
                retries = 0
                while len(failL) > 0 and retries < maxRetries:
                    retries += 1
                    logger.info("Re-attempting fetch (retry %d) for %d model files: %r", retries, len(failL), failL)
                    time.sleep(60)  # Give server a minute before refeteching
                    async with aiohttp.ClientSession() as session:
                        tasks = []
                        for modelUrl in failL:
                            tasks.append(fetchFile(modelUrl, session))
                        resL = await asyncio.gather(*tasks)
                        failL = [i for i in resL if i is not True]
                    if len(failL) > 0:
                        logger.info("Re-fetch attempt %d failed for %d model files: %r", retries, len(failL), failL)
                    else:
                        logger.info("Re-fetch succeeded for all model files")
                failList += [i for i in failL]
        #
        ok = len(failList) == 0 and len(resultList) > 0
        numModelsDownloaded = len([i for i in resultList if i is True])
        return ok, numModelsDownloaded

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
                modelFiles = glob.glob(os.path.join(modelDir, "*.cif.gz"))  # may need to be ".cif" for bulk downloads, but need to check
                modelFileList = [os.path.abspath(f) for f in modelFiles]
            except Exception as e:
                logger.exception("Failing with %s", str(e))

        return modelFileList

    def getArchiveDataDownloadDate(self):
        return self.__createdDate

    def getBaseDataPath(self):
        return self.__workPath

    def getArchiveDataCacheFilePath(self):
        return self.__dataSetCacheFile

    def getModelReorganizer(self, cachePath=None, useCache=True, workPath=None, **kwargs):
        cachePath = cachePath if cachePath else self.__cachePath
        workPath = workPath if workPath else self.__workPath
        cacheFile = kwargs.get("cacheFile", self.__dataSetHoldingsFileName)
        cacheFormat = kwargs.get("cacheFormat", "json")
        return ModelReorganizer(cachePath=cachePath, useCache=useCache, workPath=workPath, cacheFile=cacheFile, cacheFormat=cacheFormat, **kwargs)

    def getComputedModelsDataPath(self):
        return self.__cachePath

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
            sourceArchiveReleaseDate = kwargs.get("sourceArchiveReleaseDate", None)  # Use for ModelArchive files until revision date information is available in mmCIF files
            #
            mR = self.getModelReorganizer(cachePath=cachePath, useCache=useCache, **kwargs)
            #
            if inputModelList:  # Only reorganize given list of model files
                _, ok = mR.reorganize(
                    inputModelList=inputModelList,
                    modelSource="ModelArchive",
                    destBaseDir=self.__cachePath,
                    useCache=useCache,
                    sourceArchiveReleaseDate=sourceArchiveReleaseDate,
                )
                if not ok:
                    logger.error("Reorganization of model files failed for inputModelList starting with item, %s", inputModelList[0])
            #
            else:  # Reorganize ALL model files for ALL available model datasets
                archiveDirList = self.getArchiveDirList()
                for archiveDir in archiveDirList:
                    #
                    # Get release date of the dataset archive (TEMPORARY workaround until revision history is included in ModelArchive mmCIF files)
                    archiveName = os.path.basename(os.path.abspath(archiveDir))
                    archiveSummaryPageApiUrl = os.path.join(self.__modelArchiveSummaryPageBaseApiUrl, archiveName)
                    response = requests.get(archiveSummaryPageApiUrl, timeout=600)
                    try:
                        sourceArchiveReleaseDate = response.json()["release_date"]  # e.g., '2022-09-28'
                        logger.info("Dataset archive %s: release date %s", archiveName, sourceArchiveReleaseDate)
                    except Exception as e:
                        logger.exception(
                            "Failing to get release date for archive dataset %s from ModelArchive site (returned status code %r), with exception %s",
                            archiveName,
                            response.status_code,
                            str(e)
                        )
                        raise ValueError("Failed to get release date for archive dataset.")
                    #
                    inputModelList = self.getModelFileList(inputPathList=[archiveDir])
                    _, ok = mR.reorganize(
                        inputModelList=inputModelList,
                        modelSource="ModelArchive",
                        destBaseDir=self.__cachePath,
                        useCache=useCache,
                        sourceArchiveReleaseDate=sourceArchiveReleaseDate,
                    )
                    if not ok:
                        logger.error("Reorganization of model files failed for dataset archive %s", archiveDir)
                        break
            return ok
        #
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            return False
