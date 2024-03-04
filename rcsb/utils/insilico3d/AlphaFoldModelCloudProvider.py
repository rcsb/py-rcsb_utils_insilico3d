##
# File:    AlphaFoldModelCloudProvider.py
# Author:  Dennis Piehl
# Date:    7-Dec-2023
#
# Updates:
##
"""
Accessors for AlphaFold 3D Models (mmCIF) from public Google Cloud datasets.

Details about downloading these data are described here: https://github.com/google-deepmind/alphafold/blob/main/afdb/README.md

View metadata here: https://console.cloud.google.com/marketplace/product/bigquery-public-data/deepmind-alphafold?pli=1&organizationId=835909686548

View datasets here: https://console.cloud.google.com/storage/browser/public-datasets-deepmind-alphafold-v4
Example by Taxonomy ID: https://console.cloud.google.com/storage/browser/public-datasets-deepmind-alphafold-v4/proteomes;tab=objects?prefix=proteome-tax_id-38727
"""

__docformat__ = "google en"
__author__ = "Dennis Piehl"
__email__ = "dennis.piehl@rcsb.org"
__license__ = "Apache 2.0"

import datetime
import logging
import os
import time
import copy
import glob
import tarfile
from google.cloud import storage

from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.insilico3d.ModelReorganizer import ModelReorganizer

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger(__name__)


class AlphaFoldModelCloudProvider:
    """Accessors for AlphaFold models (mmCIF)."""

    def __init__(self, cachePath, baseWorkPath=None, useCache=True, redownloadBulkData=False, **kwargs):
        """Initialize AlphaFoldModelCloudProvider class.

        Args:
            cachePath (str): Path to directory where model files will be reorganized and stored permanently (also contains computed-model-cache.json file)  (i.e., 'computed-models'),
                             and will serve as the parent directory of where model files will be downloaded and where AF-specific cache file will sit (under 'work-dir/AlphaFold').
            useCache (bool, optional): Start from cache of already downloaded and/or reorganized model files. Defaults to True.
                                       When True, checks if last downloaded set of files is up-to-date and downloads any newly available models.
                                       When False (default), redownloads all model files.
            redownloadBulkData (bool, optional): Peform full data reload (i.e., download/update) upon instantiation. Defaults to False.
        """
        # Use the same root cachePath for all types of insilico3D model sources, but with unique workPath names (sub-directory)
        self.__cachePath = cachePath  # Directory where model files will be reorganized and stored permanently (also contains computed-model-cache.json file) (i.e., 'computed-models')
        self.__baseWorkPath = baseWorkPath if baseWorkPath else self.__cachePath
        self.__workPath = os.path.join(self.__baseWorkPath, "work-dir", "AlphaFoldCloud")  # Directory where model files will be downloaded (also contains AF-specific cache file)
        self.__aFCTaxIdDataCacheFile = os.path.join(self.__workPath, "model-download-cache.json")

        self.__bucketName = "public-datasets-deepmind-alphafold-v4"

        self.__mU = MarshalUtil(workPath=self.__workPath)
        self.__fU = FileUtil(workPath=self.__workPath)

        self.__oD, self.__createdDate = self.__reload(useCache=useCache, redownloadBulkData=redownloadBulkData, **kwargs)

    def testCache(self, minCount=0):  # Increase minCount once we are consistently downloading more than one species data set
        if self.__oD and len(self.__oD) > minCount:
            return True
        else:
            return False

    def reload(self, useCache, redownloadBulkData, **kwargs):
        self.__oD, self.__createdDate = self.__reload(useCache=useCache, redownloadBulkData=redownloadBulkData, **kwargs)

    def __reload(self, useCache=True, redownloadBulkData=False, **kwargs):
        """Reload cached list of species-specific AlphaFold model data files and check FTP server for updated data sets,
        or re-download latest versions of all species-specific model data sets from FTP server.

        Returns:
            oD (dict): dictionary of cached/downloaded species model data sets with general metadata about each data set (local data path, size, # files, ...)
            createdDate (str): timestamp in isoformat of when the data cache was last created
        """
        try:
            oD = None
            createdDate = None
            ok = False
            startDateTime = datetime.datetime.now().isoformat()

            alphaFoldRequestedTaxIdPrefixList = kwargs.get("alphaFoldRequestedTaxIdPrefixList", [])

            self.__fU.mkdir(self.__workPath)

            logger.info("redownloadBulkData %r useCache %r self.__aFCTaxIdDataCacheFile %r", redownloadBulkData, useCache, self.__aFCTaxIdDataCacheFile)
            if useCache and self.__mU.exists(self.__aFCTaxIdDataCacheFile):
                logger.info("Loading data cache, %s.", self.__aFCTaxIdDataCacheFile)
                cacheD = self.__mU.doImport(self.__aFCTaxIdDataCacheFile, fmt="json")
                createdDate = cacheD["created"]
                oD = cacheD["data"]
            else:
                if redownloadBulkData:
                    logger.info("Refetching all files from server.")
                    cacheD = {}
                    cacheD.update({"created": startDateTime, "data": {}})
                    for taxIdPrefix in alphaFoldRequestedTaxIdPrefixList:
                        cacheD = self.fetchTaxIdArchive(taxIdPrefix, cacheD)
                    createdDate = cacheD["created"]
                    oD = cacheD["data"]
                    ok = self.__mU.doExport(self.__aFCTaxIdDataCacheFile, cacheD, fmt="json", indent=4)
                    logger.info("Redownloaded and exported AlphaFold taxonomy archive data (%d) status %r", len(oD), ok)
                else:  # Just rebuild the cache file
                    if self.__mU.exists(self.__aFCTaxIdDataCacheFile):
                        logger.error("Cache file already exists! Program will not overwrite it. Re-name to a new file before re-running script: %r", self.__aFCTaxIdDataCacheFile)
                        return None, None
                    oD = self.rebuildArchiveDirCacheFile()
                    createdDate = startDateTime
                    cacheD = {"data": oD, "created": createdDate}
                    ok = self.__mU.doExport(self.__aFCTaxIdDataCacheFile, cacheD, fmt="json", indent=4)
                    logger.info("Rebuilt AlphaFold taxonomy archive data cache file (%d) status %r: %r", len(oD), ok, self.__aFCTaxIdDataCacheFile)

        except Exception as e:
            logger.exception("Failing with %s", str(e))

        return oD, createdDate

    def rebuildArchiveDirCacheFile(self):
        # Get a dictionary of all subdirectories and files in the given directory
        absBaseDir = os.path.abspath(self.__workPath)
        archiveDirFileD = {}
        sizeLimit = 68719476736  # 64 GB - max size of any single archive prefix subset
        with os.scandir(absBaseDir) as fObjs:
            for fObj in fObjs:
                if fObj.is_dir() and str(fObj.name).isdigit():
                    archiveDir = os.path.join(absBaseDir, fObj.name)
                    archiveDirFileD[archiveDir] = {}
                    archiveFileD = {}
                    tmpSize = 0
                    archiveDirSubsetIdx = 0
                    with os.scandir(archiveDir) as archiveObjs:
                        for archiveObj in archiveObjs:
                            if archiveObj.is_file() and archiveObj.name.endswith(".tar"):
                                fSize = archiveObj.stat().st_size
                                archiveFileD[archiveObj.name] = fSize
                                tmpSize += fSize
                                if tmpSize > sizeLimit:
                                    archiveDirFileD[archiveDir][archiveDirSubsetIdx] = {"archive_files": archiveFileD}
                                    archiveDirSubsetIdx += 1
                                    archiveFileD = {}
                                    tmpSize = 0
                    archiveDirFileD[archiveDir][archiveDirSubsetIdx] = {"archive_files": archiveFileD}
        return archiveDirFileD

    def fetchTaxIdArchive(self, taxIdPrefix, cacheD):
        try:
            startTime = time.time()
            client = storage.Client()
            bucket = client.bucket(self.__bucketName)
            #
            taxIdPrefixDataDumpDir = os.path.join(self.__workPath, taxIdPrefix)
            self.__fU.mkdir(taxIdPrefixDataDumpDir)
            blobs = bucket.list_blobs(prefix="proteomes/proteome-tax_id-" + taxIdPrefix)
            #
            for blob in blobs:
                archiveFile = blob.name.split("/")[-1]
                taxId = archiveFile.split("tax_id-")[-1].split("-")[0]
                taxIdAndIdx = archiveFile.split("tax_id-")[-1].split("_v")[0]
                logger.info("archiveFile: %r, taxId: %r, taxIdAndIdx: %r", archiveFile, taxId, taxIdAndIdx)
                # tD = {"taxId": taxId}
                if not archiveFile.endswith(".tar"):
                    logger.info("archiveFile %r not a tar file - skipping...", archiveFile)
                    continue
                archiveFileDumpPath = os.path.join(taxIdPrefixDataDumpDir, archiveFile)
                # tD.update({"data_directory": taxIdPrefixDataDumpDir})
                logger.info("Fetching file %s from Google Cloud bucket %s to local path %s", archiveFile, self.__bucketName, archiveFileDumpPath)
                blob.download_to_filename(archiveFileDumpPath)
                logger.info("Completed fetch at %s (%.4f seconds)", time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime)
                if taxIdPrefixDataDumpDir not in cacheD["data"]:
                    cacheD["data"].update({taxIdPrefixDataDumpDir: {"0": {"archive_files": {}}}})
                cacheD["data"][taxIdPrefixDataDumpDir]["0"]["archive_files"][archiveFile] = os.path.getsize(archiveFileDumpPath)

        except Exception as e:
            logger.exception("Failing on fetching of taxIdPrefix %s from FTP server, with message:\n%s", taxIdPrefix, str(e))

        return cacheD

    def getArchiveDirList(self):
        archiveDirList = [k for k in self.__oD]

        return archiveDirList

    def getArchiveDataDict(self):
        archivDataDict = copy.deepcopy(self.__oD)
        return archivDataDict

    def getArchiveFileList(self, inputPathList=None):
        """Return a list of filepaths for all tar files under the provided set of directories.

        Args:
            inputPathList (list): Optional (but recommended) list of paths to directories containing model mmCIF files to return.
                                  If not provided, method will retrieve all cached species data paths, so all model files from all species will be returned in a single list;
                                  thus, it's recommended to provide a list of specific species directory to break the returned model list down into more manageable parts.

        Returns:
            (list): list of model mmCIF file paths (only matches compressed ".cif.gz" files, to ensure only one file per model)
        """

        if not inputPathList:
            inputPathList = self.getArchiveDirList()

        archiveFileList = []

        for archiveDir in inputPathList:
            try:
                archiveFiles = glob.glob(os.path.join(archiveDir, "*.tar"))
                archiveFileList = [os.path.abspath(f) for f in archiveFiles]
            except Exception as e:
                logger.exception("Failing with %s", str(e))

        return archiveFileList

    def getSpeciesDataDownloadDate(self):
        return self.__createdDate

    def getBaseDataPath(self):
        return self.__workPath

    def getAFCloudTaxIdDataCacheFilePath(self):
        return self.__aFCTaxIdDataCacheFile

    def getComputedModelsDataPath(self):
        return self.__cachePath

    def getModelReorganizer(self, cachePath=None, useCache=True, workPath=None, **kwargs):
        cachePath = cachePath if cachePath else self.__cachePath
        workPath = workPath if workPath else self.__workPath
        return ModelReorganizer(cachePath=cachePath, useCache=useCache, workPath=workPath, **kwargs)

    def reorganizeModelFiles(self, cachePath=None, useCache=True, inputTaxIdPrefixList=None, **kwargs):
        """Reorganize model files from organism-wide model listing to hashed directory structure and rename files
        to follow internal identifier naming convention.

        The "model-download-cache.json" file (in the .../source-models/work-dir/...) used here has the following structure:
        {
          "data": {
            "/mnt/vdb1/source-models/work-dir/AlphaFoldCloud/100": {    <--- directory of tar files
              "0": {                                                    <--- subset of tar files (such that no single subset is > 64 GB)
                "archive_files": {
                  "proteome-tax_id-1000965-0_v4.tar": 95232,            <--- tar file name and size
                  "proteome-tax_id-1000960-0_v4.tar": 88576,
                  "proteome-tax_id-1000888-0_v4.tar": 128512,
                  "proteome-tax_id-1000974-0_v4.tar": 96256,
                  "proteome-tax_id-1007383-0_v4.tar": 40960,
                  ...
                "reorganized": False                                    <--- whether the subset has been reorganized yet; key is absent if not
        ...}}

        NOTE: This method does NOT delete the source archive TAR files--that task should be left up to the user to do manually,
        and only when everything is assured to be reorganized correctly and there is absolutely no need to keep the source tar files.
        Once they are deleted (e.g., all 100-999 taxId directories), the only way to get them back is to redownload them all again from GoogleCloud.
        The "keepSource" argument has a different meaning in this provider class (see description below).

        Args:
            cachePath (str): Path to cache directory.
            useCache (bool): Whether to use the existing data cache or re-run entire model reorganization process.
            inputTaxIdPrefixList (list, optional): List of input model filepaths to reorganize; defaults to all models for all species model sets.
            **kwargs (optional):
                numProc (int): number of processes to use; default 2.
                chunkSize (int): incremental chunk size used for distributed work processes; default 20.
                keepSource (bool): whether to copy files to new directory (instead of moving them); default False.
                                   Note that this only impacts the CIF files extracted from the TAR files--the TAR files are
                                   retained regardless. So, this should be set to False to prevent accumulation of many thousands of CIF files.
                cacheFilePath (str): full filepath and name for cache file containing a dictionary of all reorganized models.
                dictFilePathL (str, optional): List of dictionary files to use for BCIF encoding.
                smallFileSizeCutoff (int): size in bytes to use a file size cutoff for separating out "small" vs. "big" tar files

        Returns:
            (bool): True if successful; False otherwise.
        """
        try:
            ok = False
            #
            logger.info("Beginning reorganization with cachePath %r, useCache %r, inputTaxIdPrefixList %r, kwargs %r", cachePath, useCache, inputTaxIdPrefixList, kwargs)
            #
            smallFileSizeCutoff = kwargs.get("smallFileSizeCutoff", 8388608)  # 8mb
            smallFileSizeCutoffMb = int(smallFileSizeCutoff / (1024 * 1024))
            #
            cacheD = self.__mU.doImport(self.__aFCTaxIdDataCacheFile, fmt="json")
            cacheDataD = cacheD["data"]
            #
            reorgDataD = {}
            for archiveDir, archiveD in cacheDataD.items():
                taxIdGroup = archiveDir.split("/")[-1]
                #
                if inputTaxIdPrefixList:  # Only reorganize the archive files in the given list of taxId prefix groups
                    if taxIdGroup in inputTaxIdPrefixList:
                        reorgDataD[archiveDir] = archiveD
                else:
                    # Check if cache was already reorganized
                    for archiveSubsetIdx, archiveSubsetD in archiveD.items():
                        if archiveSubsetD.get("reorganized", False):
                            logger.info("TaxID archive data group %s already reorganized", archiveDir)
                            continue
                        reorgDataD.setdefault(archiveDir, {}).update({archiveSubsetIdx: archiveSubsetD})

            for archiveDir, archiveD in reorgDataD.items():
                for archiveSubsetIdx, archiveSubsetD in archiveD.items():
                    outputModelD = {}
                    smallArchiveFileL = [fn for fn, size in archiveSubsetD["archive_files"].items() if size <= smallFileSizeCutoff]
                    bigArchiveFileL = [fn for fn, size in archiveSubsetD["archive_files"].items() if size > smallFileSizeCutoff]
                    smallArchiveFilePathL = [os.path.join(archiveDir, archiveFile) for archiveFile in smallArchiveFileL]
                    bigArchiveFilePathL = [os.path.join(archiveDir, archiveFile) for archiveFile in bigArchiveFileL]
                    #
                    logger.info(
                        "archiveDir %s subset %r - %d small files (<= %rmb), %d large files (> %rmb), (total number of tar files %d)",
                        archiveDir,
                        archiveSubsetIdx,
                        len(smallArchiveFilePathL),
                        smallFileSizeCutoffMb,
                        len(bigArchiveFilePathL),
                        smallFileSizeCutoffMb,
                        len(archiveSubsetD["archive_files"]),
                    )

                    taxIdGroup = archiveDir.split("/")[-1]
                    holdingsFileName = "alphafold-holdings-" + taxIdGroup + "-" + str(archiveSubsetIdx) + ".json.gz"
                    mR = self.getModelReorganizer(cachePath=cachePath, useCache=useCache, workPath=archiveDir, cacheFile=holdingsFileName, cacheFormat="json", **kwargs)

                    # First, reorganize small archive files (<= 16mb) -- more efficient with multiple workers acting on separate tar files
                    if smallArchiveFilePathL:
                        logger.info("Reorganizing small archive files (<= %rmb) - %d files", smallFileSizeCutoffMb, len(smallArchiveFilePathL))
                        outputModelD, ok = mR.reorganize(
                            inputModelList=smallArchiveFilePathL,
                            modelSource="AlphaFoldCloud",
                            destBaseDir=self.__cachePath,
                            useCache=useCache,
                            inputModelD=outputModelD,
                            writeCache=False
                        )
                        if not ok:
                            logger.error("Reorganization of model files failed for archive %s subset %r", archiveDir, archiveSubsetIdx)

                    # Second, reorganize large archive files (> 16mb) -- more efficient with multiple workers acting on the same tar file
                    if bigArchiveFilePathL:
                        logger.info("Reorganizing large archive files (> %rmb) - %d files", smallFileSizeCutoffMb, len(bigArchiveFilePathL))
                        for archiveFile in bigArchiveFilePathL:
                            inputModelList = self.__extractModelCifFiles(archiveFile)
                            if inputModelList and len(inputModelList) > 0:
                                logger.info("Working on reorganizing %s (%d models)", archiveFile, len(inputModelList))
                                # Proceed with reorganization
                                outputModelD, ok = mR.reorganize(
                                    inputModelList=inputModelList,
                                    modelSource="AlphaFold",
                                    destBaseDir=self.__cachePath,
                                    useCache=useCache,
                                    inputModelD=outputModelD,
                                    writeCache=False
                                )
                                if not ok:
                                    logger.error("Reorganization of model files failed for archive %s subset %r", archiveDir, archiveSubsetIdx)
                                    break

                    # Last, write out the holdings/cache files
                    ok = mR.writeCacheFiles(outputModelD)
                    if not ok:
                        logger.error("Exporting of holdings cache files failed for archive %s subset %r", archiveDir, archiveSubsetIdx)
                        break

                    # Update the cache file to indicate that the given species archive has been reorganized
                    if ok and not cacheD["data"][archiveDir][archiveSubsetIdx].get("reorganized", False):
                        cacheD["data"][archiveDir][archiveSubsetIdx].update({"reorganized": True})
                        logger.info("Reorganization of model files complete for archiveDir %s", archiveDir)
                    ok = self.__mU.doExport(self.__aFCTaxIdDataCacheFile, cacheD, fmt="json", indent=4)
        #
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            ok = False
        #
        return ok

    def __extractModelCifFiles(self, archiveFile):
        modelPathL = []
        try:
            with tarfile.open(archiveFile, "r") as tF:
                tfL = tF.getnames()  # Get a list of items (files and directories) in the tar file
                logger.debug("Tar members in archiveFile %s: %r", archiveFile, tfL)
                #
                cifFiles = [file for file in tfL if file.endswith(".cif.gz")]  # only extract model files (not PAE and pLDDT files)
                # numModels = len(cifFiles)
                archiveFileDirPath = os.path.dirname(os.path.abspath(archiveFile))
                for cifFile in cifFiles:
                    fOutputPath = os.path.join(archiveFileDirPath, cifFile)
                    fIn = tF.extractfile(cifFile)
                    with open(fOutputPath, "wb") as ofh:
                        ofh.write(fIn.read())
                    modelPathL.append(fOutputPath)
        #
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        #
        return modelPathL
