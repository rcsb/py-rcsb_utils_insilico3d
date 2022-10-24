##
# File:    AlphaFoldModelCloudProvider.py
# Author:  Dennis Piehl
# Date:    4-Oct-2022
#
# Updates:
##
"""
Accessors for AlphaFold 3D Models (mmCIF) from public Google Cloud datasets.

"""

__docformat__ = "google en"
__author__ = "Dennis Piehl"
__email__ = "dennis.piehl@rcsb.org"
__license__ = "Apache 2.0"

import datetime
import logging
import os
import time
from pathlib import Path
import copy
import glob
# import re
from google.cloud import storage

from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.insilico3d.ModelReorganizer import ModelReorganizer

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger(__name__)


class AlphaFoldModelCloudProvider:
    """Accessors for AlphaFold models (mmCIF)."""

    def __init__(self, cachePath, baseWorkPath=None, useCache=False, reload=True, **kwargs):
        """Initialize AlphaFoldModelCloudProvider class.

        Args:
            cachePath (str): Path to directory where model files will be reorganized and stored permanently (also contains computed-model-cache.json file)  (i.e., 'computed-models'),
                             and will serve as the parent directory of where model files will be downloaded and where AF-specific cache file will sit (under 'work-dir/AlphaFold').
            useCache (bool, optional): Start from cache of already downloaded and/or reorganized model files. Defaults to False.
                                       When True, checks if last downloaded set of files is up-to-date and downloads any newly available models.
                                       When False (default), redownloads all model files.
            reload (bool, optional): Peform full reload (i.e., download/update) upon instantiation. Defaults to True.
        """
        # Use the same root cachePath for all types of insilico3D model sources, but with unique workPath names (sub-directory)
        self.__cachePath = cachePath  # Directory where model files will be reorganized and stored permanently (also contains computed-model-cache.json file) (i.e., 'computed-models')
        self.__baseWorkPath = baseWorkPath if baseWorkPath else self.__cachePath
        self.__workPath = os.path.join(self.__baseWorkPath, "work-dir", "AlphaFoldCloud")  # Directory where model files will be downloaded (also contains AF-specific cache file)
        self.__speciesDataCacheFile = os.path.join(self.__workPath, "model-download-cache.json")

        self.__bucketName = "public-datasets-deepmind-alphafold"

        self.__mU = MarshalUtil(workPath=self.__workPath)
        self.__fU = FileUtil(workPath=self.__workPath)

        if reload:
            self.__oD, self.__createdDate = self.__reload(useCache=useCache, **kwargs)

    def testCache(self, minCount=0):  # Increase minCount once we are consistently downloading more than one species data set
        if self.__oD and len(self.__oD) > minCount:
            return True
        else:
            return False

    def getCacheDirPath(self):
        return self.__workPath

    def reload(self, useCache, **kwargs):
        self.__oD, self.__createdDate = self.__reload(useCache=useCache, **kwargs)

    def __reload(self, **kwargs):
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
            useCache = kwargs.get("useCache", True)

            alphaFoldRequestedSpeciesList = kwargs.get("alphaFoldRequestedSpeciesList", [])
            if not alphaFoldRequestedSpeciesList:
                alphaFoldRequestedSpeciesList = [
                    {"species": "Panicum virgatum", "common_name": "Switchgrass", "taxIds": ["38727", "206033"]}
                ]

            self.__fU.mkdir(self.__workPath)

            logger.info("useCache %r self.__speciesDataCacheFile %r", useCache, self.__speciesDataCacheFile)
            if useCache and self.__mU.exists(self.__speciesDataCacheFile):
                logger.info("Loading data cache, %s.", self.__speciesDataCacheFile)
                cacheD = self.__mU.doImport(self.__speciesDataCacheFile, fmt="json")
                createdDate = cacheD["created"]
                oD = cacheD["data"]
                cacheArchiveFileList = [sF for sF in oD]
                logger.info("Checking consistency of cached data with data available on FTP")
                for archiveD in alphaFoldRequestedSpeciesList:
                    try:
                        speciesName = archiveD.get("species", archiveD.get("label", None))
                        archiveFile = archiveD["archive_name"]
                        archiveFileSize = int(archiveD["size_bytes"])
                        speciesNumModels = int(archiveD["num_predicted_structures"])
                        if speciesName not in cacheArchiveFileList:
                            logger.info("Species archive data for %s not found in local cache. Will re-fetch.", speciesName)
                            cacheD = self.fetchSpeciesArchive(archiveD, cacheD)
                            oD = cacheD["data"]
                            ok = self.__mU.doExport(self.__speciesDataCacheFile, cacheD, fmt="json", indent=4)
                            continue
                        #
                        # Check if cache was already reorganized
                        reorganized = archiveD.get("reorganized", False)
                        reorganizedBaseDir = archiveD.get("reorganizedBaseDir", None)
                        if reorganized and reorganizedBaseDir is not None:
                            if self.__mU.exists(reorganizedBaseDir):
                                logger.info("Species archive data for %s already reorganized to: %s", speciesName, reorganizedBaseDir)
                        #
                        cacheArchiveDir = oD[speciesName]["data_directory"]
                        cacheArchiveFileSize = oD[speciesName]["size_bytes"]
                        cacheSpeciesNumModels = oD[speciesName]["num_predicted_structures"]
                        if not os.path.exists(cacheArchiveDir) and not reorganized:
                            logger.warning("Species archive data directory for %s not found at cached path %s", archiveFile, cacheArchiveDir)
                        if cacheArchiveFileSize != archiveFileSize:
                            logger.warning("Species archive data file %s not up-to-date with file available on FTP server.", archiveFile)
                        if cacheSpeciesNumModels != speciesNumModels:
                            logger.warning(
                                "Missing some or all of the bundled model files for species archive %s as available on FTP server (found %d / %d model files)",
                                archiveFile,
                                cacheSpeciesNumModels,
                                speciesNumModels,
                            )
                    except Exception as e:
                        logger.exception("Failing on checking of cache data for %s from FTP server, with message:\n%s", archiveD["archive_name"], str(e))

            else:
                logger.info("Refetching all files from server.")
                cacheD = {}
                cacheD.update({"created": startDateTime, "data": {}})
                for archiveD in alphaFoldRequestedSpeciesList:
                    cacheD = self.fetchSpeciesArchive(archiveD, cacheD)
                createdDate = cacheD["created"]
                oD = cacheD["data"]
                ok = self.__mU.doExport(self.__speciesDataCacheFile, cacheD, fmt="json", indent=4)
                logger.info("Export AlphaFold species model data (%d) status %r", len(oD), ok)

        except Exception as e:
            logger.exception("Failing with %s", str(e))

        return oD, createdDate

    def fetchSpeciesArchive(self, archiveD, cacheD):
        try:
            sD = copy.deepcopy(archiveD)
            startTime = time.time()
            client = storage.Client()
            bucket = client.bucket(self.__bucketName)
            #
            speciesName = sD.get("species", sD.get("label", None))
            speciesDataDumpDir = os.path.join(self.__workPath, speciesName.replace(" ", "_").replace("(", "").replace(")", ""))
            self.__fU.mkdir(speciesDataDumpDir)
            taxIdL = sD.get("taxIds")
            for taxId in taxIdL:
                blobs = bucket.list_blobs(prefix="proteomes/proteome-tax_id-" + taxId + "-")
                for blob in blobs:
                    archiveFile = blob.name.split("/")[-1]
                    if not archiveFile.endswith(".tar"):
                        logger.info("archiveFile %r not a tar file - skipping...", archiveFile)
                        continue
                    archiveFileDumpPath = os.path.join(speciesDataDumpDir, archiveFile)
                    sD.update({"archive_name": archiveFile, "data_directory": speciesDataDumpDir, "archive_file_path": archiveFileDumpPath})
                    logger.info("Fetching file %s from Google Cloud bucket %s to local path %s", archiveFile, self.__bucketName, archiveFileDumpPath)
                    blob.download_to_filename(archiveFileDumpPath)
                    ok = self.__fU.unbundleTarfile(archiveFileDumpPath, dirPath=speciesDataDumpDir)
                    logger.info("Completed fetch (%r) at %s (%.4f seconds)", ok, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime)

                    logger.info("Clearing confidence and PAE files from extracted tar bundle...")
                    for pdbDumpFile in Path(speciesDataDumpDir).glob("*.json.gz"):
                        pdbDumpFile.unlink()

                    if ok:
                        self.__fU.remove(archiveFileDumpPath)

                        numFiles = len(os.listdir(speciesDataDumpDir))
                        sD.update({"num_predicted_structures": numFiles})
                        cacheD["data"].update({speciesName: sD})

        except Exception as e:
            logger.exception("Failing on fetching and expansion of file %s from FTP server, with message:\n%s", archiveD["archive_name"], str(e))

        return cacheD

    def getArchiveDirList(self):
        archiveDirList = [self.__oD[k]["data_directory"] for k in self.__oD]

        return archiveDirList

    def getArchiveDataDict(self):
        archivDataDict = copy.deepcopy(self.__oD)
        return archivDataDict

    def getModelFileList(self, inputPathList=None):
        """Return a list of filepaths for all mmCIF models under the provided set of directories.

        Args:
            inputPathList (list): Optional (but recommended) list of paths to directories containing model mmCIF files to return.
                                  If not provided, method will retrieve all cached species data paths, so all model files from all species will be returned in a single list;
                                  thus, it's recommended to provide a list of specific species directory to break the returned model list down into more manageable parts.

        Returns:
            (list): list of model mmCIF file paths (only matches compressed ".cif.gz" files, to ensure only one file per model)
        """

        if not inputPathList:
            inputPathList = self.getArchiveDirList()

        modelFileList = []

        for modelDir in inputPathList:
            try:
                modelFiles = glob.glob(os.path.join(modelDir, "*.cif.gz"))
                modelFileList = [os.path.abspath(f) for f in modelFiles]
            except Exception as e:
                logger.exception("Failing with %s", str(e))

        return modelFileList

    def getSpeciesDataDownloadDate(self):
        return self.__createdDate

    def getBaseDataPath(self):
        return self.__workPath

    def getSpeciesDataCacheFilePath(self):
        return self.__speciesDataCacheFile

    def getModelReorganizer(self, cachePath=None, useCache=True, workPath=None, **kwargs):
        cachePath = cachePath if cachePath else self.__cachePath
        workPath = workPath if workPath else self.__workPath
        return ModelReorganizer(cachePath=cachePath, useCache=useCache, workPath=workPath, **kwargs)

    def getComputedModelsDataPath(self):
        return self.__cachePath

    def reorganizeModelFiles(self, cachePath=None, useCache=True, inputModelList=None, **kwargs):
        """Reorganize model files from organism-wide model listing to hashed directory structure and rename files
        to follow internal identifier naming convention.

        Args:
            cachePath (str): Path to cache directory.
            inputModelList (list, optional): List of input model filepaths to reorganize; defaults to all models for all species model sets.
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
                ok = mR.reorganize(inputModelList=inputModelList, modelSource="AlphaFold", destBaseDir=self.__cachePath, useCache=useCache)
                if not ok:
                    logger.error("Reorganization of model files failed for inputModelList starting with item, %s", inputModelList[0])
            #
            else:  # Reorganize ALL model files for ALL available species model sets
                cacheD = self.__mU.doImport(self.__speciesDataCacheFile, fmt="json")
                archiveDataD = self.getArchiveDataDict()
                for species, archiveD in archiveDataD.items():
                    archiveDir = archiveD["data_directory"]
                    # First check if cache was already reorganized
                    if cacheD["data"][species]["data_directory"] == archiveDir:
                        reorganized = archiveD.get("reorganized", False)
                        reorganizedBaseDir = archiveD.get("reorganizedBaseDir", None)
                        if reorganized and reorganizedBaseDir is not None:
                            if self.__mU.exists(reorganizedBaseDir) and reorganizedBaseDir == self.__cachePath:
                                logger.info("Species archive data for %s already reorganized to: %s", species, reorganizedBaseDir)
                                ok = True
                                continue
                    # Proceed with reorganization
                    inputModelList = self.getModelFileList(inputPathList=[archiveDir])
                    ok = mR.reorganize(inputModelList=inputModelList, modelSource="AlphaFold", destBaseDir=self.__cachePath, useCache=useCache)
                    if not ok:
                        logger.error("Reorganization of model files failed for species archive %s", archiveDir)
                        break
                    # Update the cache file to indicate that the given species archive has been reorganized
                    if cacheD["data"][species]["data_directory"] == archiveDir:
                        cacheD["data"][species].update({"reorganized": True, "reorganizedBaseDir": self.__cachePath})
                        logger.info("Reorganization of model files complete for species, %s, from archive %s", species, archiveDir)
                    ok = self.__mU.doExport(self.__speciesDataCacheFile, cacheD, fmt="json", indent=4)
        #
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            ok = False
        #
        return ok
