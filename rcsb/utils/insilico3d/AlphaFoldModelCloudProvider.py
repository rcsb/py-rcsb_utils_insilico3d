##
# File:    AlphaFoldModelCloudProvider.py
# Author:  Dennis Piehl
# Date:    4-Oct-2022
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
from pathlib import Path
import copy
import glob
# import re
import tarfile
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
        self.__AFCloudTaxIdDataCacheFile = os.path.join(self.__workPath, "model-download-cache.json")

        self.__bucketName = "public-datasets-deepmind-alphafold-v4"

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

            # alphaFoldRequestedSpeciesList = kwargs.get("alphaFoldRequestedSpeciesList", [])
            # if not alphaFoldRequestedSpeciesList:
            #     alphaFoldRequestedSpeciesList = [
            #         {"species": "Panicum virgatum", "common_name": "Switchgrass", "taxIds": ["38727", "206033"]},
            #     ]

            alphaFoldRequestedTaxIdPrefixList = kwargs.get("alphaFoldRequestedTaxIdPrefixList", [])
            if not alphaFoldRequestedTaxIdPrefixList:
                alphaFoldRequestedTaxIdPrefixList = [
                    "206033",
                    "100000",
                    # "408170",
                ]

            self.__fU.mkdir(self.__workPath)

            logger.info("useCache %r self.__AFCloudTaxIdDataCacheFile %r", useCache, self.__AFCloudTaxIdDataCacheFile)
            if useCache and self.__mU.exists(self.__AFCloudTaxIdDataCacheFile):
                logger.info("Loading data cache, %s.", self.__AFCloudTaxIdDataCacheFile)
                cacheD = self.__mU.doImport(self.__AFCloudTaxIdDataCacheFile, fmt="json")
                createdDate = cacheD["created"]
                oD = cacheD["data"]
            else:
                logger.info("Refetching all files from server.")
                cacheD = {}
                cacheD.update({"created": startDateTime, "data": {}})
                for taxIdPrefix in alphaFoldRequestedTaxIdPrefixList:
                    cacheD = self.fetchTaxIdArchive(taxIdPrefix, cacheD)
                createdDate = cacheD["created"]
                oD = cacheD["data"]
                ok = self.__mU.doExport(self.__AFCloudTaxIdDataCacheFile, cacheD, fmt="json", indent=4)
                logger.info("Export AlphaFold species model data (%d) status %r", len(oD), ok)

        except Exception as e:
            logger.exception("Failing with %s", str(e))

        return oD, createdDate

    def fetchTaxIdArchive(self, taxIdPrefix, cacheD):
        try:
            startTime = time.time()
            client = storage.Client()
            bucket = client.bucket(self.__bucketName)
            #
            taxIdPrefixDataDumpDir = os.path.join(self.__workPath, taxIdPrefix)
            self.__fU.mkdir(taxIdPrefixDataDumpDir)
            blobs = bucket.list_blobs(prefix="proteomes/proteome-tax_id-" + taxIdPrefix)
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
                    cacheD["data"].update({taxIdPrefixDataDumpDir: {"archive_files": []}})
                cacheD["data"][taxIdPrefixDataDumpDir]["archive_files"].append(archiveFile)

                # logger.info("unbundling archiveFileDumpPath %r", archiveFileDumpPath)
                # ok = self.__fU.unbundleTarfile(archiveFileDumpPath, dirPath=taxIdPrefixDataDumpDir)
                # logger.info("Clearing confidence and PAE files from extracted tar bundle...")
                # for pdbDumpFile in Path(taxIdPrefixDataDumpDir).glob("*.json.gz"):
                #     pdbDumpFile.unlink()

                # numModels = None
                # ok = False
                # with tarfile.open(archiveFileDumpPath, "r") as tF:
                #     # Get a list of items (files and directories) in the tar file
                #     fL = tF.getnames()
                #     logger.info("tarFile items: %r", fL)
                #     fL = [f for f in fL if f.endswith(".cif.gz")]  # only extract model files (not PAE and pLDDT files)
                #     numModels = len(fL)
                #     # DON'T EXTRACT THEM YET--ONLY DO THIS WHEN RE-ORGANIZING
                #     # for fName in fL:
                #     #     fOutputPath = os.path.join(taxIdPrefixDataDumpDir, fName)
                #     #     ok = self.extractTarMember(archiveFileDumpPath, fName, fOutputPath)
                #     #     if not ok:
                #     #         logger.error("Failed to extract member %r from archiveFileDumpPath %r", fName, archiveFileDumpPath)
                #     #         ok = False
                #     #         break
                # if numModels:
                #     ok = True
                # if ok:
                #     # self.__fU.remove(archiveFileDumpPath)
                #     # tD.update({"num_predicted_structures": numModels})
                #     cacheD["data"].setdefault(taxIdPrefixDataDumpDir, []).append(archiveFile)
                #     # cacheD["data"].update({archiveFile: tD})

        except Exception as e:
            logger.exception("Failing on fetching of taxIdPrefix %s from FTP server, with message:\n%s", taxIdPrefix, str(e))

        return cacheD

    def extractTarMember(self, tarFilePath, memberName, memberPath):
        ret = True
        try:
            with tarfile.open(tarFilePath) as tar:
                fIn = tar.extractfile(memberName)
                with open(memberPath, "wb") as ofh:
                    ofh.write(fIn.read())
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            ret = False
        return ret

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
        return self.__AFCloudTaxIdDataCacheFile

    def getModelReorganizer(self, cachePath=None, useCache=True, workPath=None, **kwargs):
        cachePath = cachePath if cachePath else self.__cachePath
        workPath = workPath if workPath else self.__workPath
        return ModelReorganizer(cachePath=cachePath, useCache=useCache, workPath=workPath, **kwargs)

    def getComputedModelsDataPath(self):
        return self.__cachePath

    # def reorganizeModelFiles(self, cachePath=None, useCache=True, inputModelList=None, **kwargs):
    #     """Reorganize model files from organism-wide model listing to hashed directory structure and rename files
    #     to follow internal identifier naming convention.

    #     Args:
    #         cachePath (str): Path to cache directory.
    #         inputModelList (list, optional): List of input model filepaths to reorganize; defaults to all models for all species model sets.
    #         **kwargs (optional):
    #             numProc (int): number of processes to use; default 2.
    #             chunkSize (int): incremental chunk size used for distribute work processes; default 20.
    #             keepSource (bool): whether to copy files to new directory (instead of moving them); default False.
    #             cacheFilePath (str): full filepath and name for cache file containing a dictionary of all reorganized models.

    #     Returns:
    #         (bool): True if successful; False otherwise.
    #     """
    #     try:
    #         ok = False
    #         #
    #         mR = self.getModelReorganizer(cachePath=cachePath, useCache=useCache, **kwargs)
    #         #
    #         if inputModelList:  # Only reorganize given list of model files
    #             ok = mR.reorganize(inputModelList=inputModelList, modelSource="AlphaFold", destBaseDir=self.__cachePath, useCache=useCache)
    #             if not ok:
    #                 logger.error("Reorganization of model files failed for inputModelList starting with item, %s", inputModelList[0])
    #         #
    #         else:  # Reorganize ALL model files for ALL available species model sets
    #             cacheD = self.__mU.doImport(self.__speciesDataCacheFile, fmt="json")
    #             archiveDataD = self.getArchiveDataDict()
    #             for species, archiveD in archiveDataD.items():
    #                 archiveDir = archiveD["data_directory"]
    #                 # First check if cache was already reorganized
    #                 if cacheD["data"][species]["data_directory"] == archiveDir:
    #                     reorganized = archiveD.get("reorganized", False)
    #                     reorganizedBaseDir = archiveD.get("reorganizedBaseDir", None)
    #                     if reorganized and reorganizedBaseDir is not None:
    #                         if self.__mU.exists(reorganizedBaseDir) and reorganizedBaseDir == self.__cachePath:
    #                             logger.info("Species archive data for %s already reorganized to: %s", species, reorganizedBaseDir)
    #                             ok = True
    #                             continue
    #                 # Proceed with reorganization
    #                 inputModelList = self.getModelFileList(inputPathList=[archiveDir])
    #                 ok = mR.reorganize(inputModelList=inputModelList, modelSource="AlphaFold", destBaseDir=self.__cachePath, useCache=useCache)
    #                 if not ok:
    #                     logger.error("Reorganization of model files failed for species archive %s", archiveDir)
    #                     break
    #                 # Update the cache file to indicate that the given species archive has been reorganized
    #                 if cacheD["data"][species]["data_directory"] == archiveDir:
    #                     cacheD["data"][species].update({"reorganized": True, "reorganizedBaseDir": self.__cachePath})
    #                     logger.info("Reorganization of model files complete for species, %s, from archive %s", species, archiveDir)
    #                 ok = self.__mU.doExport(self.__speciesDataCacheFile, cacheD, fmt="json", indent=4)
    #     #
    #     except Exception as e:
    #         logger.exception("Failing with %s", str(e))
    #         ok = False
    #     #
    #     return ok

    def extractAndReorganizeModelFiles(self, cachePath=None, useCache=True, inputArchiveList=None, **kwargs):
        """Extract and reorganize model files from organism-wide model listing to hashed directory structure and rename files
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
            # if inputArchiveList:  # Only reorganize given list of model files
            #     ok = mR.reorganize(inputModelList=inputModelList, modelSource="AlphaFoldCloud", destBaseDir=self.__cachePath, useCache=useCache)
            #     if not ok:
            #         logger.error("Reorganization of model files failed for inputModelList starting with item, %s", inputArchiveList[0])
            #
            # else:  # Reorganize ALL model files for ALL available species model sets
            #
            cacheD = self.__mU.doImport(self.__AFCloudTaxIdDataCacheFile, fmt="json")

            cacheDataD = cacheD["data"]
            logger.info("cacheD keys: %r", cacheDataD.keys())
            # for archiveDir, archiveFileL in cacheDataD.items():
            for archiveDir, archiveD in cacheDataD.items():
                inputPathList = []
                logger.info("ARCHIVE DIR: %r", archiveDir)
                archiveFileL = archiveD["archive_files"]
                #
                # First check if cache was already reorganized
                reorganized = archiveD.get("reorganized", False)
                # reorganizedBaseDir = archiveD.get("reorganizedBaseDir", None)
                if reorganized:
                    logger.info("Species archive data for %s already reorganized", archiveDir)
                    ok = True
                    continue
                #
                # Proceed with reorganization
                for archiveFile in archiveFileL:
                    logger.info("ARCHIVE FILE: %r", archiveFile)
                    with tarfile.open(os.path.join(archiveDir, archiveFile), "r") as tar:
                        cifFiles = [file for file in tar.getnames() if file.endswith(".cif.gz")]
                        for cifFile in cifFiles:
                            tar.extract(cifFile, path=archiveDir)  # Extract all cif files into a given taxIdPrefix directory (e.g., "100/")
                            inputPathList.append(os.path.join(archiveDir, cifFile))
                logger.info("inputPathList[0:3]: %r", inputPathList[0:3])
                ok = mR.reorganize(inputModelList=inputPathList, modelSource="AlphaFoldCloud", destBaseDir=self.__cachePath, useCache=useCache)
                if not ok:
                    logger.error("Reorganization of model files failed for species archive %s", archiveDir)
                    break
                # Now delete the source model files
                for inFile in inputPathList:
                    os.remove(inFile)

                # Update the cache file to indicate that the given species archive has been reorganized
                if ok and not cacheD["data"][archiveDir].get("reorganized", False):
                    cacheD["data"][archiveDir].update({"reorganized": True})
                    logger.info("Reorganization of model files complete for archive %s", archiveDir)
                ok = self.__mU.doExport(self.__AFCloudTaxIdDataCacheFile, cacheD, fmt="json", indent=4)
        #
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            ok = False
        #
        return ok

# import os
# from rcsb.utils.io.MarshalUtil import MarshalUtil

# workDir = "/mnt/csm/source-models/work-dir"

# def createArchiveDirFileDict(baseDir):
#     # Get a dictionary of all subdirectories and files in the given directory
#     archiveDirFileD = {}
#     absBaseDir = os.path.abspath(baseDir)
#     with os.scandir(absBaseDir) as fObjs:
#         for fObj in fObjs:
#             if fObj.is_dir():
#                 archiveDir = os.path.join(absBaseDir, fObj.name)
#                 archiveFileL = []
#                 with os.scandir(archiveDir) as archiveObjs:
#                     for archiveObj in archiveObjs:
#                         if archiveObj.is_file() and archiveObj.name.endswith(".tar"):
#                             archiveFileL.append(archiveObj.name)
#                 archiveDirFileD[archiveDir] = archiveFileL
#     return archiveDirFileD

# cacheD = createArchiveDirFileDict(workDir)
# mU = MarshalUtil()
# outFile = os.path.join(workDir, "model-download-cache.json")
# mU.doExport(outFile, cacheD, fmt="json")







