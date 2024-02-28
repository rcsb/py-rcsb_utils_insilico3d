##
# File:    AlphaFoldModelProvider.py
# Author:  Dennis Piehl
# Date:    30-Sep-2021
#
# Updates:
#  15-Dec-2021 dwp Re-introduce use of FTP instead of HTTP for downloading files; proved to be significantly faster
#  11-Mar-2022 dwp Move reorganizing method to separate class to enable multiprocessing and to make it common to other provider classes;
#                  During reorganizing process, also rename the file to use an internal identifier instead of the original source filename;
#                  Change cache file to record new internal identifier, original filename, and accession URL for each model file;
#                  Add usage of config file object for specifying location for storing model files
#
# To Do:
# - Add check that converted files are consistent with mmCIF dictionaries
##
"""
Accessors for AlphaFold 3D Models (mmCIF).

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
import copy
import glob
import re

from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.io.FtpUtil import FtpUtil
from rcsb.utils.insilico3d.ModelReorganizer import ModelReorganizer

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger(__name__)


class AlphaFoldModelProvider:
    """Accessors for AlphaFold models (mmCIF)."""

    def __init__(self, cachePath, baseWorkPath=None, useCache=False, reload=True, **kwargs):
        """Initialize AlphaFoldModelProvider class.

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
        self.__workPath = os.path.join(self.__baseWorkPath, "work-dir", "AlphaFold")  # Directory where model files will be downloaded (also contains AF-specific cache file)
        self.__speciesDataCacheFile = os.path.join(self.__workPath, "model-download-cache.json")
        self.__dataSetHoldingsFileName = "alphafold-ftp-holdings.json.gz"

        self.__ftpHost = kwargs.get("ftpHost", "ftp.ebi.ac.uk")
        self.__ftpDataPath = kwargs.get("ftpDataPath", "/pub/databases/alphafold/")

        self.__mU = MarshalUtil(workPath=self.__workPath)
        self.__fU = FileUtil(workPath=self.__workPath)
        self.__ftpU = FtpUtil(workPath=self.__workPath)

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

            # alphaFoldLatestDataList = "https://ftp.ebi.ac.uk/pub/databases/alphafold/download_metadata.json"
            alphaFoldLatestDataList = "pub/databases/alphafold/download_metadata.json"
            alphaFoldRequestedSpeciesList = kwargs.get("alphaFoldRequestedSpeciesList", [])
            excludeArchiveFileRegexList = ["swissprot_pdb_v[0-9]+.tar"]
            excludeArchiveFileRegexListCombined = "(?:% s)" % "|".join(excludeArchiveFileRegexList)

            self.__ftpU.connect(self.__ftpHost)
            self.__fU.mkdir(self.__workPath)

            latestDataListDumpPath = os.path.join(self.__workPath, self.__fU.getFileName(alphaFoldLatestDataList))
            ok = self.__ftpU.get(alphaFoldLatestDataList, latestDataListDumpPath)
            lDL = self.__mU.doImport(latestDataListDumpPath, fmt="json")

            # Exclude undesired archives (defined in excludeArchiveFileRegexList)
            lDL = [s for s in lDL if not re.match(excludeArchiveFileRegexListCombined, s["archive_name"])]

            # If a specific list of species files was requested, only iterate over those
            if alphaFoldRequestedSpeciesList:
                alphaFoldArchiveDataList = [s for s in lDL if s.get("species", s.get("label", None)) in alphaFoldRequestedSpeciesList]
            else:
                alphaFoldArchiveDataList = lDL

            logger.info("useCache %r self.__speciesDataCacheFile %r", useCache, self.__speciesDataCacheFile)
            if useCache and self.__mU.exists(self.__speciesDataCacheFile):
                logger.info("Loading data cache, %s.", self.__speciesDataCacheFile)
                cacheD = self.__mU.doImport(self.__speciesDataCacheFile, fmt="json")
                createdDate = cacheD["created"]
                oD = cacheD["data"]
                cacheArchiveFileList = [sF for sF in oD]
                logger.info("Checking consistency of cached data with data available on FTP")
                for archiveD in alphaFoldArchiveDataList:
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
                for archiveD in alphaFoldArchiveDataList:
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
            speciesName = sD.get("species", sD.get("label", None))
            archiveFile = sD["archive_name"]
            archiveFilePath = os.path.join(self.__ftpDataPath, "latest", archiveFile)
            speciesDataDumpDir = os.path.join(self.__workPath, speciesName.replace(" ", "_").replace("(", "").replace(")", ""))
            self.__fU.mkdir(speciesDataDumpDir)
            archiveFileDumpPath = os.path.join(speciesDataDumpDir, archiveFile)
            sD.update({"data_directory": speciesDataDumpDir, "archive_file_path": archiveFileDumpPath})

            logger.info("Fetching file %s from FTP server to local path %s", archiveFilePath, archiveFileDumpPath)
            ok = self.__ftpU.get(archiveFilePath, archiveFileDumpPath)
            ok = self.__fU.unbundleTarfile(archiveFileDumpPath, dirPath=speciesDataDumpDir)
            logger.info("Completed fetch (%r) at %s (%.4f seconds)", ok, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime)

            logger.info("Clearing PDB files from extracted tar bundle...")
            for pdbDumpFile in Path(speciesDataDumpDir).glob("*.pdb.gz"):
                pdbDumpFile.unlink()

            if ok:
                cacheD["data"].update({speciesName: sD})
                self.__fU.remove(archiveFileDumpPath)

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
                _, ok = mR.reorganize(inputModelList=inputModelList, modelSource="AlphaFold", destBaseDir=self.__cachePath, useCache=useCache)
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
                    _, ok = mR.reorganize(inputModelList=inputModelList, modelSource="AlphaFold", destBaseDir=self.__cachePath, useCache=useCache)
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
