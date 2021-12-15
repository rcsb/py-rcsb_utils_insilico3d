##
# File:    AlphaFoldModelProvider.py
# Author:  Dennis Piehl
# Date:    30-Sep-2021
#
# Update:
#  30-Sep-2021 dwp Started class
#  15-Dec-2021 dwp Re-introduce use of FTP instead of HTTP for downloading files; proved to be significantly faster
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

logger = logging.getLogger(__name__)


class AlphaFoldModelProvider:
    """Accessors for AlphaFold models (mmCIF)."""

    def __init__(self, **kwargs):
        # Use the same root cachePath for all types of insilico3D model sources, but with unique dirPath names (sub-directory)
        self.__cachePath = kwargs.get("cachePath", "./CACHE-insilico3d-models")
        self.__dirPath = os.path.join(self.__cachePath, "AlphaFold")
        self.__speciesDataCacheFile = os.path.join(self.__dirPath, "species-model-data.json")
        self.__dividedDataPath = os.path.join(self.__cachePath, "computed-models")

        self.__mU = MarshalUtil(workPath=self.__dirPath)
        self.__fU = FileUtil(workPath=self.__dirPath)
        self.__ftpU = FtpUtil(workPath=self.__dirPath)

        self.__oD, self.__createdDate = self.__reload(**kwargs)

    def testCache(self, minCount=0):  # Increase minCount once we are consistently downloading more than one species data set
        if self.__oD and len(self.__oD) > minCount:
            return True
        else:
            return False

    def getCacheDirPath(self):
        return self.__dirPath

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
            startTime = time.time()
            startDateTime = datetime.datetime.now().isoformat()
            useCache = kwargs.get("useCache", True)

            alphaFoldFtpHost = kwargs.get("alphaFoldFtpHost", "ftp.ebi.ac.uk")
            alphaFoldFtpDataPath = kwargs.get("alphaFoldFtpDataPath", "/pub/databases/alphafold/")

            alphaFoldBaseUrl = kwargs.get("alphaFoldBaseUrl", "https://ftp.ebi.ac.uk/pub/databases/alphafold/")
            alphaFoldLatestDataList = os.path.join(alphaFoldBaseUrl, "download_metadata.json")
            alphaFoldRequestedSpeciesList = kwargs.get("alphaFoldRequestedSpeciesList", [])
            excludeArchiveFileRegexList = ["swissprot_pdb_v[0-9]+.tar"]
            excludeArchiveFileRegexListCombined = "(?:% s)" % "|".join(excludeArchiveFileRegexList)

            self.__ftpU.connect(alphaFoldFtpHost)
            self.__fU.mkdir(self.__dirPath)

            latestDataListDumpPath = os.path.join(self.__dirPath, self.__fU.getFileName(alphaFoldLatestDataList))
            ok = self.__fU.get(alphaFoldLatestDataList, latestDataListDumpPath)
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
                for archiveData in alphaFoldArchiveDataList:
                    try:
                        speciesName = archiveData.get("species", archiveData.get("label", None))
                        archiveFile = archiveData["archive_name"]
                        archiveFileSize = int(archiveData["size_bytes"])
                        speciesNumModels = int(archiveData["num_predicted_structures"])
                        if speciesName not in cacheArchiveFileList:
                            logger.error("Species archive data file %s on FTP server not found in local cache.", archiveFile)
                            continue
                        cacheArchiveDir = oD[speciesName]["data_directory"]
                        cacheArchiveFileSize = oD[speciesName]["size_bytes"]
                        cacheSpeciesNumModels = oD[speciesName]["num_predicted_structures"]
                        if not os.path.exists(cacheArchiveDir):
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
                        logger.exception("Failing on checking of cache data for %s from FTP server, with message:\n%s", archiveData["archive_name"], str(e))

            else:
                logger.info("Refetching all files from server.")
                cacheD = {}
                cacheD.update({"created": startDateTime, "data": {}})
                for archiveData in alphaFoldArchiveDataList:
                    try:
                        sD = copy.deepcopy(archiveData)
                        speciesName = sD.get("species", sD.get("label", None))
                        archiveFile = sD["archive_name"]
                        archiveFilePath = os.path.join(alphaFoldFtpDataPath, "latest", archiveFile)
                        speciesDataDumpDir = os.path.join(self.__dirPath, speciesName.replace(" ", "_").replace("(", "").replace(")", ""))
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
                        logger.exception("Failing on fetching and expansion of file %s from FTP server, with message:\n%s", archiveData["archive_name"], str(e))

                createdDate = cacheD["created"]
                oD = cacheD["data"]
                ok = self.__mU.doExport(self.__speciesDataCacheFile, cacheD, fmt="json", indent=3)
                logger.info("Export AlphaFold species model data (%d) status %r", len(oD), ok)

        except Exception as e:
            logger.exception("Failing with %s", str(e))

        return oD, createdDate

    def getArchiveDirList(self):
        archiveDirList = [self.__oD[k]["data_directory"] for k in self.__oD]

        return archiveDirList

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
        return self.__dirPath

    def getSpeciesDataCacheFilePath(self):
        return self.__speciesDataCacheFile

    def reorganizeModelFiles(self):
        """Move model files from organism-wide model listing to hashed directory structure constructed
        from the 6-character UniProt ID (e.g., "P52078" will be moved to "./P5/20/78")"""

        try:
            archiveDirList = self.getArchiveDirList()
            for archiveDir in archiveDirList:
                newModelDirD = {}
                dividedDataCacheFile = os.path.join(archiveDir, "species-model-files.json")
                modelFileList = self.getModelFileList(inputPathList=[archiveDir])
                for model in modelFileList:
                    modelName = self.__fU.getFileName(model)
                    uniProtID = modelName.split(".cif.gz")[0].split("-")[1]
                    first2 = uniProtID[0:2]
                    mid2 = uniProtID[2:4]
                    last2 = uniProtID[4:6]
                    destDir = os.path.join(self.__dividedDataPath, first2, mid2, last2)
                    if not self.__fU.exists(destDir):
                        self.__fU.mkdir(destDir)
                    destModelPath = os.path.join(destDir, modelName)
                    self.__fU.replace(model, destModelPath)
                    newModelDirD[modelName] = destModelPath
            self.__mU.doExport(dividedDataCacheFile, newModelDirD, fmt="json", indent=3)
            return True
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            return False
