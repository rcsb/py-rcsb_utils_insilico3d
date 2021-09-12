##
# File:    AlphaFoldModelProvider.py
# Author:  Dennis Piehl
# Date:    23-Aug-2021
#
# Update:
#
#
##
"""
Accessors for AlphaFold 3D Models (mmCIF).

"""

import datetime
import logging
import os.path
import time
from pathlib import Path
import copy
import glob

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

        self.__mU = MarshalUtil(workPath=self.__dirPath)

        self.__oD, self.__createdDate = self.__reload(**kwargs)

    def testCache(self, minCount=0):  # Increase minCount once we are consistently downloading more than one species data set
        if self.__oD and len(self.__oD) > minCount:
            return True
        else:
            return False

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
            alphaFoldFtpLatestDataList = os.path.join(alphaFoldFtpDataPath, "download_metadata.json")
            alphaFoldRequestedSpeciesFileList = kwargs.get("alphaFoldRequestedSpeciesFileList", [])

            ftpU = FtpUtil()
            ftpU.connect(alphaFoldFtpHost)

            fU = FileUtil()
            fU.mkdir(self.__dirPath)

            latestDataListDumpPath = os.path.join(self.__dirPath, alphaFoldFtpLatestDataList.split("/")[-1])
            ok = ftpU.get(alphaFoldFtpLatestDataList, latestDataListDumpPath)
            lDL = self.__mU.doImport(latestDataListDumpPath, fmt="json")

            logger.info("useCache %r self.__speciesDataCacheFile %r", useCache, self.__speciesDataCacheFile)
            if useCache and self.__mU.exists(self.__speciesDataCacheFile):
                logger.info("Loading data cache, %s.", self.__speciesDataCacheFile)
                cacheD = self.__mU.doImport(self.__speciesDataCacheFile, fmt="json")
                createdDate = cacheD["created"]
                oD = cacheD["data"]
                cacheSpeciesFileList = [sF for sF in oD]

                logger.info("Checking consistency of cached data with data available on FTP")
                for speciesData in lDL:
                    try:
                        speciesFile = speciesData["archive_name"]
                        speciesFileSize = int(speciesData["size_bytes"])
                        speciesNumModels = int(speciesData["num_predicted_structures"])
                        # If a specific list of species files was requested, only check those (i.e., skip over any that weren't requested)
                        if alphaFoldRequestedSpeciesFileList and speciesFile not in alphaFoldRequestedSpeciesFileList:
                            logger.info("Skipping species file not specifically requested in provided arguments: %s", speciesFile)
                            continue
                        if speciesFile not in cacheSpeciesFileList:
                            logger.error("Species archive data file %s on FTP server not found in local cache.", speciesFile)
                            continue
                        cacheSpeciesFile = oD[speciesFile]["archive_file_path"]
                        cacheSpeciesDir = oD[speciesFile]["data_directory"]
                        cacheSpeciesFileExists = os.path.exists(cacheSpeciesFile)
                        cacheSpeciesFileSize = os.path.getsize(cacheSpeciesFile)
                        cacheSpeciesNumModels = len([f for f in os.listdir(cacheSpeciesDir) if f.endswith(".cif.gz")])
                        if not cacheSpeciesFileExists:
                            logger.warning("Species archive data file %s not found at cached path %s", speciesFile, cacheSpeciesFile)
                        if cacheSpeciesFileSize != speciesFileSize:
                            logger.warning("Species archive data file %s not up-to-date with file available on FTP server.", speciesFile)
                        if cacheSpeciesNumModels != speciesNumModels:
                            logger.warning(
                                "Missing some or all of the bundled model files for species archive %s as available on FTP server (found %d / %d model files)",
                                speciesFile,
                                cacheSpeciesNumModels,
                                speciesNumModels,
                            )
                    except Exception as e:
                        logger.exception("Failing on checking of cache data for %s from FTP server, with message:\n%s", speciesData["archive_name"], str(e))

            else:
                logger.info("Refetching all up-to-date files from FTP server.")
                cacheD = {}
                cacheD.update({"created": startDateTime, "data": {}})
                for speciesData in lDL:
                    try:
                        sD = copy.deepcopy(speciesData)
                        speciesFile = speciesData["archive_name"]

                        # If a specific list of species files was requested, only check those (i.e., skip over any that weren't requested)
                        if alphaFoldRequestedSpeciesFileList and speciesFile not in alphaFoldRequestedSpeciesFileList:
                            logger.info("Skipping species file not specifically requested in provided arguments: %s", speciesFile)
                            continue

                        speciesFilePath = os.path.join(alphaFoldFtpDataPath, speciesFile)
                        speciesDataDumpDir = os.path.join(self.__dirPath, speciesFile.split(".")[0])
                        fU.mkdir(speciesDataDumpDir)
                        speciesFileDumpPath = os.path.join(speciesDataDumpDir, speciesFile)

                        logger.info("speciesFileDumpPath %r", speciesFileDumpPath)
                        sD.update({"data_directory": speciesDataDumpDir, "archive_file_path": speciesFileDumpPath})

                        logger.info("Fetching file %s from FTP server to local path %s", speciesFilePath, speciesFileDumpPath)
                        ok = ftpU.get(speciesFilePath, speciesFileDumpPath)
                        ok = fU.unbundleTarfile(speciesFileDumpPath, dirPath=speciesDataDumpDir)
                        logger.info("Completed fetch (%r) at %s (%.4f seconds)", ok, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime)

                        logger.info("Clearing PDB files from extracted tar bundle...")
                        for pdbDumpFile in Path(speciesDataDumpDir).glob("*.pdb.gz"):
                            pdbDumpFile.unlink()

                        if ok:
                            cacheD["data"].update({speciesFile: sD})
                            fU.remove(speciesFileDumpPath)

                    except Exception as e:
                        logger.exception("Failing on fetching and expansion of file %s from FTP server, with message:\n%s", speciesData["archive_name"], str(e))

                createdDate = cacheD["created"]
                oD = cacheD["data"]
                ok = self.__mU.doExport(self.__speciesDataCacheFile, cacheD, fmt="json", indent=3)
                logger.info("Export AlphaFold species model data (%d) status %r", len(oD), ok)

        except Exception as e:
            logger.exception("Failing with %s", str(e))

        return oD, createdDate

    def getSpeciesDirList(self):
        speciesDirList = [self.__oD[k]["data_directory"] for k in self.__oD]

        return speciesDirList

    def getModelFileList(self, inputPathList=None):
        """Method to get the data path list of directories/collections of model mmCIF files.

        Args:
            inputPathList (list): Optional (and recommended) list of paths to directories containing model mmCIF files to return.
                                  If not provided, method will retrieve all cached species data paths, so all model files from all species will be returned in a single list;
                                  thus, it's recommended to provide a list of specific species directory to break the returned model list down into more manageable parts.

        Returns:
            (list): list of model mmCIF file paths (only matches ".cif.gz" files)
        """

        if not inputPathList:
            inputPathList = self.getSpeciesDirList()

        modelFileList = []

        for modelDir in inputPathList:
            try:
                mmCifModels = glob.glob(os.path.join(modelDir, "*.cif.gz"))  # + glob.glob(os.path.join(modelDir, "*.cif")) # ONLY MATCH .GZ FILES TO ENSURE ONLY ONE FILE PER MODEL
                modelFileList += mmCifModels
            except Exception as e:
                logger.exception("Failing with %s", str(e))

        # Add: Return as list of ABSOLUTE paths, not just relative to glob command
        return modelFileList

    def getSpeciesDataDownloadDate(self):
        return self.__createdDate

    def getBaseDataPath(self):
        return self.__dirPath

    def getSpeciesDataCacheFilePath(self):
        return self.__speciesDataCacheFile

    # def hasFeature(self, modelId):
    #     return modelId in self.__oD

    # def getFeature(self, modelId, featureKey):
    #     try:
    #         return self.__oD[modelId][featureKey]
    #     except Exception:
    #         return None

    # def __parseAlphaFoldModelData(self, filePath):
    #     """Parse AlphaFold model mmCIF file

    #     Args:
    #         filePath (str): mmCIF model data file

    #     Returns:
    #         (dict, string): AlphaFold model data dictionary, model version string
    #     """
    #     try:
    #         oD = {}
    #         version = None
    #     except Exception as e:
    #         logger.exception("Failing with %s", str(e))
    #     return True
