##
# File:    SwissModelProvider.py
# Author:  Dennis Piehl
# Date:    29-Sep-2021
#
# Update:
#
#
##
"""
Accessors for SWISS-MODEL 3D Models (PDB).

"""

__docformat__ = "google en"
__author__ = "Dennis Piehl"
__email__ = "dennis.piehl@rcsb.org"
__license__ = "Apache 2.0"

import datetime
import logging
import os.path
import time
import glob
import requests

from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil

logger = logging.getLogger(__name__)


class SwissModelProvider:
    """Accessors for SWISS-MODEL models (PDB)."""

    def __init__(self, **kwargs):
        # Use the same root cachePath for all types of insilico3D model sources, but with unique dirPath names (sub-directory)
        self.__cachePath = kwargs.get("cachePath", "./CACHE-insilico3d-models")
        self.__dirPath = os.path.join(self.__cachePath, "SWISS-MODEL")
        self.__speciesDataCacheFile = os.path.join(self.__dirPath, "species-model-data.json")

        self.__mU = MarshalUtil(workPath=self.__dirPath)
        self.__fU = FileUtil(workPath=self.__dirPath)

        self.__oD, self.__createdDate = self.__reload(**kwargs)

    def testCache(self, minCount=0):  # Increase minCount once we are consistently downloading more than one species data set
        if self.__oD and len(self.__oD) > minCount:
            return True
        else:
            return False

    def getCacheDirPath(self):
        return self.__dirPath

    def __reload(self, **kwargs):
        """Reload cached list of species-specific SWISS-MODEL model data files and check server for updated data sets,
        or re-download latest versions of all species-specific model data sets from server.

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
            swissModelBaseUrl = kwargs.get("swissModelBaseUrl", "https://swissmodel.expasy.org/repository/download/core_species")
            # Example URLs:
            #  https://swissmodel.expasy.org/repository/download/core_species/9606_coords.tar.gz
            #  https://swissmodel.expasy.org/repository/download/core_species/9606_meta.tar.gz
            swissModelServerSpeciesDataPathDict = kwargs.get("swissModelServerSpeciesDataPathDict", {
                "Homo sapiens": "9606_coords.tar.gz",
                "Arabidopsis thaliana": "3702_coords.tar.gz",
                "Staphylococcus aureus": "93061_coords.tar.gz",  # Used for tests only
                # "Caenorhabditis elegans": "6239_coords.tar.gz",
                # "Mus musculus": "10090_coords.tar.gz",
                # "Escherichia coli": "83333_coords.tar.gz",
                # "Drosophila melanogaster": "7227_coords.tar.gz",
                # "Saccharomyces cerevisiae": "559292_coords.tar.gz",
                # "Schizosaccharomyces pombe": "284812_coords.tar.gz",
                # "Caulobacter vibrioides": "190650_coords.tar.gz",
                # "Mycobacterium tuberculosis": "83332_coords.tar.gz",
                # "Pseudomonas aeruginosa": "208964_coords.tar.gz",
                # "Plasmodium falciparum": "36329_coords.tar.gz",
            })

            self.__fU.mkdir(self.__dirPath)

            logger.info("useCache %r self.__speciesDataCacheFile %r", useCache, self.__speciesDataCacheFile)
            if useCache and self.__mU.exists(self.__speciesDataCacheFile):
                logger.info("Loading data cache, %s.", self.__speciesDataCacheFile)
                cacheD = self.__mU.doImport(self.__speciesDataCacheFile, fmt="json")
                createdDate = cacheD["created"]
                oD = cacheD["data"]

                logger.info("Checking consistency of cached data with data available on server")
                for species, path in swissModelServerSpeciesDataPathDict.items():
                    try:
                        speciesFilePath = os.path.join(swissModelBaseUrl, path)
                        speciesFileSize = int(self.__fU.size(speciesFilePath))
                        speciesFileHeadResp = requests.head(speciesFilePath)
                        speciesFileModDate = datetime.datetime.strptime(speciesFileHeadResp.headers["Last-Modified"], "%a, %d %b %Y %H:%M:%S %Z").isoformat()
                        cacheSpeciesDir = oD[species]["dataDirectory"]
                        cacheSpeciesFileSize = oD[species]["archiveFileSizeBytes"]
                        cacheSpeciesFileModDate = oD[species]["lastModified"]
                        if not os.path.exists(cacheSpeciesDir):
                            logger.warning("Missing archive data for species %s from server: %s", species, path)
                        if cacheSpeciesFileSize != speciesFileSize or cacheSpeciesFileModDate != speciesFileModDate:
                            logger.warning("Cached archive data for species %s not up-to-date with data archive on server: %s", species, path)
                    except Exception as e:
                        logger.info("Failing on checking of cache data for species %s: %s", species, path)
                        logger.exception("Failing with %s", str(e))

            else:
                logger.info("Refetching all files from server.")
                cacheD = {}
                cacheD.update({"created": startDateTime, "data": {}})
                for species, path in swissModelServerSpeciesDataPathDict.items():
                    try:
                        speciesFileName = path.rsplit("/", maxsplit=1)[-1]
                        speciesFilePath = os.path.join(swissModelBaseUrl, path)
                        speciesFileSize = int(self.__fU.size(speciesFilePath))
                        speciesFileHeadResp = requests.head(speciesFilePath)
                        speciesFileModDate = datetime.datetime.strptime(speciesFileHeadResp.headers["Last-Modified"], "%a, %d %b %Y %H:%M:%S %Z").isoformat()
                        speciesDataDumpDir = os.path.join(self.__dirPath, species.replace(" ", "_"))
                        self.__fU.mkdir(speciesDataDumpDir)
                        speciesFileDumpPath = os.path.join(speciesDataDumpDir, speciesFileName)
                        sD = {
                            "speciesName": species,
                            "archiveFileName": speciesFileName,
                            "archiveFileSizeBytes": speciesFileSize,
                            "lastModified": speciesFileModDate,
                            "dataDirectory": speciesDataDumpDir,
                        }
                        logger.info("Fetching file %s from server to local path %s", speciesFilePath, speciesFileDumpPath)
                        ok = self.__fU.get(speciesFilePath, speciesFileDumpPath)
                        ok = self.__fU.unbundleTarfile(speciesFileDumpPath, dirPath=speciesDataDumpDir)
                        logger.info("Completed fetch (%r) at %s (%.4f seconds)", ok, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime)

                        if ok:
                            cacheD["data"].update({species: sD})
                            self.__fU.remove(speciesFileDumpPath)

                    except Exception as e:
                        logger.info("Failing on fetching and expansion of file for species %s: %s", species, path)
                        logger.exception("Failing with %s", str(e))

                createdDate = cacheD["created"]
                oD = cacheD["data"]
                ok = self.__mU.doExport(self.__speciesDataCacheFile, cacheD, fmt="json", indent=3)
                logger.info("Export SWISS-MODEL species model data (%d) status %r", len(oD), ok)

        except Exception as e:
            logger.exception("Failing with %s", str(e))

        return oD, createdDate

    def getSpeciesNameList(self):
        """Return a list of species names/keys for all cached species data directories."""
        speciesNameList = [k for k in self.__oD]
        return speciesNameList

    def getSpeciesDirList(self):
        """Return a list of paths to all species data directories."""
        speciesDirList = [self.__oD[k]["dataDirectory"] for k in self.__oD]
        return speciesDirList

    def getModelFileList(self, inputPathList=None):
        """Return a list of filepaths for all mmCIF models under the provided set of directories.

        Args:
            inputPathList (list): list of paths to directories containing model files to return.

        Returns:
            (list): list of model file paths (only matches "*.cif.gz" files)
        """

        inputPathList = inputPathList if inputPathList else []

        modelFileList = []

        for modelDir in inputPathList:
            try:
                modelFiles = glob.glob(os.path.join(modelDir, "*.cif.gz"))
                modelFileList = [os.path.abspath(f) for f in modelFiles]
            except Exception as e:
                logger.exception("Failing with %s", str(e))

        return modelFileList

    def getSpeciesPdbModelFileList(self, speciesDataDir=None):
        """Return a list of filepaths for all SWISS-MODEL PDB model files for a given species.

        Args:
            speciesDataDir (str): path to base-level of species data directory.

        Returns:
            (list): list of absolute model file paths (only matches "*.pdb" files)
        """

        modelFileList = []

        if speciesDataDir:
            try:
                modelFiles = glob.glob(os.path.join(speciesDataDir, "SWISS-MODEL_Repository", "*", "*", "*", "swissmodel", "*.pdb"))
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

    def getSpeciesConversionDict(self, speciesName=None):
        """Returns a dictionary containing necessary data for converting PDBs to mmCIF in SwissModelProcessor class."""

        speciesConversionDict = {}
        if speciesName:
            speciesConversionDict["speciesName"] = speciesName
            speciesDataDir = self.__oD[speciesName]["dataDirectory"]
            speciesConversionDict["speciesModelDir"] = speciesDataDir
            speciesConversionDict["lastModified"] = self.__oD[speciesName]["lastModified"]
            speciesConversionDict["speciesPdbModelFileList"] = self.getSpeciesPdbModelFileList(speciesDataDir=speciesDataDir)
        return speciesConversionDict

    def removePdbModelDir(self, speciesDataDir=None):
        """Remove the directory containing PDB model files for a given species directory."""

        ok = False
        if speciesDataDir:
            try:
                speciesPdbModelDir = os.path.join(speciesDataDir, "model")
                ok = self.__fU.remove(speciesPdbModelDir)
            except Exception as e:
                logger.exception("Failing with %s", str(e))
        return ok

    def removeAlignmentDir(self, speciesDataDir=None):
        """Remove the directory containing alignment files for a given species directory."""

        ok = False
        if speciesDataDir:
            try:
                speciesAlignmentDir = os.path.join(speciesDataDir, "alignment")
                ok = self.__fU.remove(speciesAlignmentDir)
            except Exception as e:
                logger.exception("Failing with %s", str(e))
        return ok

    def removeSpeciesDataDir(self, speciesName=None, updateCache=True):
        """"Remove an entire species data directory (and its corresponding cache file entry),
        provided the species name as stored in the cache file."""

        ok = False
        if speciesName:
            try:
                cacheD = self.__mU.doImport(self.__speciesDataCacheFile, fmt="json")
                dataD = cacheD["data"]
                speciesDataD = dataD.pop(speciesName)
                if updateCache:
                    ok = self.__mU.doExport(self.__speciesDataCacheFile, cacheD, fmt="json", indent=3)
                speciesDataDir = speciesDataD["dataDirectory"]
                ok = self.__fU.remove(speciesDataDir)
            except Exception as e:
                logger.exception("Failing with %s", str(e))
        return ok
