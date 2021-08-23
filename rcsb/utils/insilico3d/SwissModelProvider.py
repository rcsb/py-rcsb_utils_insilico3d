##
# File:    SwissModelProvider.py
# Author:  Dennis Piehl
# Date:    23-Aug-2021
#
# Update:
#
#
##
"""
Accessors for SWISS-MODEL 3D Models (PDB).

"""

import datetime
import logging
import os.path
import time
from pathlib import Path
# import copy
import glob

from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.io.FtpUtil import FtpUtil

logger = logging.getLogger(__name__)


class SwissModelProvider:
    """Accessors for SWISS-MODEL models (PDB)."""

    def __init__(self, **kwargs):
        # Use the same root cachePath for all types of insilico3D model sources, but with unique dirPath names (sub-directory)
        self.__cachePath = kwargs.get("cachePath", "./CACHE-insilico3d-models")
        self.__dirPath = os.path.join(self.__cachePath, "SWISS-MODEL")
        self.__speciesDataCacheFile = os.path.join(self.__dirPath, "species-model-data.json")

        self.__mU = MarshalUtil(workPath=self.__dirPath)

        self.__oD, self.__createdDate = self.__reload(**kwargs)

    def testCache(self, minCount=0):  # Increase minCount once we are consistently downloading more than one species data set
        if self.__oD and len(self.__oD) > minCount:
            return True
        else:
            return False

    def __reload(self, **kwargs):
        """Reload cached list of species-specific SWISS-MODEL model data files and check FTP server for updated data sets,
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
            swissModelBaseUrl = kwargs.get("swissModelBaseUrl", "https://swissmodel.expasy.org/repository/download/core_species")
            # "Homo Sapiens" coordinates : "9606_coords.tar.gz"
            #  -> Full URL: https://swissmodel.expasy.org/repository/download/core_species/9606_coords.tar.gz
            # "Homo Sapiens" metadata : "9606_meta.tar.gz"
            swissModelSpeciesDataIdDict = {
                "Homo sapiens": "9606",
                "Mus musculus": "10090",
                "Caenorhabditis elegans": "6239",
                "Escherichia coli": "83333",
                "Arabidopsis thaliana": "3702",
                "Drosophila melanogaster": "7227",
                "Saccharomyces cerevisiae": "559292",
                "Schizosaccharomyces pombe": "284812",
                "Caulobacter vibrioides": "190650",
                "Mycobacterium tuberculosis": "83332",
                "Pseudomonas aeruginosa": "208964",
                "Staphylococcus aureus": "93061",
                "Plasmodium falciparum": "36329"
            }
            ok = False
            fU = FileUtil()
            fU.mkdir(self.__dirPath)

            logger.info("useCache %r self.__speciesDataCacheFile %r", useCache, self.__speciesDataCacheFile)
            if useCache and self.__mU.exists(self.__speciesDataCacheFile):
                logger.info("Loading data cache, %s.", self.__speciesDataCacheFile)
                cacheD = self.__mU.doImport(self.__speciesDataCacheFile, fmt="json")
                createdDate = cacheD["created"]
                oD = cacheD["data"]

                logger.info("Checking consistency of cached data with data available on FTP")
                for species in swissModelSpeciesDataIdDict:
                    try:
                        speciesCoordFile = swissModelSpeciesDataIdDict[species]+"_coords.tar.gz"
                        speciesCoordFilePath = os.path.join(swissModelBaseUrl, swissModelSpeciesDataIdDict[species]+"_coords.tar.gz")
                        speciesFileSize = int(fU.size(speciesCoordFile))
                        cacheSpeciesFile = oD[speciesCoordFile]["archive_file_path"]
                        cacheSpeciesFileExists = os.path.exists(cacheSpeciesFile)
                        cacheSpeciesFileSize = os.path.getsize(cacheSpeciesFile)
                        if not cacheSpeciesFileExists:
                            logger.warning("Missing species archive data file for %s which is available for download online.", speciesCoordFile)
                        if cacheSpeciesFileSize != speciesFileSize:
                            logger.warning("Species archive data file %s not up-to-date with file available online.", speciesCoordFile)
                    except Exception as e:
                        logger.info("Failing on checking of cache data for %s", speciesCoordFile)
                        logger.exception("Failing with %s", str(e))

            else:
                logger.info("Refetching all up-to-date files over HTTP.")
                cacheD = {}
                cacheD.update({"created": startDateTime, "data": {}})
                for species in swissModelSpeciesDataIdDict:
                    try:
                        # sD = copy.deepcopy(speciesData)
                        speciesId = swissModelSpeciesDataIdDict[species]
                        speciesCoordFile = speciesId+"_coords.tar.gz"
                        speciesCoordFilePath = os.path.join(swissModelBaseUrl, speciesCoordFile)
                        speciesDataDumpDir = os.path.join(self.__dirPath, species.replace(" ", "_")+"_"+speciesId)
                        fU.mkdir(speciesDataDumpDir)
                        sD = {"species": species, "id": speciesId, "data_directory": speciesDataDumpDir, "archive_file_path": speciesCoordFilePath}

                        speciesFileDumpPath = os.path.join(speciesDataDumpDir, speciesCoordFile)
                        logger.info("speciesFileDumpPath %r", speciesFileDumpPath)

                        logger.info("Fetching file %s over HTTP to local path %s", speciesCoordFilePath, speciesFileDumpPath)
                        ok = fU.get(speciesCoordFilePath, speciesFileDumpPath)
                        ok = fU.unbundleTarfile(speciesFileDumpPath, dirPath=speciesDataDumpDir)
                        logger.info("Completed fetch (%r) at %s (%.4f seconds)", ok, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime)

                        # logger.info("Clearing PDB files from extracted tar bundle...")
                        # for pdbDumpFile in Path(speciesDataDumpDir).glob("*.pdb.gz"):
                        #     pdbDumpFile.unlink()

                        if ok:
                            cacheD["data"].update({speciesCoordFile: sD})

                    except Exception as e:
                        logger.info("Failing on fetching and expansion of file %s", speciesCoordFile)
                        logger.exception("Failing with %s", str(e))

                createdDate = cacheD["created"]
                oD = cacheD["data"]
                ok = self.__mU.doExport(self.__speciesDataCacheFile, cacheD, fmt="json", indent=3)
                logger.info("Export SWISS-MODEL species model data (%d) status %r", len(oD), ok)

        except Exception as e:
            logger.exception("Failing with %s", str(e))

        return oD, createdDate

    def getSpeciesDirList(self):
        speciesDirList = [k["data_directory"] for k in self.__oD]

        return speciesDirList

    def getModelFileList(self, inputPathList=None):
        """Method to get the data path list of directories/collections of model PDB files.

        Args:
            inputPathList (list): list of paths to directories containing model PDB files to return.

        Returns:
            (list): list of model PDB file paths (only matches ".pdb" files)
        """

        inputPathList = inputPathList if inputPathList else []

        modelFileList = []

        for modelDir in inputPathList:
            pdbModels = glob.glob(os.path.join(modelDir, "*.pdb"))
            modelFileList += pdbModels

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

    # def __parseSwissModelData(self, filePath):
    #     """Parse SWISS-MODEL model PDB file

    #     Args:
    #         filePath (str): PDB model data file

    #     Returns:
    #         (dict, string): SWISS-MODEL model data dictionary, model version string
    #     """
    #     try:
    #         oD = {}
    #         version = None
    #     except Exception as e:
    #         logger.exception("Failing with %s", str(e))
    #     return True
