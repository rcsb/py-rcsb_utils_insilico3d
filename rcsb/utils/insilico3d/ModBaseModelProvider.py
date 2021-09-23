##
# File:    ModBaseModelProvider.py
# Author:  Dennis Piehl
# Date:    15-Sep-2021
#
# Update:
#
#
##
"""
Accessors for ModBase 3D Models (PDB).

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

from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil

logger = logging.getLogger(__name__)


class ModBaseModelProvider:
    """Accessors for ModBase models (PDB)."""

    def __init__(self, **kwargs):
        # Use the same root cachePath for all types of insilico3D model sources, but with unique dirPath names (sub-directory)
        self.__cachePath = kwargs.get("cachePath", "./CACHE-insilico3d-models")
        self.__dirPath = os.path.join(self.__cachePath, "ModBase")
        self.__speciesDataCacheFile = os.path.join(self.__dirPath, "species-model-data.json")

        self.__mU = MarshalUtil(workPath=self.__dirPath)

        self.__oD, self.__createdDate = self.__reload(**kwargs)

    def testCache(self, minCount=0):  # Increase minCount once we are consistently downloading more than one species data set
        if self.__oD and len(self.__oD) > minCount:
            return True
        else:
            return False

    def getCacheDirPath(self):
        return self.__dirPath

    def __reload(self, **kwargs):
        """Reload cached list of species-specific ModBase model data files and check server for updated data sets,
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
            modBaseBaseUrl = kwargs.get("modBaseBaseUrl", "https://salilab.org/modbase-download/projects/genomes")
            # Example URLs:
            #  https://salilab.org/modbase-download/projects/genomes/H_sapiens/2020/Homo_sapiens_2020.tar
            #  https://salilab.org/modbase-download/projects/genomes/H_sapiens/2020/Homo_sapiens_2020.summary.txt
            #  https://salilab.org/modbase-download/projects/genomes/A_thaliana/2021/a_thaliana_2021.tar
            #  https://salilab.org/modbase-download/projects/genomes/A_thaliana/2021/a_thaliana_2021.summary.txt
            #  https://salilab.org/modbase-download/projects/genomes/S_aureus/2008/staph_aureus.tar
            modBaseSpeciesDataPathDict = kwargs.get("modBaseSpeciesDataPathDict", {
                "Homo sapiens": "H_sapiens/2020/Homo_sapiens_2020.tar",
                # "Homo sapiens": {
                #     "models": "H_sapiens/2020/Homo_sapiens_2020.tar",
                #     "summary": "H_sapiens/2020/Homo_sapiens_2020.summary.txt"},
                # "Mus musculus": "10090",
                # "Caenorhabditis elegans": "6239",
                # "Escherichia coli": "83333",
                "Arabidopsis thaliana": "A_thaliana/2021/a_thaliana_2021.tar",
                # "Drosophila melanogaster": "7227",
                # "Saccharomyces cerevisiae": "559292",
                # "Schizosaccharomyces pombe": "284812",
                # "Caulobacter vibrioides": "190650",
                # "Mycobacterium tuberculosis": "83332",
                # "Pseudomonas aeruginosa": "208964",
                "Staphylococcus aureus": "S_aureus/2008/staph_aureus.tar",
                # "Plasmodium falciparum": "36329",
            })

            fU = FileUtil()
            fU.mkdir(self.__dirPath)

            logger.info("useCache %r self.__speciesDataCacheFile %r", useCache, self.__speciesDataCacheFile)
            if useCache and self.__mU.exists(self.__speciesDataCacheFile):
                logger.info("Loading data cache, %s.", self.__speciesDataCacheFile)
                cacheD = self.__mU.doImport(self.__speciesDataCacheFile, fmt="json")
                createdDate = cacheD["created"]
                oD = cacheD["data"]

                logger.info("Checking consistency of cached data with data available on server")
                for species, path in modBaseSpeciesDataPathDict.items():
                    try:
                        speciesFilePath = os.path.join(modBaseBaseUrl, path)
                        speciesFileSize = int(fU.size(speciesFilePath))
                        cacheSpeciesDir = oD[species]["dataDirectory"]
                        cacheSpeciesFileSize = oD[species]["archiveFileSizeBytes"]
                        if not os.path.exists(cacheSpeciesDir):
                            logger.warning("Missing archive data for species %s from server: %s", species, path)
                        if cacheSpeciesFileSize != speciesFileSize:
                            logger.warning("Cached archive data for species %s not up-to-date with data archive on server: %s", species, path)
                    except Exception as e:
                        logger.info("Failing on checking of cache data for species %s: %s", species, path)
                        logger.exception("Failing with %s", str(e))

            else:
                logger.info("Refetching all files from server.")
                cacheD = {}
                cacheD.update({"created": startDateTime, "data": {}})
                for species, path in modBaseSpeciesDataPathDict.items():
                    try:
                        speciesFileName = path.rsplit("/", maxsplit=1)[-1]
                        speciesFilePath = os.path.join(modBaseBaseUrl, path)
                        speciesFileSize = int(fU.size(speciesFilePath))
                        speciesDataDumpDir = os.path.join(self.__dirPath, species.replace(" ", "_"))
                        fU.mkdir(speciesDataDumpDir)
                        speciesFileDumpPath = os.path.join(speciesDataDumpDir, speciesFileName)
                        sD = {
                            "speciesName": species,
                            "archiveFileName": speciesFileName,
                            "archiveFileSizeBytes": speciesFileSize,
                            "dataDirectory": speciesDataDumpDir,
                        }

                        logger.info("Fetching file %s from server to local path %s", speciesFilePath, speciesFileDumpPath)
                        ok = fU.get(speciesFilePath, speciesFileDumpPath)
                        ok = fU.unbundleTarfile(speciesFileDumpPath, dirPath=speciesDataDumpDir)
                        logger.info("Completed fetch (%r) at %s (%.4f seconds)", ok, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime)

                        if ok:
                            cacheD["data"].update({species: sD})
                            fU.remove(speciesFileDumpPath)

                    except Exception as e:
                        logger.info("Failing on fetching and expansion of file for species %s: %s", species, path)
                        logger.exception("Failing with %s", str(e))

                createdDate = cacheD["created"]
                oD = cacheD["data"]
                ok = self.__mU.doExport(self.__speciesDataCacheFile, cacheD, fmt="json", indent=3)
                logger.info("Export ModBase species model data (%d) status %r", len(oD), ok)

        except Exception as e:
            logger.exception("Failing with %s", str(e))

        return oD, createdDate

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
        """Return a list of filepaths for all ModBase PDB model files for a given species.

        Args:
            speciesDataDir (str): path to base-level of species data directory.

        Returns:
            (list): list of absolute model file paths (only matches "*.pdb.xz" files)
        """

        modelFileList = []

        if speciesDataDir:
            try:
                modelFiles = glob.glob(os.path.join(speciesDataDir, "model", "*.pdb.xz"))
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

    def removePdbModelFiles(self):
        # logger.info("Clearing PDB files from extracted tar bundle...")
        # for pdbDumpFile in Path(speciesDataDumpDir).glob("*.pdb.gz"):
        #     pdbDumpFile.unlink()
        return
