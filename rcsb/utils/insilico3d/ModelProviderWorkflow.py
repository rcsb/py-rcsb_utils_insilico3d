##
# File:    ModelProviderWorkflow.py
# Author:  Dennis Piehl
# Date:    6-Apr-2022
#
# Updates:
#
# To Do:
##

"""
Workflow wrapper for retrieving and reorganizing/storing computed models from various sources.

"""

__docformat__ = "google en"
__author__ = "Dennis Piehl"
__email__ = "dennis.piehl@rcsb.org"
__license__ = "Apache 2.0"

import logging

from rcsb.utils.insilico3d.AlphaFoldModelProvider import AlphaFoldModelProvider
from rcsb.utils.insilico3d.AlphaFoldModelCloudProvider import AlphaFoldModelCloudProvider
from rcsb.utils.insilico3d.ModelArchiveModelProvider import ModelArchiveModelProvider

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger(__name__)


class ModelProviderWorkflow:
    """Provider workflow class for retrieving and storing/reorganizing computed-model files from various sources."""

    def __init__(self, srcDir, destDir, modelProviders=None, useCache=True, numProc=4, chunkSize=40, **kwargs):
        """Initialize Model Provider Workflow instance.

        Args:
            srcDir (str): Directory path to which to download model files
            destDir (str): Path to directory where model files will be reorganized and stored permanently (also contains computed-model-cache.json file)
                           (i.e., path to 'computed-models'volume path). Note that this parameter corresponds to the "cachePath" parameter in the provider-specific classes.
            useCache (bool, optional): Start from cache of already downloaded and/or reorganized model files. Defaults to True.
                                       When True (default), checks if last downloaded set of files is up-to-date and downloads any newly available models.
                                       When False, redownloads all model files.
        """
        self.__srcDir = srcDir
        self.__destDir = destDir
        self.__useCache = useCache
        self.__numProc = numProc
        self.__chunkSize = chunkSize
        self.__modelProviders = modelProviders if modelProviders else ["AlphaFold", "ModelArchive"]
        #
        self.__smallFileSizeCutoff = kwargs.get("smallFileSizeCutoff", 8388608)
        self.__alphaFoldRequestedSpeciesList = kwargs.get("alphaFoldRequestedSpeciesList", [])
        self.__alphaFoldRequestedTaxIdPrefixList = kwargs.get("alphaFoldRequestedTaxIdPrefixList", [])
        self.__modelArchiveRequestedDatasetD = kwargs.get("modelArchiveRequestedDatasetD", {})

        if "AlphaFold" in self.__modelProviders:
            self.__aFMP = AlphaFoldModelProvider(
                cachePath=self.__destDir,
                baseWorkPath=self.__srcDir,
                useCache=self.__useCache,
                reload=False,
                alphaFoldRequestedSpeciesList=self.__alphaFoldRequestedSpeciesList,
            )

        if "AlphaFoldCloud" in self.__modelProviders:
            self.__aFMCP = AlphaFoldModelCloudProvider(
                cachePath=self.__destDir,
                baseWorkPath=self.__srcDir,
                useCache=self.__useCache,
                reload=False,
                alphaFoldRequestedTaxIdPrefixList=self.__alphaFoldRequestedTaxIdPrefixList,
            )

        if "ModelArchive" in self.__modelProviders:
            self.__mAMP = ModelArchiveModelProvider(
                cachePath=self.__destDir,
                baseWorkPath=self.__srcDir,
                useCache=self.__useCache,
                reload=False,
                modelArchiveRequestedDatasetD=self.__modelArchiveRequestedDatasetD
            )

    def download(self, modelProviders=None, **kwargs):
        """Download model files from model provider sources.

        Args:
            modelProviders (list, optional): List of model providers to download from. Defaults to ["AlphaFold", "ModelArchive"].

        Returns:
            (bool): True if successful; False otherwise.
        """
        modelProviders = modelProviders if modelProviders else self.__modelProviders
        alphaFoldRequestedSpeciesList = kwargs.get("alphaFoldRequestedSpeciesList", self.__alphaFoldRequestedSpeciesList)
        modelArchiveRequestedDatasetD = kwargs.get("modelArchiveRequestedDatasetD", self.__modelArchiveRequestedDatasetD)

        ok = False
        try:
            for provider in modelProviders:
                if provider == "AlphaFold":
                    self.__aFMP.reload(useCache=self.__useCache, alphaFoldRequestedSpeciesList=alphaFoldRequestedSpeciesList)
                    ok = self.__aFMP.testCache()
                if provider == "ModelArchive":
                    self.__mAMP.reload(useCache=self.__useCache, modelArchiveRequestedDatasetD=modelArchiveRequestedDatasetD)
                    ok = self.__mAMP.testCache()
                if not ok:
                    logger.error("Problem downloading models and/or testing cache for provider, %s", provider)
                    return False
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return ok

    def reorganize(self, modelProviders=None, keepSource=False, **kwargs):
        """Reorganize collection of downloaded model files into destination directory (e.g., in 'computed-models' volume path).

        Args:
            modelProviders (list, optional): List of model providers to update downloaded files from. Defaults to ["AlphaFold", "ModelArchive"].
            keepSource (bool, optional): whether to copy files to new directory (instead of moving them). Defaults to False.

        Returns:
            (bool): True if successful; False otherwise.
        """
        modelProviders = modelProviders if modelProviders else self.__modelProviders
        numProc = kwargs.get("numProc", self.__numProc)
        chunkSize = kwargs.get("chunkSize", self.__chunkSize)
        smallFileSizeCutoff = kwargs.get("smallFileSizeCutoff", self.__smallFileSizeCutoff)
        inputDirList = kwargs.get("inputDirList", [])
        dictFilePathL = kwargs.get("dictFilePathL", None)

        ok = False
        try:
            for provider in modelProviders:
                if provider == "AlphaFold":
                    ok = self.__aFMP.reorganizeModelFiles(useCache=self.__useCache, numProc=numProc, chunkSize=chunkSize, keepSource=keepSource)
                if provider == "AlphaFoldCloud":
                    ok = self.__aFMCP.reorganizeModelFiles(
                        cachePath=self.__destDir,
                        useCache=self.__useCache,
                        inputTaxIdPrefixList=inputDirList,
                        numProc=numProc,
                        chunkSize=chunkSize,
                        smallFileSizeCutoff=smallFileSizeCutoff,
                        keepSource=keepSource,
                        dictFilePathL=dictFilePathL,
                    )
                if provider == "ModelArchive":
                    ok = self.__mAMP.reorganizeModelFiles(useCache=self.__useCache, numProc=numProc, chunkSize=chunkSize, keepSource=keepSource)
                if not ok:
                    logger.info("Failed to reorganize models for provider, %s", provider)
                    return False
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return ok

    # def run(self, **kwargs):
    #     """Run the model provider workflow:
    #         1. Retrieve model files from all sources
    #         2. Reorganize and rename model files into computed-models data directory
    #     """
    #     return
