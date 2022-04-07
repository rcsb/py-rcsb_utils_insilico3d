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
import os.path
# import copy

from rcsb.utils.insilico3d.AlphaFoldModelProvider import AlphaFoldModelProvider
from rcsb.utils.insilico3d.ModelArchiveModelProvider import ModelArchiveModelProvider
from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.config.ConfigUtil import ConfigUtil

logger = logging.getLogger(__name__)


class ModelProviderWorkflow:
    """Provider workflow class for retrieving and storing/reorganizing computed-model files from various sources."""

    def __init__(self, cachePath, useCache=False, cfgOb=None, configName=None, **kwargs):
        self.__cachePath = cachePath
        self.__cachePath = os.path.abspath(self.__cachePath)
        self.__dirPath = os.path.join(self.__cachePath, "computed-models")  # Directory where CACHE file will be stored (i.e., json or pickle file of all reorganized models)
        self.__cacheFormat = kwargs.get("cacheFormat", "json")
        # self.__cacheFormat = kwargs.get("cacheFormat", "pickle")
        cacheExt = "pic" if self.__cacheFormat == "pickle" else "json"
        self.__cacheFile = kwargs.get("cacheFile", "new-computed-models-cache." + cacheExt)
        self.__cacheFilePath = kwargs.get("cacheFilePath", os.path.join(self.__cachePath, self.__cacheFile))

        configPath = kwargs.get("configPath", "exdb-config-example.yml")  # Config file to use in case no cfgOb provided
        self.__configName = configName if configName else "site_info_configuration"
        mockTopPath = kwargs.get("mockTopPath", None)
        self.__cfgOb = cfgOb if cfgOb else ConfigUtil(configPath=configPath, defaultSectionName=self.__configName, mockTopPath=mockTopPath)

        self.__computedModelsDataPath = self.__cfgOb.getPath("PDBX_COMP_MODEL_SANDBOX_PATH", sectionName=self.__configName)  # Directory where all computed models will be DOWNLOADED (temporarily) and STORED (indefinitely) (i.e., 'computed-models' volume path)
        self.__debugFlag = kwargs.get("debugFlag", False)
        if self.__debugFlag:
            logger.setLevel(logging.DEBUG)

        self.__mU = MarshalUtil(workPath=self.__computedModelsDataPath)
        self.__fU = FileUtil(workPath=self.__computedModelsDataPath)

        self.__cD = self.__reload(useCache=useCache, **kwargs)

    def testCache(self, minCount=0):
        if self.__cD and len(self.__cD) > minCount:
            return True
        else:
            return False

    def getCacheDirPath(self):
        return self.__dirPath

    def reload(self, useCache=True):
        self.__cD = self.__reload(useCache=useCache)
        return len(self.__cD) > 0  # Returns True or False, if reload was successful or not

    def __reload(self, useCache=True, **kwargs):
        """Reload from the current cache file."""
        try:
            cD = {}
            logger.info("useCache %r cacheFilePath %r", useCache, self.__cacheFilePath)
            if useCache and self.__mU.exists(self.__cacheFilePath):
                cD = self.__mU.doImport(self.__cacheFilePath, fmt=self.__cacheFormat)
                logger.info("Reorganized models (%d) in cacheFilePath %r", len(cD), self.__cacheFilePath)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        #
        return cD

    def run(self, useCache=False, **kwargs):
        """Run the model provider workflow:
            1. Retrieve model files from all sources
            2. Reorganize and rename model files into computed-models data directory
            3. Save cacheFile which contains all available models with associated IDs, file paths, and accession URLs

        Args:
            useCache (bool, optional): Use cache file to skip over re-downloading already retrieved models. Defaults to False.
                                       NOTE: Feature not functional yet.
        """
        ok = False
        try:
            redownloadBulkData = kwargs.get("redownloadBulkData", True)
            keepSource = kwargs.get("keepSource", False)
            numProc = int(kwargs.get("numProc", 2))
            chunkSize = int(kwargs.get("chunkSize", 20))
            # fileLimit = int(kwargs.get("fileLimit")) if "fileLimit" in kwargs else None
            # documentLimit = int(kwargs.get("documentLimit")) if "documentLimit" in kwargs else None
            # failedFilePath = kwargs.get("failFileListPath", None)
            # loadFileListPath = kwargs.get("loadFileListPath", None)
            # saveInputFileListPath = kwargs.get("saveFileListPath", None)
            
            ok = self.reload(useCache=useCache)
            if ok:
                # NOTE: Feature not functional yet.
                logger.info("Will skip retrieving models for already downloaded models (%d) in cacheFilePath %r", len(self.__cD), self.__cacheFilePath)
                # ... do something different ...
            else:
                if redownloadBulkData:
                    aFMP = AlphaFoldModelProvider(
                        cachePath=self.__cachePath,
                        useCache=False,
                        cfgOb=self.__cfgOb,
                        configName=self.__configName,
                        numProc=numProc,
                        chunkSize=chunkSize,
                        alphaFoldRequestedSpeciesList=["Helicobacter pylori"]
                        # alphaFoldRequestedSpeciesList=["Helicobacter pylori", "Escherichia coli", "Methanocaldococcus jannaschii", "Homo sapiens"]
                    )
                    if not aFMP.testCache():
                        logger.info("Failed to provide models from AlphaFold")
                        return False

                    mAMP = ModelArchiveModelProvider(
                        cachePath=self.__cachePath,
                        useCache=False,
                        cfgOb=self.__cfgOb,
                        configName=self.__configName,
                        numProc=numProc,
                        chunkSize=chunkSize,
                        # serverDataSetPathD={"ma-bak-cepc": {"urlEnd": "ma-bak-cepc?type=materials_procedures__accompanying_data_file_name", "fileName": "ma-bak-cepc.zip"}}
                    )
                    if not mAMP.testCache():
                        logger.info("Failed to provide models from ModelArchive")
                        return False
                else:
                    # Reload the cache
                    aFMP = AlphaFoldModelProvider(
                        cachePath=self.__cachePath,
                        useCache=True,
                        cfgOb=self.__cfgOb,
                        configName=self.__configName,
                        numProc=numProc,
                        chunkSize=chunkSize,
                        alphaFoldRequestedSpeciesList=["Helicobacter pylori"]
                        # alphaFoldRequestedSpeciesList=["Helicobacter pylori", "Escherichia coli", "Methanocaldococcus jannaschii", "Homo sapiens"]
                    )
                    speciesDirList = aFMP.getArchiveDirList()
                    ok = True if len(speciesDirList) > 0 else False
                    if not ok:
                        logger.info("Failed to reload AlphaFold cache")
                        return False

                    mAMP = ModelArchiveModelProvider(
                        cachePath=self.__cachePath,
                        useCache=True,
                        cfgOb=self.__cfgOb,
                        configName=self.__configName,
                        numProc=numProc,
                        chunkSize=chunkSize,
                        serverDataSetPathD={"ma-bak-cepc": {"urlEnd": "ma-bak-cepc?type=materials_procedures__accompanying_data_file_name", "fileName": "ma-bak-cepc.zip"}}
                    )
                    archiveDirList = mAMP.getArchiveDirList()
                    ok = True if len(archiveDirList) > 0 else False
                    if not ok:
                        logger.info("Failed to reload ModelArchive cache")
                        return False
                
                # Reorganize the models
                ok = aFMP.reorganizeModelFiles(useCache=True, numProc=numProc, chunkSize=chunkSize, keepSource=keepSource)
                if not ok:
                    logger.info("Failed to reorganize AlphaFold models")
                    return False

                ok = mAMP.reorganizeModelFiles(useCache=True, numProc=numProc, chunkSize=chunkSize, keepSource=keepSource)
                if not ok:
                    logger.info("Failed to reorganize ModelArchive models")
                    return False
            
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        #
        return ok

