##
# File:    ModelHoldingsProvider.py
# Author:  Dennis Piehl
# Date:    27-Apr-2022
#
# Updates:
#
# To Do:
#
# 21-Jul-2022  bv Allow for smaller number of test cases
# 05-Mar-2024 dwp Adjustments to provide support for CSM scaling (using multiple holdings file, not just one);
#                 but, still maintain current production support (which relies specifically on computed-models-holdings.json.gz
#                 and includes fragmented models)
##

"""
Provider class to retrieve the cache (dict) of downloaded model files and their associated metadata.

"""

__docformat__ = "google en"
__author__ = "Dennis Piehl"
__email__ = "dennis.piehl@rcsb.org"
__license__ = "Apache 2.0"

import logging
import os.path
import os
import time

from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.io.StashableBase import StashableBase

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger(__name__)


class ModelHoldingsProvider(StashableBase):
    """Provider workflow class for retrieving and storing/reorganizing computed-model files from various sources."""

    def __init__(self, **kwargs):
        #
        self.__cachePath = kwargs.get("cachePath", ".")
        useCache = kwargs.get("useCache", True)
        self.__dirName = "computed-models"
        super(ModelHoldingsProvider, self).__init__(self.__cachePath, [self.__dirName])
        self.__dirPath = os.path.join(self.__cachePath, self.__dirName)
        #
        self.__csmRemoteDirPath = kwargs.get("csmRemoteDirPath", "http://computed-models-internal-coast.rcsb.org/staging/")
        self.__holdingsListRemotePath = kwargs.get("holdingsListRemotePath", "http://computed-models-internal-coast.rcsb.org/staging/holdings/computed-models-holdings-list.json")
        self.__holdingsListLocalPath = os.path.join(self.__dirPath, os.path.basename(self.__holdingsListRemotePath))
        #
        # Remove these when done switching to 200 million BCIF load, since only applicable for original ~1 million AF load
        self.__importHoldingsFileToMemory = True
        self.__holdingsFileToImport = "computed-models-holdings.json.gz"
        #
        self.__mU = MarshalUtil(workPath=self.__dirPath)
        self.__fU = FileUtil(workPath=self.__dirPath)
        #
        self.__hD, self.__mD, self.__fD = {}, {}, {}
        #
        self.__hD = self.__reload(useCache)
        #

    def reload(self):
        self.__hD = self.__reload(useCache=True)
        return True

    def __reload(self, useCache):
        """Reload cached list of computed-model holdings files.

        Returns:
            (dict): dictionary of holdings files
        """
        startTime = time.time()
        ok = False
        #
        logger.info("useCache %r holdingsListPath %r", useCache, self.__holdingsListLocalPath)
        hD = {}
        try:
            if useCache and self.__mU.exists(self.__holdingsListLocalPath):
                hD = self.__mU.doImport(self.__holdingsListLocalPath, fmt="json")
                for holdingsFile in hD:
                    localHoldingsFilePath = os.path.join(self.__dirPath, os.path.basename(holdingsFile))
                    if not self.__mU.exists(localHoldingsFilePath):
                        ok = self.__fetchHoldingsFile(holdingsFile)
            else:
                logger.info("Refetching computed-models holdings list file from %r to %r", self.__holdingsListRemotePath, self.__holdingsListLocalPath)
                ok = self.__fU.get(self.__holdingsListRemotePath, self.__holdingsListLocalPath)
                logger.info("Computed-model holdings list fetch status is %r", ok)
                hD = self.__mU.doImport(self.__holdingsListLocalPath, fmt="json")
                for holdingsFile in hD:
                    ok = self.__fetchHoldingsFile(holdingsFile)
            if self.__importHoldingsFileToMemory:
                holdingsFileTmpL = [hF for hF in hD if self.__holdingsFileToImport in hF]
                holdingsFile = holdingsFileTmpL[0] if holdingsFileTmpL else None
                if holdingsFile:
                    self.__mD, self.__fD = self.__reloadHoldingsFile(holdingsFile)

        except Exception as e:
            logger.exception("Failing with %s", str(e))

        logger.info("Completed reload (%r) at %s (%.4f seconds)", ok, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime)

        return hD

    def __fetchHoldingsFile(self, holdingsFile):
        """Fetch single holdings file.

        Args:
            holdingsFile (str): relative holdings file path (e.g., "holdings/computed-models-holdings.json.gz" or "holdings/alphafold-holdings-100-0.json.gz")
        """
        ok = False
        try:
            localHoldingsFilePath = os.path.join(self.__dirPath, os.path.basename(holdingsFile))
            remoteHoldingsFilePath = os.path.join(self.__csmRemoteDirPath, holdingsFile)
            ok = self.__fU.get(remoteHoldingsFilePath, localHoldingsFilePath)
            logger.info("Fetched computed-model holdings file %s status %r", remoteHoldingsFilePath, ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return ok

    def __reloadHoldingsFile(self, holdingsFile):
        """Reload holdings file of downloaded computed-model files and associated metadata.
        Also generate a dictionary of all IDs that are part of fragmented models (in which case there are no AF download URLs or landing pages).

        Returns:
            (dict): dictionary of cached/downloaded computed-model holdings data
            (dict): dictionary of fragmented computed-model IDs
        """
        startTime = time.time()
        ok = False
        #
        localHoldingsFilePath = os.path.join(self.__dirPath, os.path.basename(holdingsFile))
        #
        logger.info("Loading localHoldingsFilePath into memory %r", localHoldingsFilePath)
        mD, fD = {}, {}
        try:
            if not self.__mU.exists(localHoldingsFilePath):
                ok = self.__fetchHoldingsFile(holdingsFile)
            mD = self.__mU.doImport(localHoldingsFilePath, fmt="json")
            if mD:
                ok = True
                logger.info("Imported computed-model holdings file into memory (as self.__mD) %r length %r", localHoldingsFilePath, len(mD))
            fD = self.__getFragmentedModelIds(mD)
            if fD:
                logger.info("Fragmented computed-model ID list length %r", len(fD))
        except Exception as e:
            logger.exception("Failing with %s", str(e))

        logger.info("Completed reload (%r) at %s (%.4f seconds)", ok, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime)

        return mD, fD

    def getFragmentedModelIds(self, modelD=None):
        self.__fD = self.__getFragmentedModelIds(modelD=modelD)
        return True

    def __getFragmentedModelIds(self, modelD=None):
        """Gather list of fragmented AF computed-model files.

        Returns:
            (dict): dictionary of fragmented computed-model IDs
        """
        fD = {}
        modelD = modelD if modelD else self.__mD
        fxL = [i for i in modelD if i.startswith("AF") and not i.endswith("F1")]
        f1L = [i[0:-2] + "F1" for i in fxL if i.endswith("F2")]
        fxL += f1L
        fxL.sort()
        fD = {k: None for k in fxL}
        return fD

    def testCache(self, minCount=1):
        logger.info("computed-model holdings files count %d", len(self.__hD))
        if self.__hD and len(self.__hD) >= minCount:
            return True
        else:
            return False

    def getHoldingsFileDict(self):
        return self.__hD

    def getModelHoldingsDict(self):
        return self.__mD

    def getModelFragmentsDict(self):
        return self.__fD

    def checkIfFragmentedModel(self, compModelInternalId):
        return compModelInternalId in self.__fD

    # def getCompModelIdMap(self, modelCacheD=None):
    #     return self.__fetchCompModelIdMap(modelCacheD=modelCacheD)

    # def __fetchCompModelIdMap(self, modelCacheD=None):
    #     """Get the ID mapping between the source model IDs and internal model identifiers for computational models.
    #     """
    #     #
    #     modelCacheD = modelCacheD if modelCacheD else self.__mD
    #     compModelIdMapD = {}
    #     try:
    #         for internalModelId, modelD in modelCacheD.items():
    #             compModelIdMapD.update({modelD["sourceId"]: internalModelId})
    #         logger.info("Computed-models mapped ID length: %d", len(compModelIdMapD))
    #         #
    #     except Exception as e:
    #         logger.exception("Failing with %s", str(e))
    #     return compModelIdMapD

    # def getInternalCompModelId(self, sourceId):
    #     """Get the mapped internal ID of computed model provided the original source ID.

    #     Args:
    #         sourceId (str): entry.id of computed model (e.g., "AF-P96541-F1")
    #     """
    #     compModelInternalId = None
    #     compModelIdMapD = self.__idMapD if self.__idMapD else self.getCompModelIdMap()
    #     #
    #     try:
    #         compModelInternalId = compModelIdMapD.get(sourceId, None)
    #         logger.debug("Computed-model sourceId (%s) mapped to internalId (%r)", sourceId, compModelInternalId)
    #     except Exception as e:
    #         logger.exception("Failing with %s", str(e))
    #     #
    #     return compModelInternalId

    def getCompModelData(self, compModelInternalId):
        """Get the associated metadata dict for a specific computed model.

        Args:
            compModelInternalId (str): internal computed-model identifier (e.g., "AF_AFP96541F1")
        """
        compModelD = {}
        try:
            compModelD = self.__mD.get(compModelInternalId, {})
            if not compModelD:
                logger.error("Unable to retrieve source URL for computed-model (%s)", compModelInternalId)
            logger.debug("Computed-model internalId (%s) modelData (%r)", compModelInternalId, compModelD)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        #
        return compModelD
