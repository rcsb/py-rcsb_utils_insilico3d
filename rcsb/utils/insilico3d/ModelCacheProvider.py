##
# File:    ModelCacheProvider.py
# Author:  Dennis Piehl
# Date:    27-Apr-2022
#
# Updates:
#
# To Do:
#   - Need to make this able to retrieve from git stash or build locker (to get most recent set of models),
#     so that both coasts use the exact same set of data and cache (including populated timestamps)
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


class ModelCacheProvider(StashableBase):
    """Provider workflow class for retrieving and storing/reorganizing computed-model files from various sources."""

    def __init__(self, **kwargs):
        #
        self.__cachePath = kwargs.get("cachePath", ".")
        useCache = kwargs.get("useCache", True)
        self.__dirName = "computed-models"
        super(ModelCacheProvider, self).__init__(self.__cachePath, [self.__dirName])
        self.__dirPath = os.path.join(self.__cachePath, self.__dirName)
        #
        self.__holdingsRemotePath = kwargs.get("holdingsRemotePath", "http://computed-models-internal-west.rcsb.org/holdings/computed-models-holdings.json.gz")
        self.__fallbackUrl = "https://raw.githubusercontent.com/rcsb/py-rcsb_exdb_assets_stash/testing/stash/computed-models/computed-models-holdings.json.gz"
        self.__holdingsLocalPath = os.path.join(self.__dirPath, "computed-models-holdings.json.gz")
        #
        self.__mU = MarshalUtil(workPath=self.__dirPath)
        self.__mD, self.__idMapD = self.__reload(self.__dirPath, useCache)
        #

    def reload(self):
        self.__mD = self.__reload(self.__dirPath, True)
        return True

    def __reload(self, dirPath, useCache):
        """Reload cached list of downloaded computed-model files and associated metadata.

        Returns:
            (dict): dictionary of cached/downloaded computed-model data
        """
        startTime = time.time()
        ok = False
        #
        logger.info("useCache %r holdingsPath %r", useCache, self.__holdingsLocalPath)
        mD, idMapD = {}, {}
        try:
            if useCache and self.__mU.exists(self.__holdingsLocalPath):
                mD = self.__mU.doImport(self.__holdingsLocalPath, fmt="json")
            else:
                logger.info("Refetching computed-models holdings cache file from %r to %r", self.__holdingsRemotePath, self.__holdingsLocalPath)
                fU = FileUtil()
                fU.mkdir(dirPath)
                ok = fU.get(self.__holdingsRemotePath, self.__holdingsLocalPath)
                logger.info("Computed-model cache fetch status is %r", ok)
                if not ok:
                    logger.info("Refetching computed-models holdings cache fallback file from %r to %r", self.__holdingsRemotePath, self.__holdingsLocalPath)
                    ok = fU.get(self.__fallbackUrl, self.__holdingsLocalPath)
                    logger.info("Computed-model cache fallback fetch status is %r", ok)
                mD = self.__mU.doImport(self.__holdingsLocalPath, fmt="json")
            idMapD = self.getCompModelIdMap(modelCacheD=mD)
            if mD and idMapD:
                ok = True
        except Exception as e:
            logger.exception("Failing with %s", str(e))

        logger.info("Completed reload (%r) at %s (%.4f seconds)", ok, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime)

        return mD, idMapD

    def testCache(self, minCount=10):
        logger.info("computed-model cache count %d", len(self.__mD))
        if self.__mD and len(self.__mD) > minCount:
            return True
        else:
            return False

    def getModelCacheDict(self):
        return self.__mD

    def getCompModelIdMap(self, modelCacheD=None):
        return self.__fetchCompModelIdMap(modelCacheD=modelCacheD)

    def __fetchCompModelIdMap(self, modelCacheD=None):
        """Get the ID mapping between the source model IDs and internal model identifiers for computational models.
        """
        #
        modelCacheD = modelCacheD if modelCacheD else self.__mD
        compModelIdMapD = {}
        try:
            for internalModelId, modelD in modelCacheD.items():
                compModelIdMapD.update({modelD["sourceId"]: internalModelId})
            logger.info("Computed-models mapped ID length: %d", len(compModelIdMapD))
            #
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return compModelIdMapD

    def getInternalCompModelId(self, sourceId):
        """Get the mapped internal ID of computed model provided the original source ID.

        Args:
            sourceId (str): entry.id of computed model (e.g., "AF-P96541-F1")
        """
        compModelInternalId = None
        compModelIdMapD = self.__idMapD if self.__idMapD else self.getCompModelIdMap()
        #
        try:
            compModelInternalId = compModelIdMapD.get(sourceId, None)
            logger.debug("Computed-model sourceId (%s) mapped to internalId (%r)", sourceId, compModelInternalId)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        #
        return compModelInternalId

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
