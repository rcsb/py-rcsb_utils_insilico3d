##
# File:    testModelCacheProvider.py
# Author:  Dennis Piehl
# Date:    28-Apr-2022
# Version: 0.001
#
# Update:
#
#
##
"""
Tests for utilities retrieving computed-model cache data.
"""

__docformat__ = "google en"
__author__ = "Dennis Piehl"
__email__ = "dennis.piehl@rcsb.org"
__license__ = "Apache 2.0"

import logging
import os
import platform
import resource
import time
import unittest

from rcsb.utils.config.ConfigUtil import ConfigUtil
from rcsb.utils.insilico3d.ModelCacheProvider import ModelCacheProvider

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()


class ModelCacheProviderTests(unittest.TestCase):
    def setUp(self):
        self.__cachePath = os.path.join(HERE, "test-output", "CACHE")
        self.__dataPath = os.path.join(HERE, "test-data")
        self.__mockTopPath = os.path.join(TOPDIR, "rcsb", "mock-data")
        self.__configPath = os.path.join(self.__mockTopPath, "config", "dbload-setup-example.yml")
        configName = "site_info_configuration"
        self.__configName = configName
        self.__cfgOb = ConfigUtil(configPath=self.__configPath, defaultSectionName=configName, mockTopPath=self.__mockTopPath)
        #
        self.__holdingsFilePath = self.__cfgOb.getPath("PDBX_COMP_MODEL_CACHE_LIST_PATH", sectionName=self.__configName, default=None)
        if self.__holdingsFilePath is None:
            self.__holdingsFilePath = os.path.join(self.__dataPath, "computed-models-holdings.json.gz")
        #
        self.__startTime = time.time()
        logger.info("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        unitS = "MB" if platform.system() == "Darwin" else "GB"
        rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum resident memory size %.4f %s", rusageMax / 10 ** 6, unitS)
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testGetModelCache(self):
        mcP = ModelCacheProvider(cachePath=self.__cachePath, useCache=False, holdingsRemotePath=self.__holdingsFilePath)
        ok = mcP.testCache()
        self.assertTrue(ok)
        mcP = ModelCacheProvider(cachePath=self.__cachePath, useCache=True, holdingsRemotePath=self.__holdingsFilePath)
        ok = mcP.testCache()
        self.assertTrue(ok)
        compModelId = mcP.getInternalCompModelId("ma-bak-cepc-1100")
        ok = compModelId is not None
        self.assertTrue(ok)
        souceUrl = mcP.getCompModelSourceUrl(compModelId)
        ok = souceUrl is not None
        self.assertTrue(ok)

    #
    # Need to set this up. Also, note that this will only be able to be run on a west-coast luigi instance
    # in order to access the file (since only available over local network)
    #
    # @unittest.skip("Bootstrap test")
    # def testModelCacheBootstrap(self):
    #     try:
    #         mcP = ModelCacheProvider(cachePath=self.__cachePath, useCache=True)
    #         ok = mcP.testCache()
    #         self.assertTrue(ok)
    #         configPath = os.path.join(self.__dataPath, "sabdab-config.yml")
    #         configName = "site_info_remote_configuration"
    #         cfgOb = ConfigUtil(configPath=configPath, defaultSectionName=configName)
    #         ok = mcP.backup(cfgOb, configName, useGit=True, useStash=True)
    #         self.assertTrue(ok)
    #     except Exception as e:
    #         logger.exception("Failing with %s", str(e))
    #         self.fail()


def getModelCacheSuite():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(ModelCacheProviderTests("testGetModelCache"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = getModelCacheSuite()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
