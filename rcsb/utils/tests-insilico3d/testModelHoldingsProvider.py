##
# File:    testModelHoldingsProvider.py
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
from rcsb.utils.insilico3d.ModelHoldingsProvider import ModelHoldingsProvider

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()


class ModelHoldingsProviderTests(unittest.TestCase):
    def setUp(self):
        self.__cachePath = os.path.join(HERE, "test-output", "CACHE")
        self.__dataPath = os.path.join(HERE, "test-data")
        self.__mockTopPath = os.path.join(TOPDIR, "rcsb", "mock-data")
        self.__configPath = os.path.join(self.__mockTopPath, "config", "dbload-setup-example.yml")
        configName = "site_info_configuration"
        self.__configName = configName
        self.__cfgOb = ConfigUtil(configPath=self.__configPath, defaultSectionName=configName, mockTopPath=self.__mockTopPath)
        #
        self.__csmRemoteDirPath = self.__cfgOb.getPath("PDBX_COMP_MODEL_REPO_PATH", sectionName=self.__configName, default=None)
        self.__holdingsListRemotePath = self.__cfgOb.getPath("PDBX_COMP_MODEL_HOLDINGS_LIST_PATH", sectionName=self.__configName, default=None)
        if self.__holdingsListRemotePath is None:
            self.__holdingsListRemotePath = os.path.join(self.__dataPath, "computed-models-holdings-list.json")
        #
        self.__startTime = time.time()
        logger.info("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        unitS = "MB" if platform.system() == "Darwin" else "GB"
        rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum resident memory size %.4f %s", rusageMax / 10 ** 6, unitS)
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testGetModelHoldings(self):
        mcP = ModelHoldingsProvider(cachePath=self.__cachePath, useCache=False, csmRemoteDirPath=self.__csmRemoteDirPath, holdingsListRemotePath=self.__holdingsListRemotePath)
        ok = mcP.testCache()
        self.assertTrue(ok)
        mcP = ModelHoldingsProvider(cachePath=self.__cachePath, useCache=True, csmRemoteDirPath=self.__csmRemoteDirPath, holdingsListRemotePath=self.__holdingsListRemotePath)
        ok = mcP.testCache()
        self.assertTrue(ok)
        mD = mcP.getModelHoldingsDict()
        ok = len(mD) > 5
        self.assertTrue(ok)


def getModelCacheSuite():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(ModelHoldingsProviderTests("testGetModelHoldings"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = getModelCacheSuite()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
