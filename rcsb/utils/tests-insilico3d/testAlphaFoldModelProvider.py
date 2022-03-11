##
# File:    testAlphaFoldModelProvider.py
# Author:  Dennis Piehl
# Date:    30-Sep-2021
#
# Updates:
#
#
##
"""
Tests for accessor utilities for AlphaFold 3D Model (mmCIF) data.

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

from rcsb.utils.insilico3d.AlphaFoldModelProvider import AlphaFoldModelProvider
# from rcsb.utils.config.ConfigUtil import ConfigUtil

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()


class AlphaFoldModelProviderTests(unittest.TestCase):
    def setUp(self):
        self.__cachePath = os.path.join(HERE, "test-output", "CACHE")
        self.__startTime = time.time()

        # mockTopPath = os.path.join(TOPDIR, "rcsb", "mock-data")
        # configPath = os.path.join(mockTopPath, "config", "dbload-setup-example.yml")
        # self.__configName = "site_info_configuration"
        # self.__cfgOb = ConfigUtil(configPath=configPath, defaultSectionName=self.__configName, mockTopPath=mockTopPath)

        logger.info("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        unitS = "MB" if platform.system() == "Darwin" else "GB"
        rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum resident memory size %.4f %s", rusageMax / 10 ** 6, unitS)
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testAlphaFoldModelProvider(self):
        # First test fetching model archive
        # aFMP = AlphaFoldModelProvider(cachePath=self.__cachePath, useCache=False, alphaFoldRequestedSpeciesList=["Staphylococcus aureus"])
        # ok = aFMP.testCache()
        # self.assertTrue(ok)
        #
        # Next test reloading the cache
        aFMP = AlphaFoldModelProvider(cachePath=self.__cachePath, useCache=True, alphaFoldRequestedSpeciesList=["Staphylococcus aureus"])
        speciesDirList = aFMP.getArchiveDirList()
        ok = True if len(speciesDirList) > 0 else False
        self.assertTrue(ok)
        #
        speciesModelFileList = aFMP.getModelFileList(inputPathList=speciesDirList)
        ok = True if len(speciesModelFileList) > 0 else False
        self.assertTrue(ok)
        #
        # Next test reorganizing model file directory structure
        aFMP = AlphaFoldModelProvider(cachePath=self.__cachePath, useCache=True, alphaFoldRequestedSpeciesList=["Staphylococcus aureus"])
        ok = aFMP.testCache()
        self.assertTrue(ok)
        ok = aFMP.reorganizeModelFiles(copyModelsToDestDir=True)
        self.assertTrue(ok)


def fetchAlphaFoldModels():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(AlphaFoldModelProviderTests("testAlphaFoldModelProvider"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = fetchAlphaFoldModels()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
