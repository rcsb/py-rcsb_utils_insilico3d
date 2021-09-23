##
# File:    testModBaseModelProvider.py
# Author:  Dennis Piehl
# Date:    15-Sep-2021
#
# Update:
#
#
##
"""
Tests for accessor utilities for ModBase 3D Model (mmCIF) data.

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

from rcsb.utils.insilico3d.ModBaseModelProvider import ModBaseModelProvider

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()


class ModBaseModelProviderTests(unittest.TestCase):
    def setUp(self):
        self.__cachePath = os.path.join(HERE, "test-output", "CACHE")
        self.__startTime = time.time()
        logger.info("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        unitS = "MB" if platform.system() == "Darwin" else "GB"
        rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum resident memory size %.4f %s", rusageMax / 10 ** 6, unitS)
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testFetchModBaseModels(self):
        mP = ModBaseModelProvider(cachePath=self.__cachePath, useCache=False, modBaseSpeciesDataPathDict={"Staphylococcus aureus": "S_aureus/2008/staph_aureus.tar"})
        ok = mP.testCache()
        self.assertTrue(ok)

    def testReloadCache(self):
        mP = ModBaseModelProvider(cachePath=self.__cachePath, useCache=True, modBaseSpeciesDataPathDict={"Staphylococcus aureus": "S_aureus/2008/staph_aureus.tar"})
        speciesDirList = mP.getSpeciesDirList()
        ok = True if len(speciesDirList) > 0 else False
        self.assertTrue(ok)

        speciesPdbModelFileList = mP.getSpeciesPdbModelFileList(speciesDataDir=speciesDirList[0])
        ok = True if len(speciesPdbModelFileList) > 0 else False
        self.assertTrue(ok)


def fetchModBaseModels():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(ModBaseModelProviderTests("testFetchModBaseModels"))
    suiteSelect.addTest(ModBaseModelProviderTests("testReloadCache"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = fetchModBaseModels()
    unittest.TextTestRunner(verbosity=2).run(mySuite)