##
# File:    testModBaseModelProvider.py
# Author:  Dennis Piehl
# Date:    27-Sep-2021
#
# Updates:
#
#
##
"""
Tests for accessor utilities for ModBase 3D Model (bulk download) data.

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

    def testModBaseModelProvider(self):
        # First test fetching model archive
        mProv = ModBaseModelProvider(
            cachePath=self.__cachePath,
            useCache=False,
            modBaseServerSpeciesDataPathDict={"Staphylococcus aureus": "S_aureus/2008/staph_aureus.tar"}
        )
        ok = mProv.testCache()
        self.assertTrue(ok)
        #
        # Next test reloading the cache
        mProv = ModBaseModelProvider(
            cachePath=self.__cachePath,
            useCache=True,
            modBaseServerSpeciesDataPathDict={"Staphylococcus aureus": "S_aureus/2008/staph_aureus.tar"}
        )
        speciesDirList = mProv.getSpeciesDirList()
        ok = True if len(speciesDirList) > 0 else False
        self.assertTrue(ok)
        #
        speciesPdbModelFileList = mProv.getSpeciesPdbModelFileList(speciesDataDir=speciesDirList[0])
        ok = True if len(speciesPdbModelFileList) > 0 else False
        self.assertTrue(ok)
        #
        # Last test deleting the cache
        mProv = ModBaseModelProvider(
            cachePath=self.__cachePath,
            useCache=True,
            modBaseServerSpeciesDataPathDict={"Staphylococcus aureus": "S_aureus/2008/staph_aureus.tar"}
        )
        speciesNameList = mProv.getSpeciesNameList()
        for species in speciesNameList:
            ok = mProv.removeSpeciesDataDir(speciesName=species)
            self.assertTrue(ok)


def fetchModBaseModels():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(ModBaseModelProviderTests("testModBaseModelProvider"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = fetchModBaseModels()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
