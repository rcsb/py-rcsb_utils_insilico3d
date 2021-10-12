##
# File:    testSwissModelProvider.py
# Author:  Dennis Piehl
# Date:    29-Sep-2021
#
# Update:
#
#
##
"""
Tests for accessor utilities for SWISS-MODEL 3D Model (PDB) data.

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

from rcsb.utils.insilico3d.SwissModelProvider import SwissModelProvider

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()


class SwissModelProviderTests(unittest.TestCase):
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

    def testSwissModelProvider(self):
        # First test fetching model archive
        mProv = SwissModelProvider(
            cachePath=self.__cachePath,
            useCache=False,
            swissModelServerSpeciesDataPathDict={"Staphylococcus aureus": "93061_coords.tar.gz"}
        )
        ok = mProv.testCache()
        self.assertTrue(ok)
        #
        # Next test reloading the cache
        mProv = SwissModelProvider(
            cachePath=self.__cachePath,
            useCache=True,
            swissModelServerSpeciesDataPathDict={"Staphylococcus aureus": "93061_coords.tar.gz"}
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
        mProv = SwissModelProvider(
            cachePath=self.__cachePath,
            useCache=True,
            swissModelServerSpeciesDataPathDict={"Staphylococcus aureus": "93061_coords.tar.gz"}
        )
        speciesNameList = mProv.getSpeciesNameList()
        for species in speciesNameList:
            ok = mProv.removeSpeciesDataDir(speciesName=species)
            self.assertTrue(ok)

    # def testFetchSwissModels(self):
    #     mProv = SwissModelProvider(
    #         cachePath=self.__cachePath,
    #         useCache=False,
    #         swissModelServerSpeciesDataPathDict={"Staphylococcus aureus": "93061_coords.tar.gz"}
    #     )
    #     ok = mProv.testCache()
    #     self.assertTrue(ok)

    # def testReloadCache(self):
    #     mProv = SwissModelProvider(
    #         cachePath=self.__cachePath,
    #         useCache=True,
    #         swissModelServerSpeciesDataPathDict={"Staphylococcus aureus": "93061_coords.tar.gz"}
    #     )
    #     speciesDirList = mProv.getSpeciesDirList()
    #     ok = True if len(speciesDirList) > 0 else False
    #     self.assertTrue(ok)
    #     #
    #     speciesPdbModelFileList = mProv.getSpeciesPdbModelFileList(speciesDataDir=speciesDirList[0])
    #     ok = True if len(speciesPdbModelFileList) > 0 else False
    #     self.assertTrue(ok)

    # def testDeleteCache(self):
    #     mProv = SwissModelProvider(
    #         cachePath=self.__cachePath,
    #         useCache=True,
    #         swissModelServerSpeciesDataPathDict={"Staphylococcus aureus": "93061_coords.tar.gz"}
    #     )
    #     speciesNameList = mProv.getSpeciesNameList()
    #     for species in speciesNameList:
    #         ok = mProv.removeSpeciesDataDir(speciesName=species)
    #         self.assertTrue(ok)


def fetchSwissModels():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(SwissModelProviderTests("testSwissModelProvider"))
    # suiteSelect.addTest(SwissModelProviderTests("testFetchSwissModels"))
    # suiteSelect.addTest(SwissModelProviderTests("testReloadCache"))
    # suiteSelect.addTest(SwissModelProviderTests("testDeleteCache"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = fetchSwissModels()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
