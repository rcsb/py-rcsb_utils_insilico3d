##
# File:    testModBaseModelProcessor.py
# Author:  Dennis Piehl
# Date:    27-Sep-2021
#
# Updates:
#
#
##
"""
Tests for processors for converting ModBase PDB models to mmCIF format.

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
from rcsb.utils.insilico3d.ModBaseModelProcessor import ModBaseModelProcessor
from rcsb.utils.config.ConfigUtil import ConfigUtil

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class ModBaseModelProcessorTests(unittest.TestCase):
    skipFlag = platform.system() != "Darwin"
    buildTestingCache = True

    def setUp(self):
        self.__dataPath = os.path.join(HERE, "test-data")
        self.__cachePath = os.path.join(HERE, "test-output", "CACHE")
        self.__startTime = time.time()
        mockTopPath = os.path.join(TOPDIR, "rcsb", "mock-data")
        configPath = os.path.join(mockTopPath, "config", "dbload-setup-example.yml")
        self.__configName = "site_info_configuration"
        self.__cfgOb = ConfigUtil(configPath=configPath, defaultSectionName=self.__configName, mockTopPath=self.__cachePath)

        logger.info("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        unitS = "MB" if platform.system() == "Darwin" else "GB"
        rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum resident memory size %.4f %s", rusageMax / 1.0e6, unitS)
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testModBaseModelProcessor(self):
        """Test case: convert ModBase PDB models to mmCIF"""
        try:
            ## INSTEAD OF redownloading data, use the test-data directory for processing

            # First download the archived data for processing
            # mProv = ModBaseModelProvider(
            #     cachePath=self.__cachePath,
            #     useCache=False,
            #     modBaseServerSpeciesDataPathDict={"Staphylococcus aureus": "S_aureus/2008/staph_aureus.tar"}
            # )
            # ok = mProv.testCache()
            # self.assertTrue(ok)
            #
            # Next test reload data for processing
            mProv = ModBaseModelProvider(
                cachePath=self.__cachePath,
                useCache=True,
                modBaseServerSpeciesDataPathDict={"Staphylococcus aureus": "S_aureus/2008/staph_aureus.tar"}
            )
            ok = mProv.testCache()
            self.assertTrue(ok)
            speciesNameList = mProv.getSpeciesNameList()
            ok = True if len(speciesNameList) > 0 else False
            self.assertTrue(ok)
            speciesConversionDict = mProv.getSpeciesConversionDict(speciesName=speciesNameList[0])
            speciesConversionDict["speciesPdbModelFileList"] = speciesConversionDict["speciesPdbModelFileList"][0:20]
            ok = True if len(speciesConversionDict["speciesPdbModelFileList"]) > 0 else False
            self.assertTrue(ok)
            mProc = ModBaseModelProcessor(
                useCache=False,
                numProc=2,
                speciesD=speciesConversionDict,
            )
            ok = mProc.generate(updateOnly=False)
            self.assertTrue(ok)
            ok = mProc.generate(updateOnly=True)
            self.assertTrue(ok)
            ok = mProc.reload()
            self.assertTrue(ok)
            ok = mProc.testCache(minCount=10)
            self.assertTrue(ok)
            #
            processedModelCachePath = mProc.getCachePath()
            mProc = ModBaseModelProcessor(
                cachePath=processedModelCachePath,
                useCache=True,
                numProc=2,
                speciesD=speciesConversionDict,
            )
            ok = mProc.testCache(minCount=10)
            self.assertTrue(ok)
            #
            # ok = mProv.removePdbModelDir(speciesDataDir=speciesConversionDict["speciesModelDir"])
            # self.assertTrue(ok)
            # ok = mProv.removeAlignmentDir(speciesDataDir=speciesConversionDict["speciesModelDir"])
            # self.assertTrue(ok)
            #
            # Test reorganize models
            ok = mProv.testCache()
            self.assertTrue(ok)
            ok = mProv.reorganizeModelFiles()
            self.assertTrue(ok)
            #
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()


def modelProcessorSuite():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(ModBaseModelProcessorTests("testModBaseModelProcessor"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = modelProcessorSuite()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
