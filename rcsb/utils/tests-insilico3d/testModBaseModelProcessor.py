##
# File:    testModBaseModelProcessor.py
# Author:  Dennis Piehl
# Date:    21-Sep-2021
#
# Update:
#
#
##
"""
Tests for generators and accessors for non-polymer instance target interactions

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


from ModBaseModelProvider import ModBaseModelProvider
from ModBaseModelProcessor import ModBaseModelProcessor
# from rcsb.utils.insilico3d.ModBaseModelProvider import ModBaseModelProvider
# from rcsb.utils.insilico3d.ModBaseModelProcessor import ModBaseModelProcessor

# from rcsb.utils.config.ConfigUtil import ConfigUtil

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

# logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logging.basicConfig(level=logging.DEBUG, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class ModBaseModelProcessorTests(unittest.TestCase):
    skipFlag = platform.system() != "Darwin"
    buildTestingCache = True

    def setUp(self):
        # mockTopPath = os.path.join(TOPDIR, "rcsb", "mock-data")
        self.__cachePath = os.path.join(HERE, "test-output", "CACHE")
        # configPath = os.path.join(mockTopPath, "config", "dbload-setup-example.yml")
        # self.__configName = "site_info_configuration"
        self.__startTime = time.time()
        logger.info("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))
        self.__fileLimit = None if platform.system() == "Darwin" else 5
        #
        # self.__cfgOb = ConfigUtil(configPath=configPath, defaultSectionName=self.__configName, mockTopPath=mockTopPath)

    def tearDown(self):
        unitS = "MB" if platform.system() == "Darwin" else "GB"
        rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum resident memory size %.4f %s", rusageMax / 1.0e6, unitS)
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testModBaseModelProcessorBootstrap(self):
        """Test case: generate and load neighbor and occupancy data"""
        try:
            # mP = ModBaseModelProvider(cachePath=self.__cachePath, useCache=False, modBaseSpeciesDataPathDict={"Staphylococcus aureus": "S_aureus/2008/staph_aureus.tar"})
            mP = ModBaseModelProvider(cachePath=self.__cachePath, useCache=True, modBaseSpeciesDataPathDict={"Staphylococcus aureus": "S_aureus/2008/staph_aureus.tar"})
            ok = mP.testCache()
            self.assertTrue(ok)
            speciesDirList = mP.getSpeciesDirList()
            ok = True if len(speciesDirList) > 0 else False
            self.assertTrue(ok)
            speciesModelDir = speciesDirList[0]
            speciesPdbModelFileList = mP.getSpeciesPdbModelFileList(speciesDataDir=speciesModelDir)[0:20]
            # print(speciesPdbModelFileList)
            # print(self.__cachePath)
            mPr = ModBaseModelProcessor(cachePath=self.__cachePath, useCache=False, numProc=2, fileLimit=self.__fileLimit,
                                        speciesModelDir=speciesModelDir, speciesPdbModelFileList=speciesPdbModelFileList)
            ok = mPr.generate(updateOnly=False)
            # self.assertTrue(ok)
            # ok = mPr.generate(updateOnly=True)
            # self.assertTrue(ok)
            # ok = mPr.reload()
            # self.assertTrue(ok)
            # ok = mPr.testCache(minCount=self.__fileLimit if self.__fileLimit else 30)
            # self.assertTrue(ok)
            #
            # mPr = ModBaseModelProcessor(self.__cachePath, useCache=True, cfgOb=self.__cfgOb, configName=self.__configName, numProc=2)
            # ok = mPr.testCache(minCount=self.__fileLimit if self.__fileLimit else 30)
            # self.assertTrue(ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    # @unittest.skipIf(skipFlag, "Long test")
    # def testStashRemote(self):
    #     try:
    #         mPr = ModBaseModelProcessor(self.__cachePath, useCache=True, cfgOb=self.__cfgOb, configName=self.__configName, numProc=2, fileLimit=self.__fileLimit)
    #         ok = mPr.testCache()
    #         self.assertTrue(ok)
    #         ok = mPr.restore(self.__cfgOb, self.__configName, remotePrefix=None, useGit=True, useStash=True)
    #         self.assertTrue(ok)
    #         ok = mPr.reload()
    #         self.assertTrue(ok)
    #         #
    #     except Exception as e:
    #         logger.exception("Failing with %s", str(e))
    #         self.fail()


def modelProcessorSuite():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(ModBaseModelProcessorTests("testModBaseModelProcessorBootstrap"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = modelProcessorSuite()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
