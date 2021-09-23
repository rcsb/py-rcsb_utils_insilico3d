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

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class ModBaseModelProcessorTests(unittest.TestCase):
    skipFlag = platform.system() != "Darwin"
    buildTestingCache = True

    def setUp(self):
        self.__cachePath = os.path.join(HERE, "test-output", "CACHE")
        self.__startTime = time.time()
        logger.info("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

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
            mPr = ModBaseModelProcessor(useCache=False, numProc=2, speciesModelDir=speciesModelDir, speciesPdbModelFileList=speciesPdbModelFileList)
            ok = mPr.generate(updateOnly=False)
            self.assertTrue(ok)
            ok = mPr.generate(updateOnly=True)
            self.assertTrue(ok)
            ok = mPr.reload()
            self.assertTrue(ok)
            ok = mPr.testCache(minCount=10)
            self.assertTrue(ok)
            #
            processedModelCachePath = mPr.getCachePath()
            mPr = ModBaseModelProcessor(cachePath=processedModelCachePath, useCache=True, numProc=2,
                                        speciesModelDir=speciesModelDir, speciesPdbModelFileList=speciesPdbModelFileList)
            ok = mPr.testCache(minCount=10)
            self.assertTrue(ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()


def modelProcessorSuite():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(ModBaseModelProcessorTests("testModBaseModelProcessorBootstrap"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = modelProcessorSuite()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
