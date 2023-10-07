##
# File:    testAlphaFoldModelCloudProvider.py
# Author:  Dennis Piehl
# Date:    30-Sep-2021
#
# Updates:
#
#
##
"""
Tests for accessor utilities for AlphaFold 3D Model (mmCIF) data on Google Cloud.

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

from rcsb.utils.insilico3d.AlphaFoldModelCloudProvider import AlphaFoldModelCloudProvider

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class AlphaFoldModelCloudProviderTests(unittest.TestCase):

    def setUp(self):
        self.__cachePath = os.path.join(HERE, "test-output", "CACHE", "computed-models")
        self.__startTime = time.time()
        logger.info("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        unitS = "MB" if platform.system() == "Darwin" else "GB"
        rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum resident memory size %.4f %s", rusageMax / 10 ** 6, unitS)
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testAlphaFoldModelCloudProvider(self):
        redownloadBulkData = True
        #
        # First test fetching model archive
        if redownloadBulkData:
            aFMCP = AlphaFoldModelCloudProvider(
                cachePath=self.__cachePath,
                useCache=False,
                numProc=4,
                chunkSize=20,
                # alphaFoldRequestedTaxIdPrefixList = ["100000"]
                # alphaFoldRequestedSpeciesList=[{"species": "Panicum virgatum", "common_name": "Switchgrass", "taxIds": ["206033"]}]
            )
            ok = aFMCP.testCache()
            self.assertTrue(ok)
        #
        # Next test reloading the cache
        aFMCP = AlphaFoldModelCloudProvider(
            cachePath=self.__cachePath,
            useCache=True,
            numProc=4,
            chunkSize=20,
            # alphaFoldRequestedSpeciesList=[{"species": "Panicum virgatum", "common_name": "Switchgrass", "taxIds": ["206033"]}]
        )
        taxIdPrefixDirList = aFMCP.getArchiveDirList()
        ok = True if len(taxIdPrefixDirList) > 0 else False
        self.assertTrue(ok)
        # #
        archiveFileList = aFMCP.getArchiveFileList(inputPathList=taxIdPrefixDirList)
        ok = True if len(archiveFileList) > 0 else False
        self.assertTrue(ok)
        ok = aFMCP.testCache()
        self.assertTrue(ok)
        #
        # Next test reorganizing model file directory structure
        ok = aFMCP.extractAndReorganizeModelFiles(useCache=True, inputModelList=None, numProc=4, chunkSize=20, keepSource=True)
        self.assertTrue(ok)

        # # Now test using the reorganizer object directly
        # aFMR = aFMP.getModelReorganizer(useCache=False, numProc=4, chunkSize=20, keepSource=True)
        # destBaseDir = aFMP.getComputedModelsDataPath()
        # ok = aFMR.reorganize(inputModelList=speciesModelFileList[10:20], modelSource="AlphaFold", destBaseDir=destBaseDir, useCache=False)
        # self.assertTrue(ok)
        # ok = aFMR.testCache()
        # self.assertFalse(ok)  # Confirm that testCache FAILED (< 20 in cache)
        # ok = aFMR.reorganize(inputModelList=speciesModelFileList[20:30], modelSource="AlphaFold", destBaseDir=destBaseDir, useCache=True)
        # self.assertTrue(ok)
        # ok = aFMR.testCache()
        # self.assertTrue(ok)  # Confirm that testCache SUCCEEDED (>= 20 in cache)


def fetchAlphaFoldModels():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(AlphaFoldModelCloudProviderTests("testAlphaFoldModelCloudProvider"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = fetchAlphaFoldModels()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
