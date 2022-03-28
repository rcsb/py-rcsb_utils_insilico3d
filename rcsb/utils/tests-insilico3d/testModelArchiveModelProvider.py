##
# File:    testModelArchiveModelProvider.py
# Author:  Dennis Piehl
# Date:    18-Mar-2022
#
# Updates:
#
#
##
"""
Tests for accessor utilities for ModelArchive 3D Model (mmCIF) data.

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

from rcsb.utils.insilico3d.ModelArchiveModelProvider import ModelArchiveModelProvider
from rcsb.utils.config.ConfigUtil import ConfigUtil

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()


class ModelArchiveModelProviderTests(unittest.TestCase):

    def setUp(self):
        # self.__dataPath = os.path.join(HERE, "test-data")
        self.__cachePath = os.path.join(HERE, "test-output", "CACHE")
        self.__startTime = time.time()

        mockTopPath = os.path.join(TOPDIR, "rcsb", "mock-data")
        configPath = os.path.join(mockTopPath, "config", "dbload-setup-example.yml")
        self.__configName = "site_info_configuration"
        # Set mockTopPath to test-output directory here, because for testing the downloading and storing of the models we don't want to actually modify the mock-data directory.
        # Only use the mock-data directory for testing the LOADING of models to ExDB, etc.
        self.__cfgOb = ConfigUtil(configPath=configPath, defaultSectionName=self.__configName, mockTopPath=self.__cachePath)

        logger.info("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        unitS = "MB" if platform.system() == "Darwin" else "GB"
        rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum resident memory size %.4f %s", rusageMax / 10 ** 6, unitS)
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testModelArchiveModelProvider(self):
        redownloadBulkData = True
        #
        # First test fetching model archive
        if redownloadBulkData:
            mAMP = ModelArchiveModelProvider(
                cachePath=self.__cachePath,
                useCache=False,
                cfgOb=self.__cfgOb,
                configName=self.__configName,
                numProc=4,
                chunkSize=20,
                serverDataSetPathD={"ma-bak-cepc": {"urlEnd": "ma-bak-cepc?type=materials_procedures__accompanying_data_file_name", "fileName": "ma-bak-cepc.zip"}}
            )
            ok = mAMP.testCache()
            self.assertTrue(ok)
        #
        # Next test reloading the cache
        mAMP = ModelArchiveModelProvider(
            cachePath=self.__cachePath,
            useCache=True,
            cfgOb=self.__cfgOb,
            configName=self.__configName,
            numProc=4,
            chunkSize=20,
            serverDataSetPathD={"ma-bak-cepc": {"urlEnd": "ma-bak-cepc?type=materials_procedures__accompanying_data_file_name", "fileName": "ma-bak-cepc.zip"}}
        )
        archiveDirList = mAMP.getArchiveDirList()
        ok = True if len(archiveDirList) > 0 else False
        self.assertTrue(ok)
        #
        archiveModelFileList = mAMP.getModelFileList(inputPathList=archiveDirList)
        ok = True if len(archiveModelFileList) > 0 else False
        self.assertTrue(ok)
        ok = mAMP.testCache()
        self.assertTrue(ok)
        #
        # Next test reorganizing model file directory structure
        ok = mAMP.reorganizeModelFiles(useCache=False, inputModelList=archiveModelFileList[0:10], numProc=4, chunkSize=20, keepSource=True)
        self.assertTrue(ok)
        # Now test using the reorganizer object directly
        mAMR = mAMP.getModelReorganizer(useCache=False, numProc=4, chunkSize=20, keepSource=True)
        destBaseDir = mAMP.getComputedModelsDataPath()
        ok = mAMR.reorganize(inputModelList=archiveModelFileList[10:20], modelSource="ModelArchive", destBaseDir=destBaseDir, useCache=False)
        self.assertTrue(ok)
        ok = mAMR.testCache()
        self.assertFalse(ok)  # Confirm that testCache FAILED (< 20 in cache)
        ok = mAMR.reorganize(inputModelList=archiveModelFileList[20:30], modelSource="ModelArchive", destBaseDir=destBaseDir, useCache=True)
        self.assertTrue(ok)
        ok = mAMR.testCache()
        self.assertTrue(ok)  # Confirm that testCache SUCCEEDED (>= 20 in cache)


def fetchModelArchiveModels():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(ModelArchiveModelProviderTests("testModelArchiveModelProvider"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = fetchModelArchiveModels()
    unittest.TextTestRunner(verbosity=2).run(mySuite)