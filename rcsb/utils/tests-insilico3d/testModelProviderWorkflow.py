##
# File:    testModelProviderWorkflow.py
# Author:  Dennis Piehl
# Date:    28-Mar-2022
#
# Updates:
#
#
##
"""
Tests for workflow for retrieving and storing all models.

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
from rcsb.utils.insilico3d.ModelArchiveModelProvider import ModelArchiveModelProvider
from rcsb.utils.insilico3d.ModelProviderWorkflow import ModelProviderWorkflow
from rcsb.utils.config.ConfigUtil import ConfigUtil

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()


class ModelProviderWorkflowTests(unittest.TestCase):

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

    def testModelProviderWorkflow(self):
        redownloadBulkData = True
        # redownloadBulkData = False
        # keepSource = True
        keepSource = False

        #
        # mPWf = ModelProviderWorkflow(
        #     cachePath=self.__cachePath,
        #     useCache=False,
        #     cfgOb=self.__cfgOb,
        #     configName=self.__configName,
        #     numProc=4,
        #     chunkSize=40,
        #     redownloadBulkData=redownloadBulkData
        # )
        # ok = mPWf.run()
        # self.assertTrue(ok)

        # First test fetching model archive
        if redownloadBulkData:
            aFMP = AlphaFoldModelProvider(
                cachePath=self.__cachePath,
                useCache=False,
                cfgOb=self.__cfgOb,
                configName=self.__configName,
                numProc=4,
                chunkSize=40,
                # alphaFoldRequestedSpeciesList=["Helicobacter pylori"]
                alphaFoldRequestedSpeciesList=["Helicobacter pylori", "Escherichia coli", "Methanocaldococcus jannaschii", "Homo sapiens"]
            )
            ok = aFMP.testCache()
            self.assertTrue(ok)

            mAMP = ModelArchiveModelProvider(
                cachePath=self.__cachePath,
                useCache=False,
                cfgOb=self.__cfgOb,
                configName=self.__configName,
                numProc=4,
                chunkSize=20,
                # serverDataSetPathD={"ma-bak-cepc": {"urlEnd": "ma-bak-cepc?type=materials_procedures__accompanying_data_file_name", "fileName": "ma-bak-cepc.zip"}}
            )
            ok = mAMP.testCache()
            self.assertTrue(ok)

        #
        # Next test reloading the cache
        aFMP = AlphaFoldModelProvider(
            cachePath=self.__cachePath,
            useCache=True,
            cfgOb=self.__cfgOb,
            configName=self.__configName,
            numProc=4,
            chunkSize=40,
            # alphaFoldRequestedSpeciesList=["Helicobacter pylori"]
            alphaFoldRequestedSpeciesList=["Helicobacter pylori", "Escherichia coli", "Methanocaldococcus jannaschii", "Homo sapiens"]
        )
        speciesDirList = aFMP.getArchiveDirList()
        ok = True if len(speciesDirList) > 0 else False
        self.assertTrue(ok)

        mAMP = ModelArchiveModelProvider(
            cachePath=self.__cachePath,
            useCache=True,
            cfgOb=self.__cfgOb,
            configName=self.__configName,
            numProc=4,
            chunkSize=20,
            # serverDataSetPathD={"ma-bak-cepc": {"urlEnd": "ma-bak-cepc?type=materials_procedures__accompanying_data_file_name", "fileName": "ma-bak-cepc.zip"}}
        )
        archiveDirList = mAMP.getArchiveDirList()
        ok = True if len(archiveDirList) > 0 else False
        self.assertTrue(ok)

        #
        # Next test reorganizing model file directory structure
        ok = aFMP.reorganizeModelFiles(useCache=True, numProc=4, chunkSize=40, keepSource=keepSource)
        self.assertTrue(ok)

        ok = mAMP.reorganizeModelFiles(useCache=True, numProc=4, chunkSize=20, keepSource=keepSource)
        self.assertTrue(ok)




def fetchAlphaFoldModels():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(ModelProviderWorkflowTests("testModelProviderWorkflow"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = fetchAlphaFoldModels()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
