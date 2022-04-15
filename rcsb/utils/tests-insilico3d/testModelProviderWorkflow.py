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

from rcsb.utils.insilico3d.ModelProviderWorkflow import ModelProviderWorkflow

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class ModelProviderWorkflowTests(unittest.TestCase):

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

    def testModelProviderWorkflow(self):
        mPWf = ModelProviderWorkflow(
            destDir=self.__cachePath,
            modelProviders=["AlphaFold"],
            useCache=True,
            numProc=4,
            chunkSize=40,
            alphaFoldRequestedSpeciesList=["Helicobacter pylori"],
        )
        ok = mPWf.download()
        self.assertTrue(ok)
        ok = mPWf.reorganize(keepSource=True)
        self.assertTrue(ok)


def fetchAlphaFoldModels():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(ModelProviderWorkflowTests("testModelProviderWorkflow"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = fetchAlphaFoldModels()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
