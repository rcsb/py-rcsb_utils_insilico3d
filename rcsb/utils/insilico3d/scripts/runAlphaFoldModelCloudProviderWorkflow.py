##
# File:    runAlphaFoldModelCloudProviderWorkflow.py
# Author:  Dennis Piehl
# Date:    3-Jan-2024
#
# Updates:
#
#
##
"""
Script for running through the reorganization workflow for AlphaFold models downloaded from Google Cloud.

Assumes tar files have already been downloaded.

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
    runEntireWorkflowForAllProviderSources = True

    def setUp(self):
        # This is where the models will be downloaded to and stored, prior to processing and reorganization
        # should stay the same regardless of where you want to reorganize the processed models
        self.__workPath = "/mnt/vdb1/source-models/"

        # This is where the models will be reorganized into after processing
        self.__cachePath = "/mnt/vdb1/computed-models/CSM1"
        #
        self.__keepSource = False

        self.__startTime = time.time()
        logger.info("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        unitS = "MB" if platform.system() == "Darwin" else "GB"
        rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum resident memory size %.4f %s", rusageMax / 10 ** 6, unitS)
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    @unittest.skipUnless(runEntireWorkflowForAllProviderSources, "Skip running the entire workflow")
    def runModelProviderWorkflow(self):
        if self.__workPath and self.__cachePath:

            alphaFoldSubDirL = ["100"]
            # alphaFoldSubDirL = [str(i) for i in range(100, 1000)]

            mPWf = ModelProviderWorkflow(
                srcDir=self.__workPath,
                destDir=self.__cachePath,
                modelProviders=["AlphaFoldCloud"],
                useCache=True,
                numProc=16,
                chunkSize=4,  # the smaller the faster, generally
                smallFileSizeCutoff=8388608,  # 8mb (~130 models), seems to be the optimal cutoff
            )
            for subDir in alphaFoldSubDirL:
                ok = mPWf.reorganize(
                    keepSource=self.__keepSource,
                    inputDirList=[subDir],
                )
                self.assertTrue(ok)


def modelProviderSuite():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(ModelProviderWorkflowTests("runModelProviderWorkflow"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = modelProviderSuite()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
