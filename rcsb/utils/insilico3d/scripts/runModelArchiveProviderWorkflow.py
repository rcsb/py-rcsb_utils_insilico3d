##
# File:    runModelArchiveProviderWorkflow.py
# Author:  Dennis Piehl
# Date:    15-Apr-2022
#
# Updates:
#
#
##
"""
Script for running through the entire workflow for retrieving and storing all ModelArchive models.

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


class ModelProviderWorkflowExec(unittest.TestCase):

    def setUp(self):
        # This is where the models will be downloaded to and stored, prior to processing and reorganization
        # should stay the same regardless of where you want to reorganize the processed models
        self.__workPath = "/mnt/vdb1/source-models/"  # "/PATH/TO/GIANT/_SOURCE_/DIRECTORY"

        # This is where the models will be reorganized into after processing
        self.__cachePath = "/mnt/vdb1/computed-models/CSM1"  # "/PATH/TO/GIANT/_ORGANIZED_/DIRECTORY"

        # This controls whether to keep the downloaded source files after processing and reorganizing them.
        # It's generally a good idea to keep them around in case the reorganization step of the workflow needs
        # to be run again. (*Note that this flag means something different for GoogleCloud-sourced AlphaFold models!)
        self.__keepSource = True

        self.__fetchAndReorganizeModelArchive = True
        self.__startTime = time.time()
        logger.info("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        unitS = "MB" if platform.system() == "Darwin" else "GB"
        rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum resident memory size %.4f %s", rusageMax / 10 ** 6, unitS)
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def runModelProviderWorkflow(self):
        if self.__workPath and self.__cachePath:
            if self.__fetchAndReorganizeModelArchive:
                mPWf = ModelProviderWorkflow(
                    srcDir=self.__workPath,
                    destDir=self.__cachePath,
                    modelProviders=["ModelArchive"],
                    useCache=True,
                    numProc=16,
                    chunkSize=16,
                    # modelArchiveRequestedDatasetD={"ma-bak-cepc": {}, "ma-ornl-sphdiv": {}, "ma-coffe-slac": {}, "ma-asfv-asfvg": {}, "ma-t3vr3": {}}
                )
                ok = mPWf.download()
                self.assertTrue(ok)
                ok = mPWf.reorganize(keepSource=self.__keepSource)
                self.assertTrue(ok)


def modelProviderSuite():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(ModelProviderWorkflowExec("runModelProviderWorkflow"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = modelProviderSuite()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
