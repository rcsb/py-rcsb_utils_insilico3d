##
# File:    runModelProviderWorkflow.py
# Author:  Dennis Piehl
# Date:    15-Apr-2022
#
# Updates:
#
#
##
"""
Script for running through the entire workflow for retrieving and storing all models.

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
    runEntireWorkflowForAllProviderSources = False

    def setUp(self):
        self.__cachePath = "/mnt/vdb1/computed-models/"
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
        mPWf = ModelProviderWorkflow(
            destDir=self.__cachePath,
            modelProviders=["ModelArchive", "AlphaFold"],
            useCache=True,
            numProc=6,
            chunkSize=40,
            # alphaFoldRequestedSpeciesList=["Swiss-Prot (CIF files)"],  # This sometimes slows down significantly after only a few GBs (out of 36 GB) are downloaded, so may not finish
            alphaFoldRequestedSpeciesList=[
                "Arabidopsis thaliana",  # 27,434
                "Caenorhabditis elegans",  # 19,694
                "Candida albicans",  # 5,974
                "Danio rerio",  # 24,664
                "Dictyostelium discoideum",  # 12,622
                "Drosophila melanogaster",  # 13,458
                "Escherichia coli",  # 4,363
                "Glycine max",  # 55,799
                "Homo sapiens",  # 23,391
                "Methanocaldococcus jannaschii",  # 1,773
                "Mus musculus",  # 21,615
                "Oryza sativa",  # 43,649
                "Rattus norvegicus",  # 21,272
                "Saccharomyces cerevisiae",  # 6,040
                "Schizosaccharomyces pombe",  # 5,128
                "Zea mays",  # 39,299
                "Ajellomyces capsulatus",  # 9,199
                "Brugia malayi",  # 8,743
                "Campylobacter jejuni",  # 1,620
                "Cladophialophora carrionii",  # 11,170
                "Dracunculus medinensis",  # 10,834
                "Enterococcus faecium",  # 2,823
                "Fonsecaea pedrosoi",  # 12,509
                "Haemophilus influenzae",  # 1,662
                "Helicobacter pylori",  # 1,538
            ],
        )
        ok = mPWf.download()
        self.assertTrue(ok)
        ok = mPWf.reorganize()
        self.assertTrue(ok)


def fetchAlphaFoldModels():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(ModelProviderWorkflowTests("runModelProviderWorkflow"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = fetchAlphaFoldModels()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
