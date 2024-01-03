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
        # This is where the models will be downloaded to and stored, prior to processing and reorganization
        # should stay the same regardless of where you want to reorganize the processed models
        # IMPORTANT: MAKE SURE THIS LOCATION HAS ENOUGH STORAGE TO HOLD ALL THE DOWNLOADED MODELS!
        self.__workPath = None  # "/PATH/TO/GIANT/_SOURCE_/DIRECTORY"

        # This is where the models will be reorganized into after processing
        # IMPORTANT: MAKE SURE THIS LOCATION HAS ENOUGH STORAGE TO HOLD ALL THE PROCESSED/ORGANIZED MODELS!
        self.__cachePath = None  # "/PATH/TO/GIANT/_ORGANIZED_/DIRECTORY"

        self.__keepSource = True
        self.__fetchModelArchive = True
        self.__fetchAlphaFold = True
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
            if self.__fetchModelArchive:
                mPWf = ModelProviderWorkflow(
                    srcDir=self.__workPath,
                    destDir=self.__cachePath,
                    modelProviders=["ModelArchive"],
                    useCache=True,
                    numProc=8,
                    chunkSize=100,
                )
                ok = mPWf.download()
                self.assertTrue(ok)
                ok = mPWf.reorganize(keepSource=self.__keepSource)
                self.assertTrue(ok)

            alphaFoldSpeciesList = [
                "Helicobacter pylori",  # 1538,
                "Staphylococcus aureus",  # 2888,
                "Arabidopsis thaliana",  # 27434
                "Caenorhabditis elegans",  # 19694
                "Candida albicans",  # 5974,
                "Danio rerio",  # 24664
                "Dictyostelium discoideum",  # 12622
                "Drosophila melanogaster",  # 13458
                "Escherichia coli",  # 4363,
                "Glycine max",  # 55799
                "Homo sapiens",  # 23391
                "Methanocaldococcus jannaschii",  # 1773,
                "Mus musculus",  # 21615
                "Oryza sativa",  # 43649
                "Rattus norvegicus",  # 21272
                "Saccharomyces cerevisiae",  # 6040,
                "Schizosaccharomyces pombe",  # 5128,
                "Zea mays",  # 39299
                "Ajellomyces capsulatus",  # 9199,
                "Brugia malayi",  # 8743,
                "Campylobacter jejuni",  # 1620,
                "Cladophialophora carrionii",  # 11170
                "Dracunculus medinensis",  # 10834
                "Enterococcus faecium",  # 2823,
                "Fonsecaea pedrosoi",  # 12509
                "Haemophilus influenzae",  # 1662,
                "Klebsiella pneumoniae",  # 5727,
                "Leishmania infantum",  # 7924,
                "Madurella mycetomatis",  # 9561,
                "Mycobacterium leprae",  # 1602,
                "Mycobacterium tuberculosis",  # 3988,
                "Mycobacterium ulcerans",  # 9033,
                "Neisseria gonorrhoeae",  # 2106,
                "Nocardia brasiliensis",  # 8372,
                "Onchocerca volvulus",  # 12047
                "Paracoccidioides lutzii",  # 8794,
                "Plasmodium falciparum",  # 5187,
                "Pseudomonas aeruginosa",  # 5556,
                "Salmonella typhimurium",  # 4526,
                "Schistosoma mansoni",  # 13865
                "Shigella dysenteriae",  # 3893,
                "Sporothrix schenckii",  # 8652,
                "Streptococcus pneumoniae",  # 2030,
                "Strongyloides stercoralis",  # 12613
                "Trichuris trichiura",  # 9564,
                "Trypanosoma brucei",  # 8491,
                "Trypanosoma cruzi",  # 19036
                "Wuchereria bancrofti",  # 12721
                "Swiss-Prot (CIF files)",  # 542380; This sometimes slows down significantly after only a few GBs (out of 36 GB) are downloaded, so may not finish
                "Overlap with MANE",  # 17334 (3,844 are unique)
            ]
            if self.__fetchAlphaFold and len(alphaFoldSpeciesList) > 0:
                for species in alphaFoldSpeciesList:
                    mPWf = ModelProviderWorkflow(
                        srcDir=self.__workPath,
                        destDir=self.__cachePath,
                        modelProviders=["AlphaFold"],
                        useCache=True,
                        numProc=8,
                        chunkSize=100,
                        alphaFoldRequestedSpeciesList=[species],
                    )
                    ok = mPWf.download()
                    self.assertTrue(ok)
                    ok = mPWf.reorganize(keepSource=self.__keepSource)
                    self.assertTrue(ok)


def modelProviderSuite():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(ModelProviderWorkflowTests("runModelProviderWorkflow"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = modelProviderSuite()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
