##
# File:    runModelProviderWorkflow.py
# Author:  Dennis Piehl
# Date:    20-Oct-2025
#
# Updates:
#
#
##
"""
Script for running through the entire workflow for retrieving and storing AlphaFold models from the FTP site (i.e., the
original ~1 million subset) and ModelArchive.

When we switch to sourcing AlphaFold from Google Cloud for 200M, should replace FTP-based AlphaFold step below with
Cloud-based method (see runAlphaFoldCloudReorganizingWorkflow.py).

To run:
    0. Clone this repository to the `fa-csm` instance, and make sure the package is installed in Python3 environment
    1. Make a backup of the current holdings files and source model file cache/holdings
       (i.e., under `/mnt/vdb1/computed-models/staging/holdings` and `/mnt/vdb1/source-models/work-dir/*/*.json`)
    2. Delete the pre-existing source model files from `/mnt/vdb1/source-models`
    3. Configure the settings in `setUp()` to desired values (note that it's recommeded to first run this
       for ModelArchive--if needed--and then AlphaFold)
    4. If you only want to download certain species or data sets from each source, comment in the relevant lines
       in `runModelProviderWorkflow()`. (Note that the list of AF species below only controls which to DOWNLOAD; the
       REORGANIZATION is controlled by the source model cache file `/mnt/vdb1/source-models/work-dir/AlphaFold/model-download-cache.json`,
       in that the workflow will always try to reorganize every species in there with `"reogranized": false`)
    5. Run the script in the background and save log output with:
           python3 runModelProviderWorkflow.py >& log.1 &
    6. Once you run the workflow for both AF and MA, it will create two separate holdings files for each
       These should be merged, which you can do with the separate script, `mergeHoldingsFiles.py`
    7. Update the symlink for `/mnt/vdb1/computed-models/staging`
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

        self.__fetchAndReorganizeAlphaFoldFtp = True
        self.__fetchAndReorganizeModelArchive = True

        self.__numProc = 8
        self.__chunkSize = 8

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
            alphaFoldSpeciesList = [
                "Helicobacter pylori",  # 1540
                "Staphylococcus aureus",  # 2888
                "Arabidopsis thaliana",  # 27402
                "Caenorhabditis elegans",  # 19700
                "Candida albicans",  # 5973
                "Danio rerio",  # 26290
                "Dictyostelium discoideum",  # 12612
                "Drosophila melanogaster",  # 13461
                "Escherichia coli",  # 4370,
                "Glycine max",  # 55796
                "Homo sapiens",  # 23586
                "Methanocaldococcus jannaschii",  # 1773
                "Mus musculus",  # 21452
                "Oryza sativa",  # 43645
                "Rattus norvegicus",  # 22152
                "Saccharomyces cerevisiae",  # 6055
                "Schizosaccharomyces pombe",  # 5196
                "Zea mays",  # 39139
                "Ajellomyces capsulatus",  # 9199
                "Brugia malayi",  # 10972
                "Campylobacter jejuni",  # 1620
                "Cladophialophora carrionii",  # 11170
                "Dracunculus medinensis",  # 10834
                "Fonsecaea pedrosoi",  # 12509
                "Haemophilus influenzae",  # 1660
                "Klebsiella pneumoniae",  # 5727
                "Leishmania infantum",  # 7924,
                "Madurella mycetomatis",  # 9561
                "Mycobacterium leprae",  # 1602
                "Mycobacterium tuberculosis",  # 3991
                "Neisseria gonorrhoeae",  # 2106
                "Nocardia brasiliensis",  # 8398
                "Onchocerca volvulus",  # 12039
                "Paracoccidioides lutzii",  # 8794
                "Plasmodium falciparum",  # 5168
                "Pseudomonas aeruginosa",  # 5555
                "Salmonella typhimurium",  # 4526
                "Schistosoma mansoni",  # 9735
                "Shigella dysenteriae",  # 3893
                "Sporothrix schenckii",  # 8652
                "Streptococcus pneumoniae",  # 2031
                "Strongyloides stercoralis",  # 15335
                "Trichuris trichiura",  # 9563
                "Trypanosoma brucei",  # 8491
                "Trypanosoma cruzi",  # 19036
                "Wuchereria bancrofti",  # 12725
                "Swiss-Prot (CIF files)",  # 550122; This sometimes slows down significantly after only a few GBs (out of 36 GB) are downloaded, so may not finish
                ###
                # REMOVED BY ALPHAFOLD:
                # "Overlap with MANE",  # 17334 (3,844 are unique)  ## REMOVED IN v5
                # "Mycobacterium ulcerans",  # 9033,                ## REMOVED IN v6
                # "Enterococcus faecium",  # 2823,                  ## REMOVED IN v6
            ]
            if self.__fetchAndReorganizeAlphaFoldFtp and len(alphaFoldSpeciesList) > 0:
                for species in alphaFoldSpeciesList:
                    mPWf = ModelProviderWorkflow(
                        srcDir=self.__workPath,
                        destDir=self.__cachePath,
                        modelProviders=["AlphaFold"],
                        useCache=True,
                        numProc=self.__numProc,
                        chunkSize=self.__chunkSize,
                        alphaFoldRequestedSpeciesList=[species],
                    )
                    ok = mPWf.download()
                    self.assertTrue(ok)
                    ok = mPWf.reorganize(keepSource=self.__keepSource)
                    self.assertTrue(ok)

            if self.__fetchAndReorganizeModelArchive:
                mPWf = ModelProviderWorkflow(
                    srcDir=self.__workPath,
                    destDir=self.__cachePath,
                    modelProviders=["ModelArchive"],
                    useCache=True,
                    numProc=self.__numProc,
                    chunkSize=self.__chunkSize,
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
