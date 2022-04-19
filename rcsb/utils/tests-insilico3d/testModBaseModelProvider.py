##
# File:    testModBaseModelProvider.py
# Author:  Dennis Piehl
# Date:    27-Sep-2021
#
# Updates:
#
#
##
"""
Tests for accessor utilities for ModBase 3D Model (bulk download) data.

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
import glob

from rcsb.utils.insilico3d.ModBaseModelProvider import ModBaseModelProvider
from rcsb.utils.config.ConfigUtil import ConfigUtil

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class ModBaseModelProviderTests(unittest.TestCase):
    def setUp(self):
        self.__dataPath = os.path.join(HERE, "test-data")
        self.__cachePath = os.path.join(HERE, "test-output", "CACHE")
        self.__startTime = time.time()
        mockTopPath = os.path.join(TOPDIR, "rcsb", "mock-data")
        configPath = os.path.join(mockTopPath, "config", "dbload-setup-example.yml")
        self.__configName = "site_info_configuration"
        self.__cfgOb = ConfigUtil(configPath=configPath, defaultSectionName=self.__configName, mockTopPath=self.__dataPath)
        self.__pdbxRepoPath = self.__cfgOb.getPath("PDBX_REPO_PATH", sectionName=self.__configName)
        self.__cacheSandboxPath = os.path.join(self.__cachePath, "MOCK_COMP_MODEL_REPO")
        logger.info("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        unitS = "MB" if platform.system() == "Darwin" else "GB"
        rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum resident memory size %.4f %s", rusageMax / 10 ** 6, unitS)
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testModBaseModelProvider(self):
        redownloadBulkData = True
        deleteCache = False
        #
        try:
            # First test fetching model archive
            if redownloadBulkData:
                mProv = ModBaseModelProvider(
                    cachePath=self.__cachePath,
                    useCache=False,
                    cfgOb=self.__cfgOb,
                    configName=self.__configName,
                    numProc=4,
                    chunkSize=20,
                    modBaseServerSpeciesDataPathDict={"Staphylococcus aureus": "S_aureus/2008/staph_aureus.tar"}
                )
                ok = mProv.testCache()
                self.assertTrue(ok)
            #
            # Next test reloading the cache
            mProv = ModBaseModelProvider(
                cachePath=self.__cachePath,
                useCache=True,
                cfgOb=self.__cfgOb,
                configName=self.__configName,
                numProc=4,
                chunkSize=20,
                modBaseServerSpeciesDataPathDict={"Staphylococcus aureus": "S_aureus/2008/staph_aureus.tar"}
            )
            speciesDirList = mProv.getSpeciesDirList()
            ok = True if len(speciesDirList) > 0 else False
            self.assertTrue(ok)
            #
            speciesPdbModelFileList = mProv.getSpeciesPdbModelFileList(speciesDataDir=speciesDirList[0])
            ok = True if len(speciesPdbModelFileList) > 0 else False
            self.assertTrue(ok)
            #
            # Test processing/converting the model & alignment files to mmCIF
            speciesModelDir = os.path.join(self.__dataPath, "ModBase", "Panicum_virgatum")
            speciesPdbModelFileList = [os.path.abspath(f) for f in glob.glob(os.path.join(speciesModelDir, "model", "*.pdb.xz"))]
            mProc = mProv.getModelProcessor(
                cachePath=os.path.join(self.__cachePath, "ModBase"),
                useCache=False,
                speciesModelDir=speciesModelDir,
                speciesName="Panicum virgatum",
                speciesPdbModelFileList=speciesPdbModelFileList,
                pdbxRepoPath=self.__pdbxRepoPath,  # Path to: /pub/pdb/data/structures/divided/mmCIF
                workPath=os.path.join(self.__cachePath, "ModBase", "Panicum_virgatum"),
                cacheFormat="json",
                numProc=2,
            )
            ok = mProc.generate(updateOnly=False)
            self.assertTrue(ok)
            ok = mProc.generate(updateOnly=True)
            self.assertTrue(ok)
            ok = mProc.reload()
            self.assertTrue(ok)
            ok = mProc.testCache(minCount=5)
            self.assertTrue(ok)
            #
            # Test reorganizing the files
            modelD = mProc.getModelD()
            speciesModelFileList = [v["model"] for k, v in modelD["modelsCif"].items()]
            self.assertTrue(len(speciesModelFileList) > 3)
            mBMR = mProv.getModelReorganizer(cachePath=self.__cachePath, useCache=False, numProc=4, chunkSize=20, keepSource=True)
            # destBaseDir = mProv.getComputedModelsDataPath()
            ok = mBMR.reorganize(inputModelList=speciesModelFileList, modelSource="ModBase", destBaseDir=self.__cacheSandboxPath, useCache=False)
            self.assertTrue(ok)
            #
            # Last test deleting the cache
            if deleteCache:
                mProv = ModBaseModelProvider(
                    cachePath=self.__cachePath,
                    useCache=True,
                    modBaseServerSpeciesDataPathDict={"Staphylococcus aureus": "S_aureus/2008/staph_aureus.tar"}
                )
                speciesNameList = mProv.getSpeciesNameList()
                for species in speciesNameList:
                    ok = mProv.removeSpeciesDataDir(speciesName=species)
                    self.assertTrue(ok)
            # ok = mProv.removePdbModelDir(speciesDataDir=speciesConversionDict["speciesModelDir"])
            # self.assertTrue(ok)
            # ok = mProv.removeAlignmentDir(speciesDataDir=speciesConversionDict["speciesModelDir"])
            # self.assertTrue(ok)
            #
        #
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()


def fetchModBaseModels():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(ModBaseModelProviderTests("testModBaseModelProvider"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = fetchModBaseModels()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
