##
# File:    AlphaFoldModelWorkflow.py
# Author:  Dennis Piehl
# Date:    23-Aug-2021
#
# Update:
#
#
##
"""
Worflow for generating and storing AlphaFold mmCIF model metadata.

"""

__docformat__ = "google en"
__author__ = "Dennis Piehl"
__email__ = "dennis.piehl@rcsb.org"
__license__ = "Apache 2.0"

import logging
# import os

from rcsb.utils.insilico3d.AlphaFoldModelProvider import AlphaFoldModelProvider
from rcsb.utils.insilico3d.AlphaFoldModelMetadataProvider import AlphaFoldModelMetadataProvider
# from rcsb.utils.config.ConfigUtil import ConfigUtil

# HERE = os.path.abspath(os.path.dirname(__file__))
# TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class AlphaFoldModelWorkflow(object):
    def __init__(self, **kwargs):
        # - edit as needed -
        # self.__configName = kwargs.get("configName", "site_info_remote_configuration")
        # self.__configPath = kwargs.get("configPath", os.path.join(HERE, "exdb-config-example.yml"))
        # self.__cachePath = kwargs.get("cachePath", os.path.join(HERE, "CACHE"))
        # self.__mockTopPath = kwargs.get("mockTopPath", None)
        self.__numProc = kwargs.get("numProc", 10)
        self.__chunkSize = kwargs.get("chunkSize", 10)
        # self.__modelGroup = kwargs.get("modelGroup", None)
        self.__speciesModelDir = kwargs.get("speciesModelDir", None)
        self.__useCache = kwargs.get("useCache", False)  # Set default to True?
        #
        # self.__cfgOb = ConfigUtil(configPath=self.__configPath, defaultSectionName=self.__configName, mockTopPath=self.__mockTopPath)
        # logger.info("Configuration file path %s", self.__configPath)
        self.__tiP = AlphaFoldModelMetadataProvider(
            # self.__cachePath,
            useCache=self.__useCache, speciesModelDir=self.__speciesModelDir, numProc=self.__numProc, chunkSize=self.__chunkSize
        )

    def update(self, incremental=True):
        ok = self.__tiP.generate(updateOnly=incremental, fmt="pickle")
        return ok

    # def backup(self):
    #     ok = self.__tiP.backup(self.__cfgOb, self.__configName, remotePrefix=None, useStash=True, useGit=True)
    #     return ok

    def restore(self, minCount=0):
        # ok1 = self.__tiP.restore(self.__cfgOb, self.__configName, remotePrefix=None, useStash=True, useGit=True)
        ok2 = self.__tiP.reload()
        ok3 = self.__tiP.testCache(minCount=minCount)
        return ok2 and ok3

    def convert(self, fmt1="json", fmt2="pickle"):
        ok = self.__tiP.convert(fmt1=fmt1, fmt2=fmt2)
        return ok


if __name__ == "__main__":
    aFMP = AlphaFoldModelProvider()
    speciesModelDirs = aFMP.getSpeciesDirList()
    for speciesModelDir in speciesModelDirs:
        tiWf = AlphaFoldModelWorkflow(speciesModelDir=speciesModelDir)
        tiWf.update(incremental=False)
