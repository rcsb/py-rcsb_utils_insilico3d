##
# File:    AlphaFoldModelMetadataProvider.py
# Author:  Dennis Piehl
# Date:    23-Aug-2021
#
# Update:
#
#
##

"""
Generators and accessors for Alpha Fold mmCIF model metadata.

"""

import logging
import os.path
import time
import datetime

from rcsb.utils.insilico3d import __version__
from rcsb.utils.insilico3d.AlphaFoldModelProvider import AlphaFoldModelProvider
# from rcsb.utils.dictionary.DictMethodCommonUtils import DictMethodCommonUtils, LigandTargetInstance
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.io.StashableBase import StashableBase
from rcsb.utils.multiproc.MultiProcUtil import MultiProcUtil

logger = logging.getLogger(__name__)


class AlphaFoldModelMetadataWorker(object):
    """A skeleton class that implements the interface expected by the multiprocessing
    for extracting metadata from mmCIF model files.
    """

    def __init__(self, **kwargs):
        self.__aFMP = AlphaFoldModelProvider()  # Is this the proper place/way to initialize this class? Or initialize in AlphaFoldModelMetadataProvider class below only?
        self.__dirPath = self.__aFMP.__dirPath
        _ = kwargs
        # self.__commonU = DictMethodCommonUtils()

        self.__mU = MarshalUtil(workPath=self.__dirPath)
        # self.__mU = MarshalUtil()  # Should workPath be defined up here, or down below for each species directory?

        # Need to define and/or create metadata cache files for each species below...or elsewhere...
        # Can loop over each species directory by retrieving each species model lists incrementally in '__extractAlphaFoldMetadata()' method below

    def extractMmCifMetadata(self, dataList, procName, optionsD, workingDir):
        """Enumerate mmCIF model files for extracting metadata.

        Args:
            dataList (list): list of mmCIF model files to extract metadata from
            procName (str): worker process name (?)
            optionsD (dict): dictionary of additional options that worker can access
            workingDir (str): path to working directory

        Returns:
            successList (list): list of input data items that were successfully processed
            retList (list): list of processed results
            diagList (list): list of unique diagnostics
        """
        _ = workingDir
        successList = []
        failList = []
        retList = []
        diagList = []

        try:
            for modelFile in dataList:
                # Note: no need to unzip files sice MarshalUtil willl unzip them automatically if necessary
                rD = self.__getMmCifMetadata(procName, modelFile)
                retList.append(rD)

            successList = sorted(set(dataList) - set(failList))
            if failList:
                logger.info("%s returns %d definitions with failures: %r", procName, len(failList), failList)

            logger.debug("%s built target interactions for %d/%d entries failures %d", procName, len(retList), len(dataList), len(failList))
        except Exception as e:
            logger.exception("Failing %s for %d data items %s", procName, len(dataList), str(e))

        return successList, retList, diagList

    def __getMmCifMetadata(self, procName, modelFile):
        """Internal method return a dictionary of mmCIF model metadata.

        Args:
            modelFile (str): path to mmCIF model file to extract metadata from (as either .cif or .cif.gz; but, avoid including the same model twice in both formats)
            procName (str): worker process name (?)

        Returns:
            (dict): dictionary of extracted metadata for the provided model file, in the following structure:
                    {dataContainerEntryName: [{categoryName: [{attributeName: ...}], ...}], ...}

        """
        # Add: Get uniprotID
        # Add: Remove new-line characters from outputted long string values (e.g., 1-letter AA codes)
        # Add: pre-defined lists of specific metadata tokens to extract for each data type (i.e., separate lists for integers, floats, and strings)
        #   - float_metadata_tokens_to_extract = ['_ma_qa_metric_local.metric_value']
        #   - string_metadata_tokens_to_extract = ['_entry.id', '_entity_poly.pdbx_seq_one_letter_code', '_entity_poly.pdbx_seq_one_letter_code_can']
        # Add: separate method to return tokens and values for common tokens between all AF models

        rD = {}
        catNameList = [
            "entry", "af_target_ref_db_details", "ma_target_ref_db_details", "entity", "entity_poly", "ma_qa_metric_global",
            "ma_qa_metric_local", "pdbx_audit_revision_details", "pdbx_audit_revision_history"]

        try:
            modelObj = self.__mU.doImport(modelFile, fmt="mmcif")
            for dataContainer in modelObj:
                eName = dataContainer.getName()
                for catName in catNameList:
                    if not dataContainer.exists(catName):
                        continue
                    dObj = dataContainer.getObj(catName)
                    for ii in range(dObj.getRowCount()):
                        dD = dObj.getRowAttributeDict(ii)
                        rD.setdefault(eName, {}).setdefault(catName, []).append(dD)

        except Exception as e:
            logger.exception("%s failing with %s", procName, str(e))

        # Should export metadata cache in same json format as used to index data to DB (and pickled if necessary)
        return rD


class AlphaFoldModelMetadataProvider(StashableBase):
    """Generators and accessors for mmCIF model metadata extraction."""

    def __init__(self, cachePath=None, useCache=False, speciesModelDir=None, **kwargs):
        # Set default useCache=True ...?

        # modelGroup: name of the species directory of model files to process

        # TODO: Determine proper way and place in code to initialize AlphaFoldModelProvider() class
        # TODO: Figure out declaration of cachePath and dirPath
        # TODO: Also determine necessity of using super(class) below
        # TODO: Determine appropriate default setting to use for "useCache"

        try:
            self.__aFMP = AlphaFoldModelProvider()  # Is this the proper place/way to initialize this class?

            self.__version = __version__
            # self.__fileLimit = kwargs.get("fileLimit", None)
            self.__numProc = kwargs.get("numProc", 2)
            self.__chunkSize = kwargs.get("chunkSize", 10)

            self.__speciesModelDir = speciesModelDir
            self.__speciesName = self.__speciesModelDir.split("/")[-1]

            if cachePath:
                self.__cachePath = cachePath
            else:
                self.__cachePath = self.__aFMP.__dirPath

            self.__dirPath = self.__speciesModelDir

            super(AlphaFoldModelMetadataProvider, self).__init__(self.__cachePath, [self.__speciesName])

            # This was commented out in previous exmaple usage...if deemed OK to set default to True in method(parameter) declaration above, then keep this here too
            useCache = kwargs.get("useCache", True)

            #  - Configuration for stash services -
            #    Local target directory name to be stashed.  (subdir of dirPath)

            self.__mU = MarshalUtil(workPath=self.__dirPath)

            self.__metadataD = self.__reload(fmt="pickle", useCache=useCache)

        except Exception as e:
            # May occur if 'speciesModelDir' is not specified in instantiation
            logger.exception("Failing with %s", str(e))

    def testCache(self, minCount=1):
        try:
            if minCount == 0:
                return True
            if self.__metadataD and minCount and "entries" in self.__metadataD and len(self.__metadataD["entries"]) >= minCount:
                logger.info("Model meatadata for (%d) entries created %r version %r", len(self.__metadataD["entries"]), self.__metadataD["created"], self.__metadataD["version"])
                return True
        except Exception:
            pass
        return False

    def generate(self, updateOnly=False, fmt="pickle", indent=0):
        """Generate and export metadata for all AlphaFold mmCIF model files.

        Args:
            updateOnly (bool): only extract and update metadata for new entries.  Defaults to False.
            fmt (str, optional): export file format. Defaults to "pickle".
            indent (int, optional): json format indent. Defaults to 0.

        Returns:
            bool: True for success or False otherwise
        """
        ok = False
        try:
            tS = time.strftime("%Y %m %d %H:%M:%S", time.localtime())
            # dtS = datetime.datetime.now().isoformat()
            tD = self.__extractAlphaFoldMetadata(numProc=self.__numProc, chunkSize=self.__chunkSize, updateOnly=updateOnly)
            self.__metadataD = {"version": self.__version, "created": tS, "entries": tD, "species": self.__speciesName}
            kwargs = {"indent": indent} if fmt == "json" else {"pickleProtocol": 4}
            targetFilePath = self.__getTargetFilePath(fmt=fmt)
            ok = self.__mU.doExport(targetFilePath, self.__metadataD, fmt=fmt, **kwargs)
            logger.info("Wrote %r status %r", targetFilePath, ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return ok

    def reload(self, fmt="pickle"):
        self.__metadataD = self.__reload(fmt=fmt, useCache=True)
        return self.__metadataD is not None

    def __reload(self, fmt="pickle", useCache=True):
        """Reload from the current cache file."""
        try:
            targetFilePath = self.__getTargetFilePath(fmt=fmt)
            tS = time.strftime("%Y %m %d %H:%M:%S", time.localtime())
            metadataD = {"version": self.__version, "created": tS, "entries": {}}
            logger.debug("useCache %r targetFilePath %r", useCache, targetFilePath)
            #
            if useCache and self.__mU.exists(targetFilePath):
                metadataD = self.__mU.doImport(targetFilePath, fmt=fmt)
                # if fmt != "pickle":
                #     for _, nD in metadataD["entries"].items():
                #         nD["nearestNeighbors"] = [LigandTargetInstance(*neighbor) for neighbor in nD["nearestNeighbors"]]
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        #
        return metadataD

    def __getTargetFilePath(self, fmt="pickle"):
        ext = "pic" if fmt == "pickle" else "json"
        pth = os.path.join(self.__dirPath, self.__speciesName + "-model-metadata." + ext)
        return pth

    def __extractAlphaFoldMetadata(self, numProc=2, chunkSize=10, updateOnly=False):
        """Prepare multiprocessor queue and workers for extracting metadata for all AlphaFold mmCIF model files.
        Called by generate() method -> __extractAlphaFoldMetadata() -> Specifies use of workerMethod="extractMmCifMetadata()"

        Args:
            numProc (int, optional): number of processes to use. Defaults to 2.
            chunkSize (int, optional): incremental chunk size used for distribute work processes. Defaults to 10.
            updateOnly (bool): only extract and update metadata for new entries.  Defaults to False.

        Returns:
            (dict): dictionary of extracted metadata for the provided model file, in the following structure:
                    {dataContainerEntryName: [{categoryName: [{attributeName: ...}], ...}], ...}
        """
        updateDate = datetime.datetime.now().isoformat()
        rD = {}
        exD = {}
        #
        # updateOnly - will reuse any existing data loaded when this is instantiated
        #              otherwise the cache context is cleared before the calculation.
        if updateOnly:
            exD = {k: True for k in self.getEntries()}
            logger.info("Reusing (%d) entries", len(exD))
            rD = self.__metadataD["entries"] if "entries" in self.__metadataD else {}
        #
        modelFileList = self.__aFMP.getModelFileList(inputPathList=[self.__speciesModelDir])

        logger.info("Starting with %d entries numProc %d updateOnly (%r)", len(modelFileList), self.__numProc, updateOnly)
        #
        rWorker = AlphaFoldModelMetadataWorker()
        mpu = MultiProcUtil(verbose=True)
        optD = {"updateDate": updateDate}
        mpu.setOptions(optD)
        mpu.set(workerObj=rWorker, workerMethod="extractMmCifMetadata")
        ok, failList, resultList, _ = mpu.runMulti(dataList=modelFileList, numProc=numProc, numResults=1, chunkSize=chunkSize)
        if failList:
            logger.info("mmCIF metadata extraction failures (%d): %r", len(failList), failList)
        #
        for (entryId, nD) in resultList[0]:
            rD[entryId] = nD
        #
        logger.info("Completed with multi-proc status %r failures %r total entries with data (%d)", ok, len(failList), len(rD))
        return rD

    def convert(self, fmt1="json", fmt2="pickle"):
        #
        targetFilePath = self.__getTargetFilePath(fmt=fmt1)
        self.__neighborD = self.__mU.doImport(targetFilePath, fmt=fmt1)
        #
        targetFilePath = self.__getTargetFilePath(fmt=fmt2)
        ok = self.__mU.doExport(targetFilePath, self.__neighborD, fmt=fmt2, pickleProtocol=4)
        return ok
