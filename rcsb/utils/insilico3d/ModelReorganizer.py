##
# File:    ModelReorganizer.py
# Author:  Dennis Piehl
# Date:    10-Mar-2022
#
# Updates:
#
# To Do:
# - pylint: disable=fixme
# - Add mkdssp calculation
# - Add check that model files are consistent with mmCIF dictionaries
##

"""
Common worker methods for processing and handling computed model files.

"""

__docformat__ = "google en"
__author__ = "Dennis Piehl"
__email__ = "dennis.piehl@rcsb.org"
__license__ = "Apache 2.0"

import logging
import os.path
import copy
from datetime import datetime
import pytz

from mmcif.api.DataCategory import DataCategory
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.multiproc.MultiProcUtil import MultiProcUtil

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger(__name__)


class ModelWorker(object):
    """A skeleton class that implements the interface expected by the multiprocessing
    for working on model files.
    """

    def __init__(self, **kwargs):
        _ = kwargs
        self.__workPath = kwargs.get("workPath", None)
        self.__fU = FileUtil(workPath=self.__workPath)
        self.__mU = MarshalUtil(workPath=self.__workPath)

    def reorganize(self, dataList, procName, optionsD, workingDir):
        """Enumerate and reorganize and rename model files.

        Args:
            dataList (list): list of model files to work on
            procName (str): worker process name
            optionsD (dict): dictionary of additional options that worker can access
            workingDir (str): path to working directory

        Returns:
            successList (list): list of input data items that were successfully processed
            retList (list): list of processed results
            diagList (list): list of unique diagnostics
        """
        _ = workingDir
        _ = optionsD
        successList = []
        failList = []
        retList = []
        diagList = []

        try:
            modelSourcePrefix = optionsD.get("modelSourcePrefix")  # e.g., "AF" or "MA"
            modelSourceDbMap = optionsD.get("modelSourceDbMap")
            destBaseDir = optionsD.get("destBaseDir")  # base path for all computed models (i.e., "computed-models"); Or will be root path at HTTP endpoint
            keepSource = optionsD.get("keepSource", False)  # whether to copy files over (instead of moving them)
            reorganizeDate = optionsD.get("reorganizeDate", False)  # whether to copy files over (instead of moving them)
            #
            for modelFileIn in dataList:
                modelD = {}
                success = False
                modelFileOut = None
                modelFileNameIn = self.__fU.getFileName(modelFileIn)
                modelSourceDb = modelSourceDbMap[modelSourcePrefix]
                #
                containerList = self.__mU.doImport(modelFileIn, fmt="mmcif")
                if len(containerList) > 1:
                    # Expecting all computed models to have one container per file. When this becomes no longer the case, update this to handle it accordingly.
                    logger.error("Skipping - model file %s has more than one container (%d)", modelFileNameIn, len(containerList))
                    continue
                #
                dataContainer = containerList[0]
                #
                # Create internal model ID using entry.id and strip away all punctuation and make ALL CAPS
                sourceModelEntryId = dataContainer.getObj("entry").getValue("id", 0)
                modelEntryId = "".join(char for char in sourceModelEntryId if char.isalnum()).upper()
                internalModelId = modelSourcePrefix + "_" + modelEntryId
                #
                # Get the revision date if it exists
                if dataContainer.exists("pdbx_audit_revision_history"):
                    lastModifiedDate = dataContainer.getObj("pdbx_audit_revision_history").getValue("revision_date", -1)
                    lastModifiedDate = datetime.strptime(lastModifiedDate, '%Y-%m-%d').replace(microsecond=0).replace(tzinfo=pytz.UTC).isoformat()
                else:
                    lastModifiedDate = reorganizeDate
                #
                # Insert default deposited pdbx_assembly information into CIF
                dataContainer = self.__addDepositedAssembly(dataContainer=dataContainer)
                #
                # Gzip the original file if not already (as the case for ModelArchive model files)
                if modelFileIn.endswith(".gz"):
                    modelFileInGzip = modelFileIn
                else:
                    modelFileInGzip = modelFileIn + ".gz"
                    logger.debug("Compressing model file %s --> %s", modelFileIn, modelFileInGzip)
                    ok = self.__fU.compress(modelFileIn, modelFileInGzip)
                    if not ok:
                        logger.warning("Failed to gzip input model file: %s", modelFileIn)
                #
                internalModelName = internalModelId + ".cif.gz"
                #
                # Use last six to last two characters for second-level hashed directory
                firstDir, secondDir = modelEntryId[-6:-4], modelEntryId[-4:-2]
                modelPathFromPrefixDir = os.path.join(modelSourcePrefix, firstDir, secondDir, internalModelName)
                destModelDir = os.path.join(destBaseDir, modelSourcePrefix, firstDir, secondDir)
                if not self.__fU.exists(destModelDir):
                    self.__fU.mkdir(destModelDir)
                modelFileOut = os.path.join(destModelDir, internalModelName)
                #
                sourceModelUrl = self.__getSourceUrl(modelSourcePrefix, modelFileNameIn, sourceModelEntryId)
                #
                modelD["modelId"] = internalModelId
                modelD["modelPath"] = modelPathFromPrefixDir  # Starts at prefix (e.g., "AF/XJ/E6/AF_AFA0A385XJE6F1.cif.gz"); needed like this by RepositoryProvider
                modelD["sourceId"] = sourceModelEntryId
                modelD["sourceDb"] = modelSourceDb
                modelD["sourceModelFileName"] = modelFileNameIn
                modelD["sourceModelUrl"] = sourceModelUrl
                modelD["lastModifiedDate"] = lastModifiedDate
                #
                try:
                    self.__mU.doExport(modelFileOut, containerList, fmt="mmcif")
                    if not keepSource:
                        self.__mU.remove(modelFileInGzip)  # Remove original file
                        if self.__mU.exists(modelFileIn):   # Remove unzipped file too if it exists
                            self.__mU.remove(modelFileIn)
                    success = True
                    successList.append(modelFileIn)
                except Exception as e:
                    logger.debug("Failing to reorganize %s --> %s, with %s", modelFileIn, modelFileOut, str(e))
                #
                retList.append((modelFileIn, modelD, success))
                #
            failList = sorted(set(dataList) - set(successList))
            if failList:
                logger.info("%s returns %d definitions with failures: %r", procName, len(failList), failList)
            #
            logger.debug("%s processed %d/%d models, failures %d", procName, len(retList), len(dataList), len(failList))
        except Exception as e:
            logger.exception("Failing %s for %d data items %s", procName, len(dataList), str(e))

        return successList, retList, diagList

    def __getSourceUrl(self, modelSourcePrefix, sourceModelFileName, sourceModelEntryId):
        """Construct model accession URL for each model source.

        Args:
            modelSourcePrefix (str): Prefix for given model source (e.g., "AF" or "MA").
            sourceModelFileName (str): filename of model as provided by the source.
            sourceModelEntryId (str): model entry.id value as in the provided mmCIF file

        Returns:
            (str): Accession URL for the specified model file and source.
                    Note that file-specific downloads aren't gzipped, unlike model files in species tarball.
                    E.g., "https://alphafold.ebi.ac.uk/files/AF-Q9RQP8-F1-model_v2.cif"
        """
        sourceModelUrl = None
        try:
            if modelSourcePrefix == "AF":
                modelFileNameInUrl = sourceModelFileName.split(".gz")[0]
                sourceModelUrl = os.path.join("https://alphafold.ebi.ac.uk/files/", modelFileNameInUrl)
            elif modelSourcePrefix == "MB":
                modbaseInternalId = sourceModelEntryId.split("model_")[-1]
                # sourceModelUrl = "https://salilab.org/modbase/searchbyid?modelID=" + modbaseInternalId + "&displaymode=moddetail"  # Directs to entry page
                sourceModelUrl = "https://salilab.org/modbase/retrieve/modbase/?modelID=" + modbaseInternalId + "&format=mmcif"  # Directs to mmCIF file displayed in web browser
                # E.g.: https://salilab.org/modbase/retrieve/modbase/?modelID=ecac68b60ee6877ccde36af05cdeac58&format=mmcif
            elif modelSourcePrefix == "MA":
                modelFileNameInUrl = sourceModelFileName.split(".cif")[0]
                sourceModelUrl = "https://www.modelarchive.org/api/projects/" + modelFileNameInUrl + "?type=basic__model_file_name"
            else:
                logger.error("URL generation process not ready yet for %s", modelSourcePrefix)
        except Exception as e:
            logger.exception("Failing with %s", str(e))

        return sourceModelUrl

    def __addDepositedAssembly(self, dataContainer):
        """Add the deposited coordinates as an additional separate assembly labeled as 'deposited'
        to categories, pdbx_struct_assembly and pdb_struct_assembly_gen.

        Method copied from rcsb.utils.dictionary.DictMethodAssemblyHelper.

        Args:
            dataContainer (object): mmcif.api.DataContainer object instance

        Returns:
            bool: True for success or False otherwise

        """
        if not dataContainer.exists("pdbx_struct_assembly"):
            dataContainer.append(
                DataCategory(
                    "pdbx_struct_assembly",
                    attributeNameList=["id", "details", "method_details", "oligomeric_details", "oligomeric_count", "rcsb_details", "rcsb_candidate_assembly"],
                )
            )
        if not dataContainer.exists("pdbx_struct_assembly_gen"):
            dataContainer.append(DataCategory("pdbx_struct_assembly_gen", attributeNameList=["assembly_id", "oper_expression", "asym_id_list", "ordinal"]))

        if not dataContainer.exists("pdbx_struct_oper_list"):
            row = [
                "1",
                "identity operation",
                "1_555",
                "x, y, z",
                "1.0000000000",
                "0.0000000000",
                "0.0000000000",
                "0.0000000000",
                "0.0000000000",
                "1.0000000000",
                "0.0000000000",
                "0.0000000000",
                "0.0000000000",
                "0.0000000000",
                "1.0000000000",
                "0.0000000000",
            ]
            atList = [
                "id",
                "type",
                "name",
                "symmetry_operation",
                "matrix[1][1]",
                "matrix[1][2]",
                "matrix[1][3]",
                "vector[1]",
                "matrix[2][1]",
                "matrix[2][2]",
                "matrix[2][3]",
                "vector[2]",
                "matrix[3][1]",
                "matrix[3][2]",
                "matrix[3][3]",
                "vector[3]",
            ]
            dataContainer.append(DataCategory("pdbx_struct_oper_list", attributeNameList=atList, rowList=[row]))
        #
        logger.debug("Add deposited assembly for %s", dataContainer.getName())
        cObj = dataContainer.getObj("struct_asym")
        asymIdL = cObj.getAttributeValueList("id")
        logger.debug("AsymIdL %r", asymIdL)
        #
        # Ordinal is added by subsequent attribure-level method.
        tObj = dataContainer.getObj("pdbx_struct_assembly_gen")
        rowIdx = tObj.getRowCount()
        tObj.setValue("deposited", "assembly_id", rowIdx)
        tObj.setValue("1", "oper_expression", rowIdx)
        tObj.setValue("1", "ordinal", rowIdx)
        tObj.setValue(",".join(asymIdL), "asym_id_list", rowIdx)
        #
        tObj = dataContainer.getObj("pdbx_struct_assembly")
        rowIdx = tObj.getRowCount()
        tObj.setValue("deposited", "id", rowIdx)
        tObj.setValue("deposited_coordinates", "details", rowIdx)
        tObj.setValue("Y", "rcsb_candidate_assembly", rowIdx)
        #
        for atName in ["oligomeric_details", "method_details", "oligomeric_count"]:
            if tObj.hasAttribute(atName):
                tObj.setValue("?", atName, rowIdx)
        #
        tObj.setValue(str(len(asymIdL)), "oligomeric_count", rowIdx)

        return dataContainer


class ModelReorganizer(object):
    """Generators and accessors for model files."""

    def __init__(self, cachePath=None, useCache=True, **kwargs):
        """Initialize ModelReorganizer object.

        Args:
            cachePath (str): directory path for storing cache file containing a dictionary of all reorganized models.
            useCache (bool): whether to use the existing data cache or re-run entire model reorganization process.
            cacheFile (str, optional): filename for cache file; default "computed-models-holdings.json.gz".
            cacheFilePath (str, optional): full filepath and name for cache file (will override "cachePath" and "cacheFile" if provided).
            numProc (int, optional): number of processes to use; default 2.
            chunkSize (int, optional): incremental chunk size used for distribute work processes; default 20.
            workPath (str, optional): directory path for workers to operate in; default is cachePath.
            keepSource (bool, optional): whether to copy model files to new directory instead of moving them; default False.
        """

        try:
            self.__cachePath = cachePath if cachePath else "."

            self.__numProc = kwargs.get("numProc", 2)
            self.__chunkSize = kwargs.get("chunkSize", 20)
            self.__workPath = kwargs.get("workPath", os.path.join(self.__cachePath, "work-dir"))
            self.__keepSource = kwargs.get("keepSource", False)

            self.__mU = MarshalUtil(workPath=self.__workPath)
            self.__fU = FileUtil(workPath=self.__workPath)

            self.__cacheFormat = kwargs.get("cacheFormat", "json")
            # self.__cacheFormat = kwargs.get("cacheFormat", "pickle")
            cacheExt = "pic" if self.__cacheFormat == "pickle" else "json"
            cacheFile = kwargs.get("cacheFile", "computed-models-holdings." + cacheExt + ".gz")
            cacheFilePath = kwargs.get("cacheFilePath", os.path.join(self.__cachePath, "holdings", cacheFile))
            if not cacheFilePath.lower().endswith(".gz"):
                logger.error("Holdings cache file must be gzipped, %s", cacheFilePath)
                raise ValueError("Error: Holdings cache file must be gzipped.")
            self.__cacheFilePath = cacheFilePath  # self.__cacheFilePath is full path to gzipped cache file
            self.__cacheFilePathUnzip = cacheFilePath[0:-3]

            logger.info("Reorganizing models using cachePath %s, cacheFilePath %s, cacheFilePathUnzip %s", self.__cachePath, self.__cacheFilePath, self.__cacheFilePathUnzip)

            self.__mD = self.__reload(cacheFilePath=self.__cacheFilePath, useCache=useCache)

        except Exception as e:
            logger.exception("Failing with %s", str(e))

    def testCache(self, minCount=20):
        try:
            if minCount == 0:
                return True
            if self.__mD and len(self.__mD) >= minCount:
                logger.info("Reorganized models in cache (%d)", len(self.__mD))
                return True
        except Exception:
            pass
        return False

    def reload(self, cacheFilePath, useCache=True):
        self.__mD = self.__reload(cacheFilePath, useCache=useCache)
        return len(self.__mD) > 0  # Returns True or False, if reload was successful or not

    def __reload(self, cacheFilePath, useCache=True):
        """Reload from the current cache file."""
        try:
            mD = {}
            logger.info("useCache %r cacheFilePath %r", useCache, cacheFilePath)
            if useCache and self.__mU.exists(cacheFilePath):
                if cacheFilePath != self.__cacheFilePath:
                    self.__cacheFilePath = cacheFilePath
                mD = self.__mU.doImport(cacheFilePath, fmt=self.__cacheFormat)
                logger.info("Reorganized models (%d) in cacheFilePath %r", len(mD), cacheFilePath)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        #
        return mD

    def getCachePath(self):
        return self.__cachePath

    def getCacheFilePath(self):
        return self.__cacheFilePath

    def reorganize(self, inputModelList, modelSource, destBaseDir, useCache=True):
        """Move model files from organism-wide model listing to hashed directory structure and rename files
        to follow internal identifier naming convention.

        Args:
            inputModelList (list): List of input model filepaths to reorganize.
            modelSource (str): Source of model files ("AlphaFold", "ModBase", "ModelArchive", or "SwissModel")
            destBaseDir (str): Base destination directory into which to reorganize model files (e.g., "computed-models")

        Returns:
            bool: True for success or False otherwise
        """
        ok = False
        try:
            mD, failD = self.__reorganizeModels(
                inputModelList=inputModelList,
                modelSource=modelSource,
                destBaseDir=destBaseDir,
                numProc=self.__numProc,
                chunkSize=self.__chunkSize
            )
            if len(failD) > 0:
                logger.error("Failed to process %d model files.", len(failD))
            kwargs = {"indent": 4} if self.__cacheFormat == "json" else {"pickleProtocol": 4}
            if useCache:
                self.__mD = self.__reload(cacheFilePath=self.__cacheFilePath, useCache=useCache)
                for modelId, modelD in mD.items():
                    self.__mD.update({modelId: modelD})
            else:
                self.__mD = copy.deepcopy(mD)
            ok = self.__mU.doExport(self.__cacheFilePathUnzip, self.__mD, fmt=self.__cacheFormat, **kwargs)
            logger.info("Wrote %r status %r", self.__cacheFilePathUnzip, ok)
            if ok:
                ok2 = self.__fU.compress(self.__cacheFilePathUnzip, self.__cacheFilePath)
                logger.info("Compressed %r status %r", self.__cacheFilePath, ok2)
                if ok2:
                    ok3 = self.__fU.remove(self.__cacheFilePathUnzip)
                    logger.info("Removing uncompressed holdings cache file %s status %r", self.__cacheFilePathUnzip, ok3)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return ok

    def __reorganizeModels(self, inputModelList, modelSource, destBaseDir, numProc=2, chunkSize=20):
        """Prepare multiprocessor queue and workers for generating input:output map for model files, to use in reorganizing
        them into a structured directory tree and naming them with internal identifiers.

        Args:
            inputModelList (list): List of input model filepaths to reorganize.
            modelSource (str): Source of model files ("AlphaFold", "ModBase", "ModelArchive", or "SwissModel")
            destBaseDir (str): Base destination directory into which to reorganize model files (e.g., "computed-models")
            numProc (int, optional): number of processes to use; default 2.
            chunkSize (int, optional): incremental chunk size used for distribute work processes; default 20.

        Returns:
            mD (dict): dictionary of successfully processed models, in the following structure:
                        { internalModelId_1 : {"sourceModelFileName": ..., "modelPath": ..., "sourceModelUrl": ...},
                          internalModelId_2 : {...}, ...}
                E.g:    { "AF_AFO25670F1" : {
                            "sourceModelFileName": "AF-O25670-F1-model_v2.cif.gz",
                            "modelPath": "./AF/56/70/AF_AFO25670F1.cif.gz",
                            "sourceModelUrl": "https://alphafold.ebi.ac.uk/files/AF-O25670-F1-model_v2.cif"},
                          "AF_AFO24939F1" : {...},
                        ...}
            failD (dict): dictionary of models for which processing failed and the associated output file path and source details
                          that were attempted, in the following structure:
                          {inputModelFilePath:  {"sourceModelFileName": ..., "modelPath": ..., "sourceModelUrl": ...}, ...}
        """
        mD = {}
        failD = {}
        tS = datetime.now().replace(microsecond=0).replace(tzinfo=pytz.UTC).isoformat()  # Desired format:  2022-04-15T12:00:00+00:00
        #
        logger.info("Starting with %d models, numProc %d at %r", len(inputModelList), numProc, tS)
        #
        # Create the base destination directory if it doesn't exist
        if not self.__fU.exists(destBaseDir):
            logger.info("Creating base destination directory for model file reorganization, %s", destBaseDir)
            self.__fU.mkdir(destBaseDir)
        #
        modelSourcePrefixD = {"AlphaFold": "AF", "ModBase": "MB", "ModelArchive": "MA", "SwissModelRepository": "SMR"}
        modelSourceDbMap = {"AF": "AlphaFoldDB", "MB": "ModBase", "MA": "ModelArchive", "SMR": "SwissModelRepository"}
        modelSourcePrefix = modelSourcePrefixD[modelSource]
        #
        rWorker = ModelWorker(workPath=self.__workPath)
        mpu = MultiProcUtil(verbose=True)
        optD = {
            "modelSourcePrefix": modelSourcePrefix,
            "modelSourceDbMap": modelSourceDbMap,
            "destBaseDir": destBaseDir,
            "keepSource": self.__keepSource,
            "reorganizeDate": tS
        }
        mpu.setOptions(optD)
        mpu.set(workerObj=rWorker, workerMethod="reorganize")
        mpu.setWorkingDir(workingDir=self.__workPath)
        ok, failList, resultList, _ = mpu.runMulti(dataList=inputModelList, numProc=numProc, numResults=1, chunkSize=chunkSize)
        if failList:
            logger.info("model file failures (%d): %r", len(failList), failList)
        #
        for (modelFileIn, modelD, success) in resultList[0]:
            modelIdD = copy.deepcopy(modelD)
            if success:
                modelId = modelIdD.pop("modelId")
                mD[modelId] = modelIdD
            else:
                failD[modelFileIn] = modelIdD
        #
        logger.info("Completed with multi-proc status %r, failures %r, total models with data (%d)", ok, len(failList), len(mD))
        return mD, failD
