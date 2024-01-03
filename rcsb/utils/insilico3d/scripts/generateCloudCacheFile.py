##
# File   generateCloudCacheFile.py
# Date:  3-Jan-2024
#
#  Generate a cache holdings file containing a list of all the tar archive files downloaded from Google Cloud
#
#  This should be run when finished manually downloading the tarballs (using the `gcloud` tool via `fetchFromGcloud.sh`)
#
#  This cache file must be generated prior to running the reorganization workflow
#
##
import os
import datetime
from rcsb.utils.io.MarshalUtil import MarshalUtil


workDir = "/mnt/vdb1/source-models/work-dir/AlphaFoldCloud"


def createArchiveDirFileDict(baseDir):
    # Get a dictionary of all subdirectories and files in the given directory
    archiveDirFileD = {}
    absBaseDir = os.path.abspath(baseDir)
    with os.scandir(absBaseDir) as fObjs:
        for fObj in fObjs:
            if fObj.is_dir() and not str(fObj.name).lower().startswith("v"):
                archiveDir = os.path.join(absBaseDir, fObj.name)
                archiveFileD = {}
                with os.scandir(archiveDir) as archiveObjs:
                    for archiveObj in archiveObjs:
                        if archiveObj.is_file() and archiveObj.name.endswith(".tar"):
                            archiveFileD[archiveObj.name] = archiveObj.stat().st_size
                archiveDirFileD[archiveDir] = {"archive_files": archiveFileD}
    return archiveDirFileD


cacheD = createArchiveDirFileDict(workDir)
createdTime = datetime.datetime.now().isoformat()

outD = {"data": cacheD, "created": createdTime}
mU = MarshalUtil()
outFile = os.path.join(workDir, "model-download-cache.json")
mU.doExport(outFile, outD, fmt="json")
