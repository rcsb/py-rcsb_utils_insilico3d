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


def main():
    workDir = "/mnt/vdb1/source-models/work-dir/AlphaFoldCloud"
    sizeLimit = 68719476736  # 64 GB - max size of any single archive prefix subset

    cacheD = createArchiveDirFileDict(workDir, sizeLimit)
    createdTime = datetime.datetime.now().isoformat()

    outD = {"data": cacheD, "created": createdTime}
    mU = MarshalUtil()
    outFile = os.path.join(workDir, "model-download-cache.json")
    mU.doExport(outFile, outD, fmt="json", indent=2)


def createArchiveDirFileDict(baseDir, sizeLimit):
    # Get a dictionary of all subdirectories and files in the given directory
    absBaseDir = os.path.abspath(baseDir)
    archiveDirFileD = {}
    with os.scandir(absBaseDir) as fObjs:
        for fObj in fObjs:
            if fObj.is_dir() and str(fObj.name).isdigit() and not str(fObj.name).lower().startswith("v"):
                archiveDir = os.path.join(absBaseDir, fObj.name)
                archiveDirFileD[archiveDir] = {}
                archiveFileD = {}
                tmpSize = 0
                archiveDirSubsetIdx = 0
                with os.scandir(archiveDir) as archiveObjs:
                    for archiveObj in archiveObjs:
                        if archiveObj.is_file() and archiveObj.name.endswith(".tar"):
                            fSize = archiveObj.stat().st_size
                            archiveFileD[archiveObj.name] = fSize
                            tmpSize += fSize
                            if tmpSize > sizeLimit:
                                archiveDirFileD[archiveDir][archiveDirSubsetIdx] = {"archive_files": archiveFileD}
                                archiveDirSubsetIdx += 1
                                archiveFileD = {}
                                tmpSize = 0
                archiveDirFileD[archiveDir][archiveDirSubsetIdx] = {"archive_files": archiveFileD}
    return archiveDirFileD


if __name__ == "__main__":
    main()
