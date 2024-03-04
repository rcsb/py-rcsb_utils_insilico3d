##
# File   removePAEFiles.py
# Date:  23-Jan-2024
#
#  Iterate over all the raw proteome tar files and delete the associated PAE and pLDDT files,
#  to leave just the .cif.gz model files (for space savings, as PAE files are not yet used).
#
##
import tarfile
import os
import fnmatch

# TO ADD:
# - Multi-processing support
# - Functionality to read in the tar holdings file and update with two extra keys per tar:
#   - Whether each individual TAR has been cleaned or not
#   - The number of CIF files in each TAR


def main():
    for i in range(100, 999):  # Should add multi-processing support (until then, can just run multiple copies of script for different ranges)
        modelCount = 0
        tarDir = f"/mnt/vdb1/source-models/work-dir/AlphaFoldCloud/{i}"

        for filename in os.listdir(tarDir):
            if filename.endswith(".tar"):
                tarFilePath = os.path.join(tarDir, filename)

                # First check if a temp file already exists (e.g., from a previous run that may have been interrupted)
                if filename.endswith(".tmp.tar"):
                    originalTarFilePath = tarFilePath.replace(".tmp", "")
                    if os.path.exists(tarFilePath):
                        if os.path.exists(originalTarFilePath):
                            # If the original file does exist, then delete the temp file and allow the function to re-process the original file
                            print(f"Deleting temporary file {tarFilePath} since original tar still exists {originalTarFilePath}")
                            os.remove(tarFilePath)
                        else:
                            # Else, if the original tar file doesn't exist, assume the file was already fully processed but that the script was
                            # interrupted before the original file was deleted, so rename the temp file to the original filename
                            print(f"Renaming temporary file {tarFilePath} to original tar file name {originalTarFilePath}")
                            os.rename(tarFilePath, originalTarFilePath)

                else:
                    modelCount += removeTarItemsByPattern(tarFilePath, "*.json.gz")

        print(f"Finished directory {tarDir} - # models {modelCount}")


def removeTarItemsByPattern(tarballPath, pattern):
    cifCount = 0
    # Create a new tarball without the specified files
    newTarballPath = tarballPath.replace(".tar", ".tmp.tar")  # Use the same name as the original tarball

    # Process the tarfile and remove the tar items matching the given pattern (".json.gz")
    with tarfile.open(tarballPath, 'r') as tar:
        #
        # # First check if tarball has already been cleaned up
        # tar_file_list = tar.getnames()
        # if all([f.endswith(".cif.gz") for f in tar_file_list]):
        #     return len(tar_file_list)
        #
        with tarfile.open(newTarballPath, 'w') as newTar:
            for member in tar.getmembers():
                # Skip the files matching the pattern
                if not fnmatch.fnmatch(member.name, pattern):
                    # Extract and add the remaining files to the new tarball
                    fileContent = tar.extractfile(member)
                    newTar.addfile(member, fileobj=fileContent)
                    cifCount += 1

    os.remove(tarballPath)
    os.rename(newTarballPath, tarballPath)

    # print(f"New tarball created at {tarballPath}")

    return cifCount


if __name__ == "__main__":
    main()
