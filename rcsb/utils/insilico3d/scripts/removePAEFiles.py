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

# ADD FUNCTIONALITY TO READ IN THE TAR HOLDINGS FILE AND UPDATE WITH TWO EXTRA KEYS PER TAR:
# - Whether each individual TAR has been cleaned or not
# - The number of CIF files in each TAR


def remove_tar_items_by_pattern(tarball_path, pattern):
    cifCount = 0
    # Create a new tarball without the specified files
    new_tarball_path = tarball_path.replace(".tar", ".tmp.tar")  # Use the same name as the original tarball

    # Process the tarfile and remove the tar items matching the given pattern (".json.gz")
    with tarfile.open(tarball_path, 'r') as tar:
        #
        # First check if tarball has already been cleaned up
        tar_file_list = tar.getnames()
        if all([f.endswith(".cif.gz") for f in tar_file_list]):
            return len(tar_file_list)
        #
        with tarfile.open(new_tarball_path, 'w') as new_tar:
            for member in tar.getmembers():
                # Skip the files matching the pattern
                if not fnmatch.fnmatch(member.name, pattern):
                    # Extract and add the remaining files to the new tarball
                    file_content = tar.extractfile(member)
                    new_tar.addfile(member, fileobj=file_content)
                    cifCount += 1

    os.remove(tarball_path)
    os.rename(new_tarball_path, tarball_path)

    # print(f"New tarball created at {tarball_path}")

    return cifCount


if __name__ == "__main__":
    #
    for i in range(712, 800):
        modelCount = 0
        tar_dir = f"/mnt/vdb1/source-models/work-dir/AlphaFoldCloud/{i}"

        for filename in os.listdir(tar_dir):
            if filename.endswith(".tar"):
                tar_filepath = os.path.join(tar_dir, filename)

                # First check if a temp file already exists (e.g., from a previous run that may have been interrupted)
                if filename.endswith(".tmp.tar"):
                    original_tar_filepath = tar_filepath.replace(".tmp", "")
                    if os.path.exists(tar_filepath):
                        if os.path.exists(original_tar_filepath):
                            # If the original file does exist, then delete the temp file and allow the function to re-process the original file
                            print(f"Deleting temporary file {tar_filepath} since original tar still exists {original_tar_filepath}")
                            os.remove(tar_filepath)
                        else:
                            # Else, if the original tar file doesn't exist, assume the file was already fully processed but that the script was
                            # interrupted before the original file was deleted, so rename the temp file to the original filename
                            print(f"Renaming temporary file {tar_filepath} to original tar file name {original_tar_filepath}")
                            os.rename(tar_filepath, original_tar_filepath)

                else:
                    modelCount += remove_tar_items_by_pattern(tar_filepath, "*.json.gz")

        print(f"Finished directory {tar_dir} - # models {modelCount}")
