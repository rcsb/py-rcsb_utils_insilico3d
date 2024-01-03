#!/bin/bash

# Script to fetch AlphaFold data in bulk from Google Cloud, and store .tar files by three digit taxonomy ID prefix
#
# !!! MAKE SURE YOU HAVE ENOUGH DISK SPACE BEFORE RUNNING THIS (~25 TB) !!!
#
# Downloaded tar balls will contain PAE files too (can be extracted and removed for space savings)


for i in {000..999}
do 
    cd /mnt/vdb1/source-models/work-dir/AlphaFoldCloud/$i
    pwd
    gcloud storage cp -n -r gs://public-datasets-deepmind-alphafold-v4/proteomes/proteome-tax_id-$i* .
done
