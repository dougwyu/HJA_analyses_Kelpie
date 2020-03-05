#!/bin/bash
set -e
set -u
set -o pipefail
#######################################################################################
#######################################################################################
# a shell script to remove the original fq.gz files after running trimgalore
#######################################################################################
#######################################################################################

# RUN THIS INTERACTIVELY

# upload _post_trimgalore_cleanup_20191203.sh *into* ~/_Oregon/
cd ~/_Oregon || exit
HOMEFOLDER=$(pwd) # _Oregon/
TARGETFOLDER="2019Sep_shotgun"
echo "Home folder is ${HOMEFOLDER}/${TARGETFOLDER}"
cd ${HOMEFOLDER}/${TARGETFOLDER} || exit

# After running TrimGalore, I have the original and the trimmed fastq files in the same sample folders
# I only want the trimmed fastq files.
# Because I also have the original fastq files in 2019Sep_shotgun_1.rawdata.tar, I do the following:
     # I remove the original files from each of the sample folders in 1.rawdata/
     # I rename the 1.rawdata/ to 2.trimmeddata/

# the find command will eventually be turned into rm.  using find to check the grobbing syntax
find 1.rawdata/*/*a_{1,2}.fq.gz -type f | wc -l # 242
find 1.rawdata/*/*a_{1,2}.fq.gz -type f # filenames look correct

# remove the original, untrimmed fastq files and any trimgalore intermediate files

# rm -f 1.rawdata/*/*a_{1,2}.fq.gz # uncomment when ready
find 1.rawdata/*/*a_{1,2}.fq.gz -type f | wc -l # should return nothing
find 1.rawdata/*/*a_{1,2}_val_{1,2}.fq.gz -type f | wc -l # 242
find 1.rawdata/*/*a_{1,2}_val_{1,2}.fq.gz -type f # filenames look correct
     # SM-08-M1-S1_BDSW190603188-1a_2_val_2.fq.gz

# rename folder
mv 1.rawdata/ 2.trimmeddata/


# END




# old code

# read in folder list and make a bash array
# find ${HOMEFOLDER}/${TARGETFOLDER}/ -maxdepth 1 -type d -exec basename {} \; | sort > folderlist.txt # find all folders and then execute basename to remove folder pathnames
     # alternative command using find
          # find ${HOMEFOLDER}/${TARGETFOLDER}/ -maxdepth 1 -type d -printf '%f\n' | sort > folderlist.txt  # find all folders and print only the foldername with "-printf '%f\n'"
# https://unix.stackexchange.com/questions/118841/use-basename-to-parse-a-list-of-paths-held-in-a-file
     # grep "1.rawdata" folderlist.txt # this is the parent folder, which needs to be removed
     # wc -l folderlist.txt # 243 lines
# sed -i '/1.rawdata/d' ./folderlist.txt # GNU sed only
     # grep "1.rawdata" folderlist.txt # should return nothing
     # wc -l folderlist.txt # 242 lines == 242 samples, in each of which there are 2 fastq files


# 6223.00 GB at start
# du -sh ~/_ArcDyn # 2.4T afterwards
# remove non-trimmed fastq.gz files in PlatesA2B2/PlatesA2B2_combined/
# cd ~/_ArcDyn/PlatesA2B2/
# ls PlatesA2B2_combined/Sample*/Sample*fastq.gz  # should list only fastq.gz files
# ls PlatesA2B2_combined/Sample*/Sample*fq.gz # should list lots of fq.gz files
# ls PlatesA2B2_combined/Sample*/Sample*fastq.gz | wc -l  # count of fastq.gz files
# ls PlatesA2B2_combined/Sample*/Sample*fq.gz | wc -l # should have same count value
# # rm -f PlatesA2B2_combined/Sample*/Sample*fastq.gz  # keep commented out until ready to run this command
# ls PlatesA2B2_combined/Sample*/Sample*fastq.gz  # should list nothing
# ls PlatesA2B2_combined/Sample*/Sample*fq.gz # should list lots of fq.gz files
# ls PlatesA2B2_combined/Sample*/Sample*fq.gz | wc -l # should have same count value
