#!/bin/bash
set -e
set -u
set -o pipefail
#######################################################################################
#######################################################################################
# a shell script to loop through a set of fastq files and run cutadapt
# customised for the BWA11 folder, which are the dilution samples only
#######################################################################################
#######################################################################################

# Usage: bash _parallel_trimgalore_20200220.sh

# SCRIPT OUTLINE:
     # make list of sample folders that hold the fastq files
     # run trim_galore in each sample folder, using GNU parallel

PIPESTART=$(date)

# upload _parallel_trimgalore_20200220.sh *into* ~/_Oregon/2019Sep_shotgun/2.trimmeddata/
HOMEFOLDER=$(pwd) # this sets working directory
TARGETFOLDER="BWA11"
echo "Home folder is ${HOMEFOLDER}/${TARGETFOLDER}"
cd ${HOMEFOLDER}/${TARGETFOLDER} || exit

# read in folder list and make a bash array
find ${HOMEFOLDER}/${TARGETFOLDER}/ -mindepth 1 -maxdepth 1 -type d -exec basename {} \; | sort > folderlist.txt # find all folders and then execute basename to remove folder pathnames. -mindepth prevents BWA11/ from being listed too
     # alternative command using find
          # find ${HOMEFOLDER}/${TARGETFOLDER}/ -maxdepth 1 -type d -printf '%f\n' | sort > folderlist.txt  # find all folders and print only the foldername with "-printf '%f\n'"
# https://unix.stackexchange.com/questions/118841/use-basename-to-parse-a-list-of-paths-held-in-a-file
#      grep "BWA11" folderlist.txt # this is the parent folder, which needs to be removed
#      wc -l folderlist.txt # 33 lines
# sed -i '/BWA11/d' ./folderlist.txt # GNU sed only
#      grep "BWA11" folderlist.txt # should return nothing
#      wc -l folderlist.txt # 33 lines == 242 samples, in each of which there are 2 fastq files

sample_info=folderlist.txt # put folderlist.txt into variable
sample_names=($(cut -f 1 "$sample_info" | uniq)) # convert variable to array this way
echo "${sample_names[@]}" # echo all array elements
echo "There are" ${#sample_names[@]} "folders that will be processed." # echo number of elements in the array, which should be 33 (the same as the number of lines in folderlist.txt)

# example filenames in each folder
     # SM-PG-M1-S2_BDSW190603195-1a/
          # SM-PG-M1-S2_BDSW190603195-1a_1.fq.gz
          # SM-PG-M1-S2_BDSW190603195-1a_2.fq.gz
# GNU PARALLEL version, since trimgalore only uses one cpu per job
parallel --progress "echo Now on species {1} of ${#sample_names[@]}; trim_galore --paired --length 100 -trim-n {1}/*_1.fq.gz {1}/*_2.fq.gz -o {1}" ::: "${sample_names[@]}"
# -o {1} # to write output fq.gz files to sample folder. output is *_val_1_fq.gz
# --paired --length 100 # remove a paired read if either is less than 100 bp after trimming (default 20). Reads in this job are 150 PE

echo "Pipeline started at $PIPESTART"
NOW=$(date)
echo "Pipeline ended at $NOW"
