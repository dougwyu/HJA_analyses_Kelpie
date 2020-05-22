#!/bin/bash
set -e
set -u
set -o pipefail
##################################################################################################
##################################################################################################
# a shell script to loop through a set of fastq files and give a fastqc result for each fastq file
##################################################################################################
##################################################################################################

# Usage: bash _parallel_fastQC_20191205.sh

# upload _parallel_fastQC_20191205.bsub and _parallel_fastQC_20191205.sh *into* ~/_Oregon/

PIPESTART=$(date)
HOMEFOLDER=$(pwd) # this sets working directory, should be ~/_Oregon/
TARGETFOLDER="2019Sep_shotgun/2.trimmeddata"
echo "Home folder is ${HOMEFOLDER}/${TARGETFOLDER}" # ~/_Oregon/2019Sep_shotgun/2.trimmeddata
cd ${HOMEFOLDER}/${TARGETFOLDER} || exit

# read in folder list and make a bash array
find ${HOMEFOLDER}/${TARGETFOLDER}/ -maxdepth 1 -type d -exec basename {} \; | sort > folderlist.txt # find all folders and then execute basename to remove folder pathnames
     # wc -l folderlist.txt # 243
sed -i '/2.trimmeddata/d' ./folderlist.txt # GNU sed only
     # grep "2.trimmeddata" folderlist.txt # should return nothing because it has been deleted
     # grep "checkSize.xls" folderlist.txt # should return nothing because it's a file
     # wc -l folderlist.txt # 242 lines == 242 samples, in each of which there are 2 fastq files

# make array of folder names
sample_info=folderlist.txt # put folderlist.txt into variable
sample_names=($(cut -f 1 "$sample_info" | uniq)) # convert variable to array this way
echo "${sample_names[@]}" # echo all array elements
echo "There are" ${#sample_names[@]} "folders that will be processed." # echo number of elements in the array, which should be 242 (the same as the number of lines in folderlist.txt)

# run fastqc in parallel (24 at a time if using short-eth)
parallel "echo Now on species {1} of ${#sample_names[@]}; fastqc {1}/*.fq.gz" ::: "${sample_names[@]}"
# running on trimmed versions, so the filenames end in fq.gz
