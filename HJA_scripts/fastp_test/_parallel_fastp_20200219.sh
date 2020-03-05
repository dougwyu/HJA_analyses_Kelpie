#!/bin/bash
set -e
set -u
set -o pipefail
#######################################################################################
#######################################################################################
# a shell script to loop through a set of fastq files and run fastp
#######################################################################################
#######################################################################################
# 
# Usage: bash _parallel_fastp_20200219.sh

# SCRIPT OUTLINE:
     # make list of sample folders that hold the fastq files
     # run trim_galore in each sample folder, using GNU parallel

fastp -i 2-1-plasmid_BDSW202056829-1a_1.fq.gz -I 2-1-plasmid_BDSW202056829-1a_2.fq.gz -o 2-1-plasmid_BDSW202056829-1a_1_out.fq.gz -O 2-1-plasmid_BDSW202056829-1a_2_out.fq.gz

fastp -i 3-1-plasmid_BDSW202056828-1a_1.fq.gz -I 3-1-plasmid_BDSW202056828-1a_2.fq.gz -o 3-1-plasmid_BDSW202056828-1a_1_out.fq.gz -O 2-1-plasmid_BDSW202056828-1a_1_out.fq.gz

# PIPESTART=$(date)
#
# # upload _parallel_trimgalore_20191203.sh *into* ~/_Oregon/
# HOMEFOLDER=$(pwd) # this sets working directory
# TARGETFOLDER="2019Sep_shotgun/1.rawdata"
# echo "Home folder is ${HOMEFOLDER}/${TARGETFOLDER}"
# cd ${HOMEFOLDER}/${TARGETFOLDER} || exit
#
# # read in folder list and make a bash array
# find ${HOMEFOLDER}/${TARGETFOLDER}/ -maxdepth 1 -type d -exec basename {} \; | sort > folderlist.txt # find all folders and then execute basename to remove folder pathnames
#      # alternative command using find
#           # find ${HOMEFOLDER}/${TARGETFOLDER}/ -maxdepth 1 -type d -printf '%f\n' | sort > folderlist.txt  # find all folders and print only the foldername with "-printf '%f\n'"
# # https://unix.stackexchange.com/questions/118841/use-basename-to-parse-a-list-of-paths-held-in-a-file
#      # grep "1.rawdata" folderlist.txt # this is the parent folder, which needs to be removed
#      # wc -l folderlist.txt # 243 lines
# sed -i '/1.rawdata/d' ./folderlist.txt # GNU sed only
#      # grep "1.rawdata" folderlist.txt # should return nothing
#      # wc -l folderlist.txt # 242 lines == 242 samples, in each of which there are 2 fastq files
#
# sample_info=folderlist.txt # put folderlist.txt into variable
# sample_names=($(cut -f 1 "$sample_info" | uniq)) # convert variable to array this way
# echo "${sample_names[@]}" # echo all array elements
# echo "There are" ${#sample_names[@]} "folders that will be processed." # echo number of elements in the array, which should be 242 (the same as the number of lines in folderlist.txt)
#
# # example filenames in each folder
#      # SM-PG-M1-S2_BDSW190603195-1a/
#           # SM-PG-M1-S2_BDSW190603195-1a_1.fq.gz
#           # SM-PG-M1-S2_BDSW190603195-1a_2.fq.gz
# # GNU PARALLEL version, since trimgalore only uses one cpu per job
# parallel --progress "echo Now on species {1} of ${#sample_names[@]}; trim_galore --paired --length 100 -trim-n {1}/*_1.fq.gz {1}/*_2.fq.gz -o {1}" ::: "${sample_names[@]}"
# # -o {1} # to write output fq.gz files to sample folder. output is *_val_1_fq.gz
# # --paired --length 100 # remove a paired read if either is less than 100 bp after trimming (default 20). Reads in this job are 150 PE
#
# echo "Pipeline started at $PIPESTART"
# NOW=$(date)
# echo "Pipeline ended at $NOW"
