#!/bin/bash
set -e
set -u
set -o pipefail
##################################################################################################
##################################################################################################
# a shell script to loop through a set of fastq files and run kelpie
##################################################################################################
##################################################################################################

# Usage: bash _parallel_kelpie_20191216.sh
# Upload _parallel_kelpie_20191216.bsub and _parallel_kelpie_20191216.sh *into* ~/_Oregon/

# go to correct folder
# REMOVE THIS CHUNK when running in 2.trimmeddata/
# HOMEFOLDER="/gpfs/home/b042/_Oregon/2019Sep_shotgun/" # should be ~/_Oregon/2019Sep_shotgun/
# TARGETFOLDER="testkelpie" #
# echo "Home folder is ${HOMEFOLDER}${TARGETFOLDER}" # should be ~/_Oregon/2019Sep_shotgun/testkelpie
# cd ${HOMEFOLDER}${TARGETFOLDER} || exit

# use _launch_kelpie.sh to move the *.sh and *.bsub scripts into each BWA folder, containing 24 or 25 samples
# read in folder list and make a bash array
find ./ -maxdepth 1 -mindepth 1 -type d -exec basename {} \; | sort > folderlist.txt
# find ${HOMEFOLDER}/${TARGETFOLDER}/ -maxdepth 1 -mindepth 1 -type d -exec basename {} \; | sort > folderlist.txt # find all folders 1 levels down only, and then execute basename to remove folder pathnames
# using -mindepth 1 prevents parent folders from being included in output
     # wc -l folderlist.txt
     # cat folderlist.txt
# sed -i '/testkelpie/d' ./folderlist.txt # GNU sed only
#      grep "testkelpie" folderlist.txt # should return nothing because it has been deleted
#      wc -l folderlist.txt # 5, in each of which there are 2 fastq files
#      cat folderlist.txt

# make array of folder names
sample_info=folderlist.txt # put folderlist.txt into variable
sample_names=($(cut -f 1 "$sample_info" | uniq)) # convert variable to array this way
echo "${sample_names[@]}" # echo all array elements
echo "There are" ${#sample_names[@]} "folders that will be processed." # echo number of elements in the array, which should be 5 (the same as the number of lines in folderlist.txt)

# run Kelpie in parallel. -j 2 means 2 samples, -k keep same order as in array, --dryrun see the commands
     # on short-ib queue
     # parallel -j 1 # 270 mins per 20 samples = 13.5 mins per sample
     # parallel -j 2 # 230 mins per 20 samples = 11.5 mins per sample
     # parallel -j 3 # 210 mins per 20 samples = 10.5 mins per sample
     # parallel -j 5 # 200 mins per 20 samples = 10 mins per sample
     # running 10 simultaneous jobs, each on 24 or 25 samples.  If using -j 2, can do each job in 4.5 hrs, so the whole thing should take less than 7 hrs, if some jobs are delayed in starting
# version with gzip -d
parallel -k -j 5 "gzip -d < {1}/{1}_1_val_1.fq.gz > {1}/{1}_1_val_1.fq; gzip -d < {1}/{1}_2_val_2.fq.gz > {1}/{1}_2_val_2.fq; singularity exec --bind /gpfs/home/b042/scripts:/scripts --bind /gpfs/scratch/b042:/mnt ~/ubuntukelpie.sif /scripts/Kelpie_v2.0.4/ubuntu-16.04/Kelpie_v2 -f CCHGAYATRGCHTTYCCHCG -r TCDGGRTGNCCRAARAAYCA -unfiltered -strict -min 400 -max 440 {1}/{1}_?_val_?.fq {1}/{1}.fas; rm {1}/{1}_1_val_1.fq; rm {1}/{1}_2_val_2.fq" ::: "${sample_names[@]}"

# tested with unpigz -k # -k keep original file, but it was even slower than gunzip
# tested -tmp /mnt/{1} and -tmp /mnt # neither works when running as parallel

# old command using the kelpie binary inside the singularity container and with tmp files written to working directory instead of the scratch disk
# parallel -k -j 16 "gzip -d < {1}/{1}_1_val_1.fq.gz > {1}/{1}_1_val_1.fq; gzip -d < {1}/{1}_2_val_2.fq.gz > {1}/{1}_2_val_2.fq; singularity exec ~/ubuntukelpie.sif Kelpie_v2 -f CCHGAYATRGCHTTYCCHCG -r TCDGGRTGNCCRAARAAYCA -unfiltered -strict -min 400 -max 440 {1}/{1}_?_val_?.fq {1}/{1}.fas; rm {1}/{1}_1_val_1.fq; rm {1}/{1}_2_val_2.fq" ::: "${sample_names[@]}"

# don't dereplicate; wait to concatenate all fasta outputs before dereplicating
# parallel "vsearch --derep_fulllength {1}/{1}.fas --sizeout --output {1}/{1}_derep.fas" ::: "${sample_names[@]}"
