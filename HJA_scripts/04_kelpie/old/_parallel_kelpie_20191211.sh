#!/bin/bash
set -e
set -u
set -o pipefail
##################################################################################################
##################################################################################################
# a shell script to loop through a set of fastq files and run kelpie
##################################################################################################
##################################################################################################

# Usage: bash _parallel_kelpie_20191211.sh
# Upload _parallel_kelpie_20191211.bsub and _parallel_kelpie_20191211.sh *into* ~/_Oregon/

# go to correct folder
HOMEFOLDER="/gpfs/home/b042/_Oregon/2019Sep_shotgun/" # should be ~/_Oregon/2019Sep_shotgun/
TARGETFOLDER="testkelpie"
echo "Home folder is ${HOMEFOLDER}${TARGETFOLDER}" # ~/_Oregon/2019Sep_shotgun/2.trimmeddata
cd ${HOMEFOLDER}${TARGETFOLDER} || exit

# read in folder list and make a bash array
find ${HOMEFOLDER}/${TARGETFOLDER}/ -maxdepth 1 -type d -exec basename {} \; | sort > folderlist.txt # find all folders and then execute basename to remove folder pathnames
     # wc -l folderlist.txt # 6
     # cat folderlist.txt
sed -i '/testkelpie/d' ./folderlist.txt # GNU sed only
     # grep "testkelpie" folderlist.txt # should return nothing because it has been deleted
     # wc -l folderlist.txt # 5, in each of which there are 2 fastq files
     # cat folderlist.txt

# make array of folder names
sample_info=folderlist.txt # put folderlist.txt into variable
sample_names=($(cut -f 1 "$sample_info" | uniq)) # convert variable to array this way
echo "${sample_names[@]}" # echo all array elements
echo "There are" ${#sample_names[@]} "folders that will be processed." # echo number of elements in the array, which should be 5 (the same as the number of lines in folderlist.txt)

# run Kelpie in parallel. -j 2 means 2 samples, -k keep same order as in array, --dryrun see the commands
     # test run on 20191210 with -j 16 used 38 GB max and 7 hrs (approx) on 20 samples
     # 200 GB / 38 = 5X
     # i want 4 cores per sample, with #BSUB -n 24, i choose -j 6
     # 242 samples / 6 samples = 40 rounds, each taking around 20 mins
     # 40 * 0.33 hrs = 13.2 hrs
# version with gzip -d
parallel -k -j 1 "gzip -d < {1}/{1}_1_val_1.fq.gz > {1}/{1}_1_val_1.fq; gzip -d < {1}/{1}_2_val_2.fq.gz > {1}/{1}_2_val_2.fq; singularity exec --bind /gpfs/home/b042/scripts:/scripts --bind /gpfs/scratch/b042:/mnt ~/ubuntukelpie.sif /scripts/Kelpie_V2/ubuntu-16.04/Kelpie_v2 -f CCHGAYATRGCHTTYCCHCG -r TCDGGRTGNCCRAARAAYCA -unfiltered -strict -min 400 -max 440 {1}/{1}_?_val_?.fq {1}/{1}.fas; rm {1}/{1}_1_val_1.fq; rm {1}/{1}_2_val_2.fq" ::: "${sample_names[@]}"

# tested with unpigz -k # -k keep original file, but it was even slower than gunzip
# tested -tmp /mnt/{1} and -tmp /mnt # neither works when running as parallel
# parallel -j 1 # 13.5 mins per sample
# parallel -j 3 # 36 mins per 3 samples = 12 mins per sample
# parallel -j 5 # 3.3 hrs per 20 samples = 10 mins per sample

# old command using the kelpie binary inside the singularity container and with tmp files written to working directory instead of the scratch disk
# parallel -k -j 16 "gzip -d < {1}/{1}_1_val_1.fq.gz > {1}/{1}_1_val_1.fq; gzip -d < {1}/{1}_2_val_2.fq.gz > {1}/{1}_2_val_2.fq; singularity exec ~/ubuntukelpie.sif Kelpie_v2 -f CCHGAYATRGCHTTYCCHCG -r TCDGGRTGNCCRAARAAYCA -unfiltered -strict -min 400 -max 440 {1}/{1}_?_val_?.fq {1}/{1}.fas; rm {1}/{1}_1_val_1.fq; rm {1}/{1}_2_val_2.fq" ::: "${sample_names[@]}"

# don't dereplicate; wait to concatenate all fasta outputs before dereplicating
# parallel "vsearch --derep_fulllength {1}/{1}.fas --sizeout --output {1}/{1}_derep.fas" ::: "${sample_names[@]}"
