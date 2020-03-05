#!/bin/bash
set -e
set -u
set -o pipefail
##################################################################################################
##################################################################################################
# a shell script to interactively test Kelpie_v2 on hpc
##################################################################################################
##################################################################################################

# on HPC
# create a large interactive job
# interactive -x -R "rusage[mem=20000]" -M 20000
interactive -R "rusage[mem=8000]" -M 8000

# create testkelpie folder with 20 samples
cd ~/_Oregon/2019Sep_shotgun/
# copy 20 sample folders to testkelpie/
# mkdir testkelpie
# cp -R 2.trimmeddata/SM-PG-M1-S2_BDSW190603195-1a/ testkelpie/
# cp -R 2.trimmeddata/SM-PG-M1-S1_BDSW190603194-2a/ testkelpie/
# cp -R 2.trimmeddata/SM-09-M2-S2_BDSW190603193-2a/ testkelpie/
# cp -R 2.trimmeddata/SM-09-M2-S1_BDSW190603192-2a/ testkelpie/
# cp -R 2.trimmeddata/SM-09-M1-S2_BDSW190603191-2a/ testkelpie/
# cp -R 2.trimmeddata/076361-M1-S1_BDSW190602952-1a/ testkelpie/
# cp -R 2.trimmeddata/076361-M1-S2_BDSW190602953-1a/ testkelpie/
# cp -R 2.trimmeddata/081090-M1-S1_BDSW190602954-1a/ testkelpie/
# cp -R 2.trimmeddata/081090-M1-S2_BDSW190602955-1a/ testkelpie/
# cp -R 2.trimmeddata/123545-M1-S1_BDSW190602956-1a/ testkelpie/
# cp -R 2.trimmeddata/123545-M1-S2_BDSW190602957-1a/ testkelpie/
# cp -R 2.trimmeddata/124031-M1-S1_BDSW190602958-1a/ testkelpie/
# cp -R 2.trimmeddata/124031-M1-S2_BDSW190602959-1a/ testkelpie/
# cp -R 2.trimmeddata/124031-M2-S1_BDSW190602960-1a/ testkelpie/
# cp -R 2.trimmeddata/124031-M2-S2_BDSW190602961-1a/ testkelpie/
# cp -R 2.trimmeddata/HOBO-014-M1-S1_BDSW190603089-1a/ testkelpie/
# cp -R 2.trimmeddata/HOBO-014-M1-S2_BDSW190603090-1a/ testkelpie/
# cp -R 2.trimmeddata/HOBO-016-M1-S2_BDSW190603091-1a/ testkelpie/
# cp -R 2.trimmeddata/HOBO-032-M1-S1_BDSW190603092-1a/ testkelpie/
# cp -R 2.trimmeddata/HOBO-032-M1-S2-1_BDSW190603094-1a/ testkelpie/

# set paths and add modules
PATH=$PATH:~/scripts/vsearch-2.14.1-linux-x86_64/bin # downloaded 10 Dec 2019 from github
PATH=$PATH:~/scripts/parallel-20170722/bin/ # GNU Parallel
PATH=$PATH:~/scripts/pigz-2.4 # downloaded 11 Dec 2019 from https://zlib.net/pigz/
module add singularity/3.3.0
singularity exec ~/ubuntukelpie.sif Kelpie_v2 -h # test if binary runs

# go to correct folder
HOMEFOLDER="/gpfs/home/b042/_Oregon/2019Sep_shotgun/" # should be ~/_Oregon/2019Sep_shotgun/
TARGETFOLDER="testkelpie"
echo "Home folder is ${HOMEFOLDER}${TARGETFOLDER}" # ~/_Oregon/2019Sep_shotgun/2.trimmeddata
cd ${HOMEFOLDER}${TARGETFOLDER} # || exit

#### run Kelpie on one sample
# unzip input files
unpigz -k SM-PG-M1-S2_BDSW190603195-1a/SM-PG-M1-S2_BDSW190603195-1a_1_val_1.fq.gz # increases 5X from 2.4 to 11.3 GB
unpigz -k SM-PG-M1-S2_BDSW190603195-1a/SM-PG-M1-S2_BDSW190603195-1a_2_val_2.fq.gz
# gzip -d < SM-PG-M1-S2_BDSW190603195-1a/SM-PG-M1-S2_BDSW190603195-1a_1_val_1.fq.gz > SM-PG-M1-S2_BDSW190603195-1a/SM-PG-M1-S2_BDSW190603195-1a_1_val_1.fq # increases 5X from 2.4 to 11.3 GB
# gzip -d < SM-PG-M1-S2_BDSW190603195-1a/SM-PG-M1-S2_BDSW190603195-1a_2_val_2.fq.gz > SM-PG-M1-S2_BDSW190603195-1a/SM-PG-M1-S2_BDSW190603195-1a_2_val_2.fq

# run Kelpie.  temp kelpie files are written to my scratch disk
singularity exec --bind /gpfs/home/b042/scripts:/scripts --bind /gpfs/scratch/b042:/mnt ~/ubuntukelpie.sif /scripts/Kelpie_V2/ubuntu-16.04/Kelpie_v2 -f CCHGAYATRGCHTTYCCHCG -r TCDGGRTGNCCRAARAAYCA -unfiltered -strict -min 400 -max 440 -tmp /mnt SM-PG-M1-S2_BDSW190603195-1a/SM-PG-M1-S2_BDSW190603195-1a_?_val_?.fq SM-PG-M1-S2_BDSW190603195-1a/SM-PG-M1-S2_BDSW190603195-1a.fas


# dereplicate output fasta file
# vsearch --derep_fulllength SM-PG-M1-S2_BDSW190603195-1a/SM-PG-M1-S2_BDSW190603195-1a.fas --sizeout --output SM-PG-M1-S2_BDSW190603195-1a/SM-PG-M1-S2_BDSW190603195-1a_derep.fas --threads 3

# remove unzipped fastq files
# rm SM-PG-M1-S2_BDSW190603195-1a/SM-PG-M1-S2_BDSW190603195-1a_1_val_1.fq
# rm SM-PG-M1-S2_BDSW190603195-1a/SM-PG-M1-S2_BDSW190603195-1a_2_val_2.fq



#### run Kelpie on multiple samples
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
parallel -k -j 2 "gzip -d < {1}/{1}_1_val_1.fq.gz > {1}/{1}_1_val_1.fq; gzip -d < {1}/{1}_2_val_2.fq.gz > {1}/{1}_2_val_2.fq; singularity exec ~/ubuntukelpie.sif Kelpie_v2 -f CCHGAYATRGCHTTYCCHCG -r TCDGGRTGNCCRAARAAYCA -unfiltered -strict -min 400 -max 440 {1}/{1}_?_val_?.fq {1}/{1}.fas; rm {1}/{1}_1_val_1.fq; rm {1}/{1}_2_val_2.fq" ::: "${sample_names[@]}"


# concatenate all the fasta output filenames
find */*.fas -maxdepth 1 -type f | sort > fastalist.txt # find all fasta files one level down
sed -i '/discards/d' ./fastalist.txt  # -i save output back to the input file. remove discards.fas
sed -i '/derep/d' ./fastalist.txt  # -i save output back to the input file. remove derep.fas
cat fastalist.txt
cat fastalist.txt | wc -l

timestamp="20191211"
# concatenate the contents of the fasta files in fastalist.txt into a large fasta file
cat $(grep -v '^#' fastalist.txt) > kelpietest_${timestamp}.fas
# https://stackoverflow.com/questions/11619500/how-to-cat-multiple-files-from-a-list-of-files-in-bash
head kelpietest_${timestamp}.fas
grep ">" kelpietest_${timestamp}.fas | wc -l # 51731 sequences
# dereplicate the fasta file
# --threads 0 # zero for all cores
vsearch --derep_fulllength kelpietest_${timestamp}.fas --sizeout --threads 0 --output kelpietest_${timestamp}_derep.fas
grep ">" kelpietest_${timestamp}_derep.fas | wc -l # 193 sequences

# cleanup, uncomment when ready to run
# find */*.fas -maxdepth 1 -type f # check that i list only the fasta files one level down
# find */*.fas -maxdepth 1 -type f | xargs rm # rm those files
#      # https://shapeshed.com/unix-xargs/
#      # alternative if i already have a file list
#      # xargs rm < fasta_all_list.txt
#      # https://askubuntu.com/questions/596489/how-to-delete-files-listed-in-a-text-file
# rm -f kelpietest_${timestamp}.fas
# rm -f folderlist.txt
# rm -f fastalist.txt
