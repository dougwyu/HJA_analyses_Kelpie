#!/bin/bash
set -e
set -u
set -o pipefail
##################################################################################################
##################################################################################################
# a shell script to interactively clean up after running _parallel_kelpie
##################################################################################################
##################################################################################################

# on HPC
# create a large interactive job
# interactive -x -R "rusage[mem=20000]" -M 20000
interactive -R "rusage[mem=8000]" -M 8000

# set paths and add modules
PATH=$PATH:~/scripts/vsearch-2.14.1-linux-x86_64/bin # downloaded 10 Dec 2019 from github
PATH=$PATH:~/scripts/parallel-20170722/bin/ # GNU Parallel
PATH=$PATH:~/scripts/seqkitdir/ # downloaded 12 Dec 2019 from github

# go to correct folder
HOMEFOLDER="/gpfs/home/b042/_Oregon/2019Sep_shotgun/" # should be ~/_Oregon/2019Sep_shotgun/
TARGETFOLDER="testkelpie" # 2.trimmeddata/
echo "Home folder is ${HOMEFOLDER}${TARGETFOLDER}" # ~/_Oregon/2019Sep_shotgun/2.trimmeddata
cd ${HOMEFOLDER}${TARGETFOLDER} # || exit

#### After running kelpie, there will be fasta files in each sample folder
# Now, i need to concatenate, deduplicate names, and clean up

# find all the fasta output filenames
find ./ -maxdepth 3 -mindepth 2 -type f -name "*.fas" | sort > fastalist.txt # find all fasta files 3 levels down
sed -i '/discards/d' ./fastalist.txt  # -i save output back to the input file. remove discards.fas
sed -i '/derep/d' ./fastalist.txt  # -i save output back to the input file. remove derep.fas
cat fastalist.txt
cat fastalist.txt | wc -l

# concatenate the fasta files into a single large fasta file
timestamp="20191214"
cat $(grep -v '^#' fastalist.txt) > kelpietest_${timestamp}.fas
# https://stackoverflow.com/questions/11619500/how-to-cat-multiple-files-from-a-list-of-files-in-bash

# filter out sequences < 400 bp
seqkit stats kelpietest_${timestamp}.fas
seqkit seq -m 400 kelpietest_${timestamp}.fas -o kelpietest_${timestamp}_min400.fas
seqkit stats kelpietest_${timestamp}_min400.fas kelpietest_${timestamp}.fas
mv kelpietest_${timestamp}_min400.fas kelpietest_${timestamp}.fas
seqkit stats kelpietest_${timestamp}.fas

# a problem is that the amplicon names are reused across multiple samples (e.g. >R1 is used in each sample fasta)
# this causes a problem after clustering if the names collide (although usually avoided because sizes are different)
# so i use seqkit rename to rename duplicated names
# deduplicate header names
grep "R1$" kelpietest_${timestamp}.fas # there are duplicate headers, because i concatenated multiple fasta files
seqkit rename kelpietest_${timestamp}.fas > kelpietest_${timestamp}_rename.fas # renames duplicate headers
mv kelpietest_${timestamp}_rename.fas kelpietest_${timestamp}.fas # replace old with new, deduplicated version
grep "R1$" kelpietest_${timestamp}.fas # duplicates have new names, plus the original name info
seqkit stats kelpietest_${timestamp}.fas # 51722 sequences

# dereplicate the fasta file
# --threads 0 # zero to use all avail cores
vsearch --derep_fulllength kelpietest_${timestamp}.fas --sizeout --threads 0 --output kelpietest_${timestamp}_derep.fas
seqkit stats kelpietest_${timestamp}_derep.fas # 675 sequences

# cleanup, uncomment when ready to run
# find ./ -maxdepth 3 -mindepth 2 -type f -name "*.fas" # check that i list only the fasta files 3 levels down
# find ./ -maxdepth 3 -mindepth 2 -type f -name "*.fas" | wc -l # should be twice the number of samples
# find ./ -maxdepth 3 -mindepth 2 -type f -name "*.fas" | xargs rm # rm those files
# find ./ -maxdepth 3 -mindepth 2 -type f -name "*.fas" # should return nothing
#      # https://shapeshed.com/unix-xargs/
#      # alternative if i already have a file list
#      # xargs rm < fastalist.txt
#      # https://askubuntu.com/questions/596489/how-to-delete-files-listed-in-a-text-file
# rm -f kelpietest_${timestamp}.fas
# rm -f folderlist.txt
# rm -f fastalist.txt
