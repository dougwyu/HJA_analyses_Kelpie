#!/bin/bash
set -e
set -u
set -o pipefail
#######################################################################################
#######################################################################################
# a shell script to loop through a set of fastq files and run cutadapt
#######################################################################################
#######################################################################################

# Usage: bash _parallel_trimgalore_20191203.sh

# Batch procedure:  bsub script and submit to mellanox-ib queue
     ######### _parallel_trimgalore_20191203.sh ##############################################################################
     #######################################################################################
     # #!/bin/sh
     # #BSUB -q long-ib     #long-ib & mellanox-ib (168 hours = 7 days), short-ib (24 hours)
     # #BSUB -J trimgal01
     # #BSUB -oo trimgal01.out
     # #BSUB -eo trimgal01.err
     # #BSUB -R "rusage[mem=5000]"  # set at 5000 for the max (n.b. short-ib, long-ib, mellanox-ib)
     # #BSUB -M 5000
     # #BSUB -x # exclusive access to node, should be in anything that is threaded
     # #BSUB -B        # sends email to me when job starts
     # #BSUB -N        # sends email to me when job finishes
     #
     # . /etc/profile
     # module purge
     # module load samtools # samtools 1.5
     # module load python/anaconda/4.2/2.7 # Python 2.7.12
     # module load java  # java/jdk1.8.0_51
     # module load gcc # needed to run bedtools
     # PATH=$PATH:~/scripts/cutadapt-1.11/build/scripts-2.7
     # PATH=$PATH:~/scripts/vsearch-2.6.2-linux-x86_64/bin/ # downloaded 22 Jan 2018 from github
     # PATH=$PATH:~/scripts/TrimGalore-0.4.5/ # downloaded 22 Jan 2018 from github
     # PATH=$PATH:~/scripts/minimap2/  # made 22 Jan 2018 from github 2.7-r659-dirty
     # PATH=$PATH:~/scripts/bedtools2/bin # made 22 Jan 2018 from github 2.27.1
     #
     #
     # bash _parallel_trimgalore_20191203.sh
     #######################################################################################
     #######################################################################################

# SCRIPT OUTLINE:
     # make list of sample folders that hold the fastq files
     # run trim_galore in each sample folder, using GNU parallel

PIPESTART=$(date)

# upload _parallel_trimgalore_20191203.sh *into* ~/_Oregon/
HOMEFOLDER=$(pwd) # this sets working directory
TARGETFOLDER="2019Sep_shotgun/1.rawdata"
echo "Home folder is ${HOMEFOLDER}/${TARGETFOLDER}"
cd ${HOMEFOLDER}/${TARGETFOLDER} || exit

# read in folder list and make a bash array
find ${HOMEFOLDER}/${TARGETFOLDER}/ -maxdepth 1 -type d -exec basename {} \; | sort > folderlist.txt # find all folders and then execute basename to remove folder pathnames
     # alternative command using find
          # find ${HOMEFOLDER}/${TARGETFOLDER}/ -maxdepth 1 -type d -printf '%f\n' | sort > folderlist.txt  # find all folders and print only the foldername with "-printf '%f\n'"
# https://unix.stackexchange.com/questions/118841/use-basename-to-parse-a-list-of-paths-held-in-a-file
     # grep "1.rawdata" folderlist.txt # this is the parent folder, which needs to be removed
     # wc -l folderlist.txt # 243 lines
sed -i '/1.rawdata/d' ./folderlist.txt # GNU sed only
     # grep "1.rawdata" folderlist.txt # should return nothing
     # wc -l folderlist.txt # 242 lines == 242 samples, in each of which there are 2 fastq files

sample_info=folderlist.txt # put folderlist.txt into variable
sample_names=($(cut -f 1 "$sample_info" | uniq)) # convert variable to array this way
echo "${sample_names[@]}" # echo all array elements
echo "There are" ${#sample_names[@]} "folders that will be processed." # echo number of elements in the array, which should be 242 (the same as the number of lines in folderlist.txt)

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
