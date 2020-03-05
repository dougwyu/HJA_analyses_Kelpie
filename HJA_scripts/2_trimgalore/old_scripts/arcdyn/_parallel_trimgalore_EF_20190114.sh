#!/bin/bash
set -e
set -u
set -o pipefail
#######################################################################################
#######################################################################################
# a shell script to loop through a set of fastq files and run cutadapt
#######################################################################################
#######################################################################################

# Usage: bash _parallel_trimgalore_EF_20190114.sh

# Interactive procedure:  first run these commands if running _loop_trimgalore_20180216.sh interactively
     # interactive -x -R "rusage[mem=20000]" -M 22000  # alternative is:  interactive -q interactive-lm
     # module load samtools # samtools 1.5
     # module load python/anaconda/4.2/2.7 # Python 2.7.12
     # module load java  # java/jdk1.8.0_51
     # module load gcc # needed to run bedtools
     # PATH=$PATH:~/scripts/cutadapt-1.11/build/scripts-2.7
     # PATH=$PATH:~/scripts/bwa_0.7.17_r1188 # made 22 Jan 2018 from github
     # PATH=$PATH:~/scripts/bfc/ # r181 # made 3 Aug 2015 from github
     # PATH=$PATH:~/scripts/vsearch-2.6.2-linux-x86_64/bin/ # downloaded 22 Jan 2018 from github
     # PATH=$PATH:~/scripts/TrimGalore-0.4.5/ # downloaded 22 Jan 2018 from github
     # PATH=$PATH:~/scripts/minimap2/  # made 22 Jan 2018 from github 2.7-r659-dirty
     # PATH=$PATH:~/scripts/bedtools2/bin # made 22 Jan 2018 from github 2.27.1

# Batch procedure:  bsub script and submit to mellanox-ib queue
     ######### _parallel_trimgalore_EF_20190114.sh ##############################################################################
     #######################################################################################
     # #!/bin/sh
     # #BSUB -q short-ib     #long-ib & mellanox-ib (168 hours = 7 days), short-ib (24 hours)
     # #BSUB -J trimgal01
     # #BSUB -oo trimgal01.out
     # #BSUB -eo trimgal01.err
     # #BSUB -R "rusage[mem=5000]"  # set at 5000 for the max (n.b. short-ib, long-ib, mellanox-ib)
     # #BSUB -M 5000
     # #BSUB -x # exclusive access to node, should be in anything that is threaded
     # #BSUB -B        #sends email to me when job starts
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
     # bash _parallel_trimgalore_20190114.sh
     #######################################################################################
     #######################################################################################

# SCRIPT OUTLINE:
     # make list of sample folders that hold the fastq files
     # run trim_galore in each sample folder, using GNU parallel

PIPESTART=$(date)

# upload _parallel_trimgalore_EF_20190114.sh *into* folder that contains the sample folders that i want to map
HOMEFOLDER=$(pwd) # this sets working directory to a given BWA folder
echo "Home folder is ${HOMEFOLDER}"

COMBINEDFASTQ=$(basename ${HOMEFOLDER})
echo "COMBINEDFASTQ folder will be ${COMBINEDFASTQ}_combined/"

# read in folder list and make a bash array
# because this is PlatesEF, I change the prefix to Plate*
find ${COMBINEDFASTQ}_combined/ -maxdepth 1 -type d -name "Plate*" > folderlist.txt  # find all folders starting with "Plate*"
     # head folderlist.txt
     # wc -l folderlist.txt
sed -i '/PlatesEF_combined\/$/d' ./folderlist.txt  # remove first line, the enclosing folder:  'PlatesEF_combined/'
     # head folderlist.txt
     # tail folderlist.txt
     # wc -l folderlist.txt

sample_info=folderlist.txt # put folderlist.txt into variable
sample_names=($(cut -f 1 "$sample_info" | uniq)) # convert variable to array this way
# echo "${sample_names[@]}" # echo all array elements
echo "There are" ${#sample_names[@]} "folders that will be processed." # echo number of elements in the array


# GNU PARALLEL version, since trimgalore only uses one cpu per job
parallel --progress "echo Now on species {1} of ${#sample_names[@]}; trim_galore --paired --length 100 -trim-n {1}/*_R1.fastq.gz {1}/*_R2.fastq.gz -o {1}" ::: "${sample_names[@]}"
# -o {1} # to write output fq.gz files to sample folder. output is *_val_1_fq.gz
# --paired --length 100 # remove a paired read if either is less than 100 bp after trimming (default 20)

echo "Pipeline started at $PIPESTART"
NOW=$(date)
echo "Pipeline ended at $NOW"
