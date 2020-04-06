#!/bin/bash
set -e
set -u
set -o pipefail
##################################################################################################
##################################################################################################
# a shell script to loop through a set of fastq files and run kelpie
##################################################################################################
##################################################################################################

# Usage: bash _parallel_kelpie_20200319.sh
# upload _parallel_kelpie_20200319.bsub and _parallel_kelpie_20200319.sh *into* ~/_Oregon/2019Sep_shotgun/2.trimmeddata/

# to run in testkelpie folder, COMMENT OUT WHEN RUNNING FOR REAL IN 2.trimmeddata/
# HOMEFOLDER="/gpfs/home/b042/_Oregon/2019Sep_shotgun/" # should be ~/_Oregon/2019Sep_shotgun/
# TARGETFOLDER="testkelpie" #
# echo "Home folder is ${HOMEFOLDER}${TARGETFOLDER}" # should be ~/_Oregon/2019Sep_shotgun/testkelpie
# cd ${HOMEFOLDER}${TARGETFOLDER} || exit

# use _launch_kelpie.sh to move the *.sh and *.bsub scripts into each BWA folder, containing 24 or 25 samples
# read in folder list and make a bash array
find ./ -maxdepth 1 -mindepth 1 -type d -exec basename {} \; | sort > folderlist.txt
sed -i '/minimap2_outputs/d' ./folderlist.txt # remove minimap2 folder from folderlist.txt
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

# run Kelpie in parallel. -j n means n samples at a time, -k keep same order as in array, --dryrun see the generated commands
     # on short-ib queue
     # parallel -j 1 # 270 mins per 20 samples = 13.5 mins per sample
     # parallel -j 2 # 230 mins per 20 samples = 11.5 mins per sample
     # parallel -j 3 # 210 mins per 20 samples = 10.5 mins per sample
     # parallel -j 5 # 200 mins per 20 samples = 10 mins per sample
     # running 10 simultaneous jobs, each on 24 or 25 samples.  If using -j 2, can do each job in 4.5 hrs, so the whole thing should take less than 7 hrs, if some jobs are delayed in starting
# version with gzip -d
# fwhF2 Fol-degen-rev, -f CCHGAYATRGCHTTYCCHCG -r TCDGGRTGNCCRAARAAYCA, -max 500 (because some OTUs seem to get insertions, which i can remove later)
parallel -k -j 3 "gzip -d < {1}/{1}_1_val_1.fq.gz > {1}/{1}_1_val_1.fq; gzip -d < {1}/{1}_2_val_2.fq.gz > {1}/{1}_2_val_2.fq; singularity exec --bind /gpfs/home/b042/scripts:/scripts --bind /gpfs/scratch/b042:/mnt ~/ubuntukelpie.sif /scripts/Kelpie_v2.0.6/ubuntu-16.04/Kelpie_v2 -f CCHGAYATRGCHTTYCCHCG -r TCDGGRTGNCCRAARAAYCA -unfiltered -min 400 -max 500 {1}/{1}_?_val_?.fq {1}/{1}.fas; rm {1}/{1}_1_val_1.fq; rm {1}/{1}_2_val_2.fq" ::: "${sample_names[@]}"

# BF3BR2, -f CCHGAYATRGCHTTYCCHCG -r TCDGGRTGNCCRAARAAYCA, -max 500 (because some OTUs seem to get insertions, which i can remove later)
# parallel -k -j 3 "gzip -d < {1}/{1}_1_val_1.fq.gz > {1}/{1}_1_val_1.fq; gzip -d < {1}/{1}_2_val_2.fq.gz > {1}/{1}_2_val_2.fq; singularity exec --bind /gpfs/home/b042/scripts:/scripts --bind /gpfs/scratch/b042:/mnt ~/ubuntukelpie.sif /scripts/Kelpie_v2.0.6/ubuntu-16.04/Kelpie_v2 -f CCHGAYATRGCHTTYCCHCG -r TCDGGRTGNCCRAARAAYCA -unfiltered -min 400 -max 500 {1}/{1}_?_val_?.fq {1}/{1}.fas; rm {1}/{1}_1_val_1.fq; rm {1}/{1}_2_val_2.fq" ::: "${sample_names[@]}"

# 12S rRNA (MT-RNR1) (leech project, 82-150 bp), -f ACTGGGATTAGATACCCC -r YRGAACAGGCTCCTCTAG
# parallel -k -j 1 "gzip -d < {1}/{1}_1_val_1.fq.gz > {1}/{1}_1_val_1.fq; gzip -d < {1}/{1}_2_val_2.fq.gz > {1}/{1}_2_val_2.fq; singularity exec --bind /gpfs/home/b042/scripts:/scripts --bind /gpfs/scratch/b042:/mnt ~/ubuntukelpie.sif /scripts/Kelpie_v2.0.6/ubuntu-16.04/Kelpie_v2 -f ACTGGGATTAGATACCCC -r YRGAACAGGCTCCTCTAG -unfiltered {1}/{1}_?_val_?.fq {1}/{1}.fas; rm {1}/{1}_1_val_1.fq; rm {1}/{1}_2_val_2.fq" ::: "${sample_names[@]}"

# 16Smam rRNA (MT-RNR2, 81-117 bp), -f CGGTTGGGGTGACCTCGGA -r GCTGTTATCCCTAGGGTAACT
# parallel -k -j 1 "gzip -d < {1}/{1}_1_val_1.fq.gz > {1}/{1}_1_val_1.fq; gzip -d < {1}/{1}_2_val_2.fq.gz > {1}/{1}_2_val_2.fq; singularity exec --bind /gpfs/home/b042/scripts:/scripts --bind /gpfs/scratch/b042:/mnt ~/ubuntukelpie.sif /scripts/Kelpie_v2.0.6/ubuntu-16.04/Kelpie_v2 -f CGGTTGGGGTGACCTCGGA -r GCTGTTATCCCTAGGGTAACT -unfiltered {1}/{1}_?_val_?.fq {1}/{1}.fas; rm {1}/{1}_1_val_1.fq; rm {1}/{1}_2_val_2.fq" ::: "${sample_names[@]}"

# trnL_P6 -f GGGCAATCCTGAGCCAA -r CCATYGAGTCTCTGCACCTATC
# parallel -k -j 1 "gzip -d < {1}/{1}_1_val_1.fq.gz > {1}/{1}_1_val_1.fq; gzip -d < {1}/{1}_2_val_2.fq.gz > {1}/{1}_2_val_2.fq; singularity exec --bind /gpfs/home/b042/scripts:/scripts --bind /gpfs/scratch/b042:/mnt ~/ubuntukelpie.sif /scripts/Kelpie_v2.0.6/ubuntu-16.04/Kelpie_v2 -f GGGCAATCCTGAGCCAA -r CCATYGAGTCTCTGCACCTATC -unfiltered {1}/{1}_?_val_?.fq {1}/{1}.fas; rm {1}/{1}_1_val_1.fq; rm {1}/{1}_2_val_2.fq" ::: "${sample_names[@]}"

# 16S rRNA V4 (Earth Microbiome Project) 515F (Parada)–806R (Apprill) -f GTGYCAGCMGCCGCGGTAA -r GGACTACNVGGGTWTCTAAT
# parallel -k -j 1 "gzip -d < {1}/{1}_1_val_1.fq.gz > {1}/{1}_1_val_1.fq; gzip -d < {1}/{1}_2_val_2.fq.gz > {1}/{1}_2_val_2.fq; singularity exec --bind /gpfs/home/b042/scripts:/scripts --bind /gpfs/scratch/b042:/mnt ~/ubuntukelpie.sif /scripts/Kelpie_v2.0.6/ubuntu-16.04/Kelpie_v2 -f GTGYCAGCMGCCGCGGTAA -r GGACTACNVGGGTWTCTAAT -unfiltered {1}/{1}_?_val_?.fq {1}/{1}.fas; rm {1}/{1}_1_val_1.fq; rm {1}/{1}_2_val_2.fq" ::: "${sample_names[@]}"

# tested with unpigz -k # -k keep original file, but it was even slower than gunzip
# tested -tmp /mnt/{1} and -tmp /mnt # neither works when running as parallel

# old command using the kelpie binary inside the singularity container and with tmp files written to working directory instead of the scratch disk
# parallel -k -j 16 "gzip -d < {1}/{1}_1_val_1.fq.gz > {1}/{1}_1_val_1.fq; gzip -d < {1}/{1}_2_val_2.fq.gz > {1}/{1}_2_val_2.fq; singularity exec ~/ubuntukelpie.sif Kelpie_v2 -f CCHGAYATRGCHTTYCCHCG -r TCDGGRTGNCCRAARAAYCA -unfiltered -strict -min 400 -max 440 {1}/{1}_?_val_?.fq {1}/{1}.fas; rm {1}/{1}_1_val_1.fq; rm {1}/{1}_2_val_2.fq" ::: "${sample_names[@]}"

# don't dereplicate; wait to concatenate all fasta outputs before dereplicating
# parallel "vsearch --derep_fulllength {1}/{1}.fas --sizeout --output {1}/{1}_derep.fas" ::: "${sample_names[@]}"
