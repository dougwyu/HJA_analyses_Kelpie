#!/bin/bash
set -e
set -u
set -o pipefail
#######################################################################################
#######################################################################################
# a shell script to loop through a set of fastq files and run minimap2, samtools, bedtools
#######################################################################################
#######################################################################################

# Usage: bash _loop_minimap2_only_20191217.sh
# this version uses the 308 mitogenome reference fasta or the COI barcode fasta

# SCRIPT OUTLINE:
     # make list of folders
     # loop through folder names for each SAMPLE
     # move trimmed fastq files from $SAMPLE name folder to a working folder
     # run minimap2, and samtools view -b, samtools sort, samtools index, and samtools flagstat
     # move the outputs to minimap2_outputs
     # rm working folder
     # next folder

PIPESTART=$(date)

# upload _loop_minimap2_only_20191217.sh *into* folder that contains the sample folders that i want to map
# when i have lots of sample files, i break it up by putting the sample files into BWA**/ folders and running these scripts inside each one
HOMEFOLDER=$(pwd) # this sets working directory to that folder
echo "Home folder is ${HOMEFOLDER}"

# set variables
INDEX=1
if [ ! -d minimap2_outputs ] # if directory minimap2_outputs does not exist.  this is within the HOMEFOLDER
then
	mkdir minimap2_outputs
fi

# read in folder list and make a bash array
# run inside each BWA folder
find ./ -maxdepth 1 -mindepth 1 -type d -exec basename {} \; | sort > folderlist.txt
sed -i '/minimap2_outputs/d' ./folderlist.txt # GNU sed only
     # grep "minimap2_outputs" folderlist.txt # should return nothing because it has been deleted
     # wc -l folderlist.txt # 20
     # cat folderlist.txt
sample_info=folderlist.txt # put folderlist.txt into variable
sample_names=($(cut -f 1 "$sample_info" | uniq)) # convert variable to array this way
echo "${sample_names[@]}" # echo all array elements
echo "There are" ${#sample_names[@]} "folders that will be processed." # echo number of elements in the array

for sample in "${sample_names[@]}"  # ${sample_names[@]} is the full bash array
do
     cd "${HOMEFOLDER}"
     echo "Now on Sample" ${INDEX} of ${#sample_names[@]}". Moved back to starting directory $(date):"
     INDEX=$((INDEX+1))
     # pwd
     FOLDER="${sample}" # this sets the folder name to the value in the bash array (which is a list of the folders)
     echo "**** Working folder is" $FOLDER

     mkdir _${sample}_working

     echo "**** copying trimmed fastq.gz files to working folder $(date)"
     # this makes it easier to clean up if the job is interrupted.  i just rm the working folder
     cp ${FOLDER}/${sample}_1_val_1.fq.gz "_${sample}_working/"
     cp ${FOLDER}/${sample}_2_val_2.fq.gz "_${sample}_working/"
     echo "**** trimmed fastq.gz files moved to working folder $(date)"
     cd _${sample}_working; # pwd

     echo "**** start of minimap2, sam to bam conversion, sorting of bam file $(date)"
     #### minimap2 ####
     # against kelpie_20191217_OTUs and 2 COI_spike barcodes
     # minimap2 -ax sr ~/_Oregon/2019Sep_shotgun/reference_seqs/kelpie_20191217_vsearch97_filtered_geneious_filtered.fas ${sample}_1_val_1.fq.gz ${sample}_2_val_2.fq.gz | samtools view -b | samtools sort -@15 - -o ${sample}_sorted.bam
	minimap2 -ax sr ~/_Oregon/2019Sep_shotgun/reference_seqs/kelpie_20191217_vsearch97_filtered_geneious_filtered.fas ${sample}_1_val_1.fq.gz ${sample}_2_val_2.fq.gz | samtools sort -@15 - -o ${sample}_sorted.bam # skip the samtools view -b step 

	#### bwa mem ####
	# bwa-mem2 did not work, as of 19 Dec 2019
	# interactive -R "rusage[mem=16000]" -M 16000
	# # PATH=$PATH:~/scripts/bwa-mem2-2.0pre1_x64-linux/ # downloaded from github 18 Dec 2019
	# PATH=$PATH:~/scripts/bwa_0.7.17_r1188/ # downloaded from github 19 Dec 2019
	# module load samtools # samtools 1.5
	# bwa index ~/_Oregon/2019Sep_shotgun/reference_seqs/kelpie_20191217_vsearch97_filtered_geneious_filtered.fas
	# cd ~/_Oregon/2019Sep_shotgun/testkelpie/BWA01/076361-M1-S1_BDSW190602952-1a/; ls
	# bwa mem ~/_Oregon/2019Sep_shotgun/reference_seqs/kelpie_20191217_vsearch97_filtered_geneious_filtered.fas 076361-M1-S1_BDSW190602952-1a_1_val_1.fq.gz 076361-M1-S1_BDSW190602952-1a_2_val_2.fq.gz | samtools sort -@8 - -o 076361-M1-S1_BDSW190602952-1a_sorted.bam

     echo "**** end of minimap2, sam to bam conversion, sorting of bam file $(date)"

     # calculate flagstats
     samtools flagstat ${sample}_sorted.bam > ${sample}_sorted.bam.flagstat.txt # basic stats

     echo "**** start of samtools filtering of unmapped reads $(date)"
     # filter out unmapped reads and delete original sorted.bam
     samtools view -b -F 0x4 ${sample}_sorted.bam > ${sample}_F0x4_sorted.bam
     rm -f ${sample}_sorted.bam
     samtools index ${sample}_F0x4_sorted.bam # creates ${sample}_F0x4_sorted.bam.bai file
     echo "**** end of samtools filtering of unmapped reads $(date)"

     echo "**** start of moving outputs to minimap2_outputs/ $(date)"
     mv ${sample}_sorted.bam.flagstat.txt ../minimap2_outputs/
     mv ${sample}_F0x4_sorted.bam.bai ../minimap2_outputs/
     mv ${sample}_F0x4_sorted.bam ../minimap2_outputs/
     echo "**** end of moving outputs to minimap2_outputs/ $(date)"

     cd "${HOMEFOLDER}"
     rm -rf "_${sample}_working" # remove the working directory to make space

done

mv folderlist.txt minimap2_outputs/

echo "Pipeline started at $PIPESTART"
NOW=$(date)
echo "Pipeline ended at   $NOW"
