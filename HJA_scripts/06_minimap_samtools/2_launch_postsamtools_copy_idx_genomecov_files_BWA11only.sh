#!/bin/bash
set -e
set -u
set -o pipefail
#######################################################################################
#######################################################################################
# a shell script to copy idx and genomcov files into a directory to be downloaded to my computer
# run after running the samtools script, does not have to be uploaded to hpc
#######################################################################################
#######################################################################################

ssh hpc
interactive
# to use parallel without a pathname in bsub scripts
PATH=$PATH:~/scripts/parallel-20170722/bin/


####### set up output folder to hold everything before running
cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA11; ls # 2.trimmeddata, testkelpie

# set filters and minimum mapping quality scores
FILTER1="F2308_f0x2" # filter 1
FILTER2="F2308" # filter 2
echo $FILTER1
echo $FILTER2
QUAL1=""; echo $QUAL1 # set to "" if i don't want to use this variable for, say, q1
QUAL2=48; echo $QUAL2

MAPDATE="20200221"
TARGET="kelpie_20200214_BF3BR2_derep_filtered_geneious_vsearch97_min2_spikes" #
mkdir outputs_${FILTER1}_q${QUAL2}_${OUTPUTFOLDER}_${MAPDATE}_${TARGET}
ls
# copy output files from minimap2_outputs/ into outputs_${FILTER1}_q${QUAL2}_${OUTPUTFOLDER}_${MAPDATE}_${TARGET}
echo outputs_${FILTER1}_q${QUAL2}_${OUTPUTFOLDER}_${MAPDATE}_${TARGET}

cp minimap2_outputs/*.bam_idxstats.txt outputs_${FILTER1}_q${QUAL2}_${OUTPUTFOLDER}_${MAPDATE}_${TARGET}
ls outputs_${FILTER1}_q${QUAL2}_${OUTPUTFOLDER}_${MAPDATE}_${TARGET}
ls outputs_${FILTER1}_q${QUAL2}_${OUTPUTFOLDER}_${MAPDATE}_${TARGET}/*.bam_idxstats.txt | wc -l # 66

cp minimap2_outputs/*_sorted_genomecov_d.txt.gz outputs_${FILTER1}_q${QUAL2}_${OUTPUTFOLDER}_${MAPDATE}_${TARGET}
ls outputs_${FILTER1}_q${QUAL2}_${OUTPUTFOLDER}_${MAPDATE}_${TARGET}
ls outputs_${FILTER1}_q${QUAL2}_${OUTPUTFOLDER}_${MAPDATE}_${TARGET}/*_sorted_genomecov_d.txt.gz | wc -l # 66

cp minimap2_outputs/*_sorted.bam.flagstat.txt outputs_${FILTER1}_q${QUAL2}_${OUTPUTFOLDER}_${MAPDATE}_${TARGET}
ls outputs_${FILTER1}_q${QUAL2}_${OUTPUTFOLDER}_${MAPDATE}_${TARGET}
ls outputs_${FILTER1}_q${QUAL2}_${OUTPUTFOLDER}_${MAPDATE}_${TARGET}/*_sorted.bam.flagstat.txt | wc -l # 33


# rename, tar, and gzip for download
du -sh minimap2_outputs/ # 169 M
# set filename to something that i can understand after download
# tar gzip for download
tar -czvf outputs_${FILTER1}_q${QUAL2}_${OUTPUTFOLDER}_${MAPDATE}_${TARGET}.tar.gz outputs_${FILTER1}_q${QUAL2}_${OUTPUTFOLDER}_${MAPDATE}_${TARGET}/
# filename format:  outputs_F2308_f0x2_q48_minimap2_outputs_20191219_kelpie2091217_vsearch97_filtered_geneious_filtered.tar.gz
ls
rm -rf outputs_${FILTER1}_q${QUAL2}_${OUTPUTFOLDER}_${MAPDATE}_${TARGET}/
ls

# use Transmit to download *.gz file to Luo_Mingjie_Oregon/Kelpie_maps/
