#!/bin/bash
set -e
set -u
set -o pipefail
#######################################################################################
#######################################################################################
# a shell script to remove the pre-trimmed fastq.gz files after running trimgalore
#######################################################################################
#######################################################################################
# 6223.00 GB at start
# du -sh ~/_ArcDyn # 2.4T afterwards
# remove non-trimmed fastq.gz files in PlatesA2B2/PlatesA2B2_combined/
cd ~/_ArcDyn/PlatesA2B2/
ls PlatesA2B2_combined/Sample*/Sample*fastq.gz  # should list only fastq.gz files
ls PlatesA2B2_combined/Sample*/Sample*fq.gz # should list lots of fq.gz files
ls PlatesA2B2_combined/Sample*/Sample*fastq.gz | wc -l  # count of fastq.gz files
ls PlatesA2B2_combined/Sample*/Sample*fq.gz | wc -l # should have same count value
# rm -f PlatesA2B2_combined/Sample*/Sample*fastq.gz  # keep commented out until ready to run this command
ls PlatesA2B2_combined/Sample*/Sample*fastq.gz  # should list nothing
ls PlatesA2B2_combined/Sample*/Sample*fq.gz # should list lots of fq.gz files
ls PlatesA2B2_combined/Sample*/Sample*fq.gz | wc -l # should have same count value

# remove non-trimmed fastq.gz files in PlatesAB/PlatesAB_combined/
cd ~/_ArcDyn/PlatesAB/
ls PlatesAB_combined/Sample*/Sample*fastq.gz  # should list only fastq.gz files
ls PlatesAB_combined/Sample*/Sample*fq.gz # should list lots of fq.gz files
ls PlatesAB_combined/Sample*/Sample*fastq.gz | wc -l  # count of fastq.gz files
ls PlatesAB_combined/Sample*/Sample*fq.gz | wc -l  # should have same count value
# rm -f PlatesAB_combined/Sample*/Sample*fastq.gz  # keep commented out until ready to run this command
ls PlatesAB_combined/Sample*/Sample*fastq.gz  # should list nothing
ls PlatesAB_combined/Sample*/Sample*fq.gz # should list lots of fq.gz files
ls PlatesAB_combined/Sample*/Sample*fq.gz | wc -l  # should have same count value

# remove non-trimmed fastq.gz files in PlatesEF/PlatesEF_combined/
cd ~/_ArcDyn/PlatesEF/
ls PlatesEF_combined/Plate*/Plate*fastq.gz  # should list only fastq.gz files
ls PlatesEF_combined/Plate*/Plate*fq.gz # should list lots of fq.gz files
ls PlatesEF_combined/Plate*/Plate*fastq.gz | wc -l  # count of fastq.gz files
ls PlatesEF_combined/Plate*/Plate*fq.gz | wc -l  # should have same count value
# rm -f PlatesEF_combined/Plate*/Plate*fastq.gz  # keep commented out until ready to run this command
ls PlatesEF_combined/Plate*/Plate*fastq.gz  # should list nothing
ls PlatesEF_combined/Plate*/Plate*fq.gz # should list lots of fq.gz files
ls PlatesEF_combined/Plate*/Plate*fq.gz | wc -l  # should have same count value

# remove non-trimmed fastq.gz files in PlatesGH/PlatesGH_combined/
cd ~/_ArcDyn/PlatesGH/
du -sh PlatesGH_combined/
ls PlatesGH_combined/Sample*/Sample*fastq.gz  # should list only fastq.gz files
ls PlatesGH_combined/Sample*/Sample*fq.gz # should list lots of fq.gz files
ls -1 PlatesGH_combined/Sample*/Sample*fastq.gz | wc -l  # count of fastq.gz files
ls -1 PlatesGH_combined/Sample*/Sample*fq.gz | wc -l # should have same count value
# rm -f PlatesGH_combined/Sample*/Sample*fastq.gz  # keep commented out until ready to run this command
ls PlatesGH_combined/Sample*/Sample*fastq.gz  # should list nothing
ls PlatesGH_combined/Sample*/Sample*fq.gz # should list lots of fq.gz files
ls -1 PlatesGH_combined/Sample*/Sample*fq.gz | wc -l # should have same count value
du -sh PlatesGH_combined/
