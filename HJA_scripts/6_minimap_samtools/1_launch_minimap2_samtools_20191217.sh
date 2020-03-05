PlatesA2B2_combined#!/bin/bash
set -e
set -u
set -o pipefail
#######################################################################################
#######################################################################################
# a shell script to launch bsub files.
#######################################################################################
#######################################################################################

# upload (scp) the new minimap and samtools sh and bsub files to _Oregon/2.trimmeddata
# upload the reference fasta file kelpie_20200214_BF3BR2_derep_filtered_geneious.fas to _Oregon/2.trimmeddata

# defunct code for future bwa-mem2
     # one time only, index reference fasta
     # PATH=$PATH:~/scripts/bwa-mem2-2.0pre1_x64-linux/ # downloaded from github 18 Dec 2019
     # cd ~/_Oregon/2019Sep_shotgun/reference_seqs
     # bwa-mem2 index kelpie_20200214_BF3BR2_derep_filtered_geneious.fas
     # this creates a bunch of files that bwa-mem2 can use

############# edit minimap2 and samtools scripts #############
# ssh hpc
# interactive
# to use parallel without a pathname in bsub scripts
PATH=$PATH:~/scripts/parallel-20170722/bin/

cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/; ls # 2.trimmeddata, testkelpie

############## by hand, copy 1/10 the sample folders into each BWA folder
# mkdir BWA{01,02,03,04,05,06,07,08,09,10}; ls # BWA is prefix because this was the original mapping software
# if there are 171 sample folders:  hand move 17 into each BWA folder (easier than writing a script)
# BWA11 holds the dilution series samples
############# copy the minimap and samtools shell and bsub scripts into each BWA folder and edit the jobIDs
MINIMAP2_BSUB="_loop_minimap2_only_20200219.bsub"; echo ${MINIMAP2_BSUB}
MINIMAP2_SH="_loop_minimap2_only_20200219.sh"; echo ${MINIMAP2_SH}
SAMTOOLS_BSUB="_loop_samtools_only_20200219.bsub"; echo ${SAMTOOLS_BSUB}
SAMTOOLS_SH="_loop_samtools_only_20200219.sh"; echo ${SAMTOOLS_SH}

cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/; ls # 2.trimmeddata
parallel cp ${MINIMAP2_BSUB} BWA{} ::: 01 02 03 04 05 06 07 08 09 10
parallel cp ${MINIMAP2_SH} BWA{} ::: 01 02 03 04 05 06 07 08 09 10
parallel cp ${SAMTOOLS_BSUB} BWA{} ::: 01 02 03 04 05 06 07 08 09 10
parallel cp ${SAMTOOLS_SH} BWA{} ::: 01 02 03 04 05 06 07 08 09 10
ls BWA{01,02,03,04,05,06,07,08,09,10}

# edit the bsub files so that the correct jobID will show up (i suppose i could have instead run a job array...)
cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/; ls # 2.trimmeddata

parallel "sed 's/mnmploop01/mnmp{}/g' BWA{}/${MINIMAP2_BSUB} > BWA{}/${MINIMAP2_BSUB}_tmp" ::: 01 02 03 04 05 06 07 08 09 10
parallel "mv BWA{}/${MINIMAP2_BSUB}_tmp BWA{}/${MINIMAP2_BSUB}" ::: 01 02 03 04 05 06 07 08 09 10
head -n7 BWA{01,02,03,04,05,06,07,08,09,10}/${MINIMAP2_BSUB} # check.  should be mnmp{01,02,03,...}
     # check if i'm using mellanox-ib or short-eth
tail -n2 BWA{01,02,03,04,05,06,07,08,09,10}/${MINIMAP2_BSUB} # check.  should be the correct shell file

parallel "sed 's/samtools01/smtl{}/g' BWA{}/${SAMTOOLS_BSUB} > BWA{}/${SAMTOOLS_BSUB}_tmp" ::: 01 02 03 04 05 06 07 08 09 10
parallel "mv BWA{}/${SAMTOOLS_BSUB}_tmp BWA{}/${SAMTOOLS_BSUB}" ::: 01 02 03 04 05 06 07 08 09 10
head -n 7 BWA{01,02,03,04,05,06,07,08,09,10}/${SAMTOOLS_BSUB} # check.  should have the correct index number
     # check if i'm using mellanox-ib or short-eth
tail -n 1 BWA{01,02,03,04,05,06,07,08,09,10}/${SAMTOOLS_BSUB} # check.  should have the correct samtools shell filename

ls # BWA* folders should now sort to bottom

####### launch the minimap2 scripts #######
# cd into each BWA folder and submit bsub job
cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA01; ls # 2.trimmeddata
echo "${MINIMAP2_BSUB}" # if incorrect or NULL, go up to top of script
bsub < ${MINIMAP2_BSUB}
bjobs
ls # should show working folder

cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA02; ls # 2.trimmeddata
bsub < ${MINIMAP2_BSUB}
bjobs

cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA03; ls # 2.trimmeddata
bsub < ${MINIMAP2_BSUB}
bjobs

cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA04; ls # 2.trimmeddata
bsub < ${MINIMAP2_BSUB}
bjobs

cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA05; ls # 2.trimmeddata
bsub < ${MINIMAP2_BSUB}
bjobs

cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA06; ls # 2.trimmeddata
bsub < ${MINIMAP2_BSUB}
bjobs

cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA07; ls # 2.trimmeddata
bsub < ${MINIMAP2_BSUB}
bjobs

cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA08; ls # 2.trimmeddata
bsub < ${MINIMAP2_BSUB}
bjobs

cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA09; ls # 2.trimmeddata
bsub < ${MINIMAP2_BSUB}
bjobs

cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA10; ls # 2.trimmeddata
bsub < ${MINIMAP2_BSUB}
bjobs

cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA11; ls # 2.trimmeddata
bsub < ${MINIMAP2_BSUB}
bjobs


bqueues
ls ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA01/ # check if there is a working folder




######  WAIT FOR THE MINIMAP2 JOBS TO FINISH BEFORE LAUNCHING THE SAMTOOLS SCRIPTS
# ssh hpc
# interactive

####### launch samtools scripts #######
cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA01; ls # 2.trimmeddata
echo "${SAMTOOLS_BSUB}" # if incorrect or NULL, go up to top of script
bsub < ${SAMTOOLS_BSUB}
bjobs

cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA02; ls # 2.trimmeddata
bsub < ${SAMTOOLS_BSUB}
bjobs

cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA03; ls # 2.trimmeddata
bsub < ${SAMTOOLS_BSUB}
bjobs

cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA04; ls # 2.trimmeddata
bsub < ${SAMTOOLS_BSUB}
bjobs

cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA05; ls # 2.trimmeddata
bsub < ${SAMTOOLS_BSUB}
bjobs

cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA06; ls # 2.trimmeddata
bsub < ${SAMTOOLS_BSUB}
bjobs

cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA07; ls # 2.trimmeddata
bsub < ${SAMTOOLS_BSUB}
bjobs

cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA08; ls # 2.trimmeddata
bsub < ${SAMTOOLS_BSUB}
bjobs

cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA09; ls # 2.trimmeddata
bsub < ${SAMTOOLS_BSUB}
bjobs

cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA10; ls # 2.trimmeddata
bsub < ${SAMTOOLS_BSUB}
bjobs

# dilution series samples
# cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA11; ls # 2.trimmeddata
# bsub < ${SAMTOOLS_BSUB}
# bjobs

ls minimap2_outputs/ # should show new bam files
