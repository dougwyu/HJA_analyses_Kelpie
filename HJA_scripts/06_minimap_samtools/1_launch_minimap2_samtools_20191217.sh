#!/bin/bash
set -e
set -u
set -o pipefail
#######################################################################################
#######################################################################################
# a shell script to launch bsub files.
#######################################################################################
#######################################################################################

# upload (scp) the new minimap and samtools sh and bsub files to _Oregon/2.trimmeddata

# defunct code for future bwa-mem2
     # one time only, index reference fasta
     # PATH=$PATH:~/scripts/bwa-mem2-2.0pre1_x64-linux/ # downloaded from github 18 Dec 2019
     # cd ~/_Oregon/2019Sep_shotgun/reference_seqs
     # bwa-mem2 index kelpie_20200214_BF3BR2_derep_filtered_geneious.fas
     # this creates a bunch of files that bwa-mem2 can use

############# edit minimap2 and samtools scripts #############
# ssh ada
# interactive
# to use parallel without a pathname in bsub scripts
PATH=$PATH:~/scripts/parallel-20170722/bin/

# cd ~/_Oregon/2019Sep_shotgun/testkelpie/; ls # test
cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/; ls # production

############## by hand, copy 1/10 the sample folders into each BWA folder
# mkdir BWA{01,02,03,04,05,06,07,08,09,10}; ls # BWA is prefix because this was the original mapping software
# if there are 171 sample folders:  hand move 17 into each BWA folder (easier than writing a script)
# BWA11 holds the dilution series samples
############# copy the minimap and samtools shell and bsub scripts into each BWA folder and edit the jobIDs
MINIMAP2_SUB="_loop_minimap2_only_20200917.sub"; echo ${MINIMAP2_SUB}
MINIMAP2_SH="_loop_minimap2_only_20200917.sh"; echo ${MINIMAP2_SH}
SAMTOOLS_SUB="_loop_samtools_only_20200917.sub"; echo ${SAMTOOLS_SUB}
SAMTOOLS_SH="_loop_samtools_only_20200917.sh"; echo ${SAMTOOLS_SH}

cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/; ls # 2.trimmeddata
parallel cp ${MINIMAP2_SUB} BWA{} ::: 01 02 03 04 05 06 07 08 09 10
parallel cp ${MINIMAP2_SH} BWA{} ::: 01 02 03 04 05 06 07 08 09 10
parallel cp ${SAMTOOLS_SUB} BWA{} ::: 01 02 03 04 05 06 07 08 09 10
parallel cp ${SAMTOOLS_SH} BWA{} ::: 01 02 03 04 05 06 07 08 09 10
ls BWA{01,02,03,04,05,06,07,08,09,10}

# edit the bsub files so that the correct jobID will show up (i suppose i could have instead run a job array...)
cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/; ls # 2.trimmeddata

parallel "sed 's/mnmp01/mnmp{}/g' BWA{}/${MINIMAP2_SUB} > BWA{}/${MINIMAP2_SUB}_tmp" ::: 01 02 03 04 05 06 07 08 09 10
parallel "mv BWA{}/${MINIMAP2_SUB}_tmp BWA{}/${MINIMAP2_SUB}" ::: 01 02 03 04 05 06 07 08 09 10
head -n7 BWA{01,02,03,04,05,06,07,08,09,10}/${MINIMAP2_SUB} # check.  should be mnmp{01,02,03,...}
     # check if i'm using mellanox-ib or short-eth
tail -n2 BWA{01,02,03,04,05,06,07,08,09,10}/${MINIMAP2_SUB} # check.  should be the correct shell file

parallel "sed 's/samtl01/samtl{}/g' BWA{}/${SAMTOOLS_SUB} > BWA{}/${SAMTOOLS_SUB}_tmp" ::: 01 02 03 04 05 06 07 08 09 10
parallel "mv BWA{}/${SAMTOOLS_SUB}_tmp BWA{}/${SAMTOOLS_SUB}" ::: 01 02 03 04 05 06 07 08 09 10
head -n 7 BWA{01,02,03,04,05,06,07,08,09,10}/${SAMTOOLS_SUB} # check.  should have the correct index number
     # check if i'm using mellanox-ib or short-eth
tail -n 1 BWA{01,02,03,04,05,06,07,08,09,10}/${SAMTOOLS_SUB} # check.  should have the correct samtools shell filename

ls # BWA* folders should now sort to bottom

####### launch the minimap2 scripts #######
# cd into each BWA folder and submit bsub job
cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA01; ls # 2.trimmeddata
echo "${MINIMAP2_SUB}" # if incorrect or NULL, go up to top of script
sbatch ${MINIMAP2_SUB}
squeue -u b042
ls # should show working folder

cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA02; ls # 2.trimmeddata
sbatch ${MINIMAP2_SUB}
squeue -u b042

cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA03; ls # 2.trimmeddata
sbatch ${MINIMAP2_SUB}
squeue -u b042

cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA04; ls # 2.trimmeddata
sbatch ${MINIMAP2_SUB}
squeue -u b042

cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA05; ls # 2.trimmeddata
sbatch ${MINIMAP2_SUB}
squeue -u b042

cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA06; ls # 2.trimmeddata
sbatch ${MINIMAP2_SUB}
squeue -u b042

cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA07; ls # 2.trimmeddata
sbatch ${MINIMAP2_SUB}
squeue -u b042

cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA08; ls # 2.trimmeddata
sbatch ${MINIMAP2_SUB}
squeue -u b042

cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA09; ls # 2.trimmeddata
sbatch ${MINIMAP2_SUB}
squeue -u b042

cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA10; ls # 2.trimmeddata
sbatch ${MINIMAP2_SUB}
squeue -u b042

# dilution series only
# cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA11; ls # 2.trimmeddata
# sbatch ${MINIMAP2_SUB}
# squeue -u b042


squeue -u b042
ls ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA01/ # check if there is a working folder
tail ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA01/mnmp01.out




######  WAIT FOR THE MINIMAP2 JOBS TO FINISH BEFORE LAUNCHING THE SAMTOOLS SCRIPTS
ssh ada
interactive

####### launch samtools scripts #######
cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA01; ls # 2.trimmeddata
echo "${SAMTOOLS_SUB}" # if incorrect or NULL, go up to top of script
sbatch ${SAMTOOLS_SUB}
squeue -u b042
ls minimap2_outputs/

cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA02; ls # 2.trimmeddata
sbatch ${SAMTOOLS_SUB}
squeue -u b042

cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA03; ls # 2.trimmeddata
sbatch ${SAMTOOLS_SUB}
squeue -u b042

cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA04; ls # 2.trimmeddata
sbatch ${SAMTOOLS_SUB}
squeue -u b042

cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA05; ls # 2.trimmeddata
sbatch ${SAMTOOLS_SUB}
squeue -u b042

cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA06; ls # 2.trimmeddata
sbatch ${SAMTOOLS_SUB}
squeue -u b042

cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA07; ls # 2.trimmeddata
sbatch ${SAMTOOLS_SUB}
squeue -u b042

cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA08; ls # 2.trimmeddata
sbatch ${SAMTOOLS_SUB}
squeue -u b042

cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA09; ls # 2.trimmeddata
sbatch ${SAMTOOLS_SUB}
squeue -u b042

cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA10; ls # 2.trimmeddata
sbatch ${SAMTOOLS_SUB}
squeue -u b042

# dilution series samples
# cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA11; ls # 2.trimmeddata
# sbatch ${SAMTOOLS_SUB}
# squeue -u b042

ls ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA01/minimap2_outputs/ # should show new genomecov files
squeue -u b042
