#!/bin/bash
set -e
set -u
set -o pipefail
#######################################################################################
#######################################################################################
# a shell script to launch bsub files.
#######################################################################################
#######################################################################################

# mkdir BWA{01,02,03,04,05,06,07,08,09,10} # make 10 folders
# by hand, I moved 24 folders into each of the 10 BWA folders

############# edit kelpie script #############
ssh ada
interactive
# path to GNU parallel
PATH=$PATH:~/scripts/parallel-20200922/bin/ # GNU Parallel

############# copy the FilterReads shell and sub scripts into each BWA folder and edit the jobIDs
# cd ~/_Oregon/2019Sep_shotgun/testkelpie/; ls
cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/; ls
FILTREAD2_SUB="_parallel_FilterReads_20200915.sub"; head -n30 ${FILTREAD2_SUB}
FILTREAD2_SH="_parallel_FilterReads_20200915.sh"; head -n60 ${FILTREAD2_SH}

parallel cp ${FILTREAD2_SUB} BWA{} ::: 01 02 03 04 05 06 07 08 09 10
parallel cp ${FILTREAD2_SH} BWA{} ::: 01 02 03 04 05 06 07 08 09 10
parallel ls -lrt BWA{} ::: 01 02 03 04 05 06 07 08 09 10
# ls BWA{01,02,03,04,05,06,07,08,09,10}

# edit the bsub files so that the correct jobID will show up (i suppose i could have instead run a job array...)
cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/; ls

parallel "sed 's/filtrd/filtrd{}/g' BWA{}/${FILTREAD2_SUB} > BWA{}/${FILTREAD2_SUB}_tmp" ::: 01 02 03 04 05 06 07 08 09 10
parallel "mv BWA{}/${FILTREAD2_SUB}_tmp BWA{}/${FILTREAD2_SUB}" ::: 01 02 03 04 05 06 07 08 09 10
head -n14 BWA{01,02,03,04,05,06,07,08,09,10}/${FILTREAD2_SUB} # check.  should be filtrd{01,02,03,...}
     # check that i'm using #SBATCH -p compute-24-96
tail -n2 BWA{01,02,03,04,05,06,07,08,09,10}/${FILTREAD2_SUB} # check.  should be the correct shell file

ls # BWA* folders should now sort to bottom

####### launch the FilterReads scripts #######
# cd into each BWA folder and submit bsub job
# each job takes about an hour to run, 4 jobs are allowed to run simultaneously on ada
cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA01; ls
echo ${FILTREAD2_SUB}
sbatch ${FILTREAD2_SUB}
squeue -u b042
ls

cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA02; ls
sbatch ${FILTREAD2_SUB}
squeue -u b042

cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA03; ls
sbatch ${FILTREAD2_SUB}
squeue -u b042

cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA04; ls
sbatch ${FILTREAD2_SUB}
squeue -u b042

cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA05; ls
sbatch ${FILTREAD2_SUB}
squeue -u b042

cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA06; ls
sbatch ${FILTREAD2_SUB}
squeue -u b042

cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA07; ls
sbatch ${FILTREAD2_SUB}
squeue -u b042

cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA08; ls
sbatch ${FILTREAD2_SUB}
squeue -u b042

cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA09; ls
sbatch ${FILTREAD2_SUB}
squeue -u b042

cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA10; ls
sbatch ${FILTREAD2_SUB}
squeue -u b042

ls ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA01/HOBO-357-M1-S2_BDSW190603169-1a/ # check
ls ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA01/HOBO-349-M1-S1_BDSW190603162-1a/
squeue -u b042
head -n14 ${FILTREAD2_SUB}
