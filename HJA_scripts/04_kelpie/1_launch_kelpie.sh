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
ssh b042@ada.uea.ac.uk
interactive
# path to GNU parallel
PATH=$PATH:~/scripts/parallel-20170722/bin/

############# copy the kelpie2 shell and bsub scripts into each BWA folder and edit the jobIDs
# cd ~/_Oregon/2019Sep_shotgun/testkelpie/; ls
cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/; ls
KELPIE2_SUB="_parallel_kelpie_20200716.sub"; head -n 20 ${KELPIE2_SUB}
KELPIE2_SH="_parallel_kelpie_20200716.sh"; head ${KELPIE2_SH}

parallel cp ${KELPIE2_SUB} BWA{} ::: 01 02 03 04 05 06 07 08 09 10
parallel cp ${KELPIE2_SH} BWA{} ::: 01 02 03 04 05 06 07 08 09 10
parallel ls BWA{} ::: 01 02 03 04 05 06 07 08 09 10
# ls BWA{01,02,03,04,05,06,07,08,09,10}

# edit the bsub files so that the correct jobID will show up (i suppose i could have instead run a job array...)
cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/; ls

parallel "sed 's/klplp/klplp{}/g' BWA{}/${KELPIE2_SUB} > BWA{}/${KELPIE2_SUB}_tmp" ::: 01 02 03 04 05 06 07 08 09 10
parallel "mv BWA{}/${KELPIE2_SUB}_tmp BWA{}/${KELPIE2_SUB}" ::: 01 02 03 04 05 06 07 08 09 10
head -n14 BWA{01,02,03,04,05,06,07,08,09,10}/${KELPIE2_SUB} # check.  should be klplp{01,02,03,...}
     # check that i'm using #SBATCH -p compute
tail -n2 BWA{01,02,03,04,05,06,07,08,09,10}/${KELPIE2_SUB} # check.  should be the correct shell file

ls # BWA* folders should now sort to bottom

####### launch the kelpie2 scripts #######
# cd into each BWA folder and submit bsub job
cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA01; ls
echo ${KELPIE2_SUB}
sbatch ${KELPIE2_SUB}
squeue -u b042
ls

cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA02; ls
sbatch ${KELPIE2_SUB}
squeue -u b042

cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA03; ls
sbatch ${KELPIE2_SUB}
squeue -u b042

cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA04; ls
sbatch ${KELPIE2_SUB}
squeue -u b042

cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA05; ls
sbatch ${KELPIE2_SUB}
squeue -u b042

cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA06; ls
sbatch ${KELPIE2_SUB}
squeue -u b042

cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA07; ls
sbatch ${KELPIE2_SUB}
squeue -u b042

cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA08; ls
sbatch ${KELPIE2_SUB}
squeue -u b042

cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA09; ls
sbatch ${KELPIE2_SUB}
squeue -u b042

cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA10; ls
sbatch ${KELPIE2_SUB}
squeue -u b042

ls ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA01/HOBO-357-M1-S2_BDSW190603169-1a/ # check
ls ~/_Oregon/2019Sep_shotgun/2.trimmeddata/BWA01/HOBO-349-M1-S1_BDSW190603162-1a
squeue -u b042
head -n14 ${KELPIE2_SUB}
