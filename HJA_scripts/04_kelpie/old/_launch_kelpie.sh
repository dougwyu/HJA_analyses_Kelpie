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
ssh hpc
interactive
# path to GNU parallel
PATH=$PATH:~/scripts/parallel-20170722/bin/

############# copy the kelpie2 shell and bsub scripts into each BWA folder and edit the jobIDs
cd ~/_Oregon/2019Sep_shotgun/testkelpie/; ls  # testkelpie/
KELPIE2_BSUB="_parallel_kelpie_20191213.bsub"; echo ${KELPIE2_BSUB}
KELPIE2_SH="_parallel_kelpie_20191213.sh"; echo ${KELPIE2_SH}

parallel cp ${KELPIE2_BSUB} BWA{} ::: 01 02 03 04 05 06 07 08 09 10
parallel cp ${KELPIE2_SH} BWA{} ::: 01 02 03 04 05 06 07 08 09 10
parallel ls BWA{} ::: 01 02 03 04 05 06 07 08 09 10
# ls BWA{01,02,03,04,05,06,07,08,09,10}

# edit the bsub files so that the correct jobID will show up (i suppose i could have instead run a job array...)
cd ~/_Oregon/2019Sep_shotgun/testkelpie/; ls

parallel "sed 's/klplp/klplp{}/g' BWA{}/${KELPIE2_BSUB} > BWA{}/${KELPIE2_BSUB}_tmp" ::: 01 02 03 04 05 06 07 08 09 10
parallel "mv BWA{}/${KELPIE2_BSUB}_tmp BWA{}/${KELPIE2_BSUB}" ::: 01 02 03 04 05 06 07 08 09 10
head -n7 BWA{01,02,03,04,05,06,07,08,09,10}/${KELPIE2_BSUB} # check.  should be klplp{01,02,03,...}
     # check if i'm using mellanox-ib or short-ib
tail -n2 BWA{01,02,03,04,05,06,07,08,09,10}/${KELPIE2_BSUB} # check.  should be the correct shell file

ls # BWA* folders should now sort to bottom

####### launch the kelpie2 scripts #######
# cd into each BWA folder and submit bsub job
cd ~/_Oregon/2019Sep_shotgun/testkelpie/BWA01; ls
bsub < ${KELPIE2_BSUB}
bjobs
ls

cd ~/_Oregon/2019Sep_shotgun/testkelpie/BWA02; ls
bsub < ${KELPIE2_BSUB}
bjobs

cd ~/_Oregon/2019Sep_shotgun/testkelpie/BWA03; ls
bsub < ${KELPIE2_BSUB}
bjobs

cd ~/_Oregon/2019Sep_shotgun/testkelpie/BWA04; ls
bsub < ${KELPIE2_BSUB}
bjobs

cd ~/_Oregon/2019Sep_shotgun/testkelpie/BWA05; ls
bsub < ${KELPIE2_BSUB}
bjobs

cd ~/_Oregon/2019Sep_shotgun/testkelpie/BWA06; ls
bsub < ${KELPIE2_BSUB}
bjobs

cd ~/_Oregon/2019Sep_shotgun/testkelpie/BWA07; ls
bsub < ${KELPIE2_BSUB}
bjobs

cd ~/_Oregon/2019Sep_shotgun/testkelpie/BWA08; ls
bsub < ${KELPIE2_BSUB}
bjobs

cd ~/_Oregon/2019Sep_shotgun/testkelpie/BWA09; ls
bsub < ${KELPIE2_BSUB}
bjobs

cd ~/_Oregon/2019Sep_shotgun/testkelpie/BWA10; ls
bsub < ${KELPIE2_BSUB}
bjobs

bjobs
bqueues
ls ~/_Oregon/2019Sep_shotgun/testkelpie/BWA01/ # check
