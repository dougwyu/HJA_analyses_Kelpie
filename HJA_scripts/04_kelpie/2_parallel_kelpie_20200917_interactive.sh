#!/bin/bash
set -e
set -u
set -o pipefail
##################################################################################################
##################################################################################################
# a shell script to loop through a set of fasta files, runs kelpie, and moves output files to a single folder
##################################################################################################
##################################################################################################

# I run this in an interactive session because it is fast, but it could probably be run
# as a batch job without needing editing, using sbatch
PATH=$PATH:~/scripts/vsearch-2.15.0-linux-x86_64/bin/ # downloaded 12 Jul 2020 from github
PATH=$PATH:~/scripts/Kelpie_v2.0.8/ubuntu-16.04/
PATH=$PATH:~/scripts/parallel-20170722/bin/ # GNU Parallel

# to run on testkelpie, upload into ~/_Oregon/2019Sep_shotgun/testkelpie/BWA01
# cd ~/_Oregon/2019Sep_shotgun/testkelpie/

# upload _parallel_kelpie_YYYYMMDD.sub and _parallel_kelpie_YYYYMMDD.sh *into* ~/_Oregon/2019Sep_shotgun/2.trimmeddata/
cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/ || exit
ls

#### copy all outputs from FilterReads run into allfilterreadsoutput/
# create folder
if [ ! -d allfilterreadsoutput ] # if directory allfilterreadsoutput does not exist.
then
     mkdir allfilterreadsoutput
fi
ls
# mv files
mv BWA*/filterreadsoutput/*_COI.fa ./allfilterreadsoutput
ls allfilterreadsoutput
find allfilterreadsoutput -type f -iname "*_COI.fa" | wc -l # 484 COI.fa files



#### run kelpie on indiv fasta files and save to kelpieoutputindiv/
# create folder to hold
if [ ! -d kelpieoutputindiv ] # if directory kelpieoutputindiv does not exist.
then
     mkdir kelpieoutputindiv
fi
ls
# assumption: all FilterReads outputs will be in one folder filterreadsoutput, one up from the BWA folders
find ./allfilterreadsoutput -type f -iname "*_COI.fa" -exec basename {} \; > fastalist.txt
	# remove _{1,2}_val_{1,2}_COI.fa from filenames 076361-M1-S1_BDSW190602952-1a_1_val_1_COI.fa
sed -i 's/_[1,2]_val_[1,2]_COI.fa//g' ./fastalist.txt
cat fastalist.txt
cat fastalist.txt | wc -l # 484 files
cat fastalist.txt | sort | uniq | wc -l # 242 unique names

# make array of fasta files
sample_info=fastalist.txt # put fastalist.txt into variable
sample_names=($(cut -f 1 "$sample_info" | sort | uniq)) # convert variable to array this way
echo "${sample_names[@]}" # echo all array elements
echo "There are" ${#sample_names[@]} "files that will be processed." # 242, echo number of elements in the array

# run kelpie on each individual COI.fa file
# run Kelpie in parallel. -j n means n samples at a time, -k keep same order as in array, --dryrun see the generated commands
# BF3BR2, -f CCHGAYATRGCHTTYCCHCG -r TCDGGRTGNCCRAARAAYCA, -max 500
	nohup parallel -k -j 1 "Kelpie_v2 -f CCHGAYATRGCHTTYCCHCG -r TCDGGRTGNCCRAARAAYCA -filtered -min 400 -max 500 allfilterreadsoutput/{1}_?_val_?_COI.fa kelpieoutputindiv/{1}_BF3BR2.fas" ::: "${sample_names[@]}" &
ls kelpieoutputindiv/*BF3BR2.fas | wc -l # 242

# Leray Fol-degen-rev, -f GGWACWGGWTGAACWGTWTAYCCYCC -r TANACYTCNGGRTGNCCRAARAAYCA, -min 300 -max 400
	nohup parallel -k -j 3 "Kelpie_v2 -f GGWACWGGWTGAACWGTWTAYCCYCC -r TANACYTCNGGRTGNCCRAARAAYCA -filtered -min 300 -max 400 allfilterreadsoutput/{1}_?_val_?_COI.fa kelpieoutputindiv/{1}_LERAY.fas" ::: "${sample_names[@]}" &
ls kelpieoutputindiv/*LERAY.fas | wc -l # 242



#### run kelpie on nearest-neighbor sets of files (each sample + five nearest neighbors)
     # read in each line of neighbors_20191204_wide.csv, index i
     # save the fields into separate sitename variables
     # run find command on all the sitename and cat the output files into a large filtered.fas
     # run kelpie on it, name the output by the index i, save to kelpieoutputneighbors folder
# http://www.compciv.org/topics/bash/loops/
# https://codefather.tech/blog/bash-loop-through-lines-file/
# https://www.cyberciti.biz/faq/unix-howto-read-line-by-line-from-file/
# https://bash.cyberciti.biz/guide/While_loop#Reading_A_Text_File_With_Separate_Fields

# start from here
# allfilterreadsoutput/ should already exist from above
cd ~/_Oregon/2019Sep_shotgun/2.trimmeddata/ || exit

if [ ! -d kelpieoutputneighbors ] # if directory kelpieoutputneighbors does not exist.
then
     mkdir kelpieoutputneighbors
fi

i=0 # set index to 0
while IFS=, read -r f1 f2 f3 f4 f5 f6
  do
      ((i=i+1)) # increment i by 1
      echo "starting i should be 1: $i"
      # find any files with these sitenames and concatenate into a single fasta file
      echo "creating kelpieinput_${i}.fa from these files:"
      find ./allfilterreadsoutput -type f -iname "$f1*" -o -iname "$f2*" -o -iname "$f3*" -o -iname "$f4*" -o -iname "$f5*" -o -iname "$f6*"
      echo "$f1, $f2, $f3, $f4, $f5, $f6"

      find ./allfilterreadsoutput -type f -iname "$f1*" -o -iname "$f2*" -o -iname "$f3*" -o -iname "$f4*" -o -iname "$f5*" -o -iname "$f6*" -exec cat {} + > kelpieinput_${i}.fa

      if [ ! -s kelpieinput_${i}.fa ] # if kelpieinput_${i}.fa has filesize==0, then delete it and exit loop
      then
           echo "deleting kelpieinput_${i}.fa"
           rm -f kelpieinput_${i}.fa || exit
       fi

      if [ -s kelpieinput_${i}.fa ] # if kelpieinput_$i.fa exists and filesize > 0
      then # check that the commands inside the then fi statement are preceded only by spaces, no tabs!
           echo "running kelpie on line ${i}"
           # Leray-FolDegenRev
               # Kelpie_v2 -f GGWACWGGWTGAACWGTWTAYCCYCC -r TANACYTCNGGRTGNCCRAARAAYCA -filtered -min 300 -max 400 kelpieinput_${i}.fa kelpieoutputneighbors/LERAY_${i}.fas
           # BF3BR2
               Kelpie_v2 -f CCHGAYATRGCHTTYCCHCG -r TCDGGRTGNCCRAARAAYCA -filtered -min 400 -max 500 kelpieinput_${i}.fa kelpieoutputneighbors/BF3BR2_${i}.fas
           echo "deleting kelpieinput_${i}.fa"
           rm -f kelpieinput_${i}.fa || exit
      fi
done < <(tail --lines=+2 neighbors_20191204_wide.csv)


#### other primers, with kelpie syntax

# fwhF2 Fol-degen-rev, -f CCHGAYATRGCHTTYCCHCG -r TCDGGRTGNCCRAARAAYCA, -max 500 (because some OTUs seem to get insertions, which i can remove later)

# 12S rRNA (MT-RNR1) (leech project, 82-150 bp), -f ACTGGGATTAGATACCCC -r YRGAACAGGCTCCTCTAG

# 16Smam rRNA (MT-RNR2, 81-117 bp), -f CGGTTGGGGTGACCTCGGA -r GCTGTTATCCCTAGGGTAACT

# trnL_P6 -f GGGCAATCCTGAGCCAA -r CCATYGAGTCTCTGCACCTATC

# 16S rRNA V4 (Earth Microbiome Project) 515F (Parada)–806R (Apprill) -f GTGYCAGCMGCCGCGGTAA -r GGACTACNVGGGTWTCTAAT



#### post-process kelpie outputs
# concatenate all kelpie outputs from individual and nearest-neighbor runs
# then dereplicate


# parallel "vsearch --derep_fulllength {1}/{1}.fas --sizeout --output {1}/{1}_derep.fas" ::: "${sample_names[@]}"
