#!/bin/sh
set -e
set -u
set -o pipefail
########a shell script to loop through my Oregon shotgun fastq files and run AdapterRemvoval

cd /gpfs/home/luomj/Oregon/shotgun/Result-X101SC19072198-Z01-F002-B2-21/1.rawdata/

if [ ! -d AdapterRemoval_output ] # if directory AdapterRemoval_output does not exist
then
	mkdir  AdapterRemoval_output
fi

ls -d */ > folderlist.txt #put all folder names into a text file
sed 's/\///' folderlist.txt > folders.txt #rm "/"

sample_info=folders.txt
sample_names=($(cut -f 1 "$sample_info" | uniq))

for sample in ${sample_names[@]}
do
  cd "${sample}"
  AdapterRemoval --file1 *1.fq.gz --file2 *2.fq.gz --basename ${sample} --threads 10 --trimn
  rm *.discarded
  mv *.truncated ../AdapterRemoval_output/
  cd ..
done
