#!/bin/bash
set -e
set -u
set -o pipefail
##################################################################################################
##################################################################################################
# a shell script to interactively clean up after running _parallel_kelpie
##################################################################################################
##################################################################################################

# on HPC
ssh b042@ada.uea.ac.uk
# create a large interactive job
srun -n4 -p compute -J interactive --time=04:00:00 --mem=6G --pty bash
# or
interactive
squeue -u b042

# set paths and add modules
PATH=$PATH:~/scripts/vsearch-2.15.0-linux-x86_64/bin/ # downloaded 12 Jul 2020 from github
PATH=$PATH:~/scripts/parallel-20170722/bin/ # GNU Parallel
PATH=$PATH:~/scripts/seqkitdir/ # downloaded 12 Dec 2019 from github
# clustering methods not used
    # PATH=$PATH:~/scripts/swarm-3.0.0-linux-x86_64/bin # downloaded 16 Dec 2019 from github
    # module add Sumaclust/1.0.34 # not avail on ada

# go to correct folder
HOMEFOLDER="/gpfs/home/b042/_Oregon/2019Sep_shotgun/" # should be ~/_Oregon/2019Sep_shotgun/
# TARGETFOLDER="testkelpie/" # testkelpie
TARGETFOLDER="2.trimmeddata/" # testkelpie
echo "Home folder is ${HOMEFOLDER}${TARGETFOLDER}" # ~/_Oregon/2019Sep_shotgun/2.trimmeddata
cd ${HOMEFOLDER}${TARGETFOLDER} # || exit

#### After running kelpie, there will be fasta files in each sample folder
# Now, i need to concatenate, deduplicate names, and clean up

# find all the fasta output filenames
find ./ -maxdepth 3 -mindepth 3 -type f -name "*.fas" | sort > fastalist.txt # find all fasta files 3 levels down
sed -i '/discards/d' ./fastalist.txt  # -i save output back to the input file. remove discards.fas
sed -i '/derep/d' ./fastalist.txt  # -i save output back to the input file. remove derep.fas
cat fastalist.txt
cat fastalist.txt | sort | uniq | wc -l # should see 242 fasta files

# concatenate the fasta files into a single large fasta file
timestamp="20200717_BF3BR2" # e.g. 20200715_LERAYFOL, 20200717_BF3BR2
echo $timestamp
cat $(grep -v '^#' fastalist.txt) > kelpie_${timestamp}.fas
ls
# https://stackoverflow.com/questions/11619500/how-to-cat-multiple-files-from-a-list-of-files-in-bash

# filter out sequences < 300 bp
seqkit stats kelpie_${timestamp}.fas # 625,646 BF3BR2 seqs 2.0.4, 590,200 BF3BR2 seqs 2.0.6, 360,852 Leray seqs 2.0.8 (w/out 1 sample), 606,517 BF3BR2 seqs 2.0.8
seqkit seq -m 400 kelpie_${timestamp}.fas -o kelpie_${timestamp}_min400.fas
seqkit stats kelpie_${timestamp}_min400.fas kelpie_${timestamp}.fas # 625,646 to 625,592 seqs 2.0.4, 590,200 to 590,166 2.0.6, 360,748 Leray seqs 2.0.8, 606,460 BF3BR2 seqs 2.0.8
mv kelpie_${timestamp}_min400.fas kelpie_${timestamp}.fas
seqkit stats kelpie_${timestamp}.fas

# a problem is that the amplicon names are reused across multiple samples (e.g. >R1 is used in each sample fasta)
# this causes a problem after clustering if the names collide (although usually avoided because sizes are different)
# so i use seqkit rename to rename duplicated names
# deduplicate header names
grep "R1$" kelpie_${timestamp}.fas # there are duplicate headers, because i concatenated multiple fasta files
grep "R1$" kelpie_${timestamp}.fas | wc -l # should equal to number of samples
seqkit rename kelpie_${timestamp}.fas > kelpie_${timestamp}_rename.fas # renames duplicate headers
mv kelpie_${timestamp}_rename.fas kelpie_${timestamp}.fas # replace old with new, deduplicated version
grep "R1$" kelpie_${timestamp}.fas # duplicates have new names, plus the original name info
seqkit stats kelpie_${timestamp}.fas # 360,748 seqs

# dereplicate the fasta file
    # orig command
    # vsearch --derep_fulllength kelpie_${timestamp}.fas --sizeout --sortbysize --threads 0 --output kelpie_${timestamp}_derep.fas
# --threads 0 # zero to use all avail cores
# --relabel_sha1 # hash each sequence as the name (to compare sequences)
# --fastq_width 0 # no line breaks in sequence
# --sizeout # include size= information
# command formatted for input to swarm 3.0.0
vsearch --derep_fulllength kelpie_${timestamp}.fas --sizeout --fasta_width 0 --threads 0 --output kelpie_${timestamp}_derep.fas # --relabel_sha1

seqkit stats kelpie_${timestamp}_derep.fas # 4,849 BF3BR2 uniq seqs 2.0.4, 5,426 BF3BR2 uniq seqs 2.0.6, 3,215 uniq Leray seqs 2.0.8, 4,062 uniq BF3BR2 seqs 2.0.8

# uncomment when ready to run
# rm -f kelpie_${timestamp}.fas
ls

# this is what i am going to use as input to GBIF, after which i will align, correct indels, and then cluster into 97% OTUs (and check if i lose any of the species)

# cleanup, uncomment when ready to run
find ./ -maxdepth 3 -mindepth 3 -type f -name "*.fas" # check that i list only the fasta files 3 levels down
find ./ -maxdepth 3 -mindepth 3 -type f -name "*.fas" | wc -l # should be twice the number of samples
find ./ -maxdepth 3 -mindepth 3 -type f -name "*_discards.fas" # should return only discards
find ./ -maxdepth 3 -mindepth 3 -type f -name "*_discards.fas" | wc -l # should be the number of samples (242)
find ./ -maxdepth 3 -mindepth 3 -type f -name "*_discards.fas" | sort > discard_fastalist.txt # make list of discard fastas
head discard_fastalist.txt; tail discard_fastalist.txt
cat discard_fastalist.txt | sort | uniq | wc -l # 242 fasta files
# https://stackoverflow.com/questions/11619500/how-to-cat-multiple-files-from-a-list-of-files-in-bash
#      # https://shapeshed.com/unix-xargs/
#      # alternative if i already have a file list
#      # xargs rm < fastalist.txt
#      # https://askubuntu.com/questions/596489/how-to-delete-files-listed-in-a-text-file


# concatenate the discard fasta files into a single large fasta file
timestamp="20200717_BF3BR2" # e.g. 20200716_LERAYFOL_discards
echo $timestamp
cat $(grep -v '^#' discard_fastalist.txt) > kelpie_${timestamp}_discards.fas
ls
head kelpie_${timestamp}_discards.fas
grep "D1" kelpie_${timestamp}_discards.fas # duplicates have new names, plus the original name info

# cleanup, uncomment when ready to run
# find ./ -maxdepth 3 -mindepth 3 -type f -name "*discards.fas" # check that i list only the fasta files 3 levels down
# find ./ -maxdepth 3 -mindepth 3 -type f -name "*discards.fas" | wc -l # should be twice the number of samples
# find ./ -maxdepth 3 -mindepth 3 -type f -name "*discards.fas" | xargs rm # rm those files
# find ./ -maxdepth 3 -mindepth 3 -type f -name "*.fas" | wc -l # should be twice the number of samples
# find ./ -maxdepth 3 -mindepth 3 -type f -name "*.fas" | xargs rm # rm those files
rm -f kelpie_${timestamp}.fas
rm -f folderlist.txt
rm -f fastalist.txt
rm -f discard_fastalist.txt
ls

# END



#### vsearch clustering code
# vsearch --cluster_fast kelpie_${timestamp}_derep.fas --sizein --sizeout --id 0.97 --centroids kelpie_${timestamp}_vsearch97.fas # --uc kelpie_${timestamp}_clusters.uc
# vsearch --sortbysize kelpie_${timestamp}_vsearch97.fas --output kelpie_${timestamp}_vsearch97_allseqs.fas
# seqkit stats kelpie_${timestamp}_vsearch97_allseqs.fas # 1,155 OTUs 2.0.4, 1,285 OTUs 2.0.6
# swarm results in multiple OTUs that get assigned to the same species in BOLD, which means that they will not receive mappings due to competing, similar sequences
# sumaclust does not have an option to output only the representative sequences
# five species only appeared as OTUs with swarm and pure dereplication (e.g. Helina troene, Hybomitra affinis, Hybomitra rhombica, Laphria columbica, Phaonia errans). I will see if they reappear when i have a larger dataset


#### swarm clustering code
# # swarm clustering
# swarm -v # Swarm 3.0.0
# # -d 1 # default and recommended
# # -f # fastidious
# # -z # use usearch size= for abundance
# # -s filename # statistics file
# # -w filename # OTU representative fasta file
# # > /dev/null # discard std output
# swarm -t 16 -d 1 -f -z -s OTU_stats_${timestamp}.txt -w kelpie_${timestamp}_swarmOTUs.fas kelpie_${timestamp}_derep.fas > /dev/null
# #
# seqkit stats kelpie_${timestamp}_swarmOTUs.fas # 439 OTUs
# # for download to use on bold
# # split representative OTU set into subsets of 100 seqs for upload to BOLDSystems (provides the BOLD consensus)
# seqkit split -s 100 kelpie_${timestamp}_swarmOTUs.fas


# split representative OTU set into subsets of 100 seqs for upload to BOLDSystems (provides the BOLD consensus)
# not needed, use https://www.gbif.org/tools/sequence-id
# seqkit split -s 100 kelpie_${timestamp}_centroids_sort.fas

# sumaclust clustering
# sumaclust -t 0.97 -e kelpie_${timestamp}.fas > kelpie_${timestamp}_sumaclust97.fas
# tail kelpie_${timestamp}_sumaclust97.fas
