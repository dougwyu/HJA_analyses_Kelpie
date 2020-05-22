## load packages
# devtools::install_github("ropensci/bold") # for the dev version
library("bold") # now at 0.9.0  use ≥0.4.1.9300 If not avail from cran, use github version
library("plyr") # always load plyr before dplyr
library("seqinr") # used to read fasta file
library("dplyr")
library("readr")
# library("here")
library(conflicted)
    conflict_prefer("filter", "dplyr")
    conflict_prefer("select", "dplyr")
    conflict_prefer("arrange", "dplyr")
    conflict_prefer("here", "here")

## bold_identify multiple sequences from a fasta file ####
setwd("~/Dropbox/Working_docs/Luo_Mingjie_Oregon/HJA_scripts/")
# set_here(path = "~/Dropbox/Working_docs/Luo_Mingjie_Oregon/HJA_scripts/")
# here() # "/Users/Negorashi2011/Dropbox/Working_docs/Luo_Mingjie_Oregon/HJA_scripts"
# here::dr_here()


#### load a fasta file
testseqsall <- read.fasta(file = file.path("kelpietest_20191216_centroids_sort.fas"), seqtype = "DNA", set.attributes = FALSE, as.string = TRUE, forceDNAtolower = FALSE)  # read in fasta file using seqinr

#### bold_identify sequences
# all sequences at once, useful for smaller fasta files (<50 seqs)
# testseqsallbold <- bold_identify(testseqsall, db = "COX1")

# one sequence at a time (for large files that take a long time and might time out partway, this saves the output of each sequence)
testseqsallbold <- vector('list')
for (i in 191:length(testseqsall)) {
    print(i)
    boldoutput <- bold_identify(testseqsall[i], db = "COX1")
    testseqsallbold <- c(testseqsallbold, boldoutput)
}

# cleanup
rm(boldoutput)

#### to remove any null elements
#### WATCH IF list size gets smaller 
testseqsallbold <- testseqsallbold[!sapply(testseqsallbold, is.null)] 


#### bold_identify_parents, add parent taxa to output, wide format 
# all sequences at once (faster than bold_identify so might work for large files)
# testseqsallbold_parents <- bold_identify_parents(testseqsallbold, wide = TRUE)

# one sequence at a time (for large files that take a long time and might time out part way, this saves the output of each sequence)
testseqsallbold_parents <- vector('list')
for (i in 135:length(testseqsallbold)) {
    print(i)
    boldoutput <- bold_identify_parents(testseqsallbold[i], wide = TRUE)
    testseqsallbold_parents <- c(testseqsallbold_parents, boldoutput)
}
# if the loop times out, you can change the 1 in 1:length(testseqsallbold) to the last sequence number that was running when the loop timed out. 

# For example, if the loop times out at 82, restart the loop at 82
# [1] 80
# [1] 81
# [1] 82
# Error in curl::curl_fetch_memory(x$url$url, handle = x$url$handle) : 
#     Timeout was reached: [v4.boldsystems.org] Connection timed out after 10004 milliseconds

# restart at 82
# for (i in 82:length(testseqsallbold)) {
#     print(i)
#     boldoutput <- bold_identify_parents(testseqsallbold[i], wide = TRUE)
#     testseqsallbold_parents <- c(testseqsallbold_parents, boldoutput)
# }

# I have also found that BOLD cuts off access after a certain amount of time (10 mins?), but i can continue the loop by signing in with a vpn from a different city. 


#### convert list to dataframe
testseqsallbold_parents.df <- ldply(testseqsallbold_parents, data.frame)

# testseqsallbold_parents_docall.df <- do.call(rbind, testseqsallbold_parents) # this does not work because "Error in rbind(deparse.level, ...) : numbers of columns of arguments do not match"

#### filter for similarity threshold (optional)
testseqsallbold_parents_97_to_100.df <- testseqsallbold_parents.df %>% 
    filter(similarity>=0.97) %>% 
    arrange(.id, desc(similarity))

#### write dataframe to disk
write_tsv(testseqsallbold_parents_97_to_100.df, "kelpietest_20191210_derep_97_to_100.tsv")
write_tsv(testseqsallbold_parents.df, "kelpietest_20191210_derep.tsv")
