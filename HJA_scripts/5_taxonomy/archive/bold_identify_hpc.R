# Usage
# Rscript --vanilla bold_identify_hpc.R "CO1_1sequence_perBIN_040915_COIspiking.fas"  # i think i need to quote this


## START OF PROGRAM
# set some general options
options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE) # if you want to pass arguments to R from the shell (bash) command line
print(args)

# test if there is one argument: if not, return an error
if (length(args) < 1) {
    stop("At least one filename argument must be supplied", call.=FALSE)
}

## load packages
# devtools::install_github("ropensci/bold") # for the dev version
library("bold") # now at 0.5.0.  use ≥0.4.1.9300 If not avail from cran, use github version
library("plyr") # always load plyr before dplyr
library("seqinr") # used to read fasta file
library("dplyr")
library("readr")

## bold_identify multiple sequences from a fasta file ####
## 
# args[1] # this is the filename passed to the Rscript

#### load a fasta file
testseqsall <- read.fasta(file = args[1], seqtype = "DNA", set.attributes = FALSE, as.string = TRUE, forceDNAtolower = FALSE)  # read in fasta file using seqinr

#### bold_identify the sequences
testseqsallbold <- bold_identify(testseqsall, db = "COX1")
testseqsallbold <- testseqsallbold[!sapply(testseqsallbold, is.null)] # to remove any null elements

#### add parent taxa to output, wide format
testseqsallbold_parents <- bold_identify_parents(testseqsallbold, wide = TRUE)

#### convert list to dataframe
# testseqsallbold_parents.df <- do.call(rbind, testseqsallbold_parents) # i think this should replace ldply
testseqsallbold_parents.df <- ldply(testseqsallbold_parents, data.frame)

#### filter for similarity threshold (optional)
testseqsallbold_parents_96_to_100.df <- dplyr::filter(testseqsallbold_parents.df, similarity>=0.96)
testseqsallbold_parents_96_to_100.df <- dplyr::arrange(testseqsallbold_parents_96_to_100.df, .id, desc(similarity))

#### write dataframe to disk
write_tsv(testseqsallbold_parents_96_to_100.df, "CO1_1sequence_perBIN_040915_COIspiking_96_to_100.tsv")
write_tsv(testseqsallbold_parents.df, "CO1_1sequence_perBIN_040915_COIspiking.tsv")


# #### bold_identify one seq at a time.  
# testseqsallbold <- vector('list')
# for (i in 1:length(testseqsall)) {
#     print(i)
#     boldoutput <- bold_identify(testseqsall[i], db = "COX1")
#     testseqsallbold <- c(testseqsallbold, boldoutput)
# }
