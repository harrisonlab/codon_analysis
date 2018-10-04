args = commandArgs(trailingOnly=TRUE)

# check for arguments
if (length(args) < 3) {
  stop("**************************************** \n
       Three arguments must be supplied \n
       Argument1: Fasta formatted CDS sequences \n
       Argument2: A tabulated FPKM values per condition \n
       Argument3: A list of codons n\
       **************************************** \n", call.=FALSE)
}

# Load libraries
library(Biostrings)
library(SADEG, lib.loc = "/home/bonthas/R/libs/")
source("functions2.R")

# Load files
mygenes = readDNAStringSet(args[1])
myfpkm <- read.table(args[2])
codonw = read.delim(args[3],header = FALSE,sep="\n",quote="",stringsAsFactors = FALSE)

# Calculate average FPKM per condition
RH <- cbind(myfpkm[,c(1,2)],
      rowMeans(myfpkm[,c(3:5)]),
      rowMeans(myfpkm[,c(6:8)]),
      rowMeans(myfpkm[,c(9:11)]),
      rowMeans(myfpkm[,c(12:14)]),
      rowMeans(myfpkm[,c(15:17)]),
      rowMeans(myfpkm[,c(18:20)]),
      rowMeans(myfpkm[,c(21:23)]),
      rowMeans(myfpkm[,c(24:26)])
)
colnames(RH) <- c("GeneID","Length","RH1_mean","RH2_mean","RH3_mean","RH4_mean","RH5_mean","RH6_mean","RH7_mean","RH8_mean")

# Extract IDs from FPKM table
ids <- extractID(RH)

# Extract CDS sequences
Seq <- extractSeq(ids,mygenes)

# Calculate RSCU values
matH <- lapply(Seq[[1]], function(x) {lapply(x,SADEG.RSCU)})

matL <- lapply(Seq[[2]], function(x) {lapply(x,SADEG.RSCU)})

# Process the extracted RSCU values and calculate deltaRSCU for both highly and weakly expressing genes.
processMat2(matH=matH,matL=matL,codonw=codonw)

# End of program
