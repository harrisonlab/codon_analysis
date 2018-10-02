args = commandArgs(trailingOnly=TRUE)

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
#print(head(mygenes))
myfpkm <- read.table(args[2])
#print(head(myfpkm))
codonw = read.delim(args[3],header = FALSE,sep="\n",quote="",stringsAsFactors = FALSE)
#print(head(codonw))

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

#print(head(RH))
print(dim(RH))

ids <- extractID(RH)

#print(head(ids))
print(length(ids))

# Extract CDS sequences
Seq <- extractSeq(ids,mygenes)
print("++++++++++++++++++++++++++++++++++++")
print(length(Seq))
print(length(Seq[[1]]))
print(length(Seq[[2]]))
print(head(Seq[[1]][[3]]))
print("++++++++++++++++++++++++++++++++++++")

matH <- lapply(Seq[[1]], function(x) {lapply(x,SADEG.RSCU)})

matL <- lapply(Seq[[2]], function(x) {lapply(x,SADEG.RSCU)})

print("======================================")

#print(head(matH))
#print(head(matL))

print(length(matH))
print(length(matL))
print("======================================")

#codon <- names(matH[[3]][[1]])

#print(codon)

processMat2(matH=matH,matL=matL,codonw=codonw)
