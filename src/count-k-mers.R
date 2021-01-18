
# clear data structures
rm(list=ls())

library(kmer)
require(kmer)

DATASET = "../data/180706"
DATASET = "../data/180903"

OUTFILESUFFIX = "k-mers"

# length of k-mers to check
kmers = c(1,2,3,4);    k = kmers[1]

# load data
data = read.table(paste(DATASET,"csv",sep="."), sep=";", header=T, row.names=1)
seq = read.table(paste(DATASET,"seq.csv",sep="."), sep=";", header=T, row.names=1,colClasses="character")
# get index vector
miRNAs = row.names(data);   miRNA = miRNAs[1]

for ( k in kmers ) {

# init data
kmerCount = c();

# iterate all
for ( miRNA in miRNAs ) {

kmerCount = rbind( kmerCount, kcount(unlist(seq[miRNA,]), k, residues="RNA") )

} # for all miRNAs

rownames(kmerCount) = miRNAs;
rowsum <- apply(kmerCount, 1, sum)

# write table with absolute counts
write.table( cbind( miRNAs, kmerCount)
 , paste(DATASET,OUTFILESUFFIX,k,"abs","csv",sep="."), sep=";", quote=F, row.names=F
 , col.names=c("ID",paste("abs",colnames(kmerCount),sep=".")))

# write table with rounded values
write.table( cbind( miRNAs, round(kmerCount/rowsum, 4))
 , paste(DATASET,OUTFILESUFFIX,k,"rel","csv",sep="."), sep=";", quote=F, row.names=F
 , col.names=c("ID",paste("rel",colnames(kmerCount),sep=".")))


} # for all kmers




