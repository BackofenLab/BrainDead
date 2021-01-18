
# clear data structures
rm(list=ls())

library(kmer)
require(kmer)

DATASET = "../data/180706"
DATASET = "../data/180903"

OUTFILESUFFIX = "k-mers"


# length of k-mers to check
kmers = c(1,2,3,4);    k = kmers[1]

maxEtoConsider = -1;
maxEtoConsider = -3;

# load data
seq = read.table(paste(DATASET,"seq.csv",sep="."), sep=";", header=T, row.names=1,colClasses="character")

# get index vector
miRNAs = row.names(seq);   miRNA = miRNAs[1]

######  counting function ###############
countFreeKmers <- function (freePos,suffix) {

for ( k in kmers ) {

# init data
kmerCount = c();

# iterate all
for ( miRNA in miRNAs ) {

s = as.vector(t(seq[miRNA,!is.na(seq[miRNA,])]))

# init counter
sKcount = kcount(unlist(s), k, residues="RNA")
sKcount[,] = 0

# check intra-mol-free subsequences
free = as.vector(t(freePos[miRNA,]))
i = 1; j=0;
while( i<=length(s) ) {
 if( free[i] ) {
  j = i+1;
  while(j<=length(s) & free[j]) { j = j+1 }
  if ( j-i >= k ) {
    sKcount = rbind(sKcount, kcount(unlist(s[i:(j-1)]), k, residues="RNA"));
  }
#  subseq = paste( unlist(s[i:(j-1)]), collapse = '' ); print(paste(i,j,subseq));
  i = j;
 }
 i = i+1;
}

# accumulate normalized results of all subsequences
kmerCount = rbind( kmerCount, apply(sKcount,2, sum) )


} # for all miRNAs

rownames(kmerCount) = miRNAs;
rowsum <- apply(kmerCount, 1, sum)

# write table with absolute counts
write.table( cbind( miRNAs, kmerCount)
 , paste(DATASET,OUTFILESUFFIX,k,suffix,"abs","csv",sep="."), sep=";", quote=F, row.names=F
 , col.names=c("ID",paste(suffix,"abs",colnames(kmerCount),sep=".")))


# write table with rounded relative values
write.table( cbind( miRNAs, round(kmerCount/rowsum, 4))
 , paste(DATASET,OUTFILESUFFIX,k,suffix,"rel","csv",sep="."), sep=";", quote=F, row.names=F
 , col.names=c("ID",paste(suffix,"rel",colnames(kmerCount),sep=".")))


} # for all kmers


} # function countFreeKmers()

#########################################


freeIntra = read.table(paste(DATASET,"free-intra",maxEtoConsider,"csv",sep="."), sep=";", header=T, row.names=1)
countFreeKmers( freeIntra, paste("freeIntra",maxEtoConsider,sep="") )

freeHomo = read.table(paste(DATASET,"free-homo",maxEtoConsider,"csv",sep="."), sep=";", header=T, row.names=1)
countFreeKmers( freeHomo, paste("freeHomo",maxEtoConsider,sep="") )
countFreeKmers( freeHomo & freeIntra, paste("free",maxEtoConsider,sep="") )

freeHomoNoED = read.table(paste(DATASET,"free-homo-noED",maxEtoConsider,"csv",sep="."), sep=";", header=T, row.names=1)
countFreeKmers( freeHomoNoED, paste("freeHomoNoED",maxEtoConsider,sep="") )
# countFreeKmers( freeHomoNoED & freeIntra, paste("free",maxEtoConsider,sep="") )




