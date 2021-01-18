
# clear data structures
rm(list=ls())

library(kmer)

OUTFILESUFFIX = "k-mers"

# length of k-mers to check
kmers = c(1,2,3,4);    k = kmers[1]

maxEtoConsider = -1;
maxEtoConsider = -3;

colSep = "\t";

# load data
seq = read.table("seq.tsv", sep="\t", header=T,colClasses="character")
seq = seq[,2:ncol(seq)]
miRNAs= 1:nrow(seq)

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
write.table( cbind( "", kmerCount)
 , paste(OUTFILESUFFIX,k,suffix,"abs","tsv",sep="."), sep=colSep, quote=F, row.names=F
 , col.names=c("",paste(suffix,"abs",colnames(kmerCount),sep=".")))


# write table with rounded relative values
write.table( cbind( "", round(kmerCount/rowsum, 4))
 , paste(OUTFILESUFFIX,k,suffix,"rel","tsv",sep="."), sep=colSep, quote=F, row.names=F
 , col.names=c("",paste(suffix,"rel",colnames(kmerCount),sep=".")))


} # for all kmers


} # function countFreeKmers()

#########################################

countFreeKmers( seq != "", "kmer" )

freeIntra = read.table(paste("free-intra",maxEtoConsider,"tsv",sep="."), sep=colSep, header=T)
countFreeKmers( freeIntra[,2:ncol(freeIntra)], paste("kmer.freeIntra",maxEtoConsider,sep="") )

freeHomo = read.table(paste("free-homo",maxEtoConsider,"tsv",sep="."), sep=colSep, header=T)
countFreeKmers( freeHomo[,2:ncol(freeIntra)], paste("kmer.freeHomo",maxEtoConsider,sep="") )
countFreeKmers( freeHomo[,2:ncol(freeIntra)] & freeIntra[,2:ncol(freeIntra)], paste("kmer.free",maxEtoConsider,sep="") )

