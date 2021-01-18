
# clear data structures
rm(list=ls())

maxEtoGrep = -1; # maximal energy of an mfe to be considered stable
maxEtoGrep = -3; # maximal energy of an mfe to be considered stable


colSep = "\t";

# read data
mfeIntra = read.table("RNAsubopt.intraMinE.tsv", sep="\t", header=T)
mfeHomo = read.table("IntaRNA.homoMinE.tsv", sep="\t", header=T)

miRNAs = c();


##############################################
# compile "free of intra-mol"
##############################################

# init data
freeIntra = mfeIntra > maxEtoGrep;
colnames(freeIntra) = c("",paste("freeIntra",1:(ncol(freeIntra)-1),sep=""))
freeIntra[,1] = ""
# write data
write.table( freeIntra, paste("free-intra",maxEtoGrep,"tsv",sep=".")
	, sep=colSep, col.names=colnames(freeIntra)
	, quote=F, row.names=F, na="NA"
)

##############################################
# compile "free of homo-duplex"
##############################################

# init data
freeHomo = mfeHomo > maxEtoGrep;
colnames(freeHomo) = c("",paste("freeHomo",1:(ncol(freeHomo)-1),sep=""))
freeHomo[,1] = ""
# write data
write.table( freeHomo, paste("free-homo",maxEtoGrep,"tsv",sep=".")
	, sep=colSep, col.names=colnames(freeHomo)
	, quote=F, row.names=F, na="NA"
)

