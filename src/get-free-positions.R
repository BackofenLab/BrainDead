
# clear data structures
rm(list=ls())

DATASET = "../data/180706"
DATASET = "../data/180903"

maxEtoGrep = -1; # maximal energy of an mfe to be considered stable
maxEtoGrep = -3; # maximal energy of an mfe to be considered stable

# read data
seq = read.table(paste(DATASET,"seq.csv",sep="."), sep=";", header=T, row.names=1,colClasses="character")
mfeIntra = read.table(paste(DATASET,"RNAsubopt.intraMinE.csv",sep="."), sep=";", header=T, row.names=1)
mfeHomo = read.table(paste(DATASET,"IntaRNA.homoMinE.csv",sep="."), sep=";", header=T, row.names=1)

miRNAs = row.names(seq);   miRNA = miRNAs[1]


##############################################
# compile "free of intra-mol"
##############################################

# init data
freeIntra = mfeIntra > maxEtoGrep;
colnames(freeIntra) = paste("freeIntra",1:ncol(freeIntra),sep="")
# write data
write.table( cbind(miRNAs,freeIntra), paste(DATASET,"free-intra",maxEtoGrep,"csv",sep=".")
	, sep=";", col.names=c("ID",colnames(freeIntra))
	, quote=F, row.names=F
)

##############################################
# compile "free of homo-duplex"
##############################################

# init data
freeHomo = mfeHomo > maxEtoGrep;
colnames(freeHomo) = paste("freeHomo",1:ncol(freeHomo),sep="")
# write data
write.table( cbind(miRNAs,freeHomo), paste(DATASET,"free-homo",maxEtoGrep,"csv",sep=".")
	, sep=";", col.names=c("ID",colnames(freeHomo))
	, quote=F, row.names=F
)





