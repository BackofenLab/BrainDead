
# clear data structures
rm(list=ls())

DATASET = "../data/180706"
DATASET = "../data/180903"

maxEtoGrep = -1; # maximal energy of an mfe to be considered stable
maxEtoGrep = -3; # maximal energy of an mfe to be considered stable

# read data
seq = read.table(paste(DATASET,"seq.csv",sep="."), sep=";", header=T, row.names=1,colClasses="character")
mfeIntra = read.table(paste(DATASET,"RNAfold.csv",sep="."), sep=";", header=T, row.names=1)
mfeHomoNoEd = read.table(paste(DATASET,"IntaRNA.noED.csv",sep="."), sep=";", header=T, row.names=1)

miRNAs = row.names(seq);   miRNA = miRNAs[1]


##############################################
# compile "free of intra-mol mfe"
##############################################

# init data
freeIntra = seq
freeIntra[!is.na(seq)] = T
# generate data
for ( miRNA in miRNAs ) {
 # check if mfe is to be considered to be blocking
 if (mfeIntra[miRNA,"intra.mfe.kcal"]<=maxEtoGrep) {
	bp = unlist(strsplit(as.character(mfeIntra[miRNA,"intra.mfe.bp"]), split=""))
	freeIntra[miRNA,1:length(bp)]  = (bp == ".")
 }
} # for all miRNAs

# write data
write.table( cbind(miRNAs,freeIntra), paste(DATASET,"free-intra",maxEtoGrep,"csv",sep=".")
	, sep=";", col.names=c("ID",colnames(freeIntra))
	, quote=F, row.names=F
)

##############################################
# compile "free of homo-mol mfe (no ED)"
##############################################

# init data
freeHomo = seq
freeHomo[!is.na(seq)] = T
# generate data
for ( miRNA in miRNAs ) {
 # check if mfe is to be considered to be blocking
 if (mfeHomoNoEd[miRNA,"homo.mfe"]<=maxEtoGrep) {
	EhomoNoEdS1 = mfeHomoNoEd[miRNA,"start1"]
	EhomoNoEdE1 = mfeHomoNoEd[miRNA,"end1"]
	EhomoNoEdBP1 = unlist(strsplit(as.character(mfeHomoNoEd[miRNA,"hybridDP"]), split=""))[1:(EhomoNoEdE1-EhomoNoEdS1+1)]
	EhomoNoEdS2 = mfeHomoNoEd[miRNA,"start2"]
	EhomoNoEdE2 = mfeHomoNoEd[miRNA,"end2"]
	EhomoNoEdBP2 = unlist(strsplit(as.character(mfeHomoNoEd[miRNA,"hybridDP"]), split=""))[(EhomoNoEdE1-EhomoNoEdS1+2)+(1:(EhomoNoEdE2-EhomoNoEdS2+1))]
	free1 = freeHomo[miRNA,];
	free1[EhomoNoEdS1:EhomoNoEdE1] = (EhomoNoEdBP1 == ".")
	free2 = freeHomo[miRNA,];
	free2[EhomoNoEdS2:EhomoNoEdE2] = (EhomoNoEdBP2 == ".")
	noNA = !is.na(free1)
	freeHomo[miRNA,noNA ]  = (as.logical(free1[noNA]) | as.logical(free2[noNA]))
 }
} # for all miRNAs

# write data
write.table( cbind(miRNAs,freeHomo), paste(DATASET,"free-homo-noED",maxEtoGrep,"csv",sep=".")
	, sep=";", col.names=c("ID",colnames(freeHomo))
	, quote=F, row.names=F
)





