
# clear data structures
rm(list=ls())

DATASET = "../data/180706"
DATASET = "../data/180903"

OUTFILESUFFIX = "structuredness"

RT = 0.61
maxEtoPlot = -3;

plotNoED = FALSE;

seq = read.table(paste(DATASET,"seq.csv",sep="."), sep=";", header=T, row.names=1,colClasses="character")
Pu = read.table(paste(DATASET,"RNAplfold.Pu.csv",sep="."), sep=";", header=T, row.names=1)
Ps = read.table(paste(DATASET,"IntaRNA.Ps.csv",sep="."), sep=";", header=T, row.names=1)
mfeIntra = read.table(paste(DATASET,"RNAfold.csv",sep="."), sep=";", header=T, row.names=1)
mfeHomo = read.table(paste(DATASET,"IntaRNA.csv",sep="."), sep=";", header=T, row.names=1)
mfeHomoNoEd = read.table(paste(DATASET,"IntaRNA.noED.csv",sep="."), sep=";", header=T, row.names=1)
meta = read.table(paste(DATASET,"csv",sep="."), sep=";", header=T, row.names=1)
minEintra = read.table(paste(DATASET,"RNAsubopt.intraMinE.csv",sep="."), sep=";", header=T, row.names=1)
minEhomo = read.table(paste(DATASET,"IntaRNA.homoMinE.csv",sep="."), sep=";", header=T, row.names=1)



# get index vectors
miRNAs = row.names(seq);   miRNA = miRNAs[1]
# sanity checks
stopifnot( miRNAs == row.names(meta) );
stopifnot( miRNAs == row.names(Pu) );
stopifnot( miRNAs == row.names(Ps) );
stopifnot( miRNAs == row.names(mfeIntra) );
stopifnot( miRNAs == row.names(mfeHomo) );
stopifnot( miRNAs == row.names(mfeHomoNoEd) );
stopifnot( miRNAs == row.names(minEintra ) );
stopifnot( miRNAs == row.names(minEhomo ) );

# create output file
pdf(paste(DATASET,OUTFILESUFFIX,maxEtoPlot,"pdf",sep="."), onefile=T, width=10, height=6);

# overall energy minimum
Ebounds = c(maxEtoPlot-2, 0);
if (plotNoED) {
 EboundsAll = c( min(cbind(mfeHomo[,"homo.mfe"], mfeIntra[,"intra.mfe.kcal"], mfeHomoNoEd[,"homo.mfe"], (maxEtoPlot-2) )), 0)
 Ebounds = EboundsAll;
}else {
 EboundsAll_excludeNoED = c( min(cbind(mfeHomo[,"homo.mfe"], mfeIntra[,"intra.mfe.kcal"], (maxEtoPlot-2) )), 0)
 Ebounds = EboundsAll_excludeNoED ;
}

# process each miRNA of the DATASET
miRNA = miRNAs[1]
for (miRNA in miRNAs) {

# get data
Eintra = mfeIntra[miRNA,"intra.mfe.kcal"]
EintraBP = unlist(strsplit(as.character(mfeIntra[miRNA,"intra.mfe.bp"]), split=""))
Ehomo = mfeHomo[miRNA,"homo.mfe"]
EhomoS1 = mfeHomo[miRNA,"start1"]
EhomoE1 = mfeHomo[miRNA,"end1"]
EhomoBP1 = unlist(strsplit(as.character(mfeHomo[miRNA,"hybridDP"]), split=""))[1:(EhomoE1-EhomoS1+1)]
EhomoS2 = mfeHomo[miRNA,"start2"]
EhomoE2 = mfeHomo[miRNA,"end2"]
EhomoBP2 = unlist(strsplit(as.character(mfeHomo[miRNA,"hybridDP"]), split=""))[(EhomoE1-EhomoS1+2)+(1:(EhomoE2-EhomoS2+1))]
if (plotNoED) {
 EhomoNoEd = mfeHomoNoEd[miRNA,"homo.mfe"]
 EhomoNoEdS1 = mfeHomoNoEd[miRNA,"start1"]
 EhomoNoEdE1 = mfeHomoNoEd[miRNA,"end1"]
 EhomoNoEdBP1 = unlist(strsplit(as.character(mfeHomoNoEd[miRNA,"hybridDP"]), split=""))[1:(EhomoNoEdE1-EhomoNoEdS1+1)]
 EhomoNoEdS2 = mfeHomoNoEd[miRNA,"start2"]
 EhomoNoEdE2 = mfeHomoNoEd[miRNA,"end2"]
 EhomoNoEdBP2 = unlist(strsplit(as.character(mfeHomoNoEd[miRNA,"hybridDP"]), split=""))[(EhomoNoEdE1-EhomoNoEdS1+2)+(1:(EhomoNoEdE2-EhomoNoEdS2+1))]
}
skipcols = !is.na(seq[miRNA,]);
EDintra = na.omit(RT * log(Pu[miRNA,skipcols ]))
EDspot = na.omit(RT * log(1-Ps[miRNA,skipcols ]))
miRNAseq = as.character(seq[miRNA,skipcols])
stopifnot( length(EDintra) == length(miRNAseq ) );
stopifnot( length(EDspot) == length(miRNAseq ) );

# plot it
#Ebounds = c( min(cbind(EDintra,EDspot,Eintra,Ehomo,EhomoNoEd, (maxEtoPlot-2) )), 0)
xval = 1:sum(skipcols)
colIntra = "red"
colHomo = colors()[97]#[31]
colHomoNoEd = "blue"; #"cyan"
plot( x=xval,y=(1-Pu[miRNA,skipcols]), col="white", type="l", lwd=1, xaxt='n', ylim=Ebounds, xlab="", ylab="energy in kcal/mol", 
	main=paste(miRNA," activation =",meta[miRNA,"TNF.alpha.release"]))
axis(1,at=xval,labels=miRNAseq);
#legend("topleft", fill=c(colIntra ,colHomo), legend=c("intra-molecular (1-Pu)","homo-duplex prob"));
if (plotNoED) {
 legend("bottomleft", col=c(colIntra ,colHomo,colHomoNoEd,"black"), lty=3, lwd=3, legend=c("mfe bp intra","mfe bp homo","mfe bp homo noED","max E to be stable"));
} else {
 legend("bottomleft", col=c(colIntra ,colHomo,"black"), lty=3, lwd=3, legend=c("mfe bp intra","mfe bp homo","max E to be stable"));
}
abline( h=0)
# annotate U/G positions
points( x=xval[ miRNAseq == "U" ], y=rep(0, sum(miRNAseq == "U")), cex=2, pch='U', bg="black")
points( x=xval[ miRNAseq == "G" ], y=rep(0, sum(miRNAseq == "G")), cex=2, pch='G', bg="black")

# max E boundary
abline( h=maxEtoPlot, lty=2, lwd=3)

# plot minE profiles
lines( x=1:ncol(minEintra), y = minEintra[miRNA,], lwd=2, type="l", col=colIntra );
lines( x=1:ncol(minEhomo), y = minEhomo[miRNA,], lwd=2, type="l", col=colHomo);

# mfe intra
#abline( h=Eintra, col=colIntra , lwd=2)
lines( x=xval[ EintraBP == "(" ], y= rep(Eintra, sum(EintraBP == "(")), lwd=2, type="b", pch='(', col=colIntra )
lines( x=xval[ EintraBP == ")" ], y= rep(Eintra, sum(EintraBP == ")")), lwd=2, cex=1, pch=')', type="b", col=colIntra )
if (Eintra <= maxEtoPlot) {
#  lines( x=xval,y=EDintra, col=colIntra , lwd=2, lty=2)
#  lines( x=xval,y=(1-Pu[miRNA,skipcols]), col=colIntra , lwd=2)
  points( x=xval[ EintraBP == "(" ]-0.1, y= rep(0, sum(EintraBP == "(")), lwd=3, pch=16, cex=1,col=colIntra )
  points( x=xval[ EintraBP == ")" ]-0.1, y= rep(0, sum(EintraBP == ")")), lwd=3, pch=16, cex=1,col=colIntra )
}

if (plotNoED) {
 # mfe homo noED
 lines( x=((EhomoNoEdS1:EhomoNoEdE1)[EhomoNoEdBP1 != "."]), y=rep(EhomoNoEd, sum(EhomoNoEdBP1 != ".")), type='p', lwd=2, col=colHomoNoEd,  pch='(')
 lines( x=((EhomoNoEdS2:EhomoNoEdE2)[EhomoNoEdBP2 != "."]), y=rep(EhomoNoEd, sum(EhomoNoEdBP2 != ".")), type='p', lwd=2, col=colHomoNoEd,  pch=')')
 if (EhomoNoEd <= maxEtoPlot) {
   posBoth = intersect(((EhomoNoEdS1:EhomoNoEdE1)[EhomoNoEdBP1 != "."]),((EhomoNoEdS2:EhomoNoEdE2)[EhomoNoEdBP2 != "."]))
   points( x=posBoth+0.1, y=rep(0, length(posBoth )), cex=1.2, col=colHomoNoEd,  pch=16)
   points( x=posBoth+0.1, y=rep(0, length(posBoth )), cex=1.2, col=colHomoNoEd,  pch=16)
 }
}

# mfe homo
lines( x=((EhomoS1:EhomoE1)[EhomoBP1 != "."]), y=rep(Ehomo, sum(EhomoBP1 != ".")), type='p', lwd=2, col=colHomo,  pch='(')
lines( x=((EhomoS2:EhomoE2)[EhomoBP2 != "."]), y=rep(Ehomo, sum(EhomoBP2 != ".")), type='p', lwd=2, col=colHomo,  pch=')')
if (Ehomo <= maxEtoPlot) {
#  lines( x=xval,y=EDspot, col=colHomo, lwd=2, lty=2)
#  lines( x=xval,y=Ps[miRNA,skipcols], col=colHomo, lwd=2)
  posBoth = intersect(((EhomoS1:EhomoE1)[EhomoBP1 != "." ]), ((EhomoS2:EhomoE2)[EhomoBP2 != "."]))
  points( x=posBoth +0.1, y=rep(0, length(posBoth )), cex=1, col=colHomo,  pch=9)
  points( x=posBoth +0.1, y=rep(0, length(posBoth )), cex=1, col=colHomo,  pch=9)
}


} # end for miRNAs
# finish output file
dev.off();



# analyses

# overlaps
for ( maxE in c(-1,-2,-3)) {
  print(paste("intersect for intra/homo mfe < ",maxE ,"=",length(intersect(  miRNAs[mfeHomo$homo.mfe < maxE ], miRNAs[mfeIntra$intra.mfe.kcal < maxE ])) / max(c(sum(mfeHomo$homo.mfe < maxE ),sum(mfeIntra$intra.mfe.kcal < maxE )))*100,"%"))
}

actYes = as.vector(na.omit(miRNAs[meta$TNF.alpha.release > 0]))
actNon = as.vector(na.omit(miRNAs[meta$TNF.alpha.release == 0]))
for ( maxE in c(-1,-2,-3)) {
  print(paste("intersect activation with intra mfe < ",maxE ,"=",length(intersect( actYes , miRNAs[mfeIntra$intra.mfe.kcal < maxE ])) / length(actYes)*100,"%" ))
  print(paste("intersect non-activation with intra mfe < ",maxE ,"=",length(intersect( actNon , miRNAs[mfeIntra$intra.mfe.kcal < maxE ])) / length(actNon )*100,"%" ))
}
for ( maxE in c(-1,-2,-3)) {
  print(paste("intersect activation with homo mfe < ",maxE ,"=",length(intersect( actYes , miRNAs[mfeHomo$homo.mfe < maxE ])) / length(actYes)*100,"%" ))
  print(paste("intersect non-activation with homo mfe < ",maxE ,"=",length(intersect( actNon , miRNAs[mfeHomo$homo.mfe < maxE ])) / length(actNon )*100,"%" ))
}

