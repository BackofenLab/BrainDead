#!/usr/bin/env bash

####################################################
# creates features for a given set of miRNA sequences 
# (read file provided as first argument)
# input file format: semicolon-separated table (with column names = row1),
# col1=ID and col2=sequence
#
# the final table is written to STDOUT
####################################################

# get miRNA sequence from first argument
inputFile=$1

# handle empty input to print column names only
if [ ! -f "$inputFile" ]; then
	echo "Error: can not access file '${inputFile}'" >&2;
	exit 1;
fi

# time stamp of computation
date

#############################################
# dependency checks
#############################################

# test if required tools are available
function checkProgram() {
	binPath=$(command -v $1)
	if ! [ -x "$binPath" ]; then
	  echo "Error: $1 is not installed (or available via PATH variable)." >&2
	  exit 1
	fi
	echo $binPath
}

checkProgram IntaRNA
checkProgram RNAfold
checkProgram RNAsubopt
checkProgram R

#############################################

outMfe=RNAfold.IntaRNA.mfe.tsv

outHomoMinE=IntaRNA.homoMinE.tsv
#outHomoMinE=

outIntraMinE=RNAsubopt.intraMinE.tsv
#outIntraMinE=

keepTempFiles=yes
keepTempFiles=

# get directory of this script (to call subscripts)
BINDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# column separator used for the output
CS='\t'

# get list of miRNA ids in order of input file (to be preserved in output)
miRNAs=$(cut -d ";" -f 1 $inputFile | awk 'NR>1 {print}' | tr "\n" " ")

cut -d ";" -f 1 $inputFile > ids.tsv
echo "ids.tsv"

######################################################################
# takes a table and expands it to equal number of columns per line
# arg1 = filename
# arg2 = col delimiter
# arg3 = new colnames to be suffixed with column number
######################################################################
function expandColumns() {
	# input arguments
	file=$1
	colsep=$2
	colname=$3
	# get number of columns
	COLNUM=$(awk -F "${colsep}" 'BEGIN{m=0}{if(NF>m){m=NF}}END{print m}' ${file})

	(# write header
	for i in `seq 1 $COLNUM`; do printf "${colsep}${colname}$i"; done
	printf "\n";
	# write data rows
	cat ${file} | awk -F "${colsep}" -v m=$COLNUM '{for (j=1; j<=NF; j++){printf FS$j}; for(j=(NF+1);j<=m;j++){printf(FS"NA");}; printf("\n");}';
	) > ${file}.expandColumns
	# overwrite original file
	mv -f ${file}.expandColumns ${file}
}


###################################################
# runs a given R script file
# arg1 = scriptFilePath
###################################################
function runRscript() {
	filePathR=$1
	# handle windows cygwin path setup
	if [ -x "$(command -v cygpath)" ]; then
		filePathR=$(cygpath -w $filePathR)
	fi
	#	Rscript --vanilla $filePathR
	Rscript --no-save $filePathR
}

if [ 1 == 2 ]; then
######################################################################
# create standard features
######################################################################

# call compileFeaturesOf.sh for each miRNA
awk -F ";" -v BINDIR=$BINDIR '{if(NR>1){system(BINDIR"/compileFeaturesOf.sh "$2" "$1" > "$1".features.tsv")}}' $inputFile;

# aggregate feature output files to table
( # create header information 
$BINDIR/compileFeaturesOf.sh;
awk -F ";" -v BINDIR=$BINDIR '{if(NR>1){system("cat "$1".features.tsv")}}' $inputFile;
) > ${outMfe}
# remove tmp files
[[ ! -z "$keepTempFiles" ]] || rm -f *.features.tsv

echo "${outMfe}"

######################################################################
# compile mfe profile for homo-dimerization if needed
######################################################################

if [ ! -z "$outHomoMinE" ]; then
	# aggregate homoMinE table
	awk -F ";" '{if(NR>1){system("tail -n 1 "$1".homoMinE")}}' $inputFile | cut -d ";" -f 2- | tr ";" "${CS}" > ${outHomoMinE}
	# ensure equal number of columns in table
	expandColumns ${outHomoMinE} "${CS}" "homoMinE"
	
	echo "${outHomoMinE}"
fi

[[ ! -z "$keepTempFiles" ]] || rm -f *.homoMinE
[[ ! -z "$keepTempFiles" ]] || rm -f *.IntaRNA.csv

######################################################################
# compile mfe profile of intra-molecular structure
######################################################################

if [ ! -z "$outIntraMinE" ]; then
	# intra-molecular structures
	( # generate data
	for i in $miRNAs; do
	  printf $i;
	  # get sequence 
	  SEQ=`grep "$i;" -m 1 $inputFile | awk -F ";" '{print $2}'`
	  # call RNAsubopt
		echo $SEQ | RNAsubopt -e 100 | perl $BINDIR/RNAsubopt-minEperPos.pl;
	done ) | cut -d ";" -f 2- | tr ";" "${CS}" > $outIntraMinE
	# ensure equal number of columns in table
	expandColumns ${outIntraMinE} "${CS}" "intraMinE"
	
	echo "${outIntraMinE}"
fi

######################################################################
# identify unstructured positions
######################################################################

runRscript $BINDIR/compile-free-positions.R

echo "free-homo.*.tsv"
echo "free-intra.*.tsv"

######################################################################
# get k-mer counts
######################################################################

fi
# decompose sequences into single-char columns
cut -d ";" -f 2 $inputFile | awk '{if(NR>1){print $0}}' | sed "s/\(.\)/;\1/g" | cut -d ";" -f 2- | tr ";" "${CS}" > seq.tsv
expandColumns seq.tsv "${CS}" "S"

#runRscript $BINDIR/compile-k-mers.mfe-free.R

echo "k-mers.*"

[[ ! -z "$keepTempFiles" ]] || rm -f seq.tsv

