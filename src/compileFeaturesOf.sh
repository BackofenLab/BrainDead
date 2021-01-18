#!/usr/bin/env bash

####################################################
# creates features for a given miRNA sequence (read from first argument)
####################################################

CS='\t'

colsRNAfold="${CS}intra-mfe-bp${CS}intra-mfe-kcal"

#colsIntaRNAcall="E,start1,end1,start2,end2,hybridDP"
#colsIntaRNA="${CS}homo-mfe${CS}start1${CS}end1${CS}start2${CS}end2${CS}hybridDP"
colsIntaRNAcall="E,hybridDBfull"
colsIntaRNA="${CS}homo-mfe${CS}hybridDBfull"

cols2print="id${colsRNAfold}${colsIntaRNA}"


#############################################
# input validation
#############################################

# get miRNA sequence from first argument
miRNA=$1
miID=$2

# handle empty input to print column names only
if [ -z "$miRNA" ]; then
	echo -e "$cols2print";
	exit 0;
fi
if [ -z "$miID" ]; then
	echo "Error: no miRNA ID/name provided (second argument)" >&2;
	exit 1;
	
fi
if [[ ! $miRNA =~ ^[ACGUTacgut]+$ ]]; then
	echo "Error: given sequence is no RNA sequence [ACGUT]+" >&2;
	exit 1;
fi 

# print ID of current miRNA
printf $miID;

# get directory of this script (to call subscripts)
BINDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

#############################################
# RNAfold-based features
#############################################

# check if features are to be computed
if [ ! -z "$colsRNAfold" ]; then

curOut=$(
printf "${CS}";
OutRNAfold=$(echo $miRNA | RNAfold --noPS | tr "\n" " " | awk '{$1="";printf $0}')
echo $OutRNAfold | awk '{printf $1}';
printf "${CS}";
echo $OutRNAfold | awk '{$1=""; printf $0}' | tr -d "() \t\n") 

printf "$curOut";

fi


#############################################
# IntaRNA-based features
#############################################

# check if features are to be computed
if [ ! -z "$colsIntaRNA" ]; then

IntaRNA -t $miRNA -q $miRNA -n 1000 --outOverlap=B --out=${miID}.IntaRNA.csv --outMode=C --outCsvCols=id1,id2,${colsIntaRNAcall} --noSeed -m M --outMaxE=10
# get mfe information
OutIntaRNA=$(awk '{if(NR==2){printf $0}}' ${miID}.IntaRNA.csv | cut -d ";" -f 3- | tr ";" "${CS}" ) 
printf "${CS}${OutIntaRNA}";

# collect mfe per position
cat $miID.IntaRNA.csv | perl $BINDIR/IntaRNA-minEperPos.pl > $miID.homoMinE

fi

#############################################

# finalize output
printf "\n";
