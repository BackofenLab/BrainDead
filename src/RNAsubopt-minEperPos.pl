#!/usr/bin/env perl

############################################################################
#
# RNAsubopt-minEperPos
#
# Given an RNAsubopt output 
# this script identifies for each position the minimal 
# energy of any structure that shows a base pairing for the respective position
# 
# Input is read from STDIN
#
# Output is written in semicolon-separated table format excluding header to 
# STDOUT
#
############################################################################


my @minE = (); # the minE profile

my $minEvalue = 0;

# read STDIN
while( my $line = <> ) {
	
	# trim line
	$line =~ s/^\s+|\s+$//g;
	
	# skip empty lines
	if ($line eq "") { next; }
	
	# parse data line
	my @data = split( /\s+/, $line );
	if ($#data != 1) {next;}
	my $E = $data[ 1 ];
	my @bps = split( //, $data[0] );
	
	# check if ids are known; otherwise initialize with 0-profile
	if ( scalar @minE == 0 ) { @minE = ($minEvalue)x(length($data[0]));}

	# update data 
	for (my $i=0; $i <= $#bps; $i++) {
		if ( (! ($bps[$i] eq '.')) and $E < $minE[$i] ) { 
			$minE[$i] = $E;
		}
	}
	
}


# print data
print ";".join(";",@minE)."\n";

