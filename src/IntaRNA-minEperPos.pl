#!/usr/bin/env perl

############################################################################
#
# IntaRNA-minEperPos
#
# Given an IntaRNA output with output configuration
#  --outMode=C --outOverlap=B --outCsvCols=id1,id2,E,hybridDBfull
# this script identifies for each position of the two sequences the minimal 
# energy of any structure that shows a base pairing for the respective position
# 
# Input is read from STDIN
#
# Output is written in semicolon-separated table format including header to 
# STDOUT
#
############################################################################


my %name2col; # mapping of required column names to their respective number
my @neededCols = ("id1", "id2", "E", "hybridDBfull");

my %id2minE; # for each parsed id the minE profile
my $maxLength = 0;

my $minEvalue = 0;

# read STDIN
while( my $line = <> ) {
	
	# trim line
	$line =~ s/^\s+|\s+$//g;
	
	# skip empty lines
	if ($line eq "") { next; }
	
	# parse header if not done yet
	if ( (scalar keys %name2col) == 0 ) {
		# split header into col names
		my @header = split( /;/, $line );
		# find column numbers for the column indices of interest
		for my $name ( @neededCols ) {
			# find column index for this name 
			my $index = 0;
			++$index until  $index > $#header or $header[$index] eq $name;
			# check if name is present in the header information
			if ($index > $#header) {die("\nERROR: column name '$name' is missing in header information!\n\n"); }
			# store column index for this name
			$name2col{$name} = $index;
		}
		# skip further processing of this line
		next;
	}
	
	# parse data line
	my @data = split( /;/, $line );
	my $id1 = $data[ $name2col{"id1"} ];
	my $id2 = $data[ $name2col{"id2"} ];
	my $E = $data[ $name2col{"E"} ];
	my $hybridDBfull = $data[ $name2col{"hybridDBfull"} ];
	$hybridDBfull =~ s/\d//g;
	my @bps = split( /&/, $hybridDBfull );
	
	# check if ids are known; otherwise initialize with 0-profile
	if ( ! exists $id2minE{$id1} ) { my @nullProfile = ($minEvalue)x(length($bps[0])); $id2minE{$id1} = [@nullProfile]; };
	if ( ! exists $id2minE{$id2} ) { my @nullProfile = ($minEvalue)x(length($bps[1])); $id2minE{$id2} = [@nullProfile]; };

	# update data for id1
	my @paired1 = split( //, $bps[0] );
	for (my $i=0; $i <= $#paired1; $i++) {
		if ( $paired1[$i] eq '|' and $E < $id2minE{$id1}[$i] ) { 
			@{$id2minE{$id1}}[$i] = $E;
		}
	}
	# update data for id2
	my @paired2 = split( //, $bps[1] );
	for (my $i=0; $i <= $#paired2; $i++) {
		if ( $paired2[$i] eq '|' and $E < $id2minE{$id2}[$i] ) { 
			@{$id2minE{$id2}}[$i] = $E;
		}
	}
	
}

# get maximal sequence length among all to enable well-formatted table output 
$maxLength = 0;
for my $id (keys %id2minE) {
	@profile = @{$id2minE{$id}};
	if ($maxLength < (1+$#profile)) {
		$maxLength = (1+$#profile);
	}
}

#print header
print "id";
for( my $i=0; $i<$maxLength; $i++) {
	print ";p".($i+1);
}
print "\n";

# print data
for my $id (keys %id2minE) {
	@profile = @{$id2minE{$id}};
    print "$id;".join(";",@profile).((";")x($maxLength-$#profile-1))."\n";
}

