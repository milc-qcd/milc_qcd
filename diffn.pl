#! /usr/local/bin/perl

# diffn
# C. DeTar 18 Oct 1997

# Compares two files line-by-line and reports lines with differences in
# numeric fields when they exceed a tolerance
# Tolerance is absolute for values less than 1 in magnitude and
# relative for values greater than 1 in magnitude.

# Usage...

#    diffn file1 file2 tol

# where file1 and file2 are to be compared
# A discrepancy is reported when abs(field1 - field2) > tol

($file1,$file2,$tol) = @ARGV;

(defined($tol) && defined($file2) && defined($file1)) || 
    die "Usage $0 <file1> <file2> <tol>\n";

open(FILE1,$file1) || die "Couldn't open $file1: $!";
open(FILE2,$file2) || die "Couldn't open $file2: $!";

print "diffn $file1 $file2 $tol\n";
$lines = 0;

 FILE1LINE:
    while($line1 = <FILE1>){
	chop($line1);
	if(!($line2 = <FILE2>))
	{
	    die "Premature end of file on $file2\n";
	}
	chop($line2);
	@fields1 = split(/[ \t\n]+/,$line1);
	@fields2 = split(/[ \t\n]+/,$line2);
	$i = 0;
	$same = 1;
	for(@fields1)
	{
	    $diff = $_ - $fields2[$i];
	    if(($_ > 1.) || ($_ < -1.)){ $diff = $diff/$_; }
	    if($diff < 0){ $diff = -$diff; }
	    if($diff > $tol)
	    {
		$same = 0;
		break;
		$field = $i+1; $line = $lines+1;
		print "Field $field Line $line\n";
	    }
	    $i++;
	}
	if(!$same)
	{
	    print "< ",$line1,"\n";
	    print "> ",$line2,"\n";
	}
	$lines++;
    }

if(<FILE2>)
{
    die "Premature end of file on $file1\n";
}
printf "Compared $lines lines with a tolerance of $tol\n";
