#! /usr/local/bin/perl

# seterrfile
# C. DeTar 26 Mar 2005

# Constructs a tentative error tolerance file 
# The file is used by diffn3.pl 

# Compares two files, one the fiducial standard and one an acceptable
# test file.  A third, error tolerance file is created with the
# numeric value replaced by whichever of the following is larger:
# (1) The absolute difference.
# (2) The specified tolerance (absolute) for magnitudes greater than 1
# (3) The specified tolerance (relative) for magnitudes less than 1
# Edit the error tolerance file so it contains only lines between
# the beginning and ending patter in the application Make_test file.
# Edit to make generous tolerances for quantities irrelevant 
# to numerical accuracy, such as times to load files.

# Usage...

#    seterrfile.pl file1 file2 tol > errfile

# where file1 and file2 are to be compared

($file1,$file2,$tol) = @ARGV;

(defined($tol) && defined($file2) && defined($file1)) || 
    die "Usage $0 <file1> <file2> <tol>\n";

open(FILE1,$file1) || die "Couldn't open $file1: $!";
open(FILE2,$file2) || die "Couldn't open $file2: $!";

print STDERR "diffn $file1 $file2 $tol\n";
$lines = 0;

 FILE1LINE:
    while($line1 = <FILE1>){
	chop($line1);
	if(!($line2 = <FILE2>))
	{
	    die "Premature end of file on $file2\n";
	}
	chop($line2);
	@fields1 = split(" ",$line1);
	@fields2 = split(" ",$line2);
	@errs = @fields1;
	$i = 0;
	for(@fields1)
	{
	    # Skip nonnumeric or zero fields
	    if($_ + .001 != .001)
	    {
		# Integers should be exact
		if(/[^\d]/){
		    $diff = abs($_ - $fields2[$i]);
		    $err = $tol;
		    if(abs($_) > 1.){ $err = abs($tol*$_); }
		    if($diff > $err){ $err = $diff; }
		    # Round error to one sig fig
		    $errs[$i] = sprintf("%.1g",$err*1.5);
		}
		else
		{
		    $errs[$i] = 0;
		}
	    }
	    $i++;
	}
	$errline = join(" ",@errs);
	print "$errline\n";
	$lines++;
    }

if(<FILE2>)
{
    die "Premature end of file on $file1\n";
}
printf STDERR "Compared $lines lines with a tolerance of $tol\n";
