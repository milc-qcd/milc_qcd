#! /usr/local/bin/perl

# seterrfile
# C. DeTar 26 Mar 2005

# Constructs a tentative error tolerance file.  The file is then used
# by diffn3.pl to determine the agreement of other test output files
# and the fiducial sample output file.

# Here is how we do regression testing.  We compare test output with a
# fiducial sample output.  The comparison is done only for an excerpt
# of these files.  See headtail.pl for the constructions of these
# excerpts files.  The diffn3.pl script works with the test excerpt
# file and fiducial excerpt file and the error tolerance file.

# This script constructs an initial error tolerance file.  The file as
# constructed is somewhat speculative, since it is based on comparing
# the fiducial sample with only one test file.  Further testing
# usually reveals that some of the error tolerances are too tight.
# The script trainerrfile.pl is used to relax tolerances in the error
# tolerance file to allow diffn3.pl to pass an acceptable test result.

# Compares the fiducial standard excerpt file and an acceptable test
# excerpt file.  The error tolerance file is created on stdout with
# the same fields as the fiducial standard, except that numeric fields
# are replaced by whichever of the following is larger: (1) The
# absolute difference.  (2) The specified tolerance (absolute) for
# magnitudes greater than 1 (3) The specified tolerance (relative) for
# magnitudes less than 1.  Any nonnumeric fields that differ are
# replaced by "XXX".

# Usage...

#    seterrfile.pl file1 file2 tol > errfile

# where file1 and file2 are to be compared and tol is a reasonable
# tolerance.

($file1,$file2,$tol) = @ARGV;

(defined($tol) && defined($file2) && defined($file1)) || 
    die "Usage $0 <file1> <file2> <tol>\n";

open(FILE1,$file1) || die "Couldn't open $file1: $!";
open(FILE2,$file2) || die "Couldn't open $file2: $!";

$lines = 0;
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
	# Crude test for a numeric field. Surely, we can do better.
	if($_ + 1e-08 != 1e-08)
	{
	    # Numeric field
	    if(/[^\d]/){
		# Compute tolerance for noninteger
		$diff = abs($_ - $fields2[$i]);
		$err = $tol;
		if(abs($_) > 1.){ $err = abs($tol*$_); }
		if($diff > $err){ $err = $diff; }
		# Round error to one sig fig
		$errs[$i] = sprintf("%.1g",$err*2.0);
	    }
	    else
	    {
		# Integers should be exact
		$errs[$i] = 0;
	    }
	}
	else
	{
	    # Nonnumeric field
	    if($_ ne $fields2[$i]){
		$errs[$i] = "XXX";
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
