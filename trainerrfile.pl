#! /usr/local/bin/perl

# trainerrfile.pl
# C. DeTar 31 Mar 2005

# Trains an existing error file by accepting all differences between
# a test output file and fiducial sample output file.

# Of course, run this procedure only after you are sure the
# differences are acceptable.

# Compares two files line-by-line and compares the differences, in
# numeric fields with the allowed tolerance in the corresponding
# position in the errtol file.  When the difference is exceeded, the
# error tolerance is increased so that such a difference is accepted.
# The new error tolerance file replaces the the old one.

# Usage...

#    trainerrfile.pl file1 file2 errfile

# where file1 and file2 are to be compared
# and errfile is the error tolerance file
# A discrepancy is reported when abs(field1 - field2) > tol

($file1,$file2,$errfile) = @ARGV;

(defined($errfile) && defined($file2) && defined($file1)) || 
    die "Usage $0 <file1> <file2> <tol>\n";

# Make a backup copy of the error file
$backup = "$errfile.bak";
`mv $errfile $backup`;

open(FILE1,$file1) || die "Couldn't open $file1: $!";
open(FILE2,$file2) || die "Couldn't open $file2: $!";
open(ERRBACK,$backup) || die "Couldn't open $backup: $!";
open(ERR,">$errfile") || die "Couldn't open $errfile for writing: $!";

while($line1 = <FILE1>){
    chop($line1);
    if(!($line2 = <FILE2>))
    {
	die "Premature end of file on $file2\n";
    }
    chop($line2);
    if(!($errline = <ERRBACK>))
    {
	die "Premature end of file on $errfile\n";
    }
    chop($errline);
    @fields1 = split(/[ \t\n]+/,$line1);
    @fields2 = split(/[ \t\n]+/,$line2);
    @errs = split(/[ \t\n]+/,$errline);
    $i = 0;
    for(@fields1)
    {
	$tol = $errs[$i];
	$diff = abs($_ - $fields2[$i]);
	# Unless the corresponding errline field is XXX
	if( (($fields2[$i] + 1e-08 == 1e-08) && 
	     ($_ ne $fields2[$i]) && $tol ne "XXX") )
	{
	    $errs[$i] = "XXX";
	}
	elsif( $diff > $tol )
	{
	    # Replace error with difference
	    # Round difference up one sig fig
	    $errs[$i] = sprintf("%.1g",$diff*2.0);
	}
	$i++;
    }
    $errline = join(" ",@errs);
    print ERR "$errline\n";
}

if(<FILE2>)
{
    die "$file2 is longer than $file1\n";
}

if(<ERRBACK>)
{
    die "$errfile is longer than $file1\n";
}

