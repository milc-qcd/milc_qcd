#! /usr/local/bin/perl

# diffn3
# C. DeTar 26 Mar 2005

# Compares two files line-by-line and reports lines with differences
# in numeric fields when they exceed a tolerance and lines in
# nonnumeric fields when they are not exactly the same.  Tolerances
# are specified by a third, error tolerance file.

# Usage...

#    diffn file1 file2 errfile

# where file1 and file2 are to be compared
# and errfile is the error tolerance file
# A discrepancy is reported when abs(field1 - field2) > tol

($file1,$file2,$errfile) = @ARGV;

(defined($errfile) && defined($file2) && defined($file1)) || 
    die "Usage $0 <file1> <file2> <tol>\n";

open(FILE1,$file1) || die "Couldn't open $file1: $!";
open(FILE2,$file2) || die "Couldn't open $file2: $!";
open(ERR,$errfile) || die "Couldn't open $errfile: $!";

print "diffn $file1 $file2 $errfile\n";
$lines = 0;

$difflines = 0;
while($line1 = <FILE1>){
    chop($line1);
    if(!($line2 = <FILE2>))
    {
	die "Premature end of file on $file2\n";
    }
    chop($line2);
    if(!($errline = <ERR>))
    {
	die "Premature end of file on $errfile\n";
    }
    chop($errline);
    @fields1 = split(/[ \t\n]+/,$line1);
    @fields2 = split(/[ \t\n]+/,$line2);
    @errs = split(/[ \t\n]+/,$errline);
    $i = 0;
    $same = 1;
    for(@fields1)
    {
	$tol = $errs[$i];
	$diff = abs($_ - $fields2[$i]);
	# Nonumeric or zero fields should match exactly
	# Unless the corresponding errline field is XXX
	if( (($_ + 1e-08 == 1e-08) && 
	     ($_ ne $fields2[$i]) && $tol ne "XXX") ||
	    $diff > $tol )
	{
	    $same = 0;
	    $field = $i+1; $line = $lines+1;
	    $diffround = sprintf("%.2g",$diff);
	    print "Field $field Line $line diff $diffround >= tol $tol\n";
	}
	$i++;
    }
    if(!$same)
    {
	print "< ",$line1,"\n";
	print "> ",$line2,"\n";
	$difflines++;
    }
    $lines++;
}

if(<FILE2>)
{
    die "$file2 is longer than $file1\n";
}

if(<ERRBACK>)
{
    die "$errfile is longer than $file1\n";
}

if($difflines == 0){
    printf STDERR "$file1 OK\n";
}
else{
    printf STDERR "$file1 NOT OK $difflines lines differ\n";
}

