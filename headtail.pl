#! /usr/local/bin/perl

# headtail
# C. DeTar 18 Oct 1997

# Copies lines from stdin to stdout starting with
# pattern1 and ending with pattern2, inclusive.

# Usage...

#    headtail.pl pattern1 pattern2

# where /pattern1/ and /pattern2/ are awk/sed-type pattern-matching strings

# Example

#    headtail.pl '^POINT' 'RUNNING COMPLETED'

# starts at the first line beginning with POINT and ends at the line
# containing the string RUNNING COMPLETED anywhere in the line.


$pattern1 = $ARGV[0];
$pattern2 = $ARGV[1];

$start = 0;
 STDINLINE:
    while(<STDIN>){
	if(/$pattern1/){$start = 1;}
	if($start){print $_;}
	if(/$pattern2/){exit 0;}
    }
