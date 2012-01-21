#! /usr/local/bin/perl

# headtail
# C. DeTar 18 Oct 1997

# Copies lines from stdin to stdout starting with
# pattern1 and ending with pattern2, inclusive.

# Usage...

#    headtail.pl patterna patternb patterna patternb ...

# where /patterna/ and /patternb/ are awk/sed-type regular expression strings

# Example

#    headtail.pl '^POINT' 'RUNNING COMPLETED'

# starts at the first line beginning with POINT and ends at the line
# containing the string RUNNING COMPLETED anywhere in the line.

# The process continues from there searching for the next patterna
# and resumes copying until patternb, etc.

# If patternb is not found, copying continues to the end of file.


$start = 0;
$patterna = shift(@ARGV);
$patternb = shift(@ARGV);

while(<STDIN>){
    if(/$patterna/){$start = 1;}
    if($start){print $_;}
    if(/$patternb/){
	if($#ARGV < 0){exit;}
	$start = 0;
	$patterna = shift(@ARGV);
	$patternb = shift(@ARGV);
    }
}
