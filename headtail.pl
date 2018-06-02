#! /usr/bin/perl

# headtail
# C. DeTar 18 Oct 1997

# Copies lines from stdin to stdout starting with
# pattern1 and ending with pattern2, inclusive.
# and filtering for matches with a selection pattern

# Usage...

#    headtail.pl patterna patternb patterna patternb ... SELECT=pat1|pat2

# where /patterna/ and /patternb/ are awk/sed-type regular expression strings
# also /pat1/ and /pat2/ are regular expressions

# Example

#    headtail.pl '^POINT' 'RUNNING COMPLETED' "SELECT=PBP|FACTION"

# starts at the first line beginning with POINT and ends at the line
# containing the string RUNNING COMPLETED anywhere in the line
# and prints lines between them containing PBP or FACTION

# The process continues from there searching for the next patterna
# and resumes copying until patternb, etc.

# If patternb is not found, copying continues to the end of file.

if(@ARGV[$#ARGV] =~ /SELECT/){
    $select = pop(@ARGV);
    $select =~ s/SELECT=//;
}
$start = 0;
$patterna = shift(@ARGV);
$patternb = shift(@ARGV);

while(<STDIN>){
    if(/$patterna/){$start = 1;}
    if($start){
	if(defined($select)){
	    print grep(/$select/,$_);
	} else {
	    print $_;
	}
	if(/$patternb/){
	    if($#ARGV < 0){exit;}
	    $start = 0;
	    $patterna = shift(@ARGV);
	    $patternb = shift(@ARGV);
	}
    }
}
