#! /usr/bin/perl

# Run regression tests

# Usage:  

#     check.pl [exec-target [input-case [prec]]] < checklist

# Format of checklist: two types of lines.  One specifies the
# executable, precision, and input case.  The other specifies
# additional, optional output coming from the same test run.  Every
# test run creates a standard log file (stdout) that is compared with
# a fiducial log file.  Some test runs create additional output that
# is to be compared as well.  These additional files are specified in
# the second type of line.  The patterns in each case are used to
# demarcate the region of output that is to be compared with the
# trusted output.

# exec executable precision add-macro input-case pattern1 pattern2
#    extra-output file pattern1 pattern2
#    extra-target target


$choose_exec = shift(@ARGV);
$choose_input_case = shift(@ARGV);
$choose_prec = shift(@ARGV);

$oldexec = "";
$oldprec = "";
$oldmacro = "";

while(<>){
    if(/\#/){next;}
    if(/exec/){
	($tag, $exec, $prec, $add_macro, $input_case, $patterns) = split(" ",$_,6);
	chop($patterns);
    }
    elsif(/extra-output/){
	($tag, $extra_output, $patterns) = split(" ",$_,3);
	chop($patterns);
    }
    elsif(/extra-target/){
	($tag, $extra_target) = split(" ",$_);
    }
    else{
	next;
    }

    if(defined($choose_exec) && $exec ne $choose_exec){next;}

    if(defined($choose_prec) && $prec ne $choose_prec){next;}

    if(defined($choose_input_case) && $input_case ne $choose_input_case){next;}

    $prefix = "$exec.$prec";
    if($input_case ne "-"){$prefix = "$exec.$input_case.$prec";}

    $macro = "";
    if($add_macro ne "-"){$macro = $add_macro;}


    if($oldexec ne $exec || $oldprec ne $prec || $oldmacro ne $macro ){
	print "===========================================================\n";
	print "Making $exec for precision $prec, case $input_case, and macro $macro\n";
	if($macro eq ""){
	    print `cd .. ; make clean ; make $exec \"PRECISION=$prec\" 2>&1`;
	} else {
	    print `cd .. ; make clean ; make $exec \"PRECISION=$prec\" \"$macro\" 2>&1`;
	}
	$oldexec = $exec;
	$oldprec = $prec;
	$oldmacro = $macro;
    }

    if($tag eq "exec"){
	$output = "$prefix.test-out";
	print "Making $output\n";
	print `make $output \"EXEC=$exec\" \"PREFIX=$prefix\" 2>&1`;
	
	$diff = "$prefix.test-diff";
#	print "Making $diff\n";
	print "Making $diff with patterns $patterns\n";
	print `make $diff \"PATTERNS=$patterns\" 2>&1`;
    }
    elsif($tag eq "extra-output"){
	$prefix = $extra_output;
	$diff = "$prefix.test-diff";
	print "Making $diff\n";
	print `make $diff \"PATTERNS=$patterns\" 2>&1`;
    }
    elsif($tag eq "extra-target"){
	print "Making $extra_target\n";
	print `make $extra_target 2>&1`;
    }
}

