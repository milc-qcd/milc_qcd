#! /usr/local/bin/perl

@FILE = @ARGV;

for(@FILE){
    $filename = $_;

    $outfile = ">$filename.new";
    open(NEWCODE,$outfile) || die "Can't open $outfile\n";

    @path = split('/',$filename);
    $strip = pop(@path);
    $head1 =  "/********************** $strip (in complex.a) **********************/\n";
    $head2 = "/* MIMD version 6 */\n";

    open(CODE,$filename) || die "Can't open $filename\n";
    $first = 1;  $second = 0;
    while(<CODE>){
	if($first){
	    $first = 0;
	    if(/\/\*./ && /\.c/){
		# If header line exists, check for correct file name
		$second = 1;
		if(!/$strip/){
		    # If incorrect, fix it
		    $_ =~ s/[^\*]*\.c/ $strip/;
		}
	    }
	    else {
		# If header line does not exist, create both lines
		print NEWCODE $head1;
		print NEWCODE $head2;
		$second = 0;
	    }
	}
	elsif($second){
	    $second = 0;
	    if(!/MIMD/){
		print NEWCODE $head2;
	    }
	}
	print NEWCODE $_;
    }
}
