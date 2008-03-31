#!/usr/bin/perl
#
# nasm2c.pl
#   Translates a nasm source program into inline gcc assembler.
#   Intended for porting Fermilab MILC SSE routines into inline
#   versions.
#
# Syntax:
#   ./nasm2c.pl source.nas [argument list]
#
#   By default, the generated macro has (aa,bb,cc) as argument.
#   You may override this by specifying [argument list]
#
# Rules:
#   - all comments are ignored  (delimited by ;)
#   - all labels are ignored (delimited by :)
#   - we assume that [eax+y] references are to be translated into
#     "m" or "=m" constrained forms.  The exact translations have to
#     be done manually.
#
# Jul 26 2006
# added support for clobber lists
#   -Jordan Soyke

$do_pr_head = 0;

if ($#ARGV > 0) {
    if ($#ARGV >1) {
	print "Usage: nasm2c.pl source.nas [arg.list]\n";
	exit(1);
    }
    $macro_args = pop(@ARGV);
}
else {
    $macro_args = "aa,bb,cc";
}


while (<>) {
    $l = $_;
    chop($l);
    
    # translate any inline-only instructions into a recognizable form
    if ($l =~ m/;!inline/o) {
	$l = $';
	$l =~ s/<(\S+)>/\[0\]/o;
	$l .= "; <$1> ";
    }
    
    # get any parameters embedded in comments
    if ($l =~ m/<(\S+)>/o) {
	$param = $1;
    }
    else {
	undef $param;
    }
    
    # get rid of comments
    if ($l =~ m/([^;]*);/o) {
	$l = $1;
    };
    
    # get rid of labels
    $l =~ s/^\S+://o;
    
    # bypass blank lines
    if ($l =~ m/^\s*$/o || $l =~ m/^\s*bits/o) {
	next;
    };
    
    # translate "global routine" into "#define _inline_routine"
    if ($l =~ m/^global\s+(\S+)/o) {
	$l = "#define _inline_sse_$1($macro_args) \\";
	print "$l\n{ \\\n__asm__ __volatile__ (";
	next;
    }
    
    if ($l =~ m/^\s*(\S+)\s+(\S+)/o) {         # get opcode and operands
	$instr = $1;
	$args  = $2;
	
	if ($instr =~ m/^(mov|push|pop|align|dd|add)$/o) {    # ignore non-inline opcodes
	    next;
	}
	
	if ($args =~ m/,\[/o) {                      # check for input (source) parameter
	    if (defined $direct && $direct == 2) {     # check for prior dest term
		pr_ops();
		$do_pr_head = 1;
#	print "; \\\n__asm__ __volatile__ (";
		undef @ins;
		undef $direct;
	    }
	    $direct = 1;                               # set direction to source
	}
	elsif ($args =~ m/\],/o) {                   # check for output (dest) parameter
	    if (defined $direct && $direct == 1) {     # check for prior source term
		pr_ops();
		$do_pr_head = 1;
#	print "; \\\n__asm__ __volatile__ (";
		undef @ins;
		undef $direct;
	    }
	    $direct = 2;                               # set direction to dest
	}
	
	if ($do_pr_head == 1) {
	    pr_head();
	}
	print "\"$instr ";                           # start outputting source line ...
	$first = 1;                                  # last operand is special, so track
	while ($args =~ m/(.+,)*(\S+)$/o) {          # ... get next operand (read backwards)
	    if ($first == 0) {
		print ", ";                              # if second, prepend a space
	    }
	    $n = $2;
	    
	    if ($n =~ m/^0[xX]/o) {                    # hexadecimal literal?
		$n = "\$" . $n;                          #  ... yes, prepend a $
	    }
	    
	    elsif ($n =~ m/(\[\S+\])/o) {              # non-register (i.e., source or destination?)
		if (defined $param) {
		    $n = $param;
		}
		if ($first == 1) {
		    $n = "\"m\" (" . $n . ")";             # first operand, so it's a source term
		    $direct = 1;
		    push(@ins, $n);                        # keep track of operands
		    $n = "%$#ins";
		}
		
		else {                                   # dest term
		    $n = "\"=m\" (" . $n . ")";
		    $direct = 2;
		    push(@ins, $n);
		    $n = "%$#ins";
		}
	    }
	    
	    else {                                     # a register reference, so prepend %%
		$n = "%" .$n;
		$clobber_list{$n} = 1;                       #but only once for the clobber list	
		$n = "%" . $n;
	    }
	    
	    print "$n";                                
	    $args =~ s/,*[^,]+$//o;
	    if ($first == 1) {
		$first = 0;
	    }
	}
	if ($#ins == 5) {
	    print "\" \\\n                      ";
	    pr_ops();
	    $do_pr_head = 1;
#      print "; \\\n__asm__ __volatile__ (";
	    undef @ins;
	    undef $direct;
	}
	else {
	    print " \\n\\t\" \\\n                      ";
	}
    }
}
if (defined @ins) {
    pr_ops();
}
print "; \\\n}\n";

sub pr_ops {
    print ": \\\n                      ";
    if ($direct == 1) {
	print ": \\\n                      ";
    }
    foreach $idx (0 .. $#ins) {
	if ($idx != 0) {
	    print ", \\\n                      ";
	}
	print "$ins[$idx]";
    }
    
    if ( defined %clobber_list ) {
	print "\\\n                      ";
	if ($direct != 1) {
	    print ": \\\n                      ";
	}
	print ": \\\n                      ";
	
	foreach $clob_reg (sort keys %clobber_list)
	{
	    print "\"".$clob_reg."\"";
	    delete ($clobber_list{$clob_reg});
	    
	    if ( %clobber_list > 0 )
	    {
    		print ", \\\n                      ";
	    }
	    else
            {
    		print ", \\\n                      \"memory\"";
	    }
	}
	undef %clobber_list;
    }
    print ")";
}

sub pr_head {
    print "; \\\n__asm__ __volatile__ (";
    $do_pr_head = 0;
}
