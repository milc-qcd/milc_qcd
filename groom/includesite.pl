#! /usr/local/bin/perl

open(LATTICE_H,"lattice.h") || die "Can't open lattice.h\n";

while(<LATTICE_H>){
    if(/^#include/ && /site\.h/){
       print "\n";
       print "/* Begin definition of site structure */\n";
       print "\n";
       open(SITE_H,"site.h") || die "Can't open site.h\n";
       while(<SITE_H>){
	   if(!/SITE_H/){
	       print $_;
	   }
       }
       print "/* End definition of site structure */\n";
       print "\n";
       print "/* Definition of globals */\n";
    }
    else{
	print $_;
    }
}
