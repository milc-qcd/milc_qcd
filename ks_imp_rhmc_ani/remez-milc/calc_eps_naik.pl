#!/usr/bin/perl
# calculate eps-correction to the Naik term using 8-th order expansion
# eps=-27/40*(am)^2+327/1120*(am)^4-15607/268800*(am)^6-73697/3942400*(am)^8
# Follana et al., PRD75, 054502 (2007)

if($#ARGV!=0) {
  print "Usage: calc_eps_naik.pl <am_c>\n";
  exit(1);
}

$mc=$ARGV[0];
$mc2=$mc*$mc;
$mc4=$mc2*$mc2;
$mc6=$mc4*$mc2;
$mc8=$mc4*$mc4;

$eps=(-27.0/40.0)*$mc2+(327.0/1120.0)*$mc4-(15607.0/268800.0)*$mc6
    -(73697.0/3942400.0)*$mc8;

print "epsilon = ",$eps,"\n";

exit(0);

