#!/usr/bin/perl

for $f (glob "sse_*.h") {
  #print $f, "\n";
    print "#include \"$f\"\n";
}
