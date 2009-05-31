#! /bin/sh

# Extract checksum records from an extended KS source file

msglist="2 3 4 5 6 7 8 9 10 11 12 13"
rec=4

for f in $*
do
  for m in $msglist
  do
    $HOME/scidac/qio-single/bin/lime_extract_record $f $m $rec foo
    cat foo
    echo ""
  done
done
