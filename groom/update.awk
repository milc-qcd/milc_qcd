/\.lastmake\.humor_irix/{print "	-/bin/rm -f .lastmake.*";next;}
/humor_irix\.o/{print "	-/bin/rm -f *.o";next;}
{if($1=="make")$1="\t${MAKE}";print}
