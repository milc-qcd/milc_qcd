# Cruding parsing of preprocessor #if's and #else's, looking for macro names

# Handle common #ifdef XXXX and #ifdef XXXX /* comment */ pattern
/#ifdef/{if(NF==2)print $2;else if($3=="/*")print $2;else print "HELP: ",$0;next;}

# Handle common #ifndef XXXX and #ifndef XXXX /* comment */ pattern
/#ifndef/{if(NF==2)print $2;else if($3=="/*")print $2;else print "HELP: ",$0;next;}

# Handle common #if defined XXXX pattern
/#if/{if($2=="defined" && NF==3)print $3;else print "HELP: ",$0;next;}

# Handle solitary #else
/#else/{next;}

# Handle common #elif XXXX
/#elif/{if(NF==2)print $2;else print "HELP: ",$0;next;}

