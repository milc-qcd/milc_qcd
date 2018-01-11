/* Shamelessly copied from the QPHIX interface, it doesn't exit if the
 * assertion fails, but gives an error code and lets the user handle it */

#define GRID_ASSERT(cond, retval) if(!(cond)) { printf("%s:%d Assertion failed\n\t%s evaluates false\n", __FILE__, __LINE__, #cond); exit(retval);}
