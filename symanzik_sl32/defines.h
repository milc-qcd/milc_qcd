#ifndef _DEFINES_H
#define _DEFINES_H

/* Compiler macros common to all targets in this application */

#define SITERAND	/* Use site-based random number generators */

/* NOTE: this version divides the sites into 32 sublattices: the lattice
   is devided into 2^4 hypercubes, and these hypercubes are then
   "checkerboarded": this allows parallel updates for extended actions
   with loops up to side 2 in up to 2 directions.
*/

#define N_SUBL32 0x20

#endif	/* _DEFINES_H */
