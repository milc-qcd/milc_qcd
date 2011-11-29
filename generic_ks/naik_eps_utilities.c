/*************************** naik_eps_utilities.c  ******************************/
/* MIMD version 7 */

#include <stdio.h>
#include "../include/complex.h"
#include "../include/fermion_links.h"

/* Routines for indexing unique Naik term epsilon values */
/* Used to scan input parameters in various applications */

/* First member of eps_naik table is always zero */
void
start_eps_naik(double eps_naik_table[], int *n){
  eps_naik_table[0] = 0;
  *n = 1;
}

/* Return the index of the Naik epsilon value. */
/* Exact match is required */
int
index_eps_naik(double eps_naik_table[], int n, double find_eps_naik){
  int i;
  /* Look for number in table */
  for(i = 0; i < n; i++){
    if(eps_naik_table[i] == find_eps_naik)
      break;
  }
  return i;
}
    
/* Add the Naik term epsilon to table if new.  Return its index. */
int
fill_eps_naik(double eps_naik_table[], int *n, double next_eps_naik){
  /* Look for number in table */
  int i = index_eps_naik(eps_naik_table, *n, next_eps_naik);
  /* Add to table if not found */
  if(i == *n){
    if(*n >= MAX_NAIK){
      printf("eps_naik table overflowed\n");
      terminate(1);
    }
    eps_naik_table[i] = next_eps_naik;
    *n = *n + 1;
  }
  return i;
}

