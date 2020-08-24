/*************************** charge_utilities.c  ******************************/
/* MIMD version 7 */

#include <stdio.h>
#include "../include/complex.h"
#include "../include/fermion_links.h"

/* Routines for indexing unique charge values */
/* Used to scan input parameters in various applications */

/* First member of charge table is always zero */
void
start_charge(double charge_table[], int *n){
  charge_table[0] = 0;
  *n = 1;
}

/* Return the index of the charge value. */
/* Exact match is required */
int
index_charge(double charge_table[], int n, double find_charge){
  int i;
  /* Look for number in table */
  for(i = 0; i < n; i++){
    if(charge_table[i] == find_charge)
      break;
  }
  return i;
}
    
/* Add the charge to table if new.  Return its index. */
int
fill_charge(double charge_table[], int *n, double next_charge){
  /* Look for number in table */
  int i = index_charge(charge_table, *n, next_charge);
  /* Add to table if not found */
  if(i == *n){
    if(*n >= MAX_CHARGE){
      printf("charge table overflowed\n");
      terminate(1);
    }
    charge_table[i] = next_charge;
    *n = *n + 1;
  }
  return i;
}

