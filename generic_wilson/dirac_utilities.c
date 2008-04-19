/***************** dirac_utilities.c *************************************/

/* Miscellaneous (mostly) Dirac utilities */
/* MIMD version 7 */

#include "generic_wilson_includes.h"
#include <string.h>

/*--------------------------------------------------------------------*/
double start_timing(void){
  double dtime;

#ifdef PRTIME
  dtime = -dclock();
#else
  dtime = 0;
#endif
  return dtime;
}

/*--------------------------------------------------------------------*/
void print_timing(double dtime, char *str){

#ifdef PRTIME
  dtime += dclock();
  node0_printf("Time for %s %e\n",str, dtime);  fflush(stdout);
#endif

}

/*--------------------------------------------------------------------*/
spin_wilson_vector *create_swv_field(void){
  spin_wilson_vector *swv;
  
  swv = (spin_wilson_vector *)
    malloc(sizeof(spin_wilson_vector)*sites_on_node);
  
  if(swv == NULL){
    printf("create_swv_field(%d): No room for temporary\n",this_node);
    terminate(1);
  }
  return swv;
}

/*--------------------------------------------------------------------*/
void clear_swv_field(spin_wilson_vector *swv){
  memset(swv,'\0',sites_on_node*sizeof(spin_wilson_vector));
}

/*--------------------------------------------------------------------*/
void clear_wv_field(wilson_vector *wv){
  memset(wv,'\0',sites_on_node*sizeof(wilson_vector));
}

/*--------------------------------------------------------------------*/
spin_wilson_vector *extract_swv_from_wp(wilson_prop_field wp, int color){
  return wp[color];
}

/*--------------------------------------------------------------------*/
wilson_prop_field create_wp_field(void){
  wilson_prop_field wp;
  int color;
  
  wp = (wilson_prop_field) malloc(4*sizeof(spin_wilson_vector *));
  if(wp == NULL){
    printf("create_wp_field(%d): No room for temporary\n",this_node);
    terminate(1);
  }

  for(color= 0; color < 3; color++)
    wp[color] = create_swv_field();
  
  return wp;
}

/*--------------------------------------------------------------------*/
void clear_wp_field(wilson_prop_field wp){
  int color;
  
  for(color= 0; color < 3; color++)
    clear_swv_field(wp[color]);
}

/*--------------------------------------------------------------------*/
void copy_wp_field(wilson_prop_field wpcopy, wilson_prop_field wp){
  int color, i;
  site *s;

  for(color = 0; color < 3; color++){
    FORALLSITES(i,s){
      wpcopy[color][i] = wp[color][i];
    }
  }
}

/*--------------------------------------------------------------------*/
wilson_prop_field create_wp_field_copy(wilson_prop_field w){
  wilson_prop_field wp;

  wp = create_wp_field();
  copy_wp_field(wp, w);

  return wp;
}

/*--------------------------------------------------------------------*/
wilson_vector *create_wv_field(void){
  wilson_vector *wv;
  
  wv = (wilson_vector *)
    malloc(sizeof(wilson_vector)*sites_on_node);
  
  if(wv == NULL){
    printf("create_wv_field(%d): No room for temporary\n",this_node);
    terminate(1);
  }
  return wv;
}

/*--------------------------------------------------------------------*/
void copy_wv_from_wp(wilson_vector *wv, wilson_prop_field wp, 
		     int color, int spin){
  int i;
  site *s;
  
  FORALLSITES(i,s){
    wv[i] = wp[color][i].d[spin];
  }
}

/*--------------------------------------------------------------------*/
void copy_wv_from_wprop(wilson_vector *wv, wilson_propagator *wprop, 
			int color, int spin){
  int i;
  site *s;
  
  FORALLSITES(i,s){
    wv[i] = wprop[i].c[color].d[spin];
  }
}

/*--------------------------------------------------------------------*/
void copy_wv_from_swv(wilson_vector *wv, spin_wilson_vector *swv, int spin){
  int i;
  site *s;
  
  FORALLSITES(i,s){
    wv[i] = swv[i].d[spin];
  }
}

/*--------------------------------------------------------------------*/
void copy_wp_from_wv(wilson_prop_field wp, wilson_vector *wv, 
		     int color, int spin){
  int i;
  site *s;
  
  FORALLSITES(i,s){
    wp[color][i].d[spin] = wv[i];
  }
}

/*--------------------------------------------------------------------*/
void destroy_swv_field(spin_wilson_vector *swv){
  free(swv);
}

/*--------------------------------------------------------------------*/
void destroy_wv_field(wilson_vector *wv){
  free(wv);
}

/*--------------------------------------------------------------------*/
void destroy_wp_field(wilson_prop_field wp){
  int color;

  if(wp == NULL)return;
  for(color = 0; color < 3; color++){
    destroy_swv_field(wp[color]);
  }
  free(wp);
}
