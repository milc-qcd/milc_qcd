/********** gauge_utilities.c ********************************************/
/* MIMD version 7 */
/*--------------------------------------------------------------------*/
#include "generic_includes.h"

/*--------------------------------------------------------------------*/
/* A "G" type field has four su3_matrices per site                    */
/*--------------------------------------------------------------------*/

su3_matrix * create_G(void){
  static su3_matrix *t_links;

  t_links =(su3_matrix *)malloc(sites_on_node*4*sizeof(su3_matrix));
  if(t_links == NULL){
      printf("node %d can't malloc t_links\n",this_node);
      terminate(1);
  }
  memset(t_links, '\0', sites_on_node*4*sizeof(su3_matrix));
  return t_links;
}

/*--------------------------------------------------------------------*/

#ifndef NO_GAUGE_FIELD

su3_matrix * create_G_from_site(void){
  int i, dir;
  site *s;
  static su3_matrix *t_links;

  t_links = create_G();

  FORALLSITES(i,s){
    FORALLUPDIR(dir){
      t_links[4*i+dir] = lattice[i].link[dir];
    }
  }
  return t_links;
}

#endif

/*--------------------------------------------------------------------*/
void copy_G(su3_matrix *dst, su3_matrix *src){
  memcpy(dst, src, 4*sizeof(su3_matrix)*sites_on_node);
}

/*--------------------------------------------------------------------*/
void destroy_G(su3_matrix *t_links){
  free(t_links) ; 
}

