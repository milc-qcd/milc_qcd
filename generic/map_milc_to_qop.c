/*********************** map_milc_to_qop.c **************************/
/* Functions for mapping MILC data layouts to raw QOP layouts       */
/* C. DeTar 10/19/2005                                              */

#include "generic_includes.h"
#include <qop.h>
#include <string.h>

static int is_qop_initialized = 0;


/* The MILC layout routines list the coordinates and assume 4D */

static int milc_node_number(const int coords[]){
  return node_number(coords[0],coords[1],coords[2],coords[3]);
}

static int milc_node_index(const int coords[]){
  return node_index(coords[0],coords[1],coords[2],coords[3]);
}

/* Initialize QOP */

QOP_status_t initialize_qop(){
  int dir;
  QOP_status_t status;
  static int latsize[4];
  QOP_layout_t layout;
  
  latsize[0] = nx;
  latsize[1] = ny;
  latsize[2] = nz;
  latsize[3] = nt;

  if(is_qop_initialized)return QOP_SUCCESS;

  layout.node_number = milc_node_number;
  layout.node_index = milc_node_index;
  layout.latdim = 4;
  layout.latsize = latsize;
  layout.machdims = 4;
  layout.gridsize = get_logical_coordinate();
  layout.this_node = this_node;
  layout.sites_on_node = sites_on_node;

  status = QOP_init(&layout);

  if(status == QOP_SUCCESS)
      is_qop_initialized = 1;

  return status;
}

void site_coords(int coords[4],site *s){
  coords[0] = s->x;
  coords[1] = s->y;
  coords[2] = s->z;
  coords[3] = s->t;
}

/* Map MILC links to raw order */
su3_matrix **create_raw_G_from_site_links(){
  int coords[4];
  int i,j,dir;
  site *s;
  su3_matrix **rawlinks = NULL;

  rawlinks = (su3_matrix **)malloc(4*sizeof(su3_matrix *));
  FORALLUPDIR(dir){
    rawlinks[dir] = (su3_matrix *)malloc(sites_on_node*sizeof(su3_matrix));
    if(rawlinks == NULL){
      printf("create_raw_G_from_site_link: No room for rawlinks\n");
      return NULL;
    }
  }
  
  FORALLSITES(i,s){
    site_coords(coords,s);
    if(QOP_node_number_raw(coords) != this_node){
      printf("create_raw_G_from_site_link: incompatible layout\n");
      return NULL;
    }
    j = QOP_node_index_raw_G(coords);
    FORALLUPDIR(dir){
      memcpy(rawlinks[dir] + j, &(s->link[dir]), sizeof(su3_matrix));
    }
  }

  return rawlinks;
}

void destroy_raw_G(su3_matrix *rawlinks[]){
  int dir;

  FORALLUPDIR(dir){
    if(rawlinks[dir] != NULL)
      free(rawlinks[dir]);
  }

  free(rawlinks);
}

/* Map MILC momentum to raw order */
su3_matrix **create_raw_F_from_site_mom(){
  int coords[4];
  int i,j,dir;
  site *s;
  su3_matrix **rawforce = NULL;
  su3_matrix tmat;

  rawforce = (su3_matrix **)malloc(4*sizeof(su3_matrix *));
  FORALLUPDIR(dir){
    rawforce[dir] = (su3_matrix *)malloc(sites_on_node*sizeof(su3_matrix));
    if(rawforce == NULL){
      printf("create_raw_F_from_site_mom: No room for rawforce\n");
      return NULL;
    }
  }
  
  FORALLSITES(i,s){
    site_coords(coords,s);
    if(QOP_node_number_raw(coords) != this_node){
      printf("create_raw_F_from_site_mom: incompatible layout\n");
      return NULL;
    }
    j = QOP_node_index_raw_F(coords);
    FORALLUPDIR(dir){
      uncompress_anti_hermitian( &(s->mom[dir]), &tmat );
      memcpy(rawforce[dir] + j, &tmat, sizeof(su3_matrix));
    }
  }

  return rawforce;
}

/* Map raw force to MILC site structure mom */
void unload_raw_F_to_site_mom(su3_matrix *rawforce[]){
  int coords[4];
  int i,j,dir;
  site *s;
  su3_matrix tmat;

  FORALLSITES(i,s){
    site_coords(coords,s);
    if(QOP_node_number_raw(coords) != this_node){
      printf("create_raw_F_from_site_mom: incompatible layout\n");
      return NULL;
    }
    j = QOP_node_index_raw_F(coords);
    FORALLUPDIR(dir){
      memcpy(&tmat, rawforce[dir] + j, sizeof(su3_matrix));
      make_anti_hermitian( &tmat, &(s->mom[dir]) ); 
    }
  }
}


void destroy_raw_F(su3_matrix *rawforce[]){
  int dir;

  FORALLUPDIR(dir){
    if(rawforce[dir] != NULL)
      free(rawforce[dir]);
  }

  free(rawforce);
}

/* Map MILC site color vector to raw order */
su3_vector *create_raw_V_from_site(field_offset x){
  int coords[4];
  int i,j,dir;
  site *s;
  su3_vector *rawsu3vec = NULL;

  rawsu3vec = (su3_vector *)malloc(sites_on_node*sizeof(su3_vector));
  if(rawsu3vec == NULL){
    printf("create_raw_V_from_site_link: No room for rawlink\n");
    return NULL;
  }
  
  FORALLSITES(i,s){
    site_coords(coords,s);
    if(QOP_node_number_raw(coords) != this_node){
      printf("create_raw_V_from_site_link: incompatible layout\n");
      return NULL;
    }
    j = QOP_node_index_raw_V(coords);
    FORALLUPDIR(dir){
      memcpy(rawsu3vec + j, F_PT(s,x), sizeof(su3_vector));
    }
  }

  return rawsu3vec;
}

void destroy_raw_V(su3_vector *rawsu3vec){
  
  if(rawsu3vec != NULL)
    free(rawsu3vec);
}


