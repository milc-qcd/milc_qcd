/************************** w_source.c *****************************/
/* MIMD version 7 */

/* Initialize a source for the inverter */
/* Initialize a sink */
/* Interpret source and sink parameters */

#include "generic_wilson_includes.h"
#include "../include/io_wprop.h"
#include <string.h>
#ifdef HAVE_QIO
#include <qio.h>
#include "../include/io_scidac.h"
#include "../include/io_scidac_w.h"
#endif

static su3_matrix *ape_links = NULL;

/* Utilities */


void init_wqs(wilson_quark_source *wqs){
  wqs->type             = UNKNOWN;
  wqs->c_src            = NULL;
  wqs->wv_src           = NULL;
  wqs->file_initialized = 0;
  wqs->d1               = 0.;
#ifdef HAVEQIO
  wqs->infile           = NULL;
  wqs->outfile          = NULL;
#endif
  wqs->x0               = 0;
  wqs->y0               = 0;
  wqs->z0               = 0;
  wqs->t0               = 0;
  wqs->descrp[0]        = '\0';
}

/* Must be called to free */
static void reset_wqs(wilson_quark_source *wqs){
  if(wqs->c_src != NULL){free(wqs->c_src); wqs->c_src = NULL;}
  if(wqs->wv_src != NULL){free(wqs->wv_src); wqs->wv_src = NULL;}
}

void clear_wqs(wilson_quark_source *wqs){
  reset_wqs(wqs);
#ifdef HAVE_QIO
  if(wqs->infile != NULL){
    if(wqs->type == COMPLEX_FIELD_FILE)
      r_close_complex_scidac_file(wqs->infile);
    else if(wqs->type == DIRAC_FIELD_FILE)
      r_close_w_vector_scidac_file(wqs->infile);
  }
#endif
}

void alloc_wqs_wv_src(wilson_quark_source *wqs)
{
  if(wqs->wv_src == NULL)
    wqs->wv_src = (wilson_vector *)malloc(sizeof(wilson_vector)*sites_on_node);
  if(wqs->wv_src == NULL){
    printf("alloc_wqs_wv_src(%d): No room for source field\n",this_node);
    terminate(1);
  }
  memset(wqs->wv_src, 0, sizeof(wilson_vector)*sites_on_node);
}

void alloc_wqs_c_src(wilson_quark_source *wqs)
{
  if(wqs->c_src == NULL)
    wqs->c_src = (complex *)malloc(sizeof(complex)*sites_on_node);
  if(wqs->c_src == NULL){
    printf("alloc_wqs_c_src(%d): No room for source field\n",this_node);
    terminate(1);
  }
  memset(wqs->c_src, 0, sizeof(complex)*sites_on_node);
}

/* Copy c-source from wv */
static void extract_C_from_D(complex *dest, wilson_vector *src, 
			     int spin, int color){
  int i; site *s;
  FORALLSITES(i,s){
    dest[i] = src[i].d[spin].c[color];
  }
}

/* Copy wv from c-source */
static void insert_D_from_C(wilson_vector *dest, complex *src, 
			    int spin, int color){
  int i; site *s;
  FORALLSITES(i,s){
    dest[i].d[spin].c[color] = src[i];
  }
}

/* Copy wv-source */
static void copy_D(wilson_vector *dest, wilson_vector *src)
{
  int i; site *s;
  FORALLSITES(i,s){
    dest[i] = src[i];
  }
}

/* Choose the specific USQCD format appropriate to the source type */
int choose_usqcd_w_file_type(int source_type){
  int file_type;

  switch(source_type){
  case COMPLEX_FIELD_FILE:
  case COMPLEX_FIELD_FM_FILE:
  case COMPLEX_FIELD_STORE:
  case GAUSSIAN:
  case POINT:
  case WAVEFUNCTION_FILE:
  case UNKNOWN:
    file_type = FILE_TYPE_W_USQCD_C1D12;
    break;
  case COVARIANT_GAUSSIAN:
  case DIRAC_FIELD_FILE:
  case DIRAC_FIELD_FM_FILE:
  case DIRAC_FIELD_STORE:
  case FAT_COVARIANT_GAUSSIAN:
  case DERIV1:
  case DERIV2_D:
  case DERIV2_B:
  case DERIV3_A:
  case FAT_COVARIANT_GAUSSIAN_DERIV1:
  case FAT_COVARIANT_GAUSSIAN_DERIV2_D:
  case FAT_COVARIANT_GAUSSIAN_DERIV2_B:
  case ROTATE_3D:
    file_type = FILE_TYPE_W_USQCD_DD_PAIRS;
    break;
  default:
    file_type = -1;
  }
  return file_type;
}

static void insert_mom_D(wilson_vector *wv, int mom[3], 
			 int x0, int y0, int z0, int t0){
  int i; site *s;
  int c,d;
  Real pi = M_PI;
  Real th;
  Real px = 2*pi*mom[0]/nx;
  Real py = 2*pi*mom[1]/ny;
  Real pz = 2*pi*mom[2]/nz;
  complex z,y;

  if(mom[0] == 0 && mom[1] == 0 && mom[2] == 0)return;

  FORALLSITES(i,s)if(t0 == ALL_T_SLICES || s->t == t0){
    th = px*(s->x-x0) + py*(s->y-y0) + pz*(s->z-z0);
    y.real = cos(th); y.imag = sin(th);
    for(c=0;c<3;c++)
      for(d=0;d<4;d++){
	z = wv[i].d[d].c[c];
	CMUL(z,y,wv[i].d[d].c[c]);
      }
  }
}

static void insert_mom_C(complex *c, int mom[3], int x0, 
			 int y0, int z0, int t0){
  int i; site *s;
  Real pi = M_PI;
  Real th;
  Real px = 2*pi*mom[0]/nx;
  Real py = 2*pi*mom[1]/ny;
  Real pz = 2*pi*mom[2]/nz;
  complex z,y;
  
  if(mom[0] == 0 && mom[1] == 0 && mom[2] == 0)return;

  FORALLSITES(i,s)if(t0 == ALL_T_SLICES || s->t == t0){
    th = px*(s->x-x0) + py*(s->y-y0) + pz*(s->z-z0);
    y.real = cos(th); y.imag = sin(th);
    z = c[i];
    CMUL(z,y,c[i]);
  }
}

#ifdef HAVE_QIO
/********************************************************************/
/* Parse the record XML to get the color and check it               */
/********************************************************************/
/* For the source file we borrow the XML encoding from the USQCD
   Wilson propagator file  */
static int check_color_spin(QIO_String *recxml, int color, int spin){
  int status;
  int input_color, input_spin;
  QIO_USQCDPropRecordInfo recinfo;
  char myname[] = "check_color_spin";

  status = QIO_decode_usqcd_proprecord_info(&recinfo, recxml);
  if(status != QIO_SUCCESS) 
    return qio_status(status);
  input_color = QIO_get_usqcd_proprecord_color(&recinfo);
  input_spin = QIO_get_usqcd_proprecord_spin(&recinfo);
  if(color != input_color || spin  != input_spin ){
    node0_printf("%s(%d): Error: expected color %d and spin %d got %d and %d\n",
		 myname, this_node, color, spin, 
		 input_color, input_spin);
    return 1;
  }
  return 0;
}

#endif

/*--------------------------------------------------------------------*/
static su3_matrix * create_G(void){
  static su3_matrix *t_links;

  t_links =(su3_matrix *)malloc(sites_on_node*4*sizeof(su3_matrix));
  if(t_links == NULL){
      printf("node %d can't malloc t_links\n",this_node);
      terminate(1);
  }
  return t_links;
}

/*--------------------------------------------------------------------*/

static su3_matrix * create_G_from_site(void){
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

/*--------------------------------------------------------------------*/
static void copy_G(su3_matrix *dst, su3_matrix *src){
  memcpy(dst, src, 4*sizeof(su3_matrix)*sites_on_node);
}

/*--------------------------------------------------------------------*/
static void destroy_G(su3_matrix *t_links){
  free(t_links) ; 
}

/*--------------------------------------------------------------------*/
/* Do the 3D Fermilab rotation on a Wilson vector field.  That is,
   compute src <- (1 + a d1 * \gamma * D/2) src where D is a 3D
   covariant difference normalized to give twice the covariant
   derivative in the continuum limit. */
static void rotate_3D_wvec(wilson_vector *src, Real d1)
{
  int i;
  site *s;
  wilson_vector *mp, *tmp;
  
  mp  = create_wv_field();
  tmp = create_wv_field();
  
  /* Do Wilson Dslash on the source field */
  dslash_w_3D_field(src, mp,  PLUS, EVENANDODD);
  dslash_w_3D_field(src, tmp, MINUS, EVENANDODD);
  
  FORALLSITES(i,s){
    /* tmp <- mp - tmp = 2*Dslash*src */
    sub_wilson_vector(mp + i, tmp + i, tmp + i);
    /* src <- d1/4 * tmp + src */
    scalar_mult_add_wvec(src + i, tmp + i, d1/4., src + i);
  }

  cleanup_dslash_w_3D_temps();
  destroy_wv_field(mp); 
  destroy_wv_field(tmp);
}

/********************************************************************/
/* Construct the source field */
/********************************************************************/

void r_open_w_source(wilson_quark_source *wqs){
  
#ifdef HAVE_QIO
  char *source_file         = wqs->source_file;
  QIO_String *xml_file;
#endif
  char myname[] = "r_open_w_source";

  if(wqs->file_initialized != 0){
    node0_printf("%s: file already opened\n",myname);
    return;
  }
#ifdef HAVE_QIO
  if(wqs->type == COMPLEX_FIELD_FILE){
    xml_file = QIO_string_create();
    wqs->infile = r_open_complex_scidac_file_xml(source_file, QIO_SERIAL,
						 xml_file);
  }
  else if(wqs->type == DIRAC_FIELD_FILE){
    xml_file = QIO_string_create();
    wqs->infile = r_open_w_vector_scidac_file_xml(source_file, QIO_SERIAL,
						  xml_file);
  }
  else{
    node0_printf("%s: a file-type source is required\n",myname);
    return;
  }
  if(wqs->infile == NULL)terminate(1);
  wqs->file_initialized = 1;
  node0_printf("Source file info\n%s\n",QIO_string_ptr(xml_file));
  QIO_string_destroy(xml_file);
#else
  node0_printf("%s: QIO compilation required for this source\n", myname);
  terminate(1);
#endif
}

/* Do covariant derivative on Wilson vector field wv_src.  
   Result in wv_dst. */
static void cov_deriv_wv(wilson_vector *wv_dst, wilson_vector *wv_src,
			 su3_matrix *t_links, int dir, int disp,
			 Real weights[]){
  int i, j;
  msg_tag *tagm, *tagp;
  site *s;
  wilson_vector *wvm,*wvp,*wv;

  /* Check for valid dir */
  if(dir != XUP && dir != YUP && dir != ZUP){
    printf("cov_deriv_wv(%d): Bad direction %d\n",this_node,dir);
    terminate(1);
  }

  wvp  = create_wv_field();
  wvm  = create_wv_field();
  wv   = create_wv_field();

  if(t_links == NULL){
    printf("cov_deriv_wv(%d): NULL t_links\n",this_node);
    terminate(1);
  }

  /* Zero the destination field */
  clear_wv_field(wv_dst);

  /* wvp <- wv_src ;  wvm <- wv_src */
  copy_wv_field(wvp, wv_src);
  copy_wv_field(wvm, wv_src);

  /* Parallel transport over distance disp in direction dir */
  for(j = 0; j < disp; j++){
    /* Start gather from positive direction */
    tagp = start_gather_field( wvp, sizeof(wilson_vector), dir, 
			       EVENANDODD, gen_pt[dir] );

    /* Prepare gather from negative direction */
    /* wv <- U wvm */
    FORALLSITES(i,s){
      mult_adj_mat_wilson_vec( &t_links[4*i+dir], wvm + i, wv + i);
    }

    /* Gather from negative direction */
    tagm = start_gather_field(wv, sizeof(wilson_vector), OPP_DIR(dir), 
			      EVENANDODD, gen_pt[OPP_DIR(dir)] );

    /* Wait for gather from negative direction */
    wait_gather(tagm);

    /* Copy gather result to wvm */
    /* wvm <- shift wv = shift U wvm */
    FORALLSITES(i,s){
      copy_wvec((wilson_vector *)gen_pt[OPP_DIR(dir)][i], wvm + i);
    }

    cleanup_gather(tagm);

    /* Complete gather from positive direction and multiply by link mat */
    wait_gather(tagp);

    /* wv <- U shift wvp */
    FORALLSITES(i,s){
      mult_mat_wilson_vec( &t_links[4*i+dir], 
			   (wilson_vector * )(gen_pt[dir][i]), 
			   wv + i); 
    }

    cleanup_gather(tagp);

    /* Copy result back to wvp */
    /* wvp <- U shift wvp */
    FORALLSITES(i,s){
      copy_wvec(wv + i, wvp + i);
    }

    /* Take the difference and accumulate with weight to complete the
       derivative */
    /* wv_dst <- wv_dst + weights[j](wvp - wvm) */

    if(weights[j] != 0){
      FORALLSITES(i,s){
	sub_wilson_vector(wvp+i, wvm+i, wv+i);
	scalar_mult_add_wvec(wv_dst+i, wv+i, weights[j], wv_dst+i);
      }
    }
  }

  destroy_wv_field(wvp);
  destroy_wv_field(wvm);
  destroy_wv_field(wv);
}

/* Covariant derivative of Wilson propagator field */
static void cov_deriv_wp(wilson_prop_field wp_dst, wilson_prop_field wp_src,
			 su3_matrix *t_links, int dir, int disp, 
			 Real weights[]){
  wilson_vector *wv_src, *wv_dst;
  int spin, color;

  wv_src = create_wv_field();
  wv_dst = create_wv_field();
  for(color = 0; color < 3; color++)
    for(spin=0;spin<4;spin++){
      copy_wv_from_wp(wv_src, wp_src, color, spin);
      cov_deriv_wv(wv_dst, wv_src, t_links, dir, disp, weights);
      copy_wp_from_wv(wp_dst, wv_dst, color, spin);
    }
  
  destroy_wv_field(wv_src);
  destroy_wv_field(wv_dst);
}

/* Do 3D APE smearing of space-like links with SU(3) projection */
void ape_smear_3D(Real staple_weight, int iters){
  int space_only = 1;
  //  int nhits = 0;   /* Turn off SU(3) projection */
  int nhits = 10;   /* Turn off SU(3) projection */
  //  Real tol = 0;    /* Used only for SU(3) projection */
  Real tol = 1e-5;    /* Used only for SU(3) projection */
  int i,dir;
  //  Real link_u0 = 1. - 4.*staple_weight;
  Real link_u0 = 1.;
  su3_matrix *s_links;
  double dtime = -dclock();

  /* Copy site gauge links to ape_links */
  ape_links = create_G_from_site();
  s_links = create_G();

  for(i = 0; i < iters; i++){
    FORALLUPDIRBUT(TUP,dir){
      ape_smear_field_dir(ape_links, dir, s_links, staple_weight, link_u0,
			  space_only, nhits, tol);
    }
    copy_G(ape_links, s_links);
  }

  destroy_G(s_links);
  dtime += dclock();
  node0_printf("Time to APE smear %e sec\n",dtime);
}

void destroy_ape_links_3D(void){
  if(ape_links != NULL)
    destroy_G(ape_links);
  ape_links = NULL;
}

int w_source_field(wilson_vector *src, wilson_quark_source *wqs)
{
  register int i;
  register site *s; 
  
  short my_x,my_y,my_z;
  Real rx,ry,rz,radius2;
  int status = 0;
  
  /* Unpack structure */
  int color                 = wqs->color;
  int spin                  = wqs->spin;
  int source_type           = wqs->type;

  Real a                    = wqs->a;
  int x0                    = wqs->x0; 
  int y0                    = wqs->y0; 
  int z0                    = wqs->z0; 
  int t0                    = wqs->t0;
  Real r0                   = wqs->r0;
  int iters                 = wqs->iters;
  char *source_file         = wqs->source_file;
  int *mom                  = wqs->mom;
#ifdef HAVE_QIO
  int file_initialized      = wqs->file_initialized;
#endif
  
  /* zero src to be safe */
  if(src != NULL){
    FORALLSITES(i,s) { clear_wvec( src+i ); }
  }
  
  /* Unless we are taking the source from storage, clear it before
     reconstructing it */
  if(source_type != COMPLEX_FIELD_STORE &&
     source_type != DIRAC_FIELD_STORE)
    reset_wqs(wqs);
  
  if(source_type == COMPLEX_FIELD_FILE){
#ifdef HAVE_QIO
    if(file_initialized == 0)
      r_open_w_source(wqs);

    alloc_wqs_c_src(wqs);
    status = qio_status(read_complex_scidac(wqs->infile, wqs->c_src, 1 ));

    /* Do momentum insertion */
    if(status == 0)insert_mom_C(wqs->c_src, mom, x0, y0, z0, t0);

    if(src!=NULL)insert_D_from_C(src, wqs->c_src, spin, color);
#else
    node0_printf("w_source_field: QIO compilation required for this source\n");
    terminate(1);
#endif
  }
  else if(source_type == COMPLEX_FIELD_FM_FILE){
    /* Load to the specified spin, color and timeslice */
    /* zero src to be safe */
    alloc_wqs_wv_src(wqs);
    r_source_w_fm_to_field(source_file, wqs->wv_src, spin, color, t0, 
			   source_type);
    /* Do momentum insertion */
    if(status == 0)insert_mom_D(wqs->wv_src, mom, x0, y0, z0, t0);

    /* Save a copy as a complex field */
    alloc_wqs_c_src(wqs);
    extract_C_from_D(wqs->c_src, wqs->wv_src, spin, color);

    /* Copy to requested location */
    if(src != NULL)copy_D(src, wqs->wv_src);
  }
  else if(source_type == COMPLEX_FIELD_STORE){
    if(wqs->c_src == NULL){
      printf("w_source: Can't copy from null field\n");
      terminate(1);
    }
    /* Load to the specified spin, color and timeslice */
    if(src != NULL){
      FORALLSITES(i,s){
	if(t0 == ALL_T_SLICES || s->t == t0)
	  src[i].d[spin].c[color] = wqs->c_src[i];
      }
    }
  }
  else if(source_type == COVARIANT_GAUSSIAN){
    static su3_matrix *t_links;
    alloc_wqs_wv_src(wqs);

    /* Set delta function source */
    if(node_number(x0,y0,z0,t0) == this_node){
      i = node_index(x0,y0,z0,t0);
      wqs->wv_src[i].d[spin].c[color].real = 1.;
    }
    /* Then smear */
    t_links = create_G_from_site();
    gauss_smear_field(wqs->wv_src, t_links, r0, iters, t0);
    destroy_G(t_links);

    /* Copy to requested location */
    if(src != NULL)copy_D(src, wqs->wv_src);
  }
  else if(source_type == DIRAC_FIELD_FILE){
#ifdef HAVE_QIO
    QIO_String *recxml;

    if(file_initialized == 0)
      r_open_w_source(wqs);

    recxml = QIO_string_create();
    /* Load the Dirac vector from the file */
    alloc_wqs_wv_src(wqs);
    status = qio_status(
       read_w_vector_scidac_xml(wqs->infile, wqs->wv_src, 1, recxml));

    /* Do momentum insertion */
    if(status == 0)insert_mom_D(wqs->wv_src, mom, x0, y0, z0, t0);

    /* Verify that the requested color matches the color encoded in
       the record XML */
    if(status == 0)
      status = check_color_spin(recxml, color, spin);

    QIO_string_destroy(recxml);

    /* Copy to requested location */
    if(src != NULL && status == 0)
      copy_D(src, wqs->wv_src);
#else
    node0_printf("w_source_field: QIO compilation required for this source\n");
    terminate(1);
#endif
  }
  else if(source_type == DIRAC_FIELD_FM_FILE){
    alloc_wqs_wv_src(wqs);
    r_source_w_fm_to_field(source_file, wqs->wv_src, spin, color, 
			   t0, source_type);
    /* Do momentum insertion */
    if(status == 0)insert_mom_D(wqs->wv_src, mom, x0, y0, z0, t0);

    /* Copy to requested location */
    if(src!=NULL)copy_D(src, wqs->wv_src);
  }
  else if(source_type == DIRAC_FIELD_STORE){
    if(wqs->wv_src == NULL){
      printf("w_source: Can't copy from null field\n");
      terminate(1);
    }
    /* Ignore spin, color and load the whole wilson vector from storage */
    if(src != NULL){
      FORALLSITES(i,s){
	if(t0 == ALL_T_SLICES || s->t == t0)
	  src[i] = wqs->wv_src[i];
      }
    }
  }
  else if(source_type == FAT_COVARIANT_GAUSSIAN ||
	  source_type == FAT_COVARIANT_GAUSSIAN_DERIV1 ||
	  source_type == FAT_COVARIANT_GAUSSIAN_DERIV2_D ||
	  source_type == FAT_COVARIANT_GAUSSIAN_DERIV2_B ||
	  source_type == DERIV1 ||
	  source_type == DERIV2_D ||
	  source_type == DERIV2_B ||
	  source_type == DERIV3_A ){

    alloc_wqs_wv_src(wqs);

    /* Set delta function source at the requested position in wqs->wv_src */
    if(node_number(x0,y0,z0,t0) == this_node){
      i = node_index(x0,y0,z0,t0);
      wqs->wv_src[i].d[spin].c[color].real = 1.;
    }

    /* Do first derivative if requested */
    if(source_type == FAT_COVARIANT_GAUSSIAN_DERIV1 ||
       source_type == DERIV1){
      wilson_vector *wv_dst = create_wv_field();
      /* wv_dst <- D_dir1 wqs->wv_src */
      cov_deriv_wv(wv_dst, wqs->wv_src, ape_links, wqs->dir1, wqs->disp,
		   wqs->weights);
      /* wqs->wv_src <- wv_dst */
      copy_wv_field(wqs->wv_src, wv_dst);
      destroy_wv_field(wv_dst);
    }
    /* Do second derivative if requested */
    else if(source_type == FAT_COVARIANT_GAUSSIAN_DERIV2_D ||
	    source_type == FAT_COVARIANT_GAUSSIAN_DERIV2_B ||
	    source_type == DERIV2_D ||
	    source_type == DERIV2_B ){
      wilson_vector *wv_dst   = create_wv_field();
      wilson_vector *wv_dst12 = create_wv_field();
      wilson_vector *wv_dst21 = create_wv_field();

      /* Take 21 derivative */
      /* wv_dst <- D_dir1 wqs->wv_src */
      cov_deriv_wv(wv_dst, wqs->wv_src, ape_links, wqs->dir1, wqs->disp,
		   wqs->weights);
      /* wv_dst21 <- D_dir2 D_dir1 wqs->wv_src */
      cov_deriv_wv(wv_dst21, wv_dst, ape_links, wqs->dir2, wqs->disp,
		   wqs->weights);
      /* Take 12 derivative */
      /* wv_dst <- D_dir2 wqs->wv_src */
      cov_deriv_wv(wv_dst, wqs->wv_src, ape_links, wqs->dir2, wqs->disp,
		   wqs->weights);
      /* wv_dst12 <- D_dir1 D_dir2 wqs->wv_src */
      cov_deriv_wv(wv_dst12, wv_dst, ape_links, wqs->dir1, wqs->disp,
		   wqs->weights);
      /* Make symmetric or antisymmetric tensor combinations --
	 result in wqs->wv_src */
      if(source_type == FAT_COVARIANT_GAUSSIAN_DERIV2_D ||
	 source_type == DERIV2_D){
	/* Symmetric tensor component */
	/* wqs->wv_src <- wv_dst12 + wv_dst21 */
	FORALLSITES(i,s){
	  add_wilson_vector(wv_dst12+i, wv_dst21+i, wqs->wv_src+i);
	}
      } else { /* FAT_COVARIANT_GAUSSIAN_DERIV2_B || DERIV2_B */
	/* Antisymmetric tensor component */
	/* wqs->wv_src <- wv_dst12 - wv_dst21 */
	FORALLSITES(i,s){
	  sub_wilson_vector(wv_dst12+i, wv_dst21+i, wqs->wv_src+i);
	}
      }

      destroy_wv_field(wv_dst);
      destroy_wv_field(wv_dst12);
      destroy_wv_field(wv_dst21);
    }

    /* Finally do covariant Gaussian smearing on wqs->wv_src */
    if(source_type == FAT_COVARIANT_GAUSSIAN ||
       source_type == FAT_COVARIANT_GAUSSIAN_DERIV1 ||
       source_type == FAT_COVARIANT_GAUSSIAN_DERIV2_D ||
       source_type == FAT_COVARIANT_GAUSSIAN_DERIV2_B )

      gauss_smear_field(wqs->wv_src, ape_links, r0, iters, t0);
    
    /* Also copy resulting source to requested location */
    if(src != NULL)copy_D(src, wqs->wv_src);
  }
  else if(source_type == GAUSSIAN) {
    /* Gaussian trial source centered on  x0,y0,z0,t0 */
    /* Save a copy of the source in wqs->c_src */

    alloc_wqs_c_src(wqs);

    FORALLSITES(i,s) {
      if(s->t != t0)continue;	/* only do this if t==t0 */
      
      my_x = ((s->x)-x0+nx) % nx;
      rx = (my_x < (nx-my_x)) ? (Real) my_x : (Real) (my_x-nx);
      my_y = ((s->y)-y0+ny) % ny;
      ry = (my_y < (ny-my_y)) ? (Real) my_y : (Real) (my_y-ny);
      my_z = ((s->z)-z0+nz) % nz;
      rz = (my_z < (nz-my_z)) ? (Real) my_z : (Real) (my_z-nz);
      
      radius2 = rx*rx + ry*ry + rz*rz;
      radius2 /= (r0*r0);

      wqs->c_src[i].real = (Real)exp((double)(- radius2));
      if(src != NULL)
	src[i].d[spin].c[color].real = wqs->c_src[i].real;
    }
  }
  else if(source_type == POINT) {
    /* load 1.0 into source at cooordinates given by source_coord */
    /* initialize src to be a delta function at point x0,y0,z0,t0 */
    /* Save a copy of the source in wqs->c_src */

    alloc_wqs_c_src(wqs);

    if(node_number(x0,y0,z0,t0) == mynode()){
      i = node_index(x0,y0,z0,t0);
      wqs->c_src[i].real = 1.0;
      if(src != NULL)src[i].d[spin].c[color].real = 1.0;
    }
  }
  else if(source_type == ROTATE_3D){
    alloc_wqs_wv_src(wqs);

    /* Set delta function source */
    if(node_number(x0,y0,z0,t0) == this_node){
      i = node_index(x0,y0,z0,t0);
      wqs->wv_src[i].d[spin].c[color].real = 1.;
    }
    /* Then do 3D rotation */
    /* The MILC sign convention for gamma matrix in Dslash is
       opposite FNAL, so we rotate with -d1 */
    rotate_3D_wvec(wqs->wv_src, -wqs->d1);

    /* Copy to requested location */
    if(src != NULL)copy_D(src, wqs->wv_src);
  }
  else if(source_type == WAVEFUNCTION_FILE){
    /* Read the file, convert to lattice units, and copy to c_src */
    alloc_wqs_c_src(wqs);
    fnal_wavefunction(wqs->c_src, x0, y0, z0, t0, a, source_file);

    /* Do momentum insertion */
    if(status == 0)insert_mom_C(wqs->c_src, mom, x0, y0, z0, t0);

    /* Load to the specified spin, color and timeslice */
    /* zero src to be safe */
    alloc_wqs_wv_src(wqs);

    insert_D_from_C(wqs->wv_src, wqs->c_src, spin, color);

    /* Copy to requested location */
    if(src != NULL)copy_D(src, wqs->wv_src);
  }
  else if(source_type == UNKNOWN) {
    /* Dummy source */
    alloc_wqs_c_src(wqs);
  }
  else {
    node0_printf("w_source: Unrecognized source type %d\n",source_type);
    terminate(1);
  }

  return status;

} /* w_source_field */

void r_close_w_source(wilson_quark_source *wqs){

  char myname[] = "r_close_w_source";

  if(wqs->file_initialized == 0){
    node0_printf("%s: file has not been opened\n",myname);
    return;
  }
#ifdef HAVE_QIO
  if(wqs->type == COMPLEX_FIELD_FILE){
    r_close_complex_scidac_file(wqs->infile);
    wqs->file_initialized = 0;
    wqs->infile = NULL;
  }
  else if(wqs->type == DIRAC_FIELD_FILE){
    r_close_w_vector_scidac_file(wqs->infile);
    wqs->file_initialized = 0;
    wqs->infile = NULL;
  }
  else
    node0_printf("%s: a file-type source is required\n",myname);
#else
  node0_printf("%s: QIO compilation required for this source\n", myname);
  terminate(1);
#endif
}


void w_open_w_source(wilson_quark_source *wqs, char *fileinfo){

  char myname[] = "w_open_w_source";
#ifdef HAVE_QIO
  char *source_file = wqs->source_file;
  int volfmt, serpar;
#endif

  if(wqs->file_initialized != 0){
    node0_printf("%s: file already opened\n",myname);
    return;
  }
#ifdef HAVE_QIO
  interpret_usqcd_w_save_flag(&volfmt, &serpar, wqs->flag);

  if(wqs->type == COMPLEX_FIELD_FILE){
    wqs->outfile = w_open_complex_scidac_file(source_file, fileinfo,
					      volfmt, serpar);
    wqs->file_initialized = 1;
  }
  if(wqs->type == DIRAC_FIELD_FILE){
    wqs->outfile = w_open_w_vector_scidac_file(source_file, fileinfo,
					       volfmt, serpar);
    wqs->file_initialized = 1;
  }
  else
    node0_printf("%s: a file-type source is required\n",myname);
#else
  node0_printf("%s: QIO compilation required for this source\n", myname);
  terminate(1);
#endif
}

/* For now we support writing only a Dirac vector source field */

int w_source_write(wilson_vector *src, wilson_quark_source *wqs)
{
  char myname[] = "w_source_write";
  
  int status = 0;
  double dtime = 0;
  
  /* Unpack structure */
  int source_type           = wqs->type;
  int file_initialized      = wqs->file_initialized;
#ifdef HAVE_QIO
  int t0                    = wqs->t0;
  int color                 = wqs->color;
  int spin                  = wqs->spin;
  QIO_String *recxml;
  QIO_USQCDPropRecordInfo *recinfo;
#endif
  
  dtime = -dclock();

  if(source_type != DIRAC_FIELD_FILE){
    node0_printf("%s: Unrecognized source type for writing\n", myname);
    return 1;
  }

  if(file_initialized == 0){
    node0_printf("%s: File has not been initialized\n", myname);
    return 1;
  }
#ifdef HAVE_QIO
  recxml = QIO_string_create();

  /* Construct the record XML - we borrow the USQCD prop record XML */
  recinfo = QIO_create_usqcd_proprecord_sc_info(spin, color, "");
  QIO_encode_usqcd_proprecord_info(recxml, recinfo);
  QIO_destroy_usqcd_proprecord_info(recinfo);
  /* Write the Dirac source vector to the file */
  status = qio_status(
      write_wpropsource_D_usqcd_xml(wqs->outfile, recxml, src, t0) );
  node0_printf("Wrote source for color %d spin %d time slice %d\n", 
	       color, spin, t0);
  QIO_string_destroy(recxml);
  dtime += dclock();
  node0_printf("Time to save source spin %d color %d = %e\n",
	       wqs->spin,wqs->color,dtime);
  
#else
  node0_printf("%s: QIO compilation required for this operation\n", myname);
  terminate(1);
#endif
  return status;
}

void w_close_w_source(wilson_quark_source *wqs){

  char myname[] = "w_close_w_source";

  if(wqs->file_initialized == 0){
    node0_printf("%s: file has not been opened\n",myname);
    return;
  }
#ifdef HAVE_QIO
  if(wqs->type == COMPLEX_FIELD_FILE){
    w_close_complex_scidac_file(wqs->outfile);
    wqs->file_initialized = 0;
    wqs->outfile = NULL;
  }
  else if(wqs->type == DIRAC_FIELD_FILE){
    w_close_w_vector_scidac_file(wqs->outfile);
    wqs->file_initialized = 0;
    wqs->outfile = NULL;
  }
  else
    node0_printf("%s: a file-type source is required\n",myname);
#else
  node0_printf("%s: QIO compilation required for this source\n", myname);
  terminate(1);
#endif
}

int w_source_site(field_offset src, wilson_quark_source *wqs)
{
  int i;
  site *s;
  wilson_vector *t_src;
  int status;

#define PAD 0
  t_src  = (wilson_vector *) malloc((sites_on_node+PAD)*sizeof(wilson_vector));
  
  if(t_src == NULL){
    printf("w_source_site(%d): Can't allocate t_src\n",this_node);
    terminate(1);
  }
  
  status = w_source_field(t_src, wqs);

  /* copy src back */
  FORALLSITES(i,s) {
    *(wilson_vector *)F_PT(s,src) = t_src[i];
  }

  free(t_src);
  return status;

} /* w_source_site */


/* 
   Initialize a sink.  The complex sink wave function is replicated on
   each time slice.
 */

void w_sink_field(complex *snk, wilson_quark_source *wqs)
{
  int i;
  site *s; 
  char myname[] = "w_sink_field";
  
  int my_x,my_y,my_z;
  Real rx,ry,rz,radius2;
  
  /* Unpack structure. */
  int sink_type     = wqs->type;
  Real a            = wqs->a;
  int x0            = wqs->x0; 
  int y0            = wqs->y0; 
  int z0            = wqs->z0; 
  int t0            = ALL_T_SLICES;
  Real r0           = wqs->r0;
  char *sink_file   = wqs->source_file;
  int *mom          = wqs->mom;

  if(sink_type == COMPLEX_FIELD_FM_FILE){
    r_source_cmplx_fm_to_field(sink_file, snk, t0, sink_type);

    /* Do momentum insertion */
    insert_mom_C(snk, mom, x0, y0, z0, t0);
  }

  else if(sink_type == COMPLEX_FIELD_FILE){
  /* CAREFUL: This SciDAC file must have an identical sink wave
     function on each time slice */
#ifdef HAVE_QIO
    restore_complex_scidac_to_field(sink_file, QIO_SERIAL, snk, 1);

    /* Do momentum insertion */
    insert_mom_C(snk, mom, x0, y0, z0, t0);
#else
    node0_printf("%s QIO compilation required for this source\n",myname);
    terminate(1);
#endif
  }
  else if(sink_type == GAUSSIAN) {
    /* Gaussian trial sink centered on (x0,y0,z0) of each timeslice */
    
    FORALLSITES(i,s) {
      
      my_x = ((s->x)-x0+nx) % nx;
      rx = (my_x < (nx-my_x)) ? (Real) my_x : (Real) (my_x-nx);
      my_y = ((s->y)-y0+ny) % ny;
      ry = (my_y < (ny-my_y)) ? (Real) my_y : (Real) (my_y-ny);
      my_z = ((s->z)-z0+nz) % nz;
      rz = (my_z < (nz-my_z)) ? (Real) my_z : (Real) (my_z-nz);
      
      radius2 = rx*rx + ry*ry + rz*rz;
      radius2 /= (r0*r0);
      
      snk[i].real = (Real)exp((double)(- radius2));
      snk[i].imag = 0;
    }
  }
  else if(sink_type == POINT) {
    /* load 1.0 into sink at cooordinates given by sink_coord */
    /* initialize snk to be a delta function at point x0,y0,z0 */
    
    for(t0=0;t0<nt;t0++){
      if(node_number(x0,y0,z0,t0) == mynode()){
	i = node_index(x0,y0,z0,t0);
	snk[i].real = 1.0;
      }
    }
  }
  else if(sink_type == WAVEFUNCTION_FILE){
    /* Load the file, convert to lattice units and put in the snk field */
    fnal_wavefunction(snk, x0, y0, z0, t0, a, sink_file);

    /* Do momentum insertion */
    insert_mom_C(snk, mom, x0, y0, z0, t0);
  }
  else{
    node0_printf("%s: Unknown sink type %d\n",myname, sink_type);
    terminate(1);
  }
} /* w_sink_field */


/* 
   Initialize a sink for the inverter, using a wilson number as the
   storage container.
 */

/* snk must be type complex */
void w_sink_site(field_offset snk, wilson_quark_source *wqs)
{
  int i;
  site *s;
  complex *t_snk;

#define PAD 0
  t_snk  = (complex *) malloc((sites_on_node+PAD)*sizeof(complex));
  
  if(t_snk == NULL){
    printf("w_sink_site(%d): Can't allocate t_snk\n",this_node);
    terminate(1);
  }
  
  w_sink_field(t_snk, wqs);

  /* copy snk back */
  FORALLSITES(i,s) {
    *(complex *)F_PT(s,snk) = t_snk[i];
  }

  free(t_snk);

} /* w_sink_site */


/*--------------------------------------------------------------------*/
static void sink_smear_prop(wilson_prop_field wp, wilson_quark_source *wqs){
  
  int color;
  int ci,si,sf,cf;
  int i;
  site *s;
  complex *chi_cs, z;
  Real x;
  spin_wilson_vector *wps;
  double dtime = start_timing();
  int key[4] = {1,1,1,0};  /* 3D Fourier transform */
  
  /* Set up Fourier transform for smearing */
  setup_restrict_fourier(key, NULL);

  chi_cs = (complex *)malloc(sizeof(complex)*sites_on_node);

  /* Now convolute the quark propagator with a given wave function for
     the smeared mesons. This is done with FFT's */
  
  /* fft quark_propagator (in place) */
  for(color = 0; color < 3; color++){
    wps = extract_swv_from_wp(wp, color);
    /* wps points into wp, so wp is changed here */
    restrict_fourier_field((complex *)wps, sizeof(spin_wilson_vector), 
			   FORWARDS);
  }

  print_timing(dtime,"FFT");

  dtime = start_timing();

  /* Build sink smearing wave function as a complex field repeated on
     each time slice */
  w_sink_field(chi_cs, wqs);

  /* Normalize (for FFT) */
  x = 1./(((Real)nx)*((Real)ny)*(Real)nz);
  FORALLSITES(i,s){
    CMULREAL(chi_cs[i],x,chi_cs[i]);
  }
  
  print_timing(dtime,"build sink field");

  dtime = start_timing();
  restrict_fourier_field(chi_cs, sizeof(complex), FORWARDS);
  
  /* Now multiply quark by sink wave function */
  for(ci=0;ci<3;ci++){
    FORALLSITES(i,s)
      for(si=0;si<4;si++)
	for(sf=0;sf<4;sf++)for(cf=0;cf<3;cf++){
	    z = wp[ci][i].d[si].d[sf].c[cf];
	    CMUL(z, chi_cs[i], wp[ci][i].d[si].d[sf].c[cf]);
	  }
  }
  
  print_timing(dtime, "FFT of chi and multiply");

  /* Inverse FFT */
  dtime = start_timing();
  /* fft quark_propagator (in place) */
  for(color = 0; color < 3; color++){
    wps = extract_swv_from_wp(wp, color);
    /* wps points into wp, so wp is changed here */
    restrict_fourier_field((complex *)wps, sizeof(spin_wilson_vector), 
			   BACKWARDS);
  }
  print_timing(dtime,"FFT");
  cleanup_restrict_fourier();
  free(chi_cs);
}  

/*--------------------------------------------------------------------*/
static void rotate_prop_field(wilson_prop_field dst, wilson_prop_field src, 
			      wilson_quark_source *wqs){
  
  int spin, color;
  int i;
  site *s;
  wilson_vector *psi, *mp, *tmp;
  spin_wilson_vector *rp;
  /* The MILC sign convention for gamma matrix in Dslash is
     opposite FNAL, so we rotate with -d1 */
  Real d1 = -wqs->d1;
  
  psi = create_wv_field();
  mp  = create_wv_field();
  tmp = create_wv_field();

  for(color = 0; color < 3; color++){
    rp = dst[color];

    /* Construct propagator for "rotated" fields,
       psi_rot = Dslash psi, with Dslash the naive operator. */
    for(spin=0;spin<4;spin++){
      copy_wv_from_wp(psi, src, color, spin);
      
      /* Do Wilson Dslash on the psi field */
      dslash_w_3D_field(psi, mp,  PLUS, EVENANDODD);
      dslash_w_3D_field(psi, tmp, MINUS, EVENANDODD);
      
      FORALLSITES(i,s){
	/* From subtraction we get 2*Dslash */
	sub_wilson_vector(mp + i, tmp + i, &rp[i].d[spin]);
	/* Apply rotation */
	scalar_mult_add_wvec(psi + i, &rp[i].d[spin], d1/4., &rp[i].d[spin]);
      }
    }
  }
    
  cleanup_dslash_w_3D_temps();
  destroy_wv_field(psi); 
  destroy_wv_field(mp); 
  destroy_wv_field(tmp);
}


/* Covariant Gaussian smearing for a wilson propagator field */

static void gauss_smear_wprop(wilson_prop_field wp, su3_matrix *t_links, 
			      Real r0, int iters, int t0){
  wilson_vector *wv;
  int spin, color;

  wv = create_wv_field();
  for(color = 0; color < 3; color++)
    for(spin=0;spin<4;spin++){
      copy_wv_from_wp(wv, wp, color, spin);
      gauss_smear_field(wv, t_links, r0, iters, t0);
      copy_wp_from_wv(wp, wv, color, spin);
    }
  
  destroy_wv_field(wv);
}


/* 
   Apply a sink operation on the given propagator.
 */

wilson_prop_field w_sink_op(wilson_quark_source *wqs, wilson_prop_field src )
{
  
  /* Unpack structure. */
  int sink_type     = wqs->type;
  wilson_prop_field dst = NULL;

  /* Identity sink operator */
  if(sink_type == POINT) {
    dst = create_wp_field_copy(src);
  }
  /* Various wave-function-based smearings requiring convolution */
  else if(sink_type == GAUSSIAN || 
	  sink_type == COMPLEX_FIELD_FM_FILE ||
	  sink_type == COMPLEX_FIELD_FILE ||
	  sink_type == WAVEFUNCTION_FILE) {
    dst = create_wp_field_copy(src);
    sink_smear_prop(dst, wqs);
  }
  
  else if(sink_type == FAT_COVARIANT_GAUSSIAN ||
	  sink_type == FAT_COVARIANT_GAUSSIAN_DERIV1 ||
	  sink_type == FAT_COVARIANT_GAUSSIAN_DERIV2_D ||
	  sink_type == FAT_COVARIANT_GAUSSIAN_DERIV2_B ||
	  sink_type == DERIV1 ||
	  sink_type == DERIV2_D ||
	  sink_type == DERIV2_B ){

    /* The following operations are performed as in the source,
       but in reverse order */

    /* dst <- src */
    dst = create_wp_field_copy(src);

    /* First do covariant Gaussian smearing on all t slices of dst */
    if(sink_type == FAT_COVARIANT_GAUSSIAN ||
       sink_type == FAT_COVARIANT_GAUSSIAN_DERIV1 ||
       sink_type == FAT_COVARIANT_GAUSSIAN_DERIV2_D ||
       sink_type == FAT_COVARIANT_GAUSSIAN_DERIV2_B)

      gauss_smear_wprop(dst, ape_links, wqs->r0, wqs->iters, ALL_T_SLICES);
    
    /* Do second derivative if requested */
    if(sink_type == FAT_COVARIANT_GAUSSIAN_DERIV2_D ||
       sink_type == FAT_COVARIANT_GAUSSIAN_DERIV2_B ||
       sink_type == DERIV2_D ||
       sink_type == DERIV2_B ){
      wilson_vector *wv_src = create_wv_field();
      wilson_vector *wv_dst = create_wv_field();
      wilson_vector *wv_dst12 = create_wv_field();
      wilson_vector *wv_dst21 = create_wv_field();
      int spin, color;
      int i; site *s;

      for(color = 0; color < 3; color++)
	for(spin=0;spin<4;spin++){
	  /* Take 12 derivative */
	  /* wv_src <- dst_color_spin */
	  copy_wv_from_wp(wv_src, dst, color, spin);
	  /* wv_dst <- D_dir1 wv_src */ 
	  cov_deriv_wv(wv_dst, wv_src, ape_links, wqs->dir1, wqs->disp,
		       wqs->weights);
	  /* wv_dst12 <- D_dir2 D_dir1 wv_src */
	  cov_deriv_wv(wv_dst12, wv_dst, ape_links, wqs->dir2, wqs->disp,
		       wqs->weights);
	  /* Take 21 derivative */
	  /* wv_dst <- D_dir2 wv_src */
	  cov_deriv_wv(wv_dst, wv_src, ape_links, wqs->dir2, wqs->disp,
		       wqs->weights);
	  /* wv_dst21 <- D_dir1 D_dir2 wv_src */
	  cov_deriv_wv(wv_dst21, wv_dst, ape_links, wqs->dir1, wqs->disp,
		       wqs->weights);
	  /* Make symmetric or antisymmetric tensor combinations --
	     result in dst */
	  if(sink_type == FAT_COVARIANT_GAUSSIAN_DERIV2_D ||
	     sink_type == DERIV2_D){
	    /* Symmetric tensor component */
	    /* wv_dst <- wv_dst12 + wv_dst21 */
	    FORALLSITES(i,s){
	      add_wilson_vector(wv_dst12+i, wv_dst21+i, wv_dst+i);
	    }
	  } else { /* FAT_COVARIANT_GAUSSIAN_DERIV2_B || DERIV2_B */
	    /* Antisymmetric tensor component */
	    /* wv_dst <- wv_dst12 - wv_dst21 */
	    FORALLSITES(i,s){
	      sub_wilson_vector(wv_dst12+i, wv_dst21+i, wv_dst+i);
	    }
	  }
	  /* dst_color_spin <- wv_dst */
	  copy_wp_from_wv(dst, wv_dst, color, spin);
      }

      destroy_wv_field(wv_dst);
      destroy_wv_field(wv_dst12);
      destroy_wv_field(wv_dst21);

    } else if(sink_type == FAT_COVARIANT_GAUSSIAN_DERIV1 ||
	      sink_type == DERIV1){
      /* dst <- D_dir1 dst */
      cov_deriv_wp(dst, dst, ape_links, wqs->dir1, wqs->disp,
		   wqs->weights);
    }
  }
  /* FNAL 3D rotation */
  else if(sink_type == ROTATE_3D){
    dst = create_wp_field();
    rotate_prop_field(dst, src, wqs);
  }
  else{
    node0_printf("w_sink_op: Unknown sink type %d\n",sink_type);
    terminate(1);
  }

  return dst;
} /* w_sink_op */


int ask_w_quark_source( FILE *fp, int prompt, int *source_type, char *descrp )
{
  char *savebuf;
  char myname[] = "ask_w_quark_source";

  if (prompt!=0){
    printf("enter ");
    printf("'covariant_gaussian', ");
    printf("'complex_field', ");
    printf("\n     ");
    printf("'complex_field_fm', ");
    printf("'deriv1', ");
    printf("'deriv2_D', ");
    printf("'deriv2_B', ");
    printf("\n     ");
    printf("'deriv3_A', ");
    printf("'dirac_field', ");
    printf("'dirac_field_fm' ");
    printf("\n     ");
    printf("'fat_covariant_gaussian', ");
    printf("'fat_covariant_gaussian_deriv1', ");
    printf("'fat_covariant_gaussian_deriv2_D', ");
    printf("'fat_covariant_gaussian_deriv2_B', ");
    printf("\n     ");
    printf("'gaussian', ");
    printf("'point', ");
    printf("'rotate_3D', ");
    printf("'wavefunction', ");
    printf("\nfor source type\n");
  }

  savebuf = get_next_tag(fp, "quark source command", myname);
  if (savebuf == NULL)return 1;

  if(strcmp("covariant_gaussian",savebuf) == 0 ) {
    *source_type = COVARIANT_GAUSSIAN;
    strcpy(descrp,"covariant_gaussian");
  }
  else if(strcmp("complex_field",savebuf) == 0 ) {
    *source_type = COMPLEX_FIELD_FILE;
    strcpy(descrp,"complex_field");
  }
  else if(strcmp("complex_field_fm",savebuf) == 0 ) {
    *source_type = COMPLEX_FIELD_FM_FILE;
    strcpy(descrp,"complex_field_fm");
  }
  else if(strcmp("deriv1",savebuf) == 0 ) {
    *source_type = DERIV1;
    strcpy(descrp,"deriv1");
  }
  else if(strcmp("deriv2_D",savebuf) == 0 ) {
    *source_type = DERIV2_D;
    strcpy(descrp,"deriv2_D");
  }
  else if(strcmp("deriv2_B",savebuf) == 0 ) {
    *source_type = DERIV2_B;
    strcpy(descrp,"deriv2_B");
  }
  else if(strcmp("deriv3_A",savebuf) == 0 ) {
    *source_type = DERIV3_A;
    strcpy(descrp,"deriv3_A");
  }
  else if(strcmp("dirac_field",savebuf) == 0 ) {
    *source_type = DIRAC_FIELD_FILE;
    strcpy(descrp,"dirac_field");
  }
  else if(strcmp("dirac_field_fm",savebuf) == 0 ) {
    *source_type = DIRAC_FIELD_FM_FILE;
    strcpy(descrp,"dirac_field_fm");
  }
  else if(strcmp("fat_covariant_gaussian",savebuf) == 0 ) {
    *source_type = FAT_COVARIANT_GAUSSIAN;
    strcpy(descrp,"fat_covariant_gaussian");
  }
  else if(strcmp("fat_covariant_gaussian_deriv1",savebuf) == 0 ) {
    *source_type = FAT_COVARIANT_GAUSSIAN_DERIV1;
    strcpy(descrp,"fat_covariant_gaussian_deriv1");
  }
  else if(strcmp("fat_covariant_gaussian_deriv2_D",savebuf) == 0 ) {
    *source_type = FAT_COVARIANT_GAUSSIAN_DERIV2_D;
    strcpy(descrp,"fat_covariant_gaussian_deriv2_D");
  }
  else if(strcmp("fat_covariant_gaussian_deriv2_B",savebuf) == 0 ) {
    *source_type = FAT_COVARIANT_GAUSSIAN_DERIV2_B;
    strcpy(descrp,"fat_covariant_gaussian_deriv2_B");
  }
  else if(strcmp("gaussian",savebuf) == 0 ) {
    *source_type = GAUSSIAN;
    strcpy(descrp,"gaussian");
  }
  else if(strcmp("point",savebuf) == 0 ){
    *source_type = POINT;
    strcpy(descrp,"point");
  }
  else if(strcmp("rotate_3D",savebuf) == 0 ){
    *source_type = ROTATE_3D;
    strcpy(descrp,"rotate_3D");
  }
  else if(strcmp("wavefunction",savebuf) == 0 ){
    *source_type = WAVEFUNCTION_FILE;
    strcpy(descrp,"wavefunction");
  }
  else{
    printf("%s: ERROR IN INPUT: Wilson source command %s is invalid\n",myname,
	   savebuf); 
    return 1;
  }
  
  printf("%s\n",savebuf);
  return 0;
} /* ask_w_quark_source */

#define IF_OK if(status==0)

int ask_output_w_quark_source_file( FILE *fp, int prompt, 
				    int *flag, int *source_type,
				    int *t0, char *descrp, char *filename)
{
  char *savebuf;
  int status = 0;
  char myname[] = "ask_output_w_quark_source_file";

  if (prompt!=0){
    printf("enter 'save_serial_scidac_w_source', or");
    printf("enter 'save_multifile_scidac_w_source', or");
    printf("enter 'save_partfile_scidac_w_source'");
    printf(", for source type\n");
  }

  savebuf = get_next_tag(fp, "output quark source command", myname);
  if (savebuf == NULL)return 1;

  printf("%s ",savebuf);

  if(strcmp("save_serial_scidac_w_source",savebuf) == 0 ) {
#ifdef HAVE_QIO
    *flag=SAVE_SERIAL_SCIDAC;
    *source_type = DIRAC_FIELD_FILE;
    strcpy(descrp,"dirac_field");
#else
    node0_printf("requires QIO compilation!\n");
    terminate(1);
#endif
  }
  else if(strcmp("save_multifile_scidac_w_source",savebuf) == 0 ) {
#ifdef HAVE_QIO
    *flag=SAVE_MULTIFILE_SCIDAC;
    *source_type = DIRAC_FIELD_FILE;
    strcpy(descrp,"dirac_field");
#else
    node0_printf("requires QIO compilation!\n");
    terminate(1);
#endif
  }
  else if(strcmp("save_partfile_scidac_w_source",savebuf) == 0 ) {
#ifdef HAVE_QIO
    *flag=SAVE_PARTFILE_SCIDAC;
    *source_type = DIRAC_FIELD_FILE;
    strcpy(descrp,"dirac_field");
#else
    node0_printf("requires QIO compilation!\n");
    terminate(1);
#endif
  }
  else{
    printf("\n%s: ERROR IN INPUT: Wilson source file command %s is invalid\n",
	   myname, savebuf); 
    return 1;
  }
  
  /* Get file name */
  if( *flag != FORGET ){
    if(prompt!=0)printf("enter filename\n");
    if(scanf("%s",filename) != 1){
      printf("\n%s(%d): ERROR IN INPUT: Can't read filename\n",
	     myname, this_node); 
      status++;
    }
    else
      printf("%s\n",filename);
  }

  /* Get time slice*/
  IF_OK status += get_i(stdin, prompt, "t0", t0);

  return 0;

} /* ask_output_w_quark_source_file */

/* For parsing the derivative direction */
static int decode_dir(char c_dir[]){
  if(strcmp(c_dir,"x")==0)return XUP;
  if(strcmp(c_dir,"y")==0)return YUP;
  if(strcmp(c_dir,"z")==0)return ZUP;
  return NODIR;
}

/* Get the additional input parameters needed to specify the source */
int get_w_quark_source(FILE *fp, int prompt, wilson_quark_source *wqs){
  
  char myname[] = "get_w_quark_source";
  Real d1 = 0;
  char c_dir0[] = " ";
  char c_dir1[2] = " ";
  char *c_dir[2] = {c_dir0, c_dir1};
  int  disp;
  int  source_iters = 0;
  char source_file[MAXFILENAME] = "";
  char source_label[MAXSRCLABEL];
  int  source_loc[4] = { 0,0,0,0 };
  Real source_r0 = 0;
  int  source_type;
  int  status = 0;
  
  /* Get quark source type */
  IF_OK status += ask_w_quark_source(fp, prompt,&source_type,
				   wqs->descrp);
  IF_OK wqs->type  = source_type;
  
  /* Get source parameters */
  IF_OK {
    if ( source_type == COVARIANT_GAUSSIAN ){
      IF_OK status += get_vi(stdin, prompt, "origin", source_loc, 4);
      IF_OK status += get_f(stdin, prompt, "r0", &source_r0);
      IF_OK status += get_i(stdin, prompt, "source_iters", &source_iters);
    }
    else if ( source_type == COMPLEX_FIELD_FILE ){
      IF_OK status += get_vi(stdin, prompt, "origin", source_loc, 4);
      IF_OK status += get_s(stdin, prompt, "load_source", source_file);
      IF_OK status += get_vi(stdin, prompt, "momentum", wqs->mom, 3);
    }
    else if ( source_type == COMPLEX_FIELD_FM_FILE ){
      IF_OK status += get_vi(stdin, prompt, "origin", source_loc, 4);
      IF_OK status += get_s(stdin, prompt, "load_source", source_file);
      IF_OK status += get_vi(stdin, prompt, "momentum", wqs->mom, 3);
    }
    else if( source_type == DERIV1){
      /* Parameters for derivative */
      IF_OK status += get_vs(stdin, prompt, "dir", c_dir, 1);
      IF_OK status += get_i(stdin, prompt, "disp", &disp);
      IF_OK status += get_vf(stdin, prompt, "weights", wqs->weights, disp);
    }
    else if( source_type == DERIV2_D ||
	     source_type == DERIV2_B ){
      /* Parameters for derivatives */
      IF_OK status += get_vs(stdin, prompt, "dir", c_dir, 2);
      IF_OK status += get_i(stdin, prompt, "disp", &disp);
      IF_OK status += get_vf(stdin, prompt, "weights", wqs->weights, disp);
    }
    else if( source_type == DERIV3_A){
      /* Parameters for derivative */
      IF_OK status += get_i(stdin, prompt, "disp", &disp);
      IF_OK status += get_vf(stdin, prompt, "weights", wqs->weights, disp);
    }
    else if ( source_type == DIRAC_FIELD_FILE ){
      IF_OK status += get_i(stdin, prompt, "t0", &source_loc[3]);
      IF_OK status += get_s(stdin, prompt, "load_source", source_file);
    }
    else if ( source_type == DIRAC_FIELD_FM_FILE ){
      IF_OK status += get_i(stdin, prompt, "t0", &source_loc[3]);
      IF_OK status += get_s(stdin, prompt, "load_source", source_file);
    }
    else if( source_type == FAT_COVARIANT_GAUSSIAN){
      /* Parameters for covariant Gaussian */
      IF_OK status += get_vi(stdin, prompt, "origin", source_loc, 4);
      IF_OK status += get_f(stdin, prompt, "r0", &source_r0);
      IF_OK status += get_i(stdin, prompt, "source_iters", &source_iters);

    }
    else if( source_type == FAT_COVARIANT_GAUSSIAN_DERIV1){
      /* Parameters for covariant Gaussian */
      IF_OK status += get_vi(stdin, prompt, "origin", source_loc, 4);
      IF_OK status += get_f(stdin, prompt, "r0", &source_r0);
      IF_OK status += get_i(stdin, prompt, "source_iters", &source_iters);

      /* Parameters for derivative */
      IF_OK status += get_vs(stdin, prompt, "dir", c_dir, 1);
      IF_OK status += get_i(stdin, prompt, "disp", &disp);
      IF_OK status += get_vf(stdin, prompt, "weights", wqs->weights, disp);
    }
    else if( source_type == FAT_COVARIANT_GAUSSIAN_DERIV2_D ||
	     source_type == FAT_COVARIANT_GAUSSIAN_DERIV2_B ){
      /* Parameters for covariant Gaussian */
      IF_OK status += get_vi(stdin, prompt, "origin", source_loc, 4);
      IF_OK status += get_f(stdin, prompt, "r0", &source_r0);
      IF_OK status += get_i(stdin, prompt, "source_iters", &source_iters);

      /* Parameters for derivatives */
      IF_OK status += get_vs(stdin, prompt, "dir", c_dir, 2);
      IF_OK status += get_i(stdin, prompt, "disp", &disp);
      IF_OK status += get_vf(stdin, prompt, "weights", wqs->weights, disp);
    }
    else if ( source_type == GAUSSIAN ){
      IF_OK status += get_vi(stdin, prompt, "origin", source_loc, 4);
      /* width: psi=exp(-(r/r0)^2) */
      IF_OK status += get_f(stdin, prompt,"r0", &source_r0 );
    }
    else if ( source_type == POINT ){
      IF_OK status += get_vi(stdin, prompt, "origin", source_loc, 4);
    }
    else if ( source_type == ROTATE_3D ){
      IF_OK status += get_vi(stdin, prompt, "origin", source_loc, 4);
      IF_OK status += get_f(stdin, prompt, "d1", &d1);
    }
    else if ( source_type == WAVEFUNCTION_FILE ){
      IF_OK status += get_vi(stdin, prompt, "origin", source_loc, 4);
      IF_OK status += get_s(stdin, prompt, "load_source", source_file);
      IF_OK status += get_f(stdin, prompt, "a", &wqs->a);
      IF_OK status += get_vi(stdin, prompt, "momentum", wqs->mom, 3);
    }
    else {
      printf("(%s)Source type not supported in this application\n",myname);
      status++;
    }
  }

  IF_OK status += get_s(stdin, prompt, "source_label", source_label);
  
  wqs->d1            = d1;
  wqs->dir1          = decode_dir(c_dir[0]);
  wqs->dir2          = decode_dir(c_dir[1]);
  wqs->disp          = disp;
  wqs->iters         = source_iters;
  wqs->r0            = source_r0;
  wqs->t0            = source_loc[3];
  wqs->x0            = source_loc[0];
  wqs->y0            = source_loc[1];
  wqs->z0            = source_loc[2];
  strcpy(wqs->source_file,source_file);
  /* Source label '(null)' suppresses the label */
  if(strcmp(source_label,"(null)")==0)source_label[0]='\0';
  strncpy(wqs->label,source_label,MAXSRCLABEL);
  
  return status;
} /* get_w_quark_source */

int get_w_quark_sink(FILE *fp, int prompt, wilson_quark_source *wqs){

  char myname[] = "get_w_quark_sink";
  Real d1 = 0;
  char c_dir0[] = " "; 
  char c_dir1[] = " ";
  char *c_dir[2] = {c_dir0, c_dir1};
  int  disp;
  int  sink_loc[3] = { 0,0,0 };  /* Defaults to zero */
  char sink_file[MAXFILENAME] = "";
  int  sink_iters = 0;
  char sink_label[MAXSRCLABEL];
  Real sink_r0 = 0;
  int  sink_type;
  int  status = 0;

  /* Get antiquark source type */
  IF_OK status += ask_w_quark_source(fp, prompt,&sink_type,
				     wqs->descrp);
  IF_OK wqs->type  = sink_type;

  /* Get parameters for sink */
  IF_OK {
    if ( sink_type == COMPLEX_FIELD_FILE ||
	 sink_type == COMPLEX_FIELD_FM_FILE){
      IF_OK status += get_s(stdin, prompt, "load_sink", sink_file);
      IF_OK status += get_vi(stdin, prompt, "origin", sink_loc, 3);
      IF_OK status += get_vi(stdin, prompt, "momentum", wqs->mom, 3);
    }
    else if( sink_type == DERIV1){
      /* Parameters for derivative */
      IF_OK status += get_vs(stdin, prompt, "dir", c_dir, 1);
      IF_OK status += get_i(stdin, prompt, "disp", &disp);
      IF_OK status += get_vf(stdin, prompt, "weights", wqs->weights, disp);
    }
    else if( sink_type == DERIV2_D ||
	     sink_type == DERIV2_B ){
      /* Parameters for derivatives */
      IF_OK status += get_vs(stdin, prompt, "dir", c_dir, 2);
      IF_OK status += get_i(stdin, prompt, "disp", &disp);
      IF_OK status += get_vf(stdin, prompt, "weights", wqs->weights, disp);
    }
    else if( sink_type == DERIV3_A){
      /* Parameters for derivative */
      IF_OK status += get_i(stdin, prompt, "disp", &disp);
      IF_OK status += get_vf(stdin, prompt, "weights", wqs->weights, disp);
    }
    else if( sink_type == FAT_COVARIANT_GAUSSIAN){
      /* Parameters for covariant Gaussian */
      IF_OK status += get_f(stdin, prompt, "r0", &sink_r0);
      IF_OK status += get_i(stdin, prompt, "source_iters", &sink_iters);

    }
    else if( sink_type == FAT_COVARIANT_GAUSSIAN_DERIV1){
      /* Parameters for covariant Gaussian */
      IF_OK status += get_f(stdin, prompt, "r0", &sink_r0);
      IF_OK status += get_i(stdin, prompt, "source_iters", &sink_iters);

      /* Parameters for derivative */
      IF_OK status += get_vs(stdin, prompt, "dir", c_dir, 1);
      IF_OK status += get_i(stdin, prompt, "disp", &disp);
      IF_OK status += get_vf(stdin, prompt, "weights", wqs->weights, disp);
    }
    else if( sink_type == FAT_COVARIANT_GAUSSIAN_DERIV2_D ||
	     sink_type == FAT_COVARIANT_GAUSSIAN_DERIV2_B ){
      /* Parameters for covariant Gaussian */
      IF_OK status += get_f(stdin, prompt, "r0", &sink_r0);
      IF_OK status += get_i(stdin, prompt, "source_iters", &sink_iters);

      /* Parameters for derivatives */
      IF_OK status += get_vs(stdin, prompt, "dir", c_dir, 2);
      IF_OK status += get_i(stdin, prompt, "disp", &disp);
      IF_OK status += get_vf(stdin, prompt, "weights", wqs->weights, disp);
    }
    else if ( sink_type == GAUSSIAN ){
      /* width: psi=exp(-(r/r0)^2) */
      IF_OK status += get_f(stdin, prompt,"r0", &sink_r0 );
    }
    else if ( sink_type == POINT ){ /* No parameters needed */ }
    else if ( sink_type == ROTATE_3D ){
      IF_OK status += get_f(stdin, prompt, "d1", &d1);
    }
    else if ( sink_type == WAVEFUNCTION_FILE ){
      IF_OK status += get_vi(stdin, prompt, "origin", sink_loc, 3);
      IF_OK status += get_s(stdin, prompt, "load_sink", sink_file);
      IF_OK status += get_f(stdin, prompt, "a", &wqs->a);
      IF_OK status += get_vi(stdin, prompt, "momentum", wqs->mom, 3);
    }
    else {
      printf("(%s)Sink type not supported in this application\n",myname);
      status++;
    }
  }
    
  IF_OK status += get_s(stdin, prompt, "sink_label", sink_label);
  
  wqs->color         = 0;
  wqs->d1            = d1;
  wqs->dir1          = decode_dir(c_dir[0]);
  wqs->dir2          = decode_dir(c_dir[1]);
  wqs->disp          = disp;
  wqs->iters         = sink_iters;
  wqs->r0            = sink_r0;
  wqs->spin          = 0;
  wqs->t0            = ALL_T_SLICES;
  wqs->x0            = sink_loc[0];
  wqs->y0            = sink_loc[1];
  wqs->z0            = sink_loc[2];
  /* For a sink, "source_file" is a misnomer! */
  strcpy(wqs->source_file,sink_file);
  /* Sink label '(null)' suppresses the label */
  if(strcmp(sink_label,"(null)")==0)sink_label[0]='\0';
  strncpy(wqs->label,sink_label,MAXSRCLABEL);
  
  return status;
} /* get_w_quark_sink */


