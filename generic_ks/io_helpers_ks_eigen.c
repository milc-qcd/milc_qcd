/********************** io_helpers_ks_eigen.c *********************************/
/* MIMD version 7 */
/*
  H.Ohno: 11/20/2014 -- dreived from io_helpers_ks.c.

  High level KS eigenvector I/O routines, 
  to be used by any application that wants them.

  Eigenvalus are those of -Dslash^2.
  Eigenvectors can be saved/loaded for EVEN, ODD or EVENANDODD sites.
  The ODD (EVEN) part can be restored as follows:

  eigVec_o(e) = i/sqrt(eigVal) Dslash_oe(eo) eigVec_e(o).

  Norm of a full eigenvector is sqrt(2) since that of the EVEN (ODD) part is unity.  
*/

#include "generic_ks_includes.h"
#include <string.h>

/* Non-SciDAC file formats are deprecated. These non-QIO procedures are kept temporarily */

#ifndef HAVE_QIO

#include "../include/io_ks_eigen.h"

/* Restore the ODD (EVEN) part of KS eigenvectors from the EVEN (ODD) part */
void restore_eigVec(int Nvecs, double *eigVal, su3_vector **eigVec, int parity,
		    imp_ferm_links_t *fn){

  register int i, j;
  complex c;
  su3_vector *ttt = NULL;
  char myname[] = "restore_eigVec";

  ttt = (su3_vector *)malloc(sites_on_node*sizeof(su3_vector));
  if(ttt == NULL){
    printf("%s: Can't malloc ttt\n", myname);
    terminate(1);
  }

  for(i=0; i < Nvecs; i++){
    dslash_fn_field(eigVec[i], ttt, parity, fn);
    FORSOMEFIELDPARITY(j, parity){
      c = cmplx(0.0,1.0/sqrt(eigVal[i]));
      c_scalar_mult_su3vec(ttt + j, &c, eigVec[i] + j) ;
    }
  }

  free(ttt);
}

/*----------------------------------------------------------------*/

/* Open KS eigenvector file for reading eigenvectors */
ks_eigen_file *r_open_ks_eigen(int flag, char *filename){

  ks_eigen_file *kseigf = NULL;
  char myname[] = "r_open_ks_eigen";

  switch(flag){
  case FRESH:
    kseigf = NULL;
    break;
   case RELOAD_ASCII:
    kseigf = r_ascii_ks_eigen_i(filename);
    break;
  case RELOAD_SERIAL:
    kseigf = r_serial_ks_eigen_i(filename);
    break;
  default:
    node0_printf("%s: Unsupported read flag %d\n", myname, flag);
    kseigf = NULL;
  }

  return kseigf;
} /* r_open_ks_eigen */

/*---------------------------------------------------------------*/

/* Open KS eigenvector file for writing eigenvectors. */
ks_eigen_file *w_open_ks_eigen(int flag, char *filename, int parity) {

  ks_eigen_file *kseigf = NULL;
  char myname[] = "w_open_ks_eigen";
  
  switch(flag){
  case FORGET:
    kseigf = NULL;
    break;
  case SAVE_ASCII:
    kseigf = w_ascii_ks_eigen_i(filename, parity);
    break;
  case SAVE_SERIAL:
    kseigf = w_serial_ks_eigen_i(filename, parity);
    break;
  default:
    node0_printf("%s: Unsupported save flag %d\n", myname, flag);
    kseigf = NULL;
  }

  return kseigf;
} /* w_open_ks_eigen */

/*---------------------------------------------------------------*/

/* Close KS eigenvector file for reading eigenvectors. */
void r_close_ks_eigen(int flag, ks_eigen_file *kseigf){

  char myname[] = "r_close_ks_eigen";
  
  if(kseigf == NULL) return;

  switch(flag){
  case RELOAD_ASCII:
    r_ascii_ks_eigen_f(kseigf);
    break;
  case RELOAD_SERIAL:
    r_serial_ks_eigen_f(kseigf);
    break;
  default:
    node0_printf("%s: Unrecognized read flag %d", myname, flag);
  }
} /* r_close_ks_eigne */

/*---------------------------------------------------------------*/

/* Close KS eigenvector file for writing eigenvectors. */
void w_close_ks_eigen(int flag, ks_eigen_file *kseigf){

  char myname[] = "w_close_ks_eigen";
  
  if(kseigf == NULL) return;

  switch(flag){
  case SAVE_ASCII:
    w_ascii_ks_eigen_f(kseigf);
    break;
  case SAVE_SERIAL:
    w_serial_ks_eigen_f(kseigf); 
    break;
  default:
    node0_printf("%s: Unrecognized save flag %d\n", myname, flag);
  }
} /* w_close_ks_eigen */

#endif

#ifdef HAVE_QIO

/* QIO version */

#include "../include/io_scidac_ks.h"

/*---------------------------------------------------------------*/

/* Reload the lowest Nvecs KS eigenvectors:
   FRESH, RELOAD_ASCII, RELOAD_SERIAL
   0 is normal exit code
   >1 for seek, read error, or missing data error 
*/
int reload_ks_eigen(int flag, char *eigfile, int *Nvecs, double *eigVal,
		    su3_vector **eigVec, imp_ferm_links_t *fn, int timing){
  
  register int i, j;
  int status = 0;
  int serpar;
  int qio_status;
  int packed = 0, file_type = 0;
  QIO_Reader *infile;
  double dtime = (double)0.0;
  char myname[] = "reload_ks_eigen";
  
  if(timing && flag != FRESH) dtime = -dclock();
  
  switch(flag){
  case FRESH:
    for(i = 0; i < *Nvecs; i++){
      FORALLFIELDSITES(j){
	clearvec(eigVec[i]+j);
      }
    }
    break;
  case RELOAD_SERIAL:
  case RELOAD_PARALLEL:
    if(flag == RELOAD_SERIAL)serpar = QIO_SERIAL;
    else serpar = QIO_PARALLEL;
    
    infile = open_ks_eigen_infile(eigfile, Nvecs, &packed, &file_type, serpar);
    if(infile == NULL){
      node0_printf("ERROR: Can't open %s for reading\n", eigfile);
      status = 1;
      break;
    }
    /* Read using MILC format */
    if(file_type == 0){
      for(int i = 0; i < *Nvecs; i++){
	qio_status = read_ks_eigenvector(infile, packed, eigVec[i], &eigVal[i]);
	if(qio_status != QIO_SUCCESS){
	  if(qio_status == QIO_EOF){
	    node0_printf("WARNING: Premature EOF at %d eigenvectors\n", i);
	    *Nvecs = i;
	  } else {
	    node0_printf("ERROR: Can't read an eigenvector. Error %d\n", qio_status);
	    status = 1;
	  }
	  break;
	}
      }
    } else {
      /* Read using QUDA format */
      qio_status = read_quda_ks_eigenvectors(infile, eigVec, eigVal, Nvecs, EVEN);
      if(qio_status != QIO_SUCCESS){
	if(qio_status == QIO_EOF){
	  node0_printf("WARNING: Premature EOF at eigenvectors\n");
	  terminate(1);
	} else {
	  node0_printf("ERROR: Can't read eigenvectors. Error %d\n", qio_status);
	  status = 1;
	}
	break;
      }
    }
    close_ks_eigen_infile(infile);

    /* Reconstruct eigenvalues */
    reset_eigenvalues( eigVec, eigVal, *Nvecs, EVEN, fn);

    break;
  default:
    node0_printf("%s: Unrecognized reload flag.\n", myname);
    terminate(1);
  }
  
  if(timing && flag != FRESH){
    dtime += dclock();
    node0_printf("Time to reload %d eigenvectors = %e\n", *Nvecs, dtime);
  }

  return status;
} /* reload_ks_eigen */

#else

/* Custom version */
/*---------------------------------------------------------------*/

/* Reload the lowest Nvecs KS eigenvectors:
   FRESH, RELOAD_ASCII, RELOAD_SERIAL
   0 is normal exit code
   >1 for seek, read error, or missing data error 
*/
int reload_ks_eigen(int flag, char *eigfile, int *Nvecs, double *eigVal,
		    su3_vector **eigVec, imp_ferm_links_t *fn, int timing){

  register int i, j;
  int status = 0;
  ks_eigen_file *kseigf = NULL;
  double dtime = (double)0.0;
  char myname[] = "reload_ks_eigen";

  if(timing && flag != FRESH) dtime = -dclock();

  switch(flag){
  case FRESH:
    for(i = 0; i < *Nvecs; i++){
      FORALLFIELDSITES(j){
	clearvec(eigVec[i]+j);
      }
    }
    break;
  case RELOAD_ASCII:
    kseigf = r_open_ks_eigen(flag, eigfile);
    status = r_ascii_ks_eigen(kseigf, *Nvecs, eigVal, eigVec);
    r_close_ks_eigen(flag, kseigf);
    break;
  case RELOAD_SERIAL:
    kseigf = r_open_ks_eigen(flag, eigfile);
    status = r_serial_ks_eigen(kseigf, *Nvecs, eigVal, eigVec);
    r_close_ks_eigen(flag, kseigf);
    break;
  default:
    node0_printf("%s: Unrecognized reload flag.\n", myname);
    terminate(1);
  }
  
  if(timing && flag != FRESH){
    dtime += dclock();
    node0_printf("Time to reload %d eigenvectors = %e\n", *Nvecs, dtime);
  }

  if(kseigf != NULL && kseigf->parity == EVENANDODD){
    node0_printf("ERROR: EVENANDODD not supported in this code version\n");
    status = 1;
  }

  return status;
} /* reload_ks_eigen */
 
#endif
 
#ifdef HAVE_QIO
/* QIO version */
 
/*---------------------------------------------------------------*/
 
/* Save the lowest Nvecs KS eigenvectors:
   FORGET, SAVE_ASCII, SAVE_SERIAL
*/
 
int save_ks_eigen(int flag, char *savefile, int Nvecs, double *eigVal,
		  su3_vector **eigVec, double *resid, int timing){

  QIO_Writer *outfile;
  int status = 0;
  int serpar;
  int packed;
  double dtime = (double)0.0;
  char myname[] = "save_ks_eigen";
  
  if(timing && flag != FORGET) dtime = -dclock();

  switch(flag){
  case FORGET:
    break;
  case SAVE_SERIAL:
  case SAVE_PARALLEL:
  case SAVE_SERIAL_PACKED:
  case SAVE_PARALLEL_PACKED:

    if(flag == SAVE_SERIAL || flag == SAVE_SERIAL_PACKED)serpar = QIO_SERIAL;
    else serpar = QIO_PARALLEL;

    if(flag == SAVE_SERIAL_PACKED || flag == SAVE_PARALLEL_PACKED)packed = 1;
    else packed = 0;

    outfile = open_ks_eigen_outfile(savefile, Nvecs, QIO_SINGLEFILE, serpar, packed);
    if(outfile == NULL){
      node0_printf("ERROR: Can't open %s for writing\n", savefile);
      status = 1;
      break;
    }

    for(int i = 0; i < Nvecs; i++){
      int status = write_ks_eigenvector(outfile, packed, eigVec[i], eigVal[i], resid[i]);
      if(status != QIO_SUCCESS){
	node0_printf("ERROR: Can't write eigenvector.\n");
	status = 1;
	break;
      }
    }

    close_ks_eigen_outfile(outfile);
    break;

  case SAVE_PARTFILE_SCIDAC:
    outfile = open_ks_eigen_outfile(savefile, Nvecs, QIO_PARTFILE, QIO_SERIAL, 0);
    if(outfile == NULL){
      node0_printf("ERROR: Can't open %s for writing\n", savefile);
      status = 1;
      break;
    }
    for(int i = 0; i < Nvecs; i++){
      int status = write_ks_eigenvector(outfile, 0, eigVec[i], eigVal[i], resid[i]);
      if(status != QIO_SUCCESS){
	node0_printf("ERROR: Can't write eigenvector.\n");
	status = 1;
	break;
      }
    }

    close_ks_eigen_outfile(outfile);
    break;
    
  default:
    node0_printf("%s: Unrecognized save flag.\n", myname);
    terminate(1);
  }
  
  if(timing && flag != FORGET){
    dtime += dclock();
    node0_printf("Time to save %d eigenvectors = %e\n", Nvecs, dtime);
  }

  return status;

} /* save_ks_eigen */

#else
/* Custom version */

/*---------------------------------------------------------------*/

/* Save the lowest Nvecs KS eigenvectors:
   FORGET, SAVE_ASCII, SAVE_SERIAL
*/
int save_ks_eigen(int flag, char *savefile, int Nvecs, double *eigVal,
		   su3_vector **eigVec, double *resid, int timing){

  int status = 0;
  ks_eigen_file *kseigf;
  double dtime = (double)0.0;
  char myname[] = "save_ks_eigen";

  if(timing && flag != FORGET) dtime = -dclock();

  switch(flag){
  case FORGET:
    break;
  case SAVE_ASCII:
    kseigf = w_open_ks_eigen(flag, savefile, EVEN);
    w_ascii_ks_eigen(kseigf, Nvecs, eigVal, eigVec, resid);
    w_close_ks_eigen(flag, kseigf);
    break;
  case SAVE_SERIAL:
    kseigf = w_open_ks_eigen(flag, savefile, EVEN);
    w_serial_ks_eigen(kseigf, Nvecs, eigVal, eigVec, resid);
    w_close_ks_eigen(flag, kseigf);
    break;
  default:
    node0_printf("%s: Unrecognized save flag.\n", myname);
    terminate(1);
  }
  
  if(timing && flag != FORGET){
    dtime += dclock();
    node0_printf("Time to save %d eigenvectors = %e\n", Nvecs, dtime);
  }

  return status;

} /* save_ks_eigen */

#endif

/* We need only even-site eigenvectors, because we can easily
   reconstruct the odd-site eigenvectors from them. To cut the
   required disk storage in half, we pack the even site values into a
   lattice with half the time values.  Then we can read/write half the
   lattice */

/*------------------------------------------------------------------*/
/* Convert rank to coordinates */
static void lex_coords(int coords[], const int dim, const int size[], 
	   const size_t rank)
{
  int d;
  size_t r = rank;

  for(d = 0; d < dim; d++){
    coords[d] = r % size[d];
    r /= size[d];
  }
}

/*------------------------------------------------------------------*/
/* Parity of the coordinate */
static int coord_parity(int r[]){
  return (r[0] + r[1] + r[2] + r[3]) % 2;
}

/*------------------------------------------------------------------*/
/* Convert coordinate to linear lexicographic rank (inverse of
   lex_coords) */

static size_t lex_rank(const int coords[], int dim, int size[])
{
  int d;
  size_t rank = coords[dim-1];

  for(d = dim-2; d >= 0; d--){
    rank = rank * size[d] + coords[d];
  }
  return rank;
}

/*------------------------------------------------------------------*/
/* Map coordinate to index, even sites first */
static size_t my_index(int x, int y, int z, int t, int *latdim) {
    int i = x + latdim[0]*( y + latdim[1]*( z + latdim[2]*t ));
    if( (x+y+z+t)%2==0 ){	/* even site */
	return( i/2 );
    }
    else {
	return( (i + volume)/2 );
    }
}

/*------------------------------------------------------------------*/
/* Map index to coordinates  */
/* (The inverse of my_index) */
/* Assumes even sites come first */
static void my_coords(int coords[], size_t index, int *latdim){
  int ir, xeo;

  /* neven = the number of even sites on the whole lattice */
  size_t neven = volume/2;
  
  /* ir = the even part of the lexicographic index within the
     sublattice on node k */
  if(index < neven){
    ir = 2*index;
    xeo = 0;
  } else {
    ir = 2*(index - neven);
    xeo = 1;
  }

  /* coords = the sublattice coordinate */
  lex_coords(coords, 4, latdim, ir);

  /* Adjust coordinate according to parity */
  if( coord_parity(coords) != xeo ){
    coords[XUP]++;
    if(coords[XUP] >= latdim[XUP]){
      coords[XUP] -= latdim[XUP];
      coords[YUP]++;
      if(coords[YUP] >= latdim[YUP]){
	coords[YUP] -= latdim[YUP];
	coords[ZUP]++;
	if(coords[ZUP] >= latdim[ZUP]){
	  coords[ZUP] -= latdim[ZUP];
	  coords[TUP]++;
	}
      }
    }
  }
}

/* Define map for packing even sites into half the lattice */
static void pack_map_layouts(int x, int y, int z, int t, int *args, int fb,
			     int *xp, int *yp, int *zp, int *tp){

  int latdim[4] = {nx, ny, nz, nt};
  int latdimhalf[4] = {nx, ny, nz, nt/2};
  int coords[4];
  size_t neven = volume/2;
  size_t index;

  if(fb == FORWARDS){
    /* Map for packing even sites into a half lattice with t < nt/2 */
    if((x+y+z+t)%2 == 0){
      /* Even sites for t < nt/2 stay in place */
      if(t < nt/2){
	*xp = x; *yp = y; *zp = z; *tp = t;
      }  else {
	/* Even sites for t >= nt/2 map to odd sites with t < nt/2 */
	*xp = (x + 1)%nx; *yp = y; *zp = z; *tp = t - nt/2;
      }
    } else {
      /* Odd sites for t >= nt/2 stay in place */
      if(t >= nt/2){
	*xp = x; *yp = y; *zp = z; *tp = t;
      } else {
	/* Odd sites for t < nt/2 map to even sites with t >= nt/2 */
	*xp = (x + 1)%nx; *yp = y; *zp = z; *tp = t + nt/2;
      }
    }
  } else {  /* BACKWARDS */
    if((x+y+z+t)%2 == 0){
      /* Even sites for t < nt/2 stay in place */
      if(t < nt/2){
	*xp = x; *yp = y; *zp = z; *tp = t;
      } else {
	/* Even sites for t >= nt/2 map to odd sites for t < nt/2 */
	*xp = (x - 1 + nx)%nx; *yp = y; *zp = z; *tp = t - nt/2;
      }
    } else {
      /* Odd sites for t >= nt/2 stay in place */
      if(t >= nt/2){
	*xp = x; *yp = y; *zp = z; *tp = t;
      } else {
	/* Odd sites for t < nt/2 map to even sites for t >= nt/2 */
	*xp = (x - 1 + nx)%nx; *yp = y; *zp = z; *tp = t + nt/2;
      }
    }
  }
}

int pack_dir;
int unpack_dir;
static int pack_unpack_initialized = 0;

/* Make the packing map */
static void pack_make_gather(void){
  pack_dir =  make_gather(pack_map_layouts, NULL, WANT_INVERSE,
			  ALLOW_EVEN_ODD, SCRAMBLE_PARITY);
  unpack_dir = pack_dir + 1;  /* Convention for the inverse map */
  pack_unpack_initialized = 1;

#if 0
  /* Debug */
  int i;
  FORALLFIELDSITES(i){
    int x = lattice[i].x;
    int y = lattice[i].y;
    int z = lattice[i].z;
    int t = lattice[i].t;
    int xp, yp, zp ,tp;
    pack_map_layouts(x,y,z,t,NULL,FORWARDS,&xp,&yp,&zp,&tp);
    int xpp, ypp, zpp ,tpp;
    pack_map_layouts(xp,yp,zp,tp,NULL,BACKWARDS,&xpp,&ypp,&zpp,&tpp);
    if(xpp != x || ypp != y || zpp != z || tpp != t){
      printf("ERROR ");
      printf("%d %d %d %d -> %d %d %d %d -> %d %d %d %d\n",
	     x, y, z, t, xp, yp, zp, tp, xpp, ypp, zpp, tpp);
    }
  }
#endif
}

/* Packing and unpacking routines -- done in place */
static void pack_unpack_field(void *data, int size, int dir){
  msg_tag *mtag;
  char *temp = (char *)malloc(sites_on_node*size);
  if(temp==NULL){
    printf("pack_field: No room\n");
    terminate(1);
  }
  
  mtag = start_gather_field(data, size, dir, EVENANDODD, gen_pt[0]);
  wait_gather(mtag);

  /* First copy gathered data to temporary */
  for(int i = 0; i < sites_on_node; i++)
    memcpy(temp + size*i, gen_pt[0][i], size);

  cleanup_gather(mtag);

  /* Then copy temp back to field */
  memcpy((char *)data, temp, sites_on_node*size);

  free(temp);
}

void pack_field(void *data, int size){

  if(!pack_unpack_initialized)
    pack_make_gather();
  
  pack_unpack_field(data, size, pack_dir);
}

void unpack_field(void *data, int size){

  if(!pack_unpack_initialized)
    pack_make_gather();
  
  pack_unpack_field(data, size, unpack_dir);
}

/*---------------------------------------------------------------*/
/* Translate output flag to the appropriate input flag for restoring
   KS eigenvectors that were temporarily written to disk
*/
int convert_outflag_to_inflag_ks_eigen(int outflag){

  switch(outflag){
  case SAVE_ASCII:
    return RELOAD_ASCII;
  case SAVE_SERIAL:
    return RELOAD_SERIAL;
  default:
    return FRESH;  /* Error return */
  }
}

/*---------------------------------------------------------------*/

static void print_read_options(void){

  printf("'fresh_ks_eigen', 'reload_ascii_ks_eigen' or 'reload_serial_ks_eigen' or 'reload_parallel_ks_eigen'");
}

/*---------------------------------------------------------------*/

/* Find out what if any KS eigenvectors should be loaded.
   This routine is only called by node 0.
*/
int ask_starting_ks_eigen(FILE *fp, int prompt, int *flag, char *filename){

  char *savebuf;
  int status;
  char myname[] = "ask_starting_ks_eigen";
  
  if(prompt==1){
    printf("enter ");
    print_read_options();
    printf("\n");
  }

  savebuf = get_next_tag(fp, "read ks_eigen command", myname);
  if(savebuf == NULL) return 1;
  
  printf("%s ", savebuf);
  if(strcmp("fresh_ks_eigen", savebuf) == 0){
    *flag = FRESH;
    printf("\n");
  }
  else if(strcmp("reload_ascii_ks_eigen", savebuf) == 0)
    *flag = RELOAD_ASCII;
  else if(strcmp("reload_serial_ks_eigen", savebuf) == 0)
    *flag = RELOAD_SERIAL;
  else if(strcmp("reload_parallel_ks_eigen", savebuf) == 0)
    *flag = RELOAD_PARALLEL;
  else{
    printf("ERROR IN INPUT: ks_eigen input command %s is invalid\n", savebuf);
    return(1);
  }

  /* Read name of file and load it */
  if(*flag != FRESH){
    if(prompt==1) printf("enter name of file containing ks_eigen\n");
    status = scanf("%s", filename);
    if(status != 1){
      printf("\n%s: ERROR IN INPUT: Can't read filename\n", myname);
      return(1);
    }
    printf("%s\n", filename);
  }
  return(0);
  
} /* ask_starting_ks_eigen */

/*--------------------------------------------------------------------*/

static void print_save_options(void){

  printf("'forget_ks_eigen', 'save_ascii_ks_eigen', 'save_serial_ks_eigen', 'save_parallel_ks_eigen', or 'save_partfile_ks_eigen'");
}

/*--------------------------------------------------------------------*/

/* Find out what to do with KS eigenvectors at end.
   This routine is only called by node 0.
*/
int ask_ending_ks_eigen(FILE *fp, int prompt, int *flag, char *filename){

  char *savebuf;
  int status;
  char myname[] = "ask_ending_ks_eigen";

  if(prompt==1){
    print_save_options();
    printf("?\n");
  }

  savebuf = get_next_tag(fp, "write ks_eigen command", myname);
  if(savebuf == NULL)return 1;

  printf("%s ",savebuf);

  if(strcmp("forget_ks_eigen",savebuf) == 0){
    *flag = FORGET;
    printf("\n");
  }
  else if(strcmp("save_ascii_ks_eigen",savebuf) == 0)
    *flag = SAVE_ASCII;
  else if(strcmp("save_serial_ks_eigen",savebuf) == 0)
    *flag = SAVE_SERIAL;
  else if(strcmp("save_parallel_ks_eigen",savebuf) == 0)
    *flag = SAVE_PARALLEL;
  else if(strcmp("save_serial_packed_ks_eigen",savebuf) == 0)
    *flag = SAVE_SERIAL_PACKED;
  else if(strcmp("save_parallel_packed_ks_eigen",savebuf) == 0)
    *flag = SAVE_PARALLEL_PACKED;
  else if(strcmp("save_partfile_ks_eigen",savebuf) == 0)
    *flag = SAVE_PARTFILE_SCIDAC;
  else{
    printf("ERROR IN INPUT: ks_eigen output command %s is invalid.\n", savebuf);
    printf("Choices are ");
    print_save_options();
    printf("\n");
    return(1);
  }
  
  if(*flag != FORGET){
    if(prompt==1) printf("enter filename\n");
    status = scanf("%s", filename);
    if(status != 1){
      printf("\n%s: ERROR IN INPUT: Can't read filename\n", myname); 
      return(1);
    }
    printf("%s\n", filename);
  }
  return(0);
  
} /* ask_ending_ks_eigen */
