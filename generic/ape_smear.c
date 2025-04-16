/********** ape_smear.c *********************************************/
/* MIMD version 7 */

/* Does APE smearing */
/* Doug Toussaint 1/30/96 
   CD 11/25/01 Consolidated with other smearing routines 
   CD 11/15/01 malloc temp field

   Entry point ape_smear_dir for smearing links in one specified direction.
   Entry point ape_smear for smearing links in all four directions.

   Construct smeared links.   Smeared link is sum
   of link plus all staples:

                ---------               -------> dir1
                |       |
                |       |
                |       |               ^
                X--------               | dir2
                |       |               |
                |       |
                |       |
                ---------

   staple_weight is the weight of a staple relative to the original link

   That is we compute

   U_smeared = w_link U_link + w_staple * sum U_staple

   where w_staple/w_link = staple_weight

   w_link + w_staple is normalized so that

   w_link + nstaples * link_u0 * link_u0 * w_staple = 1

   where nstaples = 6 usually, but nstaples = 4 for spacelike
   links when called with space_only = 1.

   If an SU(3) projection is done, only the relative weight matters.

   This scheme takes care of four conventions in MILC use:

   (1) smooth_inst application

       with link_u0 = 1 and 
       staple_weight = ape_weight/(nstaples*(1 - ape_weight))
       we have
       w_link = 1 - ape_weight
       w_staple = ape_weight/nstaples

   (2) spectrum_hybrids4.c

       with link_u0 = u0
       we have
       simple_weight = 1/staple_weight
       norm_factor = 1/(nstaples * u0 * u0 + simple_weight)
       w_link = simple_weight * norm_factor
       w_staple = norm_factor

   (3) hvy_qpot and string_break application

       use 
       link_u0 = 1  (doesn't matter, because we project)
       staple_weight = 1/smear_fac
       space_only = 1

   (4) wilson_hybrids and clover_hybrids application
  
       use
       link_u0 
         = sqrt((1.0 - norm_factor*simple_weight)/(norm_factor*nstaples))
       staple_weight = 1/simple_weight
       space_only = 1
       nhits = 0
*/

#include "generic_includes.h"
#include "../include/openmp_defs.h"

/* Smear in a specified source direction. */
void ape_smear_dir(
  field_offset src,       /* field offset for su3_matrix[4] type 
			     input unsmeared links */
  int dir1,               /* link direction to smear */
  field_offset dest,      /* field offset for su3_matrix type 
			     pointing to a specific direction 
			     output smeared links */
  Real staple_weight,    /* single staple weight */
  Real link_u0,          /* single link weight - used in normalization
                             if SU(3) projection is turned off */
  int space_only,         /* = 1 (true) smear space-like links with
 			          only spacelike staples 
			     = 0 (false) smear all links with
			     all staples */
  int nhits,              /* reproject onto SU(3): number of 
			     SU(2) hits. 0 for no reprojection */
  Real tol               /* tolerance for SU(3) projection.
			     If nonzero, treat nhits as a maximum
			     number of hits.  If zero, treat nhits
			     as a prescribed number of hits. */ 
  )
{
  size_t i;
  int dir2;
  site *s;
  su3_matrix tmat1,tmat2;
  msg_tag *mtag0,*mtag1;
  Real w_link, w_staple, norm_factor;
  int nstaples;
  su3_matrix *temp;
  
  /* Allocate temporary space for staple calculation */
  temp = (su3_matrix *)malloc(sites_on_node*sizeof(su3_matrix));
  if(temp == NULL){
    printf("ape_smear_dir: No room for temp\n");
    terminate(1);
  }
  
  nstaples = (space_only==1 && dir1 != TUP) ? 4 : 6;
  norm_factor = 1.0/(staple_weight*nstaples*link_u0*link_u0 + 1.0);
  w_staple = norm_factor*staple_weight;
  w_link = norm_factor;
  
  /* dest <- src w_link */ 
  FORALLSITES_OMP(i,s,)
    {
      scalar_mult_su3_matrix( 
		     &(((su3_matrix *)F_PT(s,src))[dir1]), w_link, 
		     (su3_matrix *)F_PT(s,dest) );
    } END_LOOP_OMP;
  for(dir2=XUP;dir2<=(space_only==1?ZUP:TUP);dir2++)if(dir2!=dir1){
    
    /* Upper staple, and simple link */
    mtag0 = start_gather_site( src+dir2*sizeof(su3_matrix),
			  sizeof(su3_matrix), dir1, EVENANDODD, gen_pt[0] );
    mtag1 = start_gather_site( src+dir1*sizeof(su3_matrix),
			  sizeof(su3_matrix), dir2, EVENANDODD, gen_pt[1] );
    wait_gather(mtag0);
    wait_gather(mtag1);
    
    /* dest += w_staple * upper staple */ 
    FORALLSITES_OMP(i,s,private(tmat1,tmat2))
      {
	mult_su3_na( (su3_matrix *)gen_pt[1][i],
		     (su3_matrix *)gen_pt[0][i], &tmat1 );
	mult_su3_nn( &(((su3_matrix *)F_PT(s,src))[dir2]), 
		     &tmat1, &tmat2 );
	scalar_mult_add_su3_matrix( 
			   (su3_matrix *)F_PT(s,dest),
			   &tmat2, w_staple,
			   (su3_matrix *)F_PT(s,dest) );
      } END_LOOP_OMP;
    cleanup_gather(mtag0);
    cleanup_gather(mtag1);
    
    /* lower staple */
    mtag0 = start_gather_site( src+dir2*sizeof(su3_matrix),
			  sizeof(su3_matrix), dir1,
			  EVENANDODD, gen_pt[0] );
    wait_gather(mtag0);
    FORALLSITES_OMP(i,s,private(tmat1))
      {
	mult_su3_nn( &(((su3_matrix *)F_PT(s,src))[dir1]),
		     (su3_matrix *)gen_pt[0][i], &tmat1 );
	mult_su3_an( &(((su3_matrix *)F_PT(s,src))[dir2]),
		     &tmat1, &temp[i] );
      } END_LOOP_OMP;
    cleanup_gather(mtag0);
    mtag1 = start_gather_field( temp, sizeof(su3_matrix),
				    OPP_DIR(dir2), EVENANDODD, gen_pt[1] );
    wait_gather(mtag1);
    
    /* dest += w_staple * lower staple */ 
    FORALLSITES_OMP(i,s,){
      scalar_mult_add_su3_matrix( 
			 (su3_matrix *)F_PT(s,dest),
			 (su3_matrix *)gen_pt[1][i], w_staple, 
			 (su3_matrix *)F_PT(s,dest) );
    } END_LOOP_OMP;
    cleanup_gather(mtag1);
    
  } /* dir2 loop */
  
  /* project links onto SU(3) if nhits > 0 */
  if(nhits > 0){
    int status = 0;
    FORALLSITES_OMP(i,s,private(tmat1) reduction(+:status)){
      /* Use partially reunitarized link for guess */
      tmat1 = *((su3_matrix *)F_PT(s,dest));
      reunit_su3(&tmat1);
      status += project_su3(&tmat1,
			    (su3_matrix *)F_PT(s,dest),nhits,tol);
      /* Copy projected matrix to dest */
      *((su3_matrix *)F_PT(s,dest)) = tmat1;
    } END_LOOP_OMP;

    g_intsum(&status);
    if(this_node == 0 && status > 0)
      printf("WARNING %d sites report no convergence in project_su3\n", status);
  }

  free(temp);
  
} /* ape_smear_dir */


void ape_smear(
  field_offset src,       /* field offset for su3_matrix type 
			     input unsmeared links */
  field_offset dest,      /* field offset for su3_matrix type 
			     output smeared links */
  Real staple_weight,    /* single staple weight */
  Real link_u0,          /* single link weight - used in normalization
                             if SU(3) projection is turned off */
  int space_only,         /* = 1 (true) smear space-like links with
 			          only spacelike staples 
			     = 0 (false) smear all links with
			     all staples */
  int nhits,              /* reproject onto SU(3): number of 
			     SU(2) hits. 0 for no reprojection */
  Real tol               /* tolerance for SU(3) projection.
			     If nonzero, treat nhits as a maximum
			     number of hits.  If zero, treat nhits
			     as a prescribed number of hits. */ 
  )
{
  int dir1;
  
  for(dir1=XUP;dir1<=TUP;dir1++){
    ape_smear_dir(src,dir1,dest+dir1*sizeof(su3_matrix),
		  staple_weight,link_u0,space_only,nhits,tol);
  }
} /* ape_smear */

/* Smear in a specified source direction. 
   Gauge matrices are stored four per site. They are taken at
   stride 4 (matrices) starting from the matrix src[dir] 
   The result is stored likewise in dest */

/* NOTE: This code can be unified with the ape_smear_dir code if we
   specify the start and stride and use declare_strided_gather
   throughout - CD */

void ape_smear_field_dir(
  su3_matrix *src,        /* su3_matrix[4] type 
			     input unsmeared links */
  int dir1,               /* link direction to smear */
  su3_matrix *dest,       /* su3_matrix[4] type smeared links */
  Real staple_weight,    /* single staple weight */
  Real link_u0,          /* single link weight - used in normalization
                             if SU(3) projection is turned off */
  int space_only,         /* = 1 (true) smear space-like links with
 			          only spacelike staples 
			     = 0 (false) smear all links with
			     all staples */
  int nhits,              /* reproject onto SU(3): number of 
			     SU(2) hits. 0 for no reprojection */
  Real tol               /* tolerance for SU(3) projection.
			     If nonzero, treat nhits as a maximum
			     number of hits.  If zero, treat nhits
			     as a prescribed number of hits. */ 
  )
{
  size_t i;
  int dir2;
  site *s;
  su3_matrix tmat1,tmat2;
  msg_tag *mtag0,*mtag1;
  Real w_link, w_staple, norm_factor;
  int nstaples;
  su3_matrix *temp;
  
  /* Allocate temporary space for staple calculation */
  temp = (su3_matrix *)malloc(sites_on_node*sizeof(su3_matrix));
  if(temp == NULL){
    printf("ape_smear_field_dir: No room for temp\n");
    terminate(1);
  }
  
  nstaples = (space_only==1 && dir1 != TUP) ? 4 : 6;
  norm_factor = 1.0/(staple_weight*nstaples*link_u0*link_u0 + 1.0);
  w_staple = norm_factor*staple_weight;
  w_link = norm_factor;
  
  /* dest <- src w_link */ 
  FORALLSITES_OMP(i,s,)
    {
      scalar_mult_su3_matrix( &src[4*i+dir1], w_link, &dest[4*i+dir1] );
    } END_LOOP_OMP;
  for(dir2=XUP;dir2<=(space_only==1?ZUP:TUP);dir2++)if(dir2!=dir1){
    
    /* Upper staple, and simple link */
    mtag0 = declare_strided_gather( (char *)&src[dir2], 4*sizeof(su3_matrix),
				    sizeof(su3_matrix), dir1, 
				    EVENANDODD, gen_pt[0] );
    prepare_gather(mtag0);
    do_gather(mtag0);

    mtag1 = declare_strided_gather( (char *)&src[dir1], 4*sizeof(su3_matrix),
				    sizeof(su3_matrix), dir2, 
				    EVENANDODD, gen_pt[1] );
    prepare_gather(mtag1);
    do_gather(mtag1);

    wait_gather(mtag0);
    wait_gather(mtag1);
    
    /* dest += w_staple * upper staple */ 
    FORALLSITES_OMP(i,s,private(tmat1,tmat2))
      {
	mult_su3_na( (su3_matrix *)gen_pt[1][i],
		     (su3_matrix *)gen_pt[0][i], &tmat1 );
	mult_su3_nn( &src[4*i+dir2], &tmat1, &tmat2 );
	scalar_mult_add_su3_matrix( &dest[4*i+dir1], &tmat2, w_staple,
				    &dest[4*i+dir1] );
      } END_LOOP_OMP;
    cleanup_gather(mtag0);
    cleanup_gather(mtag1);
    
    /* lower staple */
    mtag0 = declare_strided_gather( (char *)&src[dir2], 4*sizeof(su3_matrix),
				    sizeof(su3_matrix), dir1,
				    EVENANDODD, gen_pt[0] );
    prepare_gather(mtag0);
    do_gather(mtag0);

    wait_gather(mtag0);
    FORALLSITES_OMP(i,s,private(tmat1))
      {
	mult_su3_nn( &src[4*i+dir1], (su3_matrix *)gen_pt[0][i], &tmat1 );
	mult_su3_an( &src[4*i+dir2], &tmat1, &temp[i] );
      } END_LOOP_OMP;
    cleanup_gather(mtag0);
    mtag1 = start_gather_field( temp, sizeof(su3_matrix),
				    OPP_DIR(dir2), EVENANDODD, gen_pt[1] );
    wait_gather(mtag1);
    
    /* dest += w_staple * lower staple */ 
    FORALLSITES_OMP(i,s,){
      scalar_mult_add_su3_matrix( &dest[4*i+dir1], 
			 (su3_matrix *)gen_pt[1][i], w_staple, 
			 &dest[4*i+dir1] );
    } END_LOOP_OMP;
    cleanup_gather(mtag1);
    
  } /* dir2 loop */
  
  /* project links onto SU(3) if nhits > 0 */
  if(nhits > 0){
    int status = 0;
    FORALLSITES_OMP(i,s,private(tmat1)){
      /* Use partially reunitarized link for guess */
      tmat1 = dest[4*i+dir1];
      reunit_su3(&tmat1);
      status += project_su3(&tmat1, &dest[4*i+dir1], nhits, tol);
      /* Copy projected matrix to dest */
      dest[4*i+dir1] = tmat1;
    } END_LOOP_OMP;
    g_intsum(&status);
    if(this_node == 0 && status > 0)
      printf("WARNING %d sites report no convergence in project_su3\n", status);
  }
  
  free(temp);
  
} /* ape_smear_dir */

/* Input field has four contigous SU(3) matrices per site */

void ape_smear_field(
  su3_matrix *src,       /* Gauge field input unsmeared */
  su3_matrix *dest,      /* Gauge field output smeared */
  Real staple_weight,    /* single staple weight */
  Real link_u0,          /* single link weight - used in normalization
                             if SU(3) projection is turned off */
  int space_only,         /* = 1 (true) smear space-like links with
 			          only spacelike staples 
			     = 0 (false) smear all links with
			     all staples */
  int nhits,              /* reproject onto SU(3): number of 
			     SU(2) hits. 0 for no reprojection */
  Real tol               /* tolerance for SU(3) projection.
			     If nonzero, treat nhits as a maximum
			     number of hits.  If zero, treat nhits
			     as a prescribed number of hits. */ 
  )
{
  int dir1;
  
  for(dir1=XUP;dir1<=TUP;dir1++){
    ape_smear_field_dir(src,dir1,dest,staple_weight,
			link_u0,space_only,nhits,tol);
  }
} /* ape_smear_field */

/* Unit gauge field. See also io_helpers.c:coldlat() */
static su3_matrix *create_unit_gauge_field(void){
  su3_matrix *ape_links;
  size_t i;
  int j,dir;

  /* Allocate and zero the field */
  ape_links = create_G();

  /* Set matrices to unity */
  FORALLFIELDSITES_OMP(i,private(dir)){
    FORALLUPDIRBUT(TUP,dir){
      for(j = 0; j < 3; j++)
	ape_links[4*i+dir].e[j][j].real = 1.0;
    }
  } END_LOOP_OMP;
  return ape_links;
}

/* Do 3D APE smearing of space-like links with SU(3) projection */
/* The time links are copied from the underlying gauge field */
/* NOTE A HACK: If iters < 0 we return a unit gauge field instead */

su3_matrix *ape_smear_3D(Real staple_weight, int iters){
  int space_only = 1;
  //  int nhits = 0;   /* Turn off SU(3) projection */
  int nhits = 10;
  //  Real tol = 0;    /* Used only for SU(3) projection */
  Real tol = 1e-5;    /* Used only for SU(3) projection */
  int i,dir;
  //  Real link_u0 = 1. - 4.*staple_weight;
  Real link_u0 = 1.;
  su3_matrix *s_links;
  su3_matrix *ape_links;
  double dtime = -dclock();

  /* KLUDGE! So we get noncovariant shifts in quark_source_sink_op.c */
  /* A negative value results in a unit gauge field */
  if(iters < 0){
    ape_links = create_unit_gauge_field();
    node0_printf("APE links are set to unit matrices\n");
    return ape_links;
  }

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
  return ape_links;
} /* ape_smear_3D */

/* Do 4D APE smearing of space-like links with SU(3) projection */
/* The time links are copied from the underlying gauge field */
/* NOTE A HACK: If iters < 0 we return a unit gauge field instead */

su3_matrix *ape_smear_4D(Real staple_weight, int iters){
  int space_only = 1;
  //  int nhits = 0;   /* Turn off SU(3) projection */
  int nhits = 10;
  //  Real tol = 0;    /* Used only for SU(3) projection */
  Real tol = 1e-5;    /* Used only for SU(3) projection */
  int i,dir;
  //  Real link_u0 = 1. - 4.*staple_weight;
  Real link_u0 = 1.;
  su3_matrix *s_links;
  su3_matrix *ape_links;
  double dtime = -dclock();

  /* KLUDGE! So we get noncovariant shifts in quark_source_sink_op.c */
  /* A negative value results in a unit gauge field */
  if(iters < 0){
    ape_links = create_unit_gauge_field();
    node0_printf("APE links are set to unit matrices\n");
    return ape_links;
  }

  /* Copy site gauge links to ape_links */
  ape_links = create_G_from_site();
  s_links = create_G();

  for(i = 0; i < iters; i++){
    FORALLUPDIR(dir){
      ape_smear_field_dir(ape_links, dir, s_links, staple_weight, link_u0,
			  space_only, nhits, tol);
    }
    copy_G(ape_links, s_links);
  }

  destroy_G(s_links);
  dtime += dclock();
  node0_printf("Time to APE smear %e sec\n",dtime);
  return ape_links;
} /* ape_smear_4D */

void destroy_ape_links_3D(su3_matrix *ape_links){
  if(ape_links != NULL)
    destroy_G(ape_links);
} /* destroy_ape_links_3D */

void destroy_ape_links_4D(su3_matrix *ape_links){
  if(ape_links != NULL)
    destroy_G(ape_links);
} /* destroy_ape_links_4D */


#ifdef APE_LINKS_FILE
/* find out how to get the input APE link file to use.  This routine
   is only called by node 0. */
int ask_starting_apelinks( FILE *fp, int prompt, int *flag, char *filename ){
  const char *savebuf;
  int status;
  const char myname[] = "ask_starting_apelinks";
  
  if (prompt==1) printf("enter 'fresh_ape' or 'reload_serial_ape'\n");
  
  savebuf = get_next_tag(fp, "read APE links command", myname);
  if (savebuf == NULL)return 1;
  
  printf("%s ",savebuf);
  if(strcmp("fresh_ape",savebuf) == 0 ) {
    *flag = FRESH;
    printf("\n");
  }
  else if(strcmp("reload_serial_ape",savebuf) == 0 ) {
    *flag = RELOAD_SERIAL;
  }
  else{
    printf(" is not a valid starting lattice command. INPUT ERROR.\n"); 
    return 1;
  }
  
  /*read name of file and load it */
  if( *flag != FRESH ){
    if(prompt==1)printf("enter name of file containing lattice\n");
    status=fscanf(fp," %s",filename);
    if(status !=1) {
      printf("\n%s(%d): ERROR IN INPUT: error reading file name\n",
	     myname, this_node); 
      return 1;
    }
    printf("%s\n",filename);
  }
  return 0;
}

/* find out how to save APE links.  This routine is only called by
   node 0.
*/
int ask_ending_apelinks(FILE *fp, int prompt, int *flag, char *filename ){
  const char *savebuf;
  int status;
  const char myname[] = "ask_ending_lattice";
  
  if (prompt==1) printf(
			"'forget_ape' links at end, 'save_serial_scidac_ape', 'save_parallel_scidac_ape', 'save_partfile_scidac_ape', 'save_serial_scidac_ape_dp', 'save_parallel_scidac_ape_dp', 'save_partfile_scidac_ape_dp'\n");
  
  savebuf = get_next_tag(fp, "save APE file command", myname);
  if (savebuf == NULL)return 1;

  printf("%s ",savebuf);
  if(strcmp("save_serial_ape",savebuf) == 0 ) {
    *flag=SAVE_SERIAL;
  }
  else if(strcmp("save_parallel_ape",savebuf) == 0 ) {
    *flag=SAVE_PARALLEL;
  }
  else if(strcmp("save_serial_scidac_ape",savebuf) == 0 ) {
#ifdef HAVE_QIO
    *flag=SAVE_SERIAL_SCIDAC;
#else
    node0_printf("requires QIO compilation!\n");
    terminate(1);
#endif
  }
  else if(strcmp("save_serial_scidac_ape_dp",savebuf) == 0 ) {
#ifdef HAVE_QIO
    *flag=SAVE_SERIAL_SCIDAC_DP;
#else
    node0_printf("requires QIO compilation!\n");
    terminate(1);
#endif
  }
  else if(strcmp("save_parallel_scidac_ape",savebuf) == 0 ) {
#ifdef HAVE_QIO
    *flag=SAVE_PARALLEL_SCIDAC;
#else
    node0_printf("requires QIO compilation!\n");
    terminate(1);
#endif
  }
  else if(strcmp("save_parallel_scidac_ape_dp",savebuf) == 0 ) {
#ifdef HAVE_QIO
    *flag=SAVE_PARALLEL_SCIDAC_DP;
#else
    node0_printf("requires QIO compilation!\n");
    terminate(1);
#endif
  }
  else if(strcmp("save_partfile_scidac_ape",savebuf) == 0 ) {
#ifdef HAVE_QIO
    *flag=SAVE_PARTFILE_SCIDAC;
#else
    node0_printf("requires QIO compilation!\n");
    terminate(1);
#endif
  }
  else if(strcmp("save_partfile_scidac_dp",savebuf) == 0 ) {
#ifdef HAVE_QIO
    *flag=SAVE_PARTFILE_SCIDAC_DP;
#else
    node0_printf("requires QIO compilation!\n");
    terminate(1);
#endif
  }
  else if(strcmp("forget_ape",savebuf) == 0 ) {
    *flag=FORGET;
    printf("\n");
  }
  else {
    printf("is not a save lattice command. INPUT ERROR\n");
    return 1;
  }
  
  if( *flag != FORGET ){
    if(prompt==1)printf("enter filename\n");
    status=fscanf(fp,"%s",filename);
    if(status !=1){
      printf("\nask_ending_apelinks: ERROR IN INPUT: error reading filename\n"); return 1;
    }
    printf("%s\n",filename);
    
  }
  return 0;
}

#endif
