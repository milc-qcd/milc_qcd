#include "generic_ks_includes.h"
#include "../include/openmp_defs.h"

static inline void _scalar_mult_su3_vector(su3_vector *a, Real s, su3_vector *c)
{
	c->c[0].real = s * a->c[0].real;
	c->c[0].imag = s * a->c[0].imag;
	c->c[1].real = s * a->c[1].real;
	c->c[1].imag = s * a->c[1].imag;
	c->c[2].real = s * a->c[2].real;
	c->c[2].imag = s * a->c[2].imag;
}

static inline void _su3vec_copy( su3_vector *a, su3_vector *b ){
    b->c[0].real = a->c[0].real;
    b->c[0].imag = a->c[0].imag;
    b->c[1].real = a->c[1].real;
    b->c[1].imag = a->c[1].imag;
    b->c[2].real = a->c[2].real;
    b->c[2].imag = a->c[2].imag;
}

#ifdef GB_BARYON

/** Get a Fourier factor for a site */
static complex
ff(Real theta, char parity, complex tmp)
{
  complex z = {0.,0.};

  if(parity == EVEN){
    z.real =  tmp.real*cos(theta);
    z.imag =  tmp.imag*cos(theta);
  }
  else if(parity == ODD){
    z.real = -tmp.imag*sin(theta);
    z.imag =  tmp.real*sin(theta);
  }
  else if(parity == EVENANDODD){
    z.real =  tmp.real*cos(theta)-tmp.imag*sin(theta);
    z.imag =  tmp.imag*cos(theta)+tmp.real*sin(theta);
  }
  else{
    printf("ff(%d): bad parity %d\n", this_node, parity);
    terminate(1);
  }
  return z;
} /* ff */

/*---------------------------------------------------------------------*/
/**
   Normalize a correlator according to the phase and factor given in the
   correlator line of the input.
 */
static void
norm_corr(int phase, Real fact, complex prop[])
{
  int k;
  complex z = {0.,0.};

  for(k=0; k<nt; k++){
    switch(phase){
    case 0:
      z =            prop[k];
      break;
    case 1:
      TIMESPLUSI(    prop[k], z);
      break;
    case 2:
      TIMESMINUSONE( prop[k], z);
      break;
    case 3:
      TIMESMINUSI(   prop[k], z);
    }
    CMULREAL(z,fact,prop[k]);
  }
}

/*------------------------------------------------------------------*/
// convert to sensible representation
static int gamma_bits[16] = {1,2,4,8,15,6,5,3,9,10,12,14,13,11,7,0};

/*------------------------------------------------------------------*/
/* slightly different implementation to account for 2-point function tie-ups */
static int
decode_gamma_taste_bits(int index){
  if (index == GB_2POINT_BACKPROP){ return 0; } // must act as scalar-scalar
  return gamma_bits[(index-128)%16];
}

/*------------------------------------------------------------------*/
static int
decode_gamma_spin_bits(int index){
  if (index == GB_2POINT_BACKPROP){ return 0; } // must act as scalar-scalar
  return gamma_bits[(index-128)/16];
}

/**
  Decode the sign factor for the given spin-taste current.
  Probably more complicated than it needs to be.
  */
static int
spin_taste_phase_factor_bits(int stIdx){
  int spin,taste,nlink;
  int shiftbit,linkbit,signbit;

  spin = decode_gamma_spin_bits(stIdx)%8;
  taste = decode_gamma_taste_bits(stIdx)%8;
  shiftbit = (decode_gamma_taste_bits(stIdx)/8)%2; // gamma4 only
  if(shiftbit) { signbit = 7; } // -1^(x+y+z)
  else         { signbit = 0; }

  linkbit = spin ^ taste; // bitwise xor
  nlink = (linkbit/4)%2 + (linkbit/2)%2 + linkbit%2;
  if(nlink > 1){
    // multi-link partner
    signbit = signbit ^ 2; // *-1^y
  }

  if(taste == 0 || taste == 7){
    // S/A_4 or V_4/P taste
    if(nlink == 1 || nlink == 2){
      // 3-dimensional irrep, eta factor
      switch(spin){
        case 3: case 4:
          signbit = signbit ^ 3; break; // \eta_z
        case 5: case 2:
          signbit = signbit ^ 1; break; // \eta_y
        case 6: case 1:
        default: break; // \eta_x
      }
    }
  }
  else{
    // T_jk/V_i or A_i/T_i4 taste
    if(nlink == 0 || nlink == 3){
      // 3-dimensional irrep, -1^(x_i)\eta_4\zeta_4 factor
      switch(spin){
        case 3: case 4:
          signbit = signbit ^ 4; break; // -1^z
        case 5: case 2:
          signbit = signbit ^ 2; break; // -1^y
        case 6: case 1:
          signbit = signbit ^ 1; break; // -1^x
        default: break; // 1
      }
      signbit = signbit ^ 0x7;
    }
    else{
      if(spin == 0 || spin == 7){
        // 3-dimensional irrep, zeta factor
        switch(taste){
          case 6: case 1:
            signbit = signbit ^ 6; break; // \zeta_x
          case 5: case 2:
            signbit = signbit ^ 4; break; // \zeta_y
          case 3: case 4:
          default: break; // \zeta_z
        } // switch taste
      }
      else {
        // 6-dimensional irrep, eta_j * -1^(x_i)\eta_4\zeta_4
        switch(spin){
          case 3: case 4:
            signbit = signbit ^ 3; break; // \eta_z
          case 5: case 2:
            signbit = signbit ^ 1; break; // \eta_y
          case 6: case 1:
          default: break; // \eta_x
        }
        switch(taste){
          case 3: case 4:
            signbit = signbit ^ 4; break; // -1^z
          case 5: case 2:
            signbit = signbit ^ 2; break; // -1^y
          case 6: case 1:
            signbit = signbit ^ 1; break; // -1^x
          default: break; // 1
        }
        signbit = signbit ^ 7;

      } // if spin = 0,7
    } // else nlink != 0,3
  } //else taste != 0,7
  return signbit;
}

/**
  Return the +-1 sign associated with the phase factor and cube corner
  */
static Real
spin_taste_sign(int signbit,int corner){
  int phsbit = signbit & corner; //bitwise and

  //node0_printf("\nMy sign bit is %d\n", signbit);
  //node0_printf("My phase bit is %d\n", phsbit);
  int nbits = (phsbit/4)%2 + (phsbit/2)%2 + phsbit%2;
  //node0_printf("My nbits is %d\n", nbits);
  if(nbits %2) {
      return -1.;
  }
  else {
      return  1.;
  }
}

/*------------------------------------------------------------------*/
/**
   Symmetric shift a quark annihilation sink.
 */
/* Apply the symmetric shift operator in direction "dir" *
 * This is the explicit version                          *
 * The KS phases MUST BE in the links                    */
static void
sym_shift_sink(int dir, su3_vector *dest, su3_vector *src, su3_matrix *links)
{
  int i;
  msg_tag *tag[2];
  su3_vector *tvec = create_v_field();
  node0_printf("Entering sym_shift_sink\n");


  tag[0] = start_gather_field( src, sizeof(su3_vector), dir, EVENANDODD, gen_pt[0] );
  /* With ONE_SIDED_SHIFT_GB defined, the shift is asymmetric */
#ifndef ONE_SIDED_SHIFT_GB
#ifdef NO_SINK_LINKS
  FORALLFIELDSITES_OMP(i,) { _su3vec_copy(src+i, tvec+i); } END_LOOP_OMP

#else
  FORALLFIELDSITES_OMP(i,)
    {
      /* Link conjugacy should be flipped so quark (creation op) *
       * has backward link in forward dir                        */
      mult_adj_su3_mat_vec( links+4*i+dir, src+i, tvec+i );
    } END_LOOP_OMP
#endif // NO_SINK_LINKS
  tag[1] = start_gather_field( tvec, sizeof(su3_vector), OPP_DIR(dir), EVENANDODD, gen_pt[1] );
#endif // ONE_SIDED_SHIFT_GB
  wait_gather(tag[0]);
#ifdef NO_SINK_LINKS
  FORALLFIELDSITES_OMP(i,) { _su3vec_copy((su3_vector *)gen_pt[0][i],dest+i); } END_LOOP_OMP

#else
  FORALLFIELDSITES_OMP(i,)
    {
      mult_su3_mat_vec( links+4*i+dir, (su3_vector *)gen_pt[0][i], dest+i );
    } END_LOOP_OMP
#endif // NO_SINK_LINKS
#ifndef ONE_SIDED_SHIFT_GB
  wait_gather(tag[1]);
  FORALLFIELDSITES_OMP(i,)
    {
      add_su3_vector(dest+i, (su3_vector*)gen_pt[1][i], dest+i );
   
      /* Now divide by 2 eq. (4.2b) of Golterman's Meson paper*/
      _scalar_mult_su3_vector( dest+i, .5, dest+i ) ;
    } END_LOOP_OMP
  cleanup_gather(tag[1]);
#endif
  cleanup_gather(tag[0]);
  destroy_v_field(tvec);
}

/*------------------------------------------------------------------*/
static void
apply_sym_shift_snk_v(int n, int *d, int *r0, ks_prop_field *dest,
                      ks_prop_field *src, su3_matrix *links){
  /* Apply the symmetric shifts listed in d             *
   * Similar to zeta_shift_field in spin_taste_ops.c    *
   * but with different sign factors                    */
  int i,c;
  ks_prop_field *tmp0 = create_ksp_field(3);

  for(i=0;i<n;i++)
  {
    /* Do the shift in d[i] */
    if(i==0)
      /* first time from source */
      for(c=0;c<3;c++) sym_shift_sink(d[i], dest->v[c], src->v[c], links);
    else {
      /* other times from dest */
      for(c=0;c<3;c++) sym_shift_sink(d[i], tmp0->v[c], dest->v[c], links);
      /* tmp0->dest */
      copy_ksp_field(dest,tmp0);
    }
  }

  destroy_ksp_field(tmp0);
}

/*------------------------------------------------------------------*/
void
apply_par_xport_snk_v(ks_prop_field *dest, ks_prop_field *src,
                      int n, int dir[], int r0[], su3_matrix *links){
  ks_prop_field *tvec0 = create_ksp_field(3);
  ks_prop_field *tvec1 = create_ksp_field(3);
  ks_prop_field *tsrc  = create_ksp_field(3);
  int d[6][3] =
    {{XUP,YUP,ZUP},
     {YUP,ZUP,XUP},
     {ZUP,XUP,YUP},
     {XUP,ZUP,YUP},
     {YUP,XUP,ZUP},
     {ZUP,YUP,XUP}};
  copy_ksp_field(tsrc,src);

  if(n == 0){
    copy_ksp_field(dest,src);
  }
  else if(n == 1){
    /* one link */
    d[0][0] = dir[0];
    apply_sym_shift_snk_v(n,d[0],r0,dest,tsrc,links);
  }
  else if (n == 2){
    /* two link */
    d[0][0] = dir[0]; d[0][1] = dir[1];
    d[1][1] = dir[0]; d[1][0] = dir[1];
    apply_sym_shift_snk_v(n,d[0],r0,tvec0,tsrc,links);
    apply_sym_shift_snk_v(n,d[1],r0,tvec1,tsrc,links);

    #pragma omp parallel for collapse(2)
    for(int c=0;c<3;c++){
      for(int i=0; i<sites_on_node; i++){
        add_su3_vector( &tvec0->v[c][i], &tvec1->v[c][i], &dest->v[c][i] );
        _scalar_mult_su3_vector( &dest->v[c][i], 0.5, &dest->v[c][i] );
      }
    } // END OMP LOOPS
  }
  else if (n == 3){
    /* three link */
    /* use the d given */
    //node0_printf("Entering n==3");
    for(int j=0;j<6;j++){
      apply_sym_shift_snk_v(n,d[j],r0,tvec0,tsrc,links);
      if(j==0) {
	copy_ksp_field(tvec1,tvec0);
      } else {
        #pragma omp parallel for collapse(2)
	for(int c=0;c<3;c++){
	  for(int i=0; i<sites_on_node; i++){
	    add_su3_vector(&tvec1->v[c][i],&tvec0->v[c][i],&tvec1->v[c][i]);
	  }
	} // END OMP LOOPS
      }
    }
    #pragma omp parallel for collapse(2)
    for(int c=0;c<3;c++){
      for(int i=0; i<sites_on_node; i++){
	_scalar_mult_su3_vector( &tvec1->v[c][i], 1./6., &dest->v[c][i] );
      }
    } // END OMP LOOPS
  }
  destroy_ksp_field(tvec1);
  destroy_ksp_field(tvec0);
  destroy_ksp_field(tsrc);
}

/**
   Create the ks_prop_field if mmap is not being used, otherwise fetch existing.
   Set remap to true to reuse an existing pointer.
   */
static void
map_ksp_field(ks_prop_field **dest, ks_prop_field **src, int qknum,
              int siSrc, int siSnk, int r0[], su3_matrix *links, short remap){
#ifdef GB_BARYON_MMAP
  if(remap){ toss_ksp_from_cache(dest); }
  /* get object from mmap index qknum */
  fetch_ksp_from_cache(dest,qknum,siSrc,siSnk);
#else
  int dir[3] = {0,0,0};
  if(!remap) {
  *dest = create_ksp_field(3);
  }

  int n = singlet_index_to_disp(siSnk);
  singlet_index_to_dir(siSnk,dir);
  apply_par_xport_snk_v(*dest,src[siSrc],n,dir,r0,links);
#endif
}

/**
   Destroy the ks_prop_field if mmap is not being used, otherwise just delete the pointer
   */
static void
unmap_ksp_field(ks_prop_field **ksp){
#ifdef GB_BARYON_MMAP
  toss_ksp_from_cache(ksp);
#else
  destroy_ksp_field(*ksp);
#endif
}

/*------------------------------------------------------------------*/
/**
   Antisymmetrize over color for a single site.
   This subroutine specifically handles the baryon sink antisymmetry.
   Result is saved to dt, overwriting previous
*/
static inline void
baryon_color_asym_v(su3_vector *qk0, su3_vector *qk1, su3_vector *qk2, complex *dt){
  /* TODO: make this faster? */
  {
    complex cc,dd,sum;
    CMUL(qk0->c[0],qk1->c[1],cc);
    CMUL(qk0->c[1],qk1->c[0],dd);
    CSUB(cc,dd,cc);
    CMUL(cc,qk2->c[2],sum);

    CMUL(qk0->c[0],qk1->c[2],cc);
    CMUL(qk0->c[2],qk1->c[0],dd);
    CSUB(cc,dd,cc);
    CMUL(cc,qk2->c[1],dd);
    CSUB(sum,dd,sum);

    CMUL(qk0->c[1],qk1->c[2],cc);
    CMUL(qk0->c[2],qk1->c[1],dd);
    CSUB(cc,dd,cc);
    CMUL(cc,qk2->c[0],dd);
    CADD(sum,dd,sum);

    *dt = sum;
  } // end inlined det_su3
}

/*------------------------------------------------------------------*/
/**
   Antisymmetrize over color for a single site.
   This subroutine specifically handles the baryon source antisymmetry
   and calls the baryon sink antisymmetry.
   Point splittings are handled by parent subroutines.
   Result is saved to dt, overwriting previous.
*/
static void
baryon_color_asym_mat(ks_prop_field *qk0, ks_prop_field *qk1, ks_prop_field *qk2,
                      int i, complex *dt){
  complex cc;
  baryon_color_asym_v(&qk0->v[0][i], &qk1->v[1][i], &qk2->v[2][i],dt);
  baryon_color_asym_v(&qk0->v[1][i], &qk1->v[2][i], &qk2->v[0][i],&cc);
  CSUM(*dt,cc);
  baryon_color_asym_v(&qk0->v[2][i], &qk1->v[0][i], &qk2->v[1][i],&cc);
  CSUM(*dt,cc);
  baryon_color_asym_v(&qk0->v[0][i], &qk1->v[2][i], &qk2->v[1][i],&cc);
  CSUB(*dt,cc,*dt);
  baryon_color_asym_v(&qk0->v[1][i], &qk1->v[0][i], &qk2->v[2][i],&cc);
  CSUB(*dt,cc,*dt);
  baryon_color_asym_v(&qk0->v[2][i], &qk1->v[1][i], &qk2->v[0][i],&cc);
  CSUB(*dt,cc,*dt);
}
/** same as above, but multiplies in momentum too */
static void
baryon_color_asym_mat_mom(ks_prop_field *qk0, ks_prop_field *qk1, ks_prop_field *qk2,
                      complex *mom, int i, complex *dt){
  complex cc;
  baryon_color_asym_v(&qk0->v[0][i], &qk1->v[1][i], &qk2->v[2][i],dt);
  baryon_color_asym_v(&qk0->v[1][i], &qk1->v[2][i], &qk2->v[0][i],&cc);
  CSUM(*dt,cc);
  baryon_color_asym_v(&qk0->v[2][i], &qk1->v[0][i], &qk2->v[1][i],&cc);
  CSUM(*dt,cc);
  baryon_color_asym_v(&qk0->v[0][i], &qk1->v[2][i], &qk2->v[1][i],&cc);
  CSUB(*dt,cc,*dt);
  baryon_color_asym_v(&qk0->v[1][i], &qk1->v[0][i], &qk2->v[2][i],&cc);
  CSUB(*dt,cc,*dt);
  baryon_color_asym_v(&qk0->v[2][i], &qk1->v[1][i], &qk2->v[0][i],&cc);
  CSUB(*dt,cc,*dt);
  set_complex_equal(dt,&cc);
  CMUL(cc,mom[i],*dt);
}

/** 
 * Antisymmetrizing for wall sink  
 * WARNING: This does not work with openMP enabled because
 * we need global reduction within the loop. Need a new implementation
 * for wall tieup in openACC maybe for GPU implementation??
 * 
*/
static void
baryon_color_asym_mat_wall(su3_vector *vqk0, su3_vector *vqk1, su3_vector *vqk2, complex *dt){
  complex cc;
  baryon_color_asym_v(vqk0+0, vqk1+1, vqk2+2, dt);
  baryon_color_asym_v(vqk0+1, vqk1+2, vqk2+0, &cc);
  CSUM(*dt,cc);
  baryon_color_asym_v(vqk0+2, vqk1+0, vqk2+1, &cc);
  CSUM(*dt,cc);
  baryon_color_asym_v(vqk0+0, vqk1+2, vqk2+1, &cc);
  CSUB(*dt,cc,*dt);
  baryon_color_asym_v(vqk0+1, vqk1+0, vqk2+2, &cc);
  CSUB(*dt,cc,*dt);
  baryon_color_asym_v(vqk0+2, vqk1+1, vqk2+0,&cc);
  CSUB(*dt,cc,*dt);
}

/*------------------------------------------------------------------*/
/**
   Accumulate the color antisymmetrization for a specific choice of
   point-splitting for all sites. Sink point-splitting is handled
   here, source point-splitting is implicit in the choice of
   ks_prop_fields.
   Result is summed with dt.
 */
static void
accum_baryon_color_asym(ks_prop_field *qk0, ks_prop_field *qk1, ks_prop_field *qk2,
                        short domom, complex *mom, int flip_snk, int orig, 
                        Real pfi, short dowall, complex *dt){
  //complex testsum;
  Real pf = pfi/36.;
  complex csum_new[nt] __attribute__((aligned (64))); // 64-byte
  #pragma omp simd aligned(csum_new:64)
  for(int t=0;t<nt;t++){
    (csum_new[t]).real = 0.;
    (csum_new[t]).imag = 0.;
  }

   
  if (dowall) {  // wall sink tieups    
    if (domom){
      node0_printf("accum_baryon_color_asym error: not yet implemented for wall sink with momenta\n");
      terminate(1);
    } else {
      su3_vector *vqk0, *vqk1, *vqk2;
      int nc = 3;

      //  Initilize zero vectors
      vqk0 = (su3_vector *) calloc(nt*nc,sizeof(su3_vector)); // TODO: posix_memalign(vqk0,64,nt*nc,sizeof(su3_vector)) & zero
      vqk1 = (su3_vector *) calloc(nt*nc,sizeof(su3_vector));
      vqk2 = (su3_vector *) calloc(nt*nc,sizeof(su3_vector));

      // Now do the global thing sum for wall source
      int disp_x = ((int) orig % 2     ) ^ flip_snk;
      int disp_y = ((int)(orig / 2) % 2) ^ flip_snk;
      int disp_z = ((int) orig / 4     ) ^ flip_snk;
      int i;
      site *s1;
      FORALLSITES_OMP(i,s1,){
        if (((s1->x+disp_x)%2==0) & ((s1->y+disp_y)%2==0) & ((s1->z+disp_z)%2==0)){
          for (int icolor=0; icolor<nc; icolor++){ 
            add_su3_vector(vqk0+(s1->t*nc)+icolor, &qk0->v[icolor][i], vqk0+s1->t*nc+icolor);
            add_su3_vector(vqk1+(s1->t*nc)+icolor, &qk1->v[icolor][i], vqk1+s1->t*nc+icolor);
            add_su3_vector(vqk2+(s1->t*nc)+icolor, &qk2->v[icolor][i], vqk2+s1->t*nc+icolor);
          }
        }
      } END_LOOP_OMP

      // Aggregate results on node0
      for (int t=0; t<nt; t++){
        for (int icolor=0; icolor<nc; icolor++){
          g_veccomplexsum((vqk0+t*nc+icolor)->c, nc);
          g_veccomplexsum((vqk1+t*nc+icolor)->c, nc);
          g_veccomplexsum((vqk2+t*nc+icolor)->c, nc);
        }
	complex cc;
        baryon_color_asym_mat_wall(vqk0+t*nc, 
                                   vqk1+t*nc, 
                                   vqk2+t*nc, &cc);
        CSUM(csum_new[t], cc);
        CMULREAL(csum_new[t], pf, csum_new[t]);
        CSUM(dt[t], csum_new[t]);
      }
      if (this_node != 0){
	#pragma omp simd
        for (int t=0; t<nt; t++){
          // Safety feature, because dt will be reduced again 
          // later and we really want to make sure it is zero for all other nodes
          // All the results are aggregated on dt of node0
          dt[t].real = 0.0;
          dt[t].imag = 0.0;
        }
      } // end aggregating results

      // Clean up
      free(vqk0);
      free(vqk1);
      free(vqk2);
    }  // no domom
  } else { // point sink tieups

    typedef struct {
      int t;
      int node_index;
    } node_index_t_tuple;

    // static variable initialization is done exactly once by first procedure call
    static bool initialized = false; // Don't initilize it everytime! Make it global variable
    static int corner_indx[8] = {0,0,0,0,0,0,0,0}; // Keep track of index for each corner_sites
    static node_index_t_tuple * corner_sites[8]; // use a pointer to malloc'd memory instead of a fixed size
    if (initialized==false) {
      int volume = (nx * ny * nz * nt);
      // parallelize over cube corners
      #pragma omp parallel for
      for(int temp_corner=0; temp_corner<8; temp_corner++){
          // malloc the memory (so we don't need size ahead of time)
          // Note: This is never free'd, so this is technically a memory leak
          // But if initialization is only called once, we only allocate once
          // and then we will use all this memory for the program duration
          corner_sites[temp_corner] = malloc(sizeof(node_index_t_tuple) * volume/8 + 1);
          int disp_x = ((int) temp_corner % 2     ) ^ flip_snk;
          int disp_y = ((int)(temp_corner / 2) % 2) ^ flip_snk;
          int disp_z = ((int) temp_corner / 4     ) ^ flip_snk;
          int i;
	  site *s1;
	  FORALLSITES(i,s1) {
            if (((s1->x+disp_x)%2==0) & ((s1->y+disp_y)%2==0) & ((s1->z+disp_z)%2==0)){
               corner_sites[temp_corner][corner_indx[temp_corner]] = (node_index_t_tuple){.t=s1->t, .node_index=i};
               corner_indx[temp_corner]++;
            }
          } // sites
      } // END OMP for
      initialized = true;
    } // end initilization
  
    #pragma omp parallel
      {
	complex th_csum_new[nt] __attribute__((aligned (64))); // thread private space
        #pragma omp simd aligned(th_csum_new:64)
	for(int t=0; t<nt; t++){ // zero
	  th_csum_new[t].real = 0.;
	  th_csum_new[t].imag = 0.;
	}
        #pragma omp for
	for(int i=0; i<corner_indx[orig]; i++){
	  complex cc;
	  node_index_t_tuple temp_tuple = corner_sites[orig][i];
	  if (domom){
	    baryon_color_asym_mat_mom(qk0, qk1, qk2, mom, temp_tuple.node_index, &cc);
	  } else {
	    baryon_color_asym_mat(qk0, qk1, qk2, temp_tuple.node_index, &cc);
	  }
	  CSUM(th_csum_new[temp_tuple.t],cc); // accumulate in thread private
	} // END OMP FOR
	#pragma omp critical
	{
          #pragma omp simd aligned(csum_new,th_csum_new:64)
	  for(int t=0; t<nt; t++) CSUM(csum_new[t],th_csum_new[t]); // accumulate
	}
      } // END OMP PARALLEL

    #pragma omp simd aligned(csum_new:64)
    for(int t=0;t<nt;t++){
      CMULREAL(csum_new[t],pf,csum_new[t]);
      CSUM(dt[t],csum_new[t]);
    }
  } // end point sink
}

/*------------------------------------------------------------------*/
/**
   General Golterman-Bailey baryon -type SYMMETRIC tie-up term.
   For tying up Golterman-Bailey baryon correlators.
   Source point-splitting is implicit in the prop_fields,
   choice is made in subroutine calling this one.
   Sink point-splitting explicitly determined by snkiN.
   Values stored to complex field dt which contains one term of
   the propagator.
 */
static void
gb_symm_sink_term(ks_prop_field **qk0, ks_prop_field **qk1, ks_prop_field **qk2,
                 su3_matrix *links, int tscIdx, int tskIdx, int r0[], int stIdx, short dowall, 
                 short docube, short domom, complex *mom, int flip_snk, Real pfi, complex *dt){
  int c;
  int cubestart,cubeend;
  int si_src[3] = {0};
  int si_snk[3] = {0};
  int orig,k_disp,s_disp;
  int stphs=0;
  short remap = 0x0;
  Real stsign;

  ks_prop_field *ksp0 = NULL;
  ks_prop_field *ksp1 = NULL;
  ks_prop_field *ksp2 = NULL;
  ks_prop_field *ksp3 = NULL;
  ks_prop_field *ksp4 = NULL;
  ks_prop_field *ksp5 = NULL;

  /** sink displacement - mod8 removes time direction
   number of links displaced is just the sum mod2
    of the spin and taste gamma singlet indices
   this uses bitwise xor operator ^
   */

  s_disp = (decode_gamma_spin_bits(stIdx)
         ^  decode_gamma_taste_bits(stIdx))%8;


  stphs = spin_taste_phase_factor_bits(stIdx);

  /* point-split the sinks appropriately */
  triplet_to_singlet_index(tscIdx,si_src);
  triplet_to_singlet_index(tskIdx,si_snk);
  if(flip_snk){
   for(c=0;c<3;c++){ si_snk[c] = offset_singlet_index(si_snk[c],7); }
  }

  if (docube){ cubestart = 0; cubeend = 8; }
  else {
#ifdef UNIT_TEST_ORIGIN_TRANSLATE
   cubestart = UNIT_TEST_ORIGIN_TRANSLATE;
#else
   cubestart = 0;
#endif
   cubeend = cubestart + 1;
  }
  for(orig=cubestart;orig<cubeend;orig++){

  /* keep source/sink offset fixed throughout unit cube */
  /* this is where the color epsilon is placed */
  k_disp = offset_singlet_index(s_disp,orig);
  /* spin-taste sign might need to be applied too */
  /* also correct for the basis of the sequential inversions! */
  //stsign = spin_taste_sign(stphs,offset_singlet_index(orig,si_src[0]));
  stsign = spin_taste_sign(stphs, orig);
  

  /** apply sink point splitting and calculate
     applying orig to source index automatically applies to sink too */
  map_ksp_field(&ksp0, qk0, 0, offset_singlet_index(si_src[0],orig),
    si_snk[0], r0, links, remap);


  map_ksp_field(&ksp1, qk1, 1, offset_singlet_index(si_src[1],orig),
    si_snk[1], r0, links, remap);

  map_ksp_field(&ksp2, qk1, 1, offset_singlet_index(si_src[1],orig),
    si_snk[2], r0, links, remap);
  map_ksp_field(&ksp3, qk2, 2, offset_singlet_index(si_src[2],orig),
    si_snk[0], r0, links, remap);
  map_ksp_field(&ksp4, qk2, 2, offset_singlet_index(si_src[2],orig),
    si_snk[1], r0, links, remap);
  map_ksp_field(&ksp5, qk2, 2, offset_singlet_index(si_src[2],orig),
    si_snk[2], r0, links, remap);

  remap = 0x1;

  accum_baryon_color_asym(ksp0,ksp1,ksp5,domom,mom,flip_snk,k_disp,stsign*pfi, dowall, dt);

  } // orig

  unmap_ksp_field(&ksp0);
  unmap_ksp_field(&ksp1);
  unmap_ksp_field(&ksp2);
  unmap_ksp_field(&ksp3);
  unmap_ksp_field(&ksp4);
  unmap_ksp_field(&ksp5);
}

/*------------------------------------------------------------------*/
/**
   General Golterman-Bailey baryon -type MIXED symmetry tie-up term.
   For tying up Golterman-Bailey baryon correlators.
   Source point-splitting is implicit in the prop_fields,
   choice is made in subroutine calling this one.
   Sink point-splitting explicitly determined by snkiN.
   Values stored to complex field dt which contains one term of
   the propagator.
 */
static void
gb_mixed_sink_term(ks_prop_field **qk0, ks_prop_field **qk1, ks_prop_field **qk2,
                  su3_matrix *links, int tscIdx, int tskIdx, int r0[], int stIdx, short dowall,
                  short docube, short domom, complex *mom, int flip_snk, Real pfi, complex *dt){
  int c;
  int cubestart,cubeend;
  int si_src[3] = {0};
  int si_snk[3] = {0};
  ks_prop_field *ksp0 = NULL;
  ks_prop_field *ksp1 = NULL;
  ks_prop_field *ksp2 = NULL;
  int orig,k_disp,s_disp;
  int stphs=0;
  short remap = 0x0;
  Real stsign;

  /** sink displacement - mod8 removes time direction
   number of links displaced is just the sum mod2
    of the spin and taste gamma singlet indices
   this uses bitwise xor operator ^
   */
  s_disp = (decode_gamma_spin_bits(stIdx)
         ^  decode_gamma_taste_bits(stIdx))%8;
  stphs = spin_taste_phase_factor_bits(stIdx);

  /* point-split the sinks appropriately */
  triplet_to_singlet_index(tscIdx,si_src);
  triplet_to_singlet_index(tskIdx,si_snk);
  if(flip_snk){
   for(c=0;c<3;c++){ si_snk[c] = offset_singlet_index(si_snk[c],7); }
  }


  if (docube){ cubestart = 0; cubeend = 8; }
  else {
#ifdef UNIT_TEST_ORIGIN_TRANSLATE
   cubestart = UNIT_TEST_ORIGIN_TRANSLATE;
#else
   cubestart = 0;
#endif
   cubeend = cubestart + 1;
  }
  for(orig=cubestart;orig<cubeend;orig++){

  /* keep source/sink offset fixed throughout unit cube */
  /* this is where the color epsilon is placed */
  k_disp = offset_singlet_index(s_disp,orig);
  /* spin-taste sign might need to be applied too */
  stsign = spin_taste_sign(stphs,orig);

  // three partial sums
  complex isum1[nt];
  int t;
  for(t=0;t<nt;t++){
    (isum1[t]).real = 0.0;
    (isum1[t]).imag = 0.0;
  }

  /* apply sink point splitting */
  map_ksp_field(&ksp0, qk0, 0, offset_singlet_index(si_src[0], k_disp),
    si_snk[0], r0, links, remap);
  map_ksp_field(&ksp1, qk1, 1, offset_singlet_index(si_src[1], k_disp),
    si_snk[1], r0, links, remap);
  map_ksp_field(&ksp2, qk2, 2, offset_singlet_index(si_src[2], k_disp),
    si_snk[2], r0, links, remap);
  accum_baryon_color_asym(ksp0,ksp1,ksp2,domom,mom,flip_snk,orig,stsign*3.*pfi, dowall, isum1);

  // Sum up partial sums
  for(t=0;t<nt;t++){
    CSUM(dt[t],isum1[t]);
  }
  } // orig

  unmap_ksp_field(&ksp0);
  unmap_ksp_field(&ksp1);
  unmap_ksp_field(&ksp2);
}

/*------------------------------------------------------------------*/
/**
   General Golterman-Bailey baryon -type ANTIsymmetric tie-up term.
   For tying up Golterman-Bailey baryon correlators.
   Source point-splitting is implicit in the prop_fields,
   choice is made in subroutine calling this one.
   Sink point-splitting explicitly determined by snkiN.
   Values stored to complex field dt which contains one term of
   the propagator.
 */
static void
gb_asymm_sink_term(ks_prop_field **qk0, ks_prop_field **qk1, ks_prop_field **qk2,
                  su3_matrix *links, int tscIdx, int tskIdx, int r0[], int stIdx, short dowall, 
                  short docube, short domom, complex *mom, int flip_snk, Real pfi, complex *dt){
  int c;
  int cubestart,cubeend;
  int si_src[3] = {0};
  int si_snk[3] = {0};
  ks_prop_field *ksp0 = NULL;
  ks_prop_field *ksp1 = NULL;
  ks_prop_field *ksp2 = NULL;
  ks_prop_field *ksp3 = NULL;
  ks_prop_field *ksp4 = NULL;
  ks_prop_field *ksp5 = NULL;
  int orig,k_disp,s_disp;
  int stphs=0;
  short remap = 0x0;
  Real stsign;

  /** sink displacement - mod8 removes time direction
   number of links displaced is just the sum mod2
    of the spin and taste gamma singlet indices
   this uses bitwise xor operator ^
   */
  s_disp = (decode_gamma_spin_bits(stIdx)
         ^  decode_gamma_taste_bits(stIdx))%8;
  stphs = spin_taste_phase_factor_bits(stIdx);

  /* point-split the sinks appropriately */
  triplet_to_singlet_index(tscIdx,si_src);
  triplet_to_singlet_index(tskIdx,si_snk);
  if(flip_snk){
   for(c=0;c<3;c++){ si_snk[c] = offset_singlet_index(si_snk[c],7); }
  }

  if (docube){ cubestart = 0; cubeend = 8; }
  else {
#ifdef UNIT_TEST_ORIGIN_TRANSLATE
   cubestart = UNIT_TEST_ORIGIN_TRANSLATE;
#else
   cubestart = 0;
#endif
   cubeend = cubestart + 1;
  }
  for(orig=cubestart;orig<cubeend;orig++){

  /* apply sink point splitting and calculate*/

  /* keep source/sink offset fixed throughout unit cube */
  /* this is where the color epsilon is placed */
  k_disp = offset_singlet_index(s_disp,orig);
  /* spin-taste sign might need to be applied too */
  stsign = spin_taste_sign(stphs,orig);

  map_ksp_field(&ksp0, qk0, 0, offset_singlet_index(si_src[0], orig),
    si_snk[0], r0, links, remap);
  map_ksp_field(&ksp1, qk1, 1, offset_singlet_index(si_src[1], orig),
    si_snk[1], r0, links, remap);
  map_ksp_field(&ksp2, qk1, 1, offset_singlet_index(si_src[1], orig),
    si_snk[2], r0, links, remap);
  map_ksp_field(&ksp3, qk2, 2, offset_singlet_index(si_src[2], orig),
    si_snk[0], r0, links, remap);
  map_ksp_field(&ksp4, qk2, 2, offset_singlet_index(si_src[2], orig),
    si_snk[1], r0, links, remap);
  map_ksp_field(&ksp5, qk2, 2, offset_singlet_index(si_src[2], orig),
    si_snk[2], r0, links, remap);

  remap = 0x1;
  accum_baryon_color_asym(ksp0,ksp1,ksp5,domom,mom,flip_snk,k_disp, stsign*pfi/6., dowall, dt);
  accum_baryon_color_asym(ksp0,ksp2,ksp4,domom,mom,flip_snk,k_disp,-stsign*pfi/6., dowall, dt);
  map_ksp_field(&ksp0, qk0, 0, offset_singlet_index(si_src[0], orig),
    si_snk[1], r0, links, remap);
  accum_baryon_color_asym(ksp0,ksp2,ksp3,domom,mom,flip_snk,k_disp, stsign*pfi/6., dowall, dt);
  map_ksp_field(&ksp2, qk1, 1, offset_singlet_index(si_src[1], orig),
    si_snk[0], r0, links, remap);
  accum_baryon_color_asym(ksp0,ksp2,ksp5,domom,mom,flip_snk,k_disp,-stsign*pfi/6., dowall, dt);
  map_ksp_field(&ksp0, qk0, 0, offset_singlet_index(si_src[0], orig),
    si_snk[2], r0, links, remap);
  accum_baryon_color_asym(ksp0,ksp2,ksp4,domom,mom,flip_snk,k_disp, stsign*pfi/6., dowall, dt);
  accum_baryon_color_asym(ksp0,ksp1,ksp3,domom,mom,flip_snk,k_disp,-stsign*pfi/6., dowall, dt);

  } // orig

  unmap_ksp_field(&ksp0);
  unmap_ksp_field(&ksp1);
  unmap_ksp_field(&ksp2);
  unmap_ksp_field(&ksp3);
  unmap_ksp_field(&ksp4);
  unmap_ksp_field(&ksp5);
}

/*------------------------------------------------------------------*/
static void
gb_sink_term_loop(ks_prop_field **qk0, ks_prop_field **qk1, ks_prop_field **qk2,
                 su3_matrix *links, int tscIdx, enum gb_baryon_op snk_op,
                 int num_d, int num_s, int r0[], int stIdx,
                 short dowall, short docube, short domom, complex *mom,
                 int flip_snk, Real pfi, complex *dt){
  int j,k;
  int nperm,nterm;              /* number of permutations/terms */
  int pfsum;                    /* prefactor sum */
  int pf[4];                    /* prefactors */
  int tskIdx[4];                /* triplet indices */
  int s_idx[3];                 /* singlet indices of quarks */
  int dflt[3];                  /* default taste index order, for mixed symmetry */
  int perm[3] = {0};            /* permuted order of taste indices, for mixed symmetry */
  enum gb_sym_type stype;   /* symmetry type of the operator */
  Real psign = 1.;              /* sign of Bailey operator mixed symmetry permutation */

  nperm = gb_get_num_permutations(snk_op, num_d, num_s);
  nterm = gb_get_num_terms(snk_op);
  stype = gb_get_sym_type(snk_op);
  gb_get_prefactors(snk_op,pf);
  gb_get_triplet_index(snk_op,tskIdx);
  pfsum=0;
  for(j=0;j<nterm;j++){
    pfsum+=pf[j]*pf[j];
  }

  for(j=0;j<nterm;j++){
    triplet_to_singlet_index(tskIdx[j],s_idx);
  }

  if(stype == GBSYM_SYM) {
    /* simple implementation */
    for(j=0;j<nterm;j++){
      gb_symm_sink_term(qk0, qk1, qk2, links, tscIdx, tskIdx[j],
       r0, stIdx, dowall, docube, domom, mom, flip_snk, pfi*(Real)(pf[j]/sqrt(pfsum)), dt);
    }
  } /* symmetric sink */
  else if(stype == GBSYM_ANTISYM) {
    /* simple implementation */
    for(j=0;j<nterm;j++){
      gb_asymm_sink_term(qk0, qk1, qk2, links, tscIdx, tskIdx[j],
       r0, stIdx, dowall, docube, domom, mom, flip_snk, pfi*(Real)pf[j]/sqrt(pfsum), dt);
    }
  } /* antisymmetric sink */
  else if(stype == GBSYM_MIXED){

    for(j=0;j<nterm;j++){
      if(pf[j] == 0) { continue; } // don't compute trivial terms
      triplet_to_singlet_index(tskIdx[j],dflt);
      for(k=0;k<nperm;k++){
        gb_get_permutation_order(snk_op,k,num_d,num_s,dflt,perm); // permute tastes
        psign = gb_get_permutation_sign(snk_op,k,num_d,num_s);
        gb_mixed_sink_term(qk0, qk1, qk2, links, tscIdx, singlet_to_triplet_index(perm),
          r0, stIdx, dowall, docube, domom, mom, flip_snk, pfi*(Real)pf[j]*psign/sqrt(nperm*pfsum), dt);

        /* Need to also tieup the other way around */
        int temp_perm[3];
        temp_perm[0] = perm[1];
        temp_perm[1] = perm[0];
        temp_perm[2] = perm[2];
        gb_mixed_sink_term(qk0, qk1, qk2, links, tscIdx, singlet_to_triplet_index(temp_perm),
          r0, stIdx, dowall, docube, domom, mom, flip_snk, pfi*(Real)pf[j]*psign/sqrt(nperm*pfsum), dt);

      }
    }
  }
  else{
    node0_printf("Unknown baryon operator!\n");
    node0_printf("Failure for sink operator %s!\n",gb_baryon_label(snk_op));
    terminate(1);
  }
}

/*------------------------------------------------------------------*/
/**
   Simple symmetrization over a single quark triplet in source term.
 */
static void
gb_symm_source_term_loop(ks_prop_field **qko0, ks_prop_field **qko1, ks_prop_field **qko2,
                   su3_matrix *links, int t_idx, enum gb_baryon_op snk_op,
                   int num_d, int num_s, int r0[],
                   int stIdx, short dowall, short docube, short domom, complex *mom,
                   int flip_snk, Real pfi, complex *dt){
  int c;
  int s_idx[3] = {0};
  int s_perm[3] = {0};
  triplet_to_singlet_index(t_idx,s_idx);

  for(c=0;c<3;c++){
    if(qko0[s_idx[c]] == NULL ||
       qko1[s_idx[c]] == NULL ||
       qko2[s_idx[c]] == NULL){
       node0_printf("gb_symm_source_term_loop called with non-existant quark!\n");
       node0_printf("triplet: %i %i %i\n",s_idx[0],s_idx[1],s_idx[2]);
       terminate(1);
    }
  }
 /* all unique quarks and tastes */

	gb_sink_term_loop(qko0, qko1, qko2, links, t_idx, snk_op, num_d, num_s,
	  r0, stIdx, dowall, docube, domom, mom, flip_snk, pfi/6., dt);
	s_perm[0] = s_idx[1]; s_perm[1] = s_idx[2]; s_perm[2] = s_idx[0];
	gb_sink_term_loop(qko0, qko1, qko2, links, singlet_to_triplet_index(s_perm),
	  snk_op, num_d, num_s, r0, stIdx, dowall, docube, domom, mom, flip_snk, pfi/6., dt);
	s_perm[0] = s_idx[2]; s_perm[1] = s_idx[0]; s_perm[2] = s_idx[1];
	gb_sink_term_loop(qko0, qko1, qko2, links, singlet_to_triplet_index(s_perm),
	  snk_op, num_d, num_s, r0, stIdx, dowall, docube, domom, mom, flip_snk, pfi/6., dt);
	s_perm[0] = s_idx[0]; s_perm[1] = s_idx[2]; s_perm[2] = s_idx[1];
	gb_sink_term_loop(qko0, qko1, qko2, links, singlet_to_triplet_index(s_perm),
	  snk_op, num_d, num_s, r0, stIdx, dowall, docube, domom, mom, flip_snk, pfi/6., dt);
	s_perm[0] = s_idx[1]; s_perm[1] = s_idx[0]; s_perm[2] = s_idx[2];
	gb_sink_term_loop(qko0, qko1, qko2, links, singlet_to_triplet_index(s_perm),
	  snk_op, num_d, num_s, r0, stIdx, dowall, docube, domom, mom, flip_snk, pfi/6., dt);
	s_perm[0] = s_idx[2]; s_perm[1] = s_idx[1]; s_perm[2] = s_idx[0];
	gb_sink_term_loop(qko0, qko1, qko2, links, singlet_to_triplet_index(s_perm),
	  snk_op, num_d, num_s, r0, stIdx, dowall, docube, domom, mom, flip_snk, pfi/6., dt);

}

/*------------------------------------------------------------------*/
/**
   Simple symmetrization over a single quark triplet in source term.
 */
static void
gb_mixed_source_term_loop(ks_prop_field **qko0, ks_prop_field **qko1, ks_prop_field **qko2,
                          su3_matrix *links, int t_idx, enum gb_baryon_op snk_op,
                          int num_d, int num_s, char *qkLoc, int r0[],
                          int stIdx, short dowall, short docube, short domom, complex *mom,
                          int flip_snk, Real pfi, complex *dt){
  int c;
  int s_idx[3] = {0};
  int s_perm[3] = {0};
  triplet_to_singlet_index(t_idx, s_idx);

  for(c=0;c<3;c++){
    if(qko0[s_idx[c]] == NULL ||
       qko1[s_idx[c]] == NULL ||
       qko2[s_idx[c]] == NULL){
       node0_printf("gb_mixed_source_term_loop called with non-existant quark!\n");
       node0_printf("triplet: %i %i %i\n",s_idx[0],s_idx[1],s_idx[2]);
       terminate(1);
    }
  }

  /* assume all combos are unique */
  if (strcmp(qkLoc, "d'ud") == 0) {
	  s_perm[0] = s_idx[0];
	  s_perm[1] = s_idx[2];
	  s_perm[2] = s_idx[1];
  } else if (strcmp(qkLoc, "ud'd") == 0) {
	  s_perm[0] = s_idx[1];
	  s_perm[1] = s_idx[2];
	  s_perm[2] = s_idx[0];
  } else if (strcmp(qkLoc, "uud") == 0) {
	  s_perm[0] = s_idx[0];
	  s_perm[1] = s_idx[1];
	  s_perm[2] = s_idx[2];
  } else {
	  node0_printf("Unknown quark composition!\n");
  }
  gb_sink_term_loop(qko0, qko1, qko2, links, singlet_to_triplet_index(s_perm),
	  snk_op, num_d, num_s, r0, stIdx, dowall, docube, domom, mom, flip_snk, 2.*pfi, dt);

  int tempperm[3] = {0};
  tempperm[0] = s_idx[1]; tempperm[1] = s_idx[2]; tempperm[2] = s_idx[0];
  /* assume all combos are unique */
  if (strcmp(qkLoc, "d'ud") == 0) {
	  s_perm[0] = tempperm[0];
	  s_perm[1] = tempperm[2];
	  s_perm[2] = tempperm[1];
  } else if (strcmp(qkLoc, "ud'd") == 0) {
	  s_perm[0] = tempperm[1];
	  s_perm[1] = tempperm[2];
	  s_perm[2] = tempperm[0];
  } else if (strcmp(qkLoc, "uud") == 0) {
	  s_perm[0] = tempperm[0];
	  s_perm[1] = tempperm[1];
	  s_perm[2] = tempperm[2];
  } else {
	  node0_printf("Unknown quark composition!\n");
  }
  gb_sink_term_loop(qko0, qko1, qko2, links, singlet_to_triplet_index(s_perm),
      snk_op, num_d, num_s, r0, stIdx, dowall, docube, domom, mom, flip_snk, -pfi, dt);

  tempperm[0] = s_idx[2]; tempperm[1] = s_idx[0]; tempperm[2] = s_idx[1];
  /* assume all combos are unique */
  if (strcmp(qkLoc, "d'ud") == 0) {
	  s_perm[0] = tempperm[0];
	  s_perm[1] = tempperm[2];
	  s_perm[2] = tempperm[1];
  } else if (strcmp(qkLoc, "ud'd") == 0) {
	  s_perm[0] = tempperm[1];
	  s_perm[1] = tempperm[2];
	  s_perm[2] = tempperm[0];
  } else if (strcmp(qkLoc, "uud") == 0) {
	  s_perm[0] = tempperm[0];
	  s_perm[1] = tempperm[1];
	  s_perm[2] = tempperm[2];
  } else {
	  node0_printf("Unknown quark composition!\n");
  }
  gb_sink_term_loop(qko0, qko1, qko2, links, singlet_to_triplet_index(s_perm),
      snk_op, num_d, num_s, r0, stIdx, dowall, docube, domom, mom, flip_snk, -pfi, dt);

}

/*------------------------------------------------------------------*/
/**
   Simple antisymmetrization over a single quark triplet in source term.
 */
static void
gb_asymm_source_term_loop(ks_prop_field **qko0, ks_prop_field **qko1, ks_prop_field **qko2,
                   su3_matrix *links, int t_idx, enum gb_baryon_op snk_op,
                   int num_d, int num_s, int r0[],
                   int stIdx, short dowall, short docube, short domom, complex *mom,
                   int flip_snk, Real pfi, complex *dt){
  int c;
  int s_idx[3] = {0};
  int s_perm[3] = {0};
  triplet_to_singlet_index(t_idx,s_idx);

  for(c=0;c<3;c++){
    if(qko0[s_idx[c]] == NULL ||
       qko1[s_idx[c]] == NULL ||
       qko2[s_idx[c]] == NULL){
       node0_printf("gb_asymm_source_term_loop called with non-existant quark!\n");
       node0_printf("triplet: %i %i %i\n",s_idx[0],s_idx[1],s_idx[2]);
       terminate(1);
    }
  }

  /* all must be unique */
  gb_sink_term_loop(qko0, qko1, qko2, links, t_idx, snk_op, num_d, num_s,
      r0, stIdx, dowall, docube, domom, mom, flip_snk,  pfi/6., dt);
  s_perm[0] = s_idx[1]; s_perm[1] = s_idx[2]; s_perm[2] = s_idx[0];
  gb_sink_term_loop(qko0, qko1, qko2, links, singlet_to_triplet_index(s_perm),
      snk_op, num_d, num_s, r0, stIdx, dowall, docube, domom, mom, flip_snk,  pfi/6., dt);
  s_perm[0] = s_idx[2]; s_perm[1] = s_idx[0]; s_perm[2] = s_idx[1];
  gb_sink_term_loop(qko0, qko1, qko2, links, singlet_to_triplet_index(s_perm),
      snk_op, num_d, num_s, r0, stIdx, dowall, docube, domom, mom, flip_snk,  pfi/6., dt);
  s_perm[0] = s_idx[0]; s_perm[1] = s_idx[2]; s_perm[2] = s_idx[1];
  gb_sink_term_loop(qko0, qko1, qko2, links, singlet_to_triplet_index(s_perm),
      snk_op, num_d, num_s, r0, stIdx, dowall, docube, domom, mom, flip_snk, -pfi/6., dt);
  s_perm[0] = s_idx[1]; s_perm[1] = s_idx[0]; s_perm[2] = s_idx[2];
  gb_sink_term_loop(qko0, qko1, qko2, links, singlet_to_triplet_index(s_perm),
      snk_op, num_d, num_s, r0, stIdx, dowall, docube, domom, mom, flip_snk, -pfi/6., dt);
  s_perm[0] = s_idx[2]; s_perm[1] = s_idx[1]; s_perm[2] = s_idx[0];
  gb_sink_term_loop(qko0, qko1, qko2, links, singlet_to_triplet_index(s_perm),
      snk_op, num_d, num_s, r0, stIdx, dowall, docube, domom, mom, flip_snk, -pfi/6., dt);
}

/*------------------------------------------------------------------*/
static void
gb_source_term_loop(ks_prop_field **qko0, ks_prop_field **qko1, ks_prop_field **qko2,
                   su3_matrix *links, enum gb_baryon_op src_op,
                   enum gb_baryon_op snk_op, int stIdx, short dowall, short docube,
                   int num_d, int num_s, int r0[], short domom,
                   complex *mom, int flip_snk, complex *dt){
  int j,k;
  int nperm,nterm;              /* number of permutations/terms */
  int pfsum;                    /* prefactor sum */
  int pf[4];                    /* prefactors */
  int t_idx[4];                 /* triplet indices */
  int dflt[3];                  /* default taste index order, for mixed symmetry */
  int perm[3] = {0};            /* permuted order of taste indices, for mixed symmetry */
  enum gb_sym_type stype;       /* symmetry type of the operator */
  Real psign = 1.;              /* sign of Bailey operator mixed symmetry permutation */

  nperm = gb_get_num_permutations(src_op, num_d, num_s);
  nterm = gb_get_num_terms(src_op);
  stype = gb_get_sym_type(src_op);
  gb_get_prefactors(src_op,pf);
  gb_get_triplet_index(src_op,t_idx);
  pfsum=0;
  for(j=0;j<nterm;j++){
    pfsum+=pf[j]*pf[j];
  }

  if(stype == GBSYM_SYM || stype == GBSYM_ANTISYM) {
    /* simple implementation */
    for(j=0;j<nterm;j++){
      /* choose source point splitting here */
      //triplet_to_singlet_index(t_idx[j],s_idx);
      if(stype == GBSYM_SYM){
		  gb_symm_source_term_loop(qko0, qko1, qko2,
			  links, t_idx[j], snk_op, num_d, num_s, r0, stIdx, dowall, docube, domom, mom,
			  flip_snk, (Real)pf[j]/sqrt(pfsum), dt);
      }

      else if(stype == GBSYM_ANTISYM)
        gb_asymm_source_term_loop(qko0, qko1, qko2,
          links, t_idx[j], snk_op, num_d, num_s, r0, stIdx, dowall, docube, domom, mom,
          flip_snk, (Real)pf[j]/sqrt(pfsum), dt);
    } /* terms in sink */
  } /* (anti)symmetric sink */
  else if(stype == GBSYM_MIXED){
    /* swap taste indices to get other terms in mixed symmetry operators */
    /* assumes strange quarks are always last */
    for(j=0;j<nterm;j++){
      if(pf[j] == 0) { continue; } // don't compute trivial terms
      triplet_to_singlet_index(t_idx[j],dflt);
      for(k=0;k<nperm;k++){
        // permute tastes rather than flavors
        gb_get_permutation_order(src_op,k,num_d,num_s,dflt,perm);
        psign = gb_get_permutation_sign(src_op,k,num_d,num_s);

        if (stIdx != GB_2POINT_BACKPROP) {
        	/* For 3pt only, the goal is to make uud with ud insertion into ddu
        	 * You need to permute the last down quark around a bit */
        	char qkLoc[] = "d'ud";
            gb_mixed_source_term_loop(qko0, qko1, qko2, links, singlet_to_triplet_index(perm),
              snk_op, num_d, num_s, qkLoc, r0, stIdx, dowall, docube, domom, mom,
              flip_snk, (Real)pf[j]*psign/sqrt(nperm*pfsum), dt);

            strcpy(qkLoc, "ud'd");
			gb_mixed_source_term_loop(qko0, qko1, qko2, links, singlet_to_triplet_index(perm),
			   snk_op, num_d, num_s, qkLoc, r0, stIdx, dowall, docube, domom, mom,
			   flip_snk, (Real)pf[j]*psign/sqrt(nperm*pfsum), dt);

        } else {
        	/* For 2pt only */
        	char qkLoc[] = "uud";

            gb_mixed_source_term_loop(qko0, qko1, qko2, links, singlet_to_triplet_index(perm),
              snk_op, num_d, num_s, qkLoc, r0, stIdx, dowall, docube, domom, mom,
              flip_snk, (Real)pf[j]*psign/sqrt(nperm*pfsum), dt);
        }

      }
    }
  }
  else{
    node0_printf("Unknown baryon operator!\n");
    node0_printf("Failure for source operator %s!\n",gb_baryon_label(src_op));
    terminate(1);
  }
}

/*------------------------------------------------------------------*/
void
gb_baryon(ks_prop_field *qko0[], ks_prop_field *qko1[], ks_prop_field *qko2[],
              su3_matrix *links, enum gb_baryon_op src_op[],
              enum gb_baryon_op snk_op[],
              int stIdx, short dowall[], short docube[], int num_d, int num_s, int r0[],
              int mom[], char par[], complex *momfld, int flip_snk[],
              int num_corr_gb, int phase[], Real fact[], complex *prop[]){
  int i,px,py,pz;
  char ex,ey,ez;
  double factx = 2.0*PI/(1.0*nx);
  double facty = 2.0*PI/(1.0*ny);
  double factz = 2.0*PI/(1.0*nz);
  site *s;
  short domom = 0x0;
  if (mom[0] != 0 || mom[1] != 0 || mom[2] != 0){
    // if nonzero, compute Fourier phase once and reuse many times
    // negative sign is necessary to get momentum conservation
    px = mom[0]; py = mom[1]; pz = mom[2];
    ex = par[0]; ey = par[1]; ez = par[2];
    FORALLSITES_OMP(i,s,) {
      complex tmp = {1.,0.};
      tmp = ff(factx*(s->x-r0[0])*px, ex, tmp);
      tmp = ff(facty*(s->y-r0[1])*py, ey, tmp);
      tmp = ff(factz*(s->z-r0[2])*pz, ez, tmp);
      set_complex_equal(&tmp,&(momfld[i]));
    } END_LOOP_OMP
   domom = 0x1;
  } else { // just specify that momentum is not needed
   domom = 0x0;
  }

  for(i=0;i<num_corr_gb;i++){
    node0_printf("gb baryon: %s %s\n", gb_baryon_label(src_op[i]),gb_baryon_label(snk_op[i]));
    double dtime = start_timing();
    gb_source_term_loop(qko0,qko1,qko2,links,src_op[i],snk_op[i],stIdx,dowall[i],docube[i],
     num_d,num_s,r0,domom,momfld,flip_snk[i],prop[i]);
    print_timing(dtime, "computing correlator");
    norm_corr( phase[i], fact[i], prop[i] );
  }
}

#endif //GB_BARYON
