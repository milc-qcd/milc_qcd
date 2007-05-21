/********************** calc_matrix.c *********************************/
/* MIMD version 7 */
/*
 *  This function calculates the static-light smearing matrix.
 *  WARNING::: this routine overwrites the gauge configuration to save
 *  workspace.
 *
 *  The variational matrix is constructed from the light quark 
 *  propagator and  the Wilson line.
 *
 *  T( \pm t )_AB =  1  \sum_{r'} f(r')_A  \sum_{p} 
 *              2 V
 *
 *                 trace S(-p,\pm t) (1 \pm \gamma_4 ) W(p,\pm t)  e^{ i p .r'}
 *
 *  Nov 13 1996 : To speed up this code, I am trying to FFT larger
 *                data arrays.
 */


#include "w_static_includes.h"

/* 
 *  Function to construct the variational matrix of static
 *  correlators.
 */

void calc_vary_matrix()
{
  int dim = nt*nosmear*nosmear  ;
  int bloop ;
  register site *s;

  Real *vary_matrix ; /* Variational matrix *******/
  int i ;
  Real ts,te ;
  const int block_size = 4 ; 
  int block_pt = 0 ;
  int block_smear[4] ; 

  /****..................................................**/

  ts = dclock();

  /* Reserve space for the variational matrix on each node ****/
  if( (vary_matrix = (Real *) calloc( (size_t) nt*nosmear*nosmear, sizeof(Real) )  ) == NULL )
  {
    printf("ERROR: could not reserve buffer space for the variational matrix");
    terminate(1);
  }

  /***** zero the variational matrix  *****/
  for(i=0 ; i < dim ;++i)
    *(vary_matrix + i ) = 0.0 ;

  /* FFT the "spinless" light quark propagator ***/
  fft_negmom_quark();

  /* ... calculate the variational matrix ****/

  for(block_pt = 0 ; block_pt < block_size ; ++block_pt)
    block_smear[ block_pt ] =  0 ;

  block_pt = 0 ; 
  for(bloop =0 ; bloop < nosmear ; ++bloop)
  {
    block_smear[ block_pt ] =  bloop ;
  
    /* Smear the Wilson line [multiply by the smearing function]  ***/
    FORALLSITES(i,s)
    {
      c_scalar_mult_su3mat(  
			   &(s->w_line),
			   &(s->smear_func[bloop]),
			   &(s->smear_w_line[ block_pt ]) )  ;
    }


    if( (block_pt + 1 == block_size ) || ( bloop + 1 == nosmear) )
    {
      /* Convolve with the A sources  ***/
      ++block_pt ;

      contract_a_src(vary_matrix,block_pt , block_smear);
      for(block_pt = 0 ; block_pt < block_size ; ++block_pt)
	block_smear[ block_pt ] =  0 ;

      block_pt = 0 ; 
    }
    else
    {
      ++block_pt ; 
    }


  }  /* end the loop over the B-source ***/


  /** sum the variational matrix over the nodes of the lattice ****/
  for(i=0 ; i < dim ;++i)
    g_floatsum(vary_matrix + i ) ;


  /* Save the variational correlation matrix to disk **/
  IF_MASTER
    write_vary_matrix(vary_matrix) ;

  /* Free up the memory used in the calculation ***/
  free(vary_matrix);


 te = dclock() - ts ;
 IF_MASTER
   printf("Time to calculate the variational smearing matrix = %g sec\n",te);

  fflush(stdout);


}



/*
 *   Contract the A - sources into the convolution of the light quark
 *   and static propagator.
 *
 */


void contract_a_src(Real *vary_matrix, int block_end, int block_smear[4] )
{
  int aloop ;
  int t,pt,i ;
  complex z ;
  register site *s;
  int bloop ;
  int block_pt ;
  su3_matrix prod;

  /*
     The routines that trace (1 + gamma_{4}) into the light
     quark propagator do not multiply by the factor of 1/2;
     this factor of 1/2 is included in this routine.
  */
  Real norm = 0.5/((Real) nx*ny*nz) ;


  /* FFT the smeared Wilson line, stored in a block structure ***/
  restrict_fourier_site(F_OFFSET(smear_w_line), 4*sizeof(su3_matrix), 1);
  
  for( block_pt = 0 ; block_pt < block_end ; ++block_pt)
  {
    /* Convolve the smeared Wilson line with the light quark propagator **/
    FORALLSITES(i,s)
    {
      mult_su3_nn(&(s->smear_w_line[block_pt]),&(s->strip_quark),&prod )  ;
      s->chi.d[block_pt].c[0]  = trace_su3(&prod);
    }
    
  }  /*** end of the loop over the block ***/
	
  /* FFT the convolution back to real space ****/
  restrict_fourier_site(F_OFFSET(chi), sizeof(wilson_vector), -1);

  for( block_pt = 0 ; block_pt < block_end ; ++block_pt)
  {
    bloop =  block_smear[block_pt]  ; 
    for(aloop=0 ; aloop < nosmear ; ++aloop )
    {
      FORALLSITES(i,s)
      {
	t = s->t ;
	pt = VM_PT ;
      
	z = cmul(&( s->chi.d[block_pt].c[0] ), &(s->smear_func[aloop])   ) ;
	*(vary_matrix +pt) += z.real*norm ;
	
      }
    } /* end the loop over the A-sources ***/
  }  /**** end of the loop over the block size ***/



}



