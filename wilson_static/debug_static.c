/*  $Header: /lqcdproj/detar/cvsroot/milc_qcd/wilson_static/debug_static.c,v 1.3 2007/10/09 21:02:32 detar Exp $
 *  This file contains a number of functions, thought
 *  to be useful to help debug the static code.
 *
 */


#include "w_static_includes.h"


/*
 *  Set the smearing function equal to a 
 *  set of random numbers.
 *
 *
 *
 */

void random_smear_func()
{
  int i ;
  register site *s;
  int spt ;
  long seed = -1020 ;
  Real ran1(long *idum);
  /***********-----------*************************/



  FORALLSITES(i,s)
  {
    for(spt = 0 ; spt < nosmear ; ++spt)
    {
      s->smear_func[spt].real = ran1(&seed) ;
      s->smear_func[spt].imag = ran1(&seed) ;
    }

  }


}




/*
 *  Dump to the screen the gauge configurations
 *  This version only dumps the configuration in the 
 *  time direction.
 */

void dump_time_gauge()
{
  int i ;
  register site *s;
  int ic ;
  int dir  ; 


  for( dir = 0 ; dir < 4 ; ++dir)
  {
    printf("******************************* \n");
    printf("**** direction = %d  *********\n",dir);
    printf("******************************* \n");

    FORALLSITES(i,s)
    {
      printf("x= [ %d %d %d ] t= %d \n",s->x,s->y,s->z ,s->t   );
      for(ic=0 ;ic < 3 ; ++ic )
      {
	printf(" ( %g %g    %g %g    %g  %g ) \n",
	     s->link[dir].e[ic][ 0 ].real , s->link[dir].e[ic][ 0 ].imag ,
	     s->link[dir].e[ic][ 1 ].real , s->link[dir].e[ic][ 1 ].imag ,
	     s->link[dir].e[ic][ 2 ].real , s->link[dir].e[ic][ 2 ].imag );
	     
      }
    
    }  /*** end the loop over sites on the lattice ***/

  }  /*** end the loop over the directions ****/



}









/*
 *  Dump to the screen the gauge configurations
 *  This version only dumps the configuration in the 
 *  time direction.
 */

void dump_gauge()
{
  int i ;
  register site *s;
  Real g_re, g_im ;
  int ic, jc ;


  printf("Here is the gauge configuration in the time direction\n");

  FORALLSITES(i,s)
  {
    printf("x= [ %d %d %d ] t= %d \n",s->x,s->y,s->z ,s->t   );
    for(ic=0 ;ic < 3 ; ++ic )
      for(jc=0 ;jc < 3 ; ++jc )
      {
	g_re = s->link[3].e[ic][jc].real ;
	g_im = s->link[3].e[ic][jc].imag ;
	printf("%g  %g\n",g_re,g_im);
      }

  }


}


/*
 *  Create a special gauge configuration in the time
 *  direction -- useful for debugging the code that calculates the 
 *  Wilson line.
 */

void test_gauge_config()
{
  int i ;
  register site *s;
  int ic, jc ;


  printf("I am creating a test gauge configuration\n");


  FORALLSITES(i,s)
  {
    for(ic=0 ;ic < 3 ; ++ic )
    {
      for(jc=ic+1 ;jc < 3 ; ++jc )
      {
	s->link[3].e[ic][jc].real = 0.0 ;
	s->link[3].e[ic][jc].imag = 0.0 ;

	s->link[3].e[jc][ic].real = 0.0 ;
	s->link[3].e[jc][ic].imag = 0.0 ;

      }
      s->link[3].e[ic][ic].real = 0.5 + s->t  ;
      s->link[3].e[ic][ic].imag = 0.0 ;

    }  


  }  /** end of the loop over the lattice **/


}





/*
 *  Dump to the screen the spinless light quark
 *  propagator.
 */

void dump_strip_quark()
{
  int i ;
  register site *s;
  Real g_re, g_im ;
  int ic, jc ;


  printf("Here is the spinless quark propagator\n");

  FORALLSITES(i,s)
  {
    printf("x= [ %d %d %d ] t= %d \n",s->x,s->y,s->z ,s->t   );
    for(ic=0 ;ic < 3 ; ++ic )
      for(jc=0 ;jc < 3 ; ++jc )
      {
	g_re = s->strip_quark.e[ic][jc].real ;
	g_im = s->strip_quark.e[ic][jc].imag ;
	printf("%g  %g\n",g_re,g_im);
      }

  }


}






/*
 *  Dump the Wilson vector loaded in from disk
 *  propagator.
 */

void dump_psi()
{
  int i ;
  register site *s;
  Real g_re, g_im ;
  int ic=1, ispin=2 ;
  int where ;
  /******----------**************************************/


  node0_printf("QUARK PROPAGATOR for colour = %d and spin = %d\n",ic,ispin);

  FORALLSITES(i,s)
  {
    where = s->x + nx*( s->y + ny*( s->z + nz*s->t)) ;

    g_re = s->psi.d[ispin].c[ic].real ;
    g_im = s->psi.d[ispin].c[ic].imag ;
    printf("%d  %g  %g\n",where,g_re,g_im);

  }


/** wait for all the nodes to dump their data ***/
      g_sync();


/******************************
  printf("Here is the quark propagator\n");
  FORALLSITES(i,s)
  {
    printf("x= [ %d %d %d ] t= %d \n",s->x,s->y,s->z ,s->t   );
    for(ic=0 ;ic < 3 ; ++ic )
      for(ispin=0 ;ispin < 4 ; ++ispin )
      {
	g_re = s->psi.d[ispin].c[ic].real ;
	g_im = s->psi.d[ispin].c[ic].imag ;
	printf("%g  %g\n",g_re,g_im);
      }

  }
****************************************/



}


/*
 *  Dump the Wilson vector loaded in fomr disk
 *  propagator.
 */

void dump_psi_smear()
{
  int i ;
  register site *s;
  Real g_re, g_im ;
  int ic, ispin ;
  /***********************************************************/

  printf("Here is the quark propagator [original and smeared]\n");

  FORALLSITES(i,s)
  {
    printf("x= [ %d %d %d ] t= %d \n",s->x,s->y,s->z ,s->t   );
    for(ic=0 ;ic < 3 ; ++ic )
      for(ispin=0 ;ispin < 4 ; ++ispin )
      {
	g_re = s->psi.d[ispin].c[ic].real ;
	g_im = s->psi.d[ispin].c[ic].imag ;

/**	g_s_re = s->psi_smear.d[ispin].c[ic].real ;
	g_s_im = s->psi_smear.d[ispin].c[ic].imag ;
****/

	/**	printf("%g  %g  smear = %g  %g  \n",g_re,g_im, g_s_re,g_s_im );**/
	printf("%g  %g   \n",g_re,g_im );
      }

  }


}


/*  Numerical recipies code for random numbers
 *  FOR debugging purposes only !!!!!!!!!!!!
 */


#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS_NR 1.2e-7
#define RNMX (1.0-EPS_NR)

Real ran1(long *idum)
{
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	Real temp;

	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum < 0) *idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = *idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS_NR
#undef RNMX
/* (C) Copr. 1986-92 Numerical Recipes Software 5.){2p491&+K45#+9. */


/*
 *  Write the variational matrix to the screen
 *
 */

void dump_vary_matrix(Real *vary_matrix)
{
  int i ;
  int dim = nt*nosmear*nosmear ;


  printf("Terminal copy of the variotional matrix \n");
  for(i=0 ; i < dim ;++i)
    printf("vary_matrix[ %d ] = %g \n",i,*(vary_matrix + i));

}




/*
 *  Dump to the screen the static propagator
 *  The gauge part anyway
 */

void dump_static_prop()
{
  int i ;
  register site *s;
  Real g_re, g_im ;
  int ic, jc ;


  printf("Here is the static quark propagator\n");

  FORALLSITES(i,s)
  {
    printf("x= [ %d %d %d ] t= %d \n",s->x,s->y,s->z ,s->t   );
    for(ic=0 ;ic < 3 ; ++ic )
      for(jc=0 ;jc < 3 ; ++jc )
      {
	g_re = s->w_line.e[ic][jc].real ;
	g_im = s->w_line.e[ic][jc].imag ;
	printf("%g  %g\n",g_re,g_im);
      }

  }


}




/*
 *  Dump to the screen the smeared Wilson line
 *
 */

void dump_smear_line()
{
  int i ;
  register site *s;
  Real g_re, g_im ;
  int ic, jc ;


  printf("Here is the smeared Wilson line \n");

  FORALLSITES(i,s)
  {
    printf("x= [ %d %d %d ] t= %d \n",s->x,s->y,s->z ,s->t   );
    for(ic=0 ;ic < 3 ; ++ic )
      for(jc=0 ;jc < 3 ; ++jc )
      {
	g_re = s->smear_w_line[0].e[ic][jc].real ;
	g_im = s->smear_w_line[0].e[ic][jc].imag ;
	printf("%g  %g\n",g_re,g_im);
      }

  }


}



/*
 *  Dump the smearing function to the screen
 *  I have modified the code so that the FFT of the smearing
 *  functions is also dumped out to disk.
 *
 *
 */

void dump_smear_func()
{
  int i ;
  register site *s;
  int spt ;
  Real s_re,s_im ;
  /***********-----------*************************/


  for(spt = 0 ; spt < nosmear ; ++spt)
  {
    printf("Here is the smearing function number %d \n",spt);
    FORALLSITES(i,s)
    {
      s_re = s->smear_func[spt].real  ;
      s_im = s->smear_func[spt].imag  ;

/***      s_fft_re = s->smear_func_fft[spt].real  ;
      s_fft_im = s->smear_func_fft[spt].imag  ;  ***/

      printf("xyzt = %d %d %d %d ",s->x, s->y, s->z, s->t );
      /**      printf("Smear = %g %g fft = %g %g \n",s_re,s_im, s_fft_re, s_fft_im  );**/
      printf("Smear = %g %g \n",s_re,s_im );
    }

  }


}







/*
 *  
 *  Write the smearing functions read in
 *  to the screen
 *
 *
 */

void write_smear_func()
{
  int i ;
  register site *s;
  int spt ;
  double s_re, s_im ;
  /***********-----------*************************/



  for(spt = 0 ; spt < nosmear ; ++spt)
  {
    printf(">> Smearing function number %d\n",spt);
    FORALLSITES(i,s)
    {
      s_re = s->smear_func[spt].real  ;
      s_im = s->smear_func[spt].imag  ;

      printf("x = [%d %d %d %d]  smear = %f %f\n",s->x,s->y,s->z,s->t, s_re, s_im  );
    }



  }




}





/*  Test the gamma matrix routine which multiplies 
 *   a Wilson vector with \gamma_4
 *
 *
 */



void test_gamma_left()
{
  wilson_matrix rrr, g4_rrr ;
  long seed = -1020 ;
  Real ran1(long *idum);
  Real diff_re , diff_im ;

  int ispin , jspin ;
  int ic, jc ;
  int ispin_perm ;
  int perm[4] ;
  /***----------..........----------..........  ****/

  perm[0] = 2 ;   perm[1] = 3 ;   perm[2] = 0 ;   perm[3] = 1 ; 

  printf("Test some of the gamma matrix routines \n");
  printf("----------------------------------------\n");

  /** create the random initial matrix ***/

  for(ic=0 ;ic < 3 ; ++ic )
    for(jc=0 ;jc < 3 ; ++jc )
      for(ispin = 0 ; ispin < 4 ; ++ispin)
	for(jspin = 0 ; jspin < 4 ; ++jspin)
	{
	  rrr.d[ispin].c[ic].d[jspin].c[jc].real = ran1(&seed) ;
	  rrr.d[ispin].c[ic].d[jspin].c[jc].imag = ran1(&seed) ;
	}

  /*** multiply with the gamma_5 matrix ****/
  mult_by_gamma_left(&rrr , &g4_rrr , TUP) ;

  /*** check the multiplication   ********/

  for(ic=0 ;ic < 3 ; ++ic )
    for(jc=0 ;jc < 3 ; ++jc )
    {
      printf("colour = %d %d \n",ic,jc);

      for(ispin = 0 ; ispin < 4 ; ++ispin)
	for(jspin = 0 ; jspin < 4 ; ++jspin)
	{
	  ispin_perm = perm[ ispin ] ;

	  diff_re = rrr.d[ispin_perm].c[ic].d[jspin].c[jc].real - 
	    g4_rrr.d[ispin].c[ic].d[jspin].c[jc].real ;

	  diff_im = rrr.d[ispin_perm].c[ic].d[jspin].c[jc].imag - 
	    g4_rrr.d[ispin].c[ic].d[jspin].c[jc].imag ;

	  printf("%g  %g  --- %g %g = %g %g \n",
		 rrr.d[ispin_perm].c[ic].d[jspin].c[jc].real,
		 rrr.d[ispin_perm].c[ic].d[jspin].c[jc].imag ,
		 g4_rrr.d[ispin].c[ic].d[jspin].c[jc].real,
		 g4_rrr.d[ispin].c[ic].d[jspin].c[jc].imag,
		 diff_re, diff_im);

	}
    } /** end the loop over colour ***/


}





/*
 *   Copy a wilson matrix into a wilson vector
 *
 *
 */

void copy_w_right_vec_mat(wilson_vector *v, wilson_matrix *m, int spin, int color)
{
  int ispin, ic ;


  for(ispin= 0 ; ispin < 4 ;++ispin )
    for(ic = 0 ; ic < 3 ; ++ic )
      v->d[ispin].c[ic] = m->d[ispin].c[ic].d[spin].c[color]  ;



}


/*
 *   Copy a wilson vector into a wilson matrix 
 *
 *
 */

void copy_w_right_mat_vec(wilson_matrix *m, wilson_vector *v, int spin, int color)
{
  int ispin, ic ;


  for(ispin= 0 ; ispin < 4 ;++ispin )
    for(ic = 0 ; ic < 3 ; ++ic )
      m->d[ispin].c[ic].d[spin].c[color]  =  v->d[ispin].c[ic] ;



}



/*
 *  Test the fourier transform routines
 *  Write out a one dimensional complex array
 *
 *
 */

void dump_complex_array(complex *data, int nodata)
{
  int i ;

  printf("Here is the smeared meson array");

  for(i=0 ; i < nodata ; ++i)
    printf("%d   (%g , %g) \n",i, data[i].real, data[i].imag  );

}






/*
 *  Dump out the smeared meson array
 *
 *
 */

void dump_smearedmeson(complex *data)
{
  int ismear,t ;
  int ic,colour ;
  int ispin,spin ;

  printf("Here is the smeared meson array \n");

  for( ismear = 0 ; ismear < nosmear ; ++ismear)
  {
    printf(">>>>>> smearing function %d \n",ismear);
    for(t=0 ; t < nt ; ++t)
    {
      printf("time slice = %d \n",t);

      for(colour=0 ; colour < 3 ; ++colour)
	for(spin=0 ; spin < 4 ; ++spin)
	  for(ic=0 ; ic < 3 ;++ic)
	    for(ispin=0 ; ispin < 4 ;++ispin)
	      printf("src_c_d %d %d   sink_c_d %d %d   (%g , %g) \n",
		     colour,spin,ic,ispin,
		     (*(data + MESON_WHERE)).real , (*(data + MESON_WHERE)).imag   );


    } /** end of the loop over time ****/
  }/** end th loop over the smearing functions ***/

}



/*
 * Time the two different ways of writing a gauge configuration
 * to disk.
 *
 *
 */


void time_gauge_write()
{
  char file_one[] = "/pfs/mcneile/unit_one" ;
  /** dummmy paramters used in the gauge configuration write *****/
  /** timing information *******/
  Real ts , te ;

  printf("Time the varous routines which write to /pfs \n");


  /*** use the standard MILC gauge configuration writing routine *****/
  ts=dclock();
   save_parallel(file_one);
  te=dclock() - ts ;
  node0_printf("Time for writing %s using >>save_parallel<< is %g sec \n",file_one,te);

  fflush(stdout);


#ifdef CUT_THIS_TIMING_STUFF_OUT
  /**** use Claude's faster gauge configuration writing routine ******/
  ts=dclock();
   save_checkpoint1(file_two,c1,c2);
  te=dclock() - ts ;
  node0_printf("Time for writing %s using >>save_checkpoint1<< is %g sec \n",file_one,te);

  fflush(stdout);

#endif


  node0_printf("End of the debug CODE premature stop \n");


    terminate(0);
}



/*
 *   Calculate the pion propagator using the wilson vector
 *   from each source colour and spin component.
 *
 *   flag = 0  :: call the set up part of the code
 *   flag = 1  :: calculate the contribution to the pion correlators
 *                from the wilson vector
 *   flag = 2  :: write to standard output the pion correlator
 */

#define MAX_PION 110

void light_quark_pion(int flag)
{
  int i ;
  register site *s;
  int ic, ispin ;
  static double pion[MAX_PION] ;
  int t ;
  /***********************************************************/

  if( nt > MAX_PION ) 
  {
    node0_printf("light_quark_pion:::  not enough room for the pion correlator\n");
    terminate(1);
  }


  if( flag == 0 )
  {
    /**       setup    *****/
    for(t=0 ; t < nt ; ++t)
      pion[t] = 0.0 ;
    node0_printf("PIon correlator calculation is set up");

    return ;
  }
  else if ( flag == 2 )
  {
    /** sum the correlators over all the nodes ****/
    g_vecdoublesum(&pion[0] ,nt );

    /**** output the pion correlators to the screen ****/
    node0_printf("Pion correlator \n");
    for(t=0 ; t < nt ; ++t)
      node0_printf("%d   %e\n",t,pion[t]);

    return ;
  }


  FORALLSITES(i,s)
  {
    for(ic=0 ;ic < 3 ; ++ic )
      for(ispin=0 ;ispin < 4 ; ++ispin )
      {
	pion[s->t] += pow(s->psi.d[ispin].c[ic].real , 2.0) ;
	pion[s->t] += pow(s->psi.d[ispin].c[ic].imag , 2.0) ;

      }

  }  /*** end the loop over points on the lattice *******/



}


/*
 *  Create a random quark wilson_vector
 *
 */

void random_wilson_vec(int colour, int spin)
{
  int i ;
  register site *s;
  int ic, ispin ;
  long seed ;
  /***********************************************************/

  /** this is a very bad trick to get a different random propagator
    for each colour spin compoenent ***/

  seed = (long) -107 + 5*(colour + 3*spin) ;



  FORALLSITES(i,s)
  {
    for(ic=0 ;ic < 3 ; ++ic )
      for(ispin=0 ;ispin < 4 ; ++ispin )
      {
	s->psi.d[ispin].c[ic].real  = ran1(&seed) ;
	s->psi.d[ispin].c[ic].imag  = ran1(&seed) ;
      }
  }


  node0_printf("Random wilson vector created for spin = %d and colour = %d\n",spin,colour);



}





/**
 **  Check an mpl routine which sums up a vector
 **  of numbers over all the nodes
 **/


void check_g_vecdoublesum()
{
  const int total = 3 ;
  double *dpt ;
  int n = mynode() ;
  int i , n_pt ; 
  double ans ; 

  dpt = (double *)calloc(total, sizeof(double)) ;

  /***** set up the data *****/
  for(i=0 ; i < total ; ++i)
    *(dpt + i ) = (n+1)*sin( (double) (i+1)) ;

  g_vecdoublesum(dpt, total );

  for(i=0 ; i < total ; ++i)
    printf("node %d v[%d] = %g \n",n,i,*(dpt + i ));

  free(dpt) ;

  /***** now print out the check of the numbers ******/
  if( mynode() == 0 )
  {
    for(i=0 ; i < total ; ++i)
    {
      ans = 0.0 ; 
      for( n_pt = 0 ; n_pt < number_of_nodes ; ++n_pt)
	ans += (n_pt + 1)*sin( (double) (i+1)) ;
      printf("SERIAL total[%d] = %g\n",i, ans) ;
    }

  }

  terminate(0);

}






/**
 **  Check an mpl routine which sums up a vector
 **  of complex numbers over all the nodes
 **/


void check_g_veccomplexsum()
{
  const int total = 3 ;
  complex *dpt ;
  int n = mynode() ;
  int i ; 
  double ans_re , ans_im ; 
  int n_pt ; 

  dpt = (complex *)calloc(total, sizeof(complex)) ;

  /***** set up the data *****/
  for(i=0 ; i < total ; ++i)
  {
    (dpt + i )->real =  0.5 * (n+1)*sin( (double) (i+1)) ;
    (dpt + i )->imag =  -0.5 * (n+1)*cos( (double) (i+1)) ;
  }
  g_veccomplexsum(dpt, total );

  for(i=0 ; i < total ; ++i)
    printf("node %d v[%d] = %g %g \n",n,i,(dpt + i )->real, (dpt + i )->imag);


  /***** now print out the check of the numbers ******/
  if( mynode() == 0 )
  {
    for(i=0 ; i < total ; ++i)
    {
      ans_re = ans_im = 0.0 ; 
      for( n_pt = 0 ; n_pt < number_of_nodes ; ++n_pt)
      {
	ans_re += 0.5 * (n_pt+1)*sin( (double) (i+1)) ;
	ans_im += -0.5 * (n_pt+1)*cos( (double) (i+1)) ;
      }
      printf("SERIAL total[%d] = (%g , %g) \n",i, ans_re, ans_im) ;
    }

  }



  free(dpt) ;



  terminate(0);

}


/*
 *  Create a fake Wilson vector which contains the node 
 *  number in each component of the wilson vector.
 *
 */

void fake_node_psi(int spin, int colour)
{
  int i ;
  register site *s;
  int j,k ;
  Real fact = spin + 1 + 4*colour ; 
  

  FORALLSITES(i,s)
  {
    for(j=0 ; j < 4 ; ++j)
      for(k=0 ; k < 3 ; ++k)
      {
	s->psi.d[j].c[k].real  =  (this_node + 0.5 ) * fact ;
	s->psi.d[j].c[k].imag  = -(this_node + 0.5) * fact ;
	
      }
  }


}


/*
 *   Calculate the pion correlator on each node
 *
 */

void node_pion(int spin, int colour)
{
  int i ;
  register site *s;
  int j,k ;
  Real g_re ,g_im , pion = 0.0 ;

  FORALLSITES(i,s)
  {
    for(j=0 ; j < 4 ; ++j)
      for(k=0 ; k < 3 ; ++k)
      {
	g_re = s->psi.d[j].c[k].real  ;
	g_im = s->psi.d[j].c[k].imag  ;

	pion = g_re*g_re + g_im*g_im ;
      }
  }


  printf("pion (spin=%d, color= %d) for node %d = %e\n",spin,colour,this_node,pion);

}









/*
 *  Zero the Wilson line
 */

void zero_w_line()
{
  int i ;
  register site *s;

  int ic, jc ;
  /***** --------------------------------------------------  *****/


  FORALLSITES(i,s)
  {
    for(ic=0 ;ic < 3 ; ++ic )
      for(jc=0 ;jc < 3 ; ++jc )
      {
	s->w_line.e[ic][jc].real = 0.0 ;
	s->w_line.e[ic][jc].imag = 0.0  ;
      }
  }




}




/*
 *  Wrapper function to check the calc_matrix
 *
 */




void check_calc_matrix()
{
  void calc_vary_matrix() ;

  zero_strip_quark() ;
  zero_w_line(); 

  node0_printf("I am calculating the static variational code with ZERO input\n");

  calc_vary_matrix() ;  /** calculate the static variational matrix **/

  terminate(0) ;

}
