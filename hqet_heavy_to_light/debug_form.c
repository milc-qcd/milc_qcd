/*  $Header: /lqcdproj/detar/cvsroot/milc_qcd/hqet_heavy_to_light/debug_form.c,v 1.1 2005/02/23 00:05:07 detar Exp $
 *  This file contains a number of functions, thought
 *  to be useful to help debug the HQET code.
 *
 */

#include "hqet_light_includes.h"
#ifdef DEBUGDEF
#include DEBUGDEF

/*
 *   Calculate the HQET pion correlator 
 *   This has no physical meaning but it is
 *   a useful debug test.
 */

#ifdef  PION_DIM
 error: The macro pion_dim has already been defined
#endif

#define PION_DIM 100

void hqet_pion()
{
  int i ;
  register site *s;
  int j,k ;
  double g_re ,g_im ;
  double pion[PION_DIM  ] ;
  int t ;
  /*****------------------------------**********/

  if( nt > PION_DIM) 
  {
    printf("Not enough space of the pion propagator") ;
    exit(1);
  }


  for( t=0 ; t < nt  ; ++t)
    pion[t] = 0.0 ;

  FORALLSITES(i,s)
  {
    for(j=0 ; j < 3 ; ++j)
      for(k=0 ; k < 3 ; ++k)
      {
	g_re = s->heavy_prop.e[j][k].real  ;
	g_im = s->heavy_prop.e[j][k].imag  ;


	pion[s->t] += pow(g_re, 2.0) + pow(g_im, 2.0) ;

      }
  }


  for( t=0 ; t < nt ; ++t)
    {
      g_doublesum( &pion[ t ]  ) ;
      if( this_node == 0 )       printf("hqet_pion[ %d ] = %g\n",t,pion[ t ] );
    }

}

#undef PION_DIM





/*
 *   Calculate the "fake" pion correlator  for 
 *   a single Wilson vector
 *   This has no physical meaning but it is
 *   a useful debug test.
 */

#ifdef  PION_DIM
  #error  The macro pion_dim has already been defined
#endif

#define PION_DIM 100

void wilson_vector_pion(field_offset prop)
{
  int i ;
  register site *s;
  int j,k ;
  double g_re ,g_im ;
  double pion[PION_DIM  ] ;
  int t ;
  complex z ;

  /*****------------------------------**********/

  if( nt > PION_DIM) 
  {
    printf("Not enough space for the pion propagator\n") ;
    exit(1);
  }


  for( t=0 ; t < nt  ; ++t)
    pion[t] = 0.0 ;

  FORALLSITES(i,s)
  {
    for(j=0 ; j < 4 ; ++j)
      for(k=0 ; k < 3 ; ++k)
      {

	z = ((wilson_vector *)  F_PT(s,prop) ) ->d[j].c[k] ; 

	g_re = z.real  ;
	g_im = z.imag  ;


	pion[s->t] += pow(g_re, 2.0) + pow(g_im, 2.0) ;

      }
  }


  for( t=0 ; t < nt ; ++t)
    printf("wilson_vector_pion[ %d ] = %g\n",t,pion[ t ] );


}



void dump_wilson_vector(field_offset prop)
{
  int i ;
  register site *s;
  int j,k ;
  double g_re ,g_im ;
  int t ;
  complex z ;

  /*****------------------------------**********/

  FORALLSITES(i,s)
  {
    printf("xyzt = [%d %d %d %d]\n",s->x,s->y,s->z,s->t); 

    for(j=0 ; j < 4 ; ++j)
      for(k=0 ; k < 3 ; ++k)
      {

	z = ((wilson_vector *)  F_PT(s,prop) ) ->d[j].c[k] ; 

	g_re = z.real  ;
	g_im = z.imag  ;

	printf("w_vec[%d][%d] = (%g,%g)\n",j,k,g_re,g_im) ;
      }

  }



}




void zero_wilson_vector(field_offset prop)
{
  int i ;
  register site *s;
  int j,k ;
  double g_re ,g_im ;
  int t ;
  complex z ;

  /*****------------------------------**********/

  FORALLSITES(i,s)
  {
    for(j=0 ; j < 4 ; ++j)
      for(k=0 ; k < 3 ; ++k)
      {

	((wilson_vector *)  F_PT(s,prop) ) ->d[j].c[k].real = 0.0  ; 
	((wilson_vector *)  F_PT(s,prop) ) ->d[j].c[k].imag = 0.0  ; 
      }

  }

  if( this_node == 0 ) printf("wilson_vector is set equal to ZERO\n") ; 

}




/*  
 *  Load a Wilson vector into a spin 
 * 
 * This is not a very general function.
 */

void loadup_debug(int spin)
{
  int i ;
  register site *s;
  int i_spin, ic ; 
  /***************************************************************/

  FORALLSITES(i,s)
  {

    for(i_spin = 0 ; i_spin < 4 ; ++i_spin )
      for(ic =0 ; ic < 3 ; ++ic)
      {
	s->quark_sequential.d[spin].d[i_spin].c[ic]  = s->psi.d[i_spin].c[ic]  ;
      }


  }  /*** end of the loop over lattice sites***/


}


/*
 *
 *
 *
 */

void  dump_spin_wilson_vector(spin_wilson_vector *in ) 
{
  int pt = 0 ;
  int i , j , k=0 ;

  for(i = 0 ; i < 4 ; ++i)
    for(j = 0 ; j < 4 ; ++j)
      for(k = 0 ; k < 3 ; ++k) 
      {
	printf("color = %d spin = [%d ,  %d] ==== %f %f\n",k,i,j,in->d[i].d[j].c[k].real,
	       in->d[i].d[j].c[k].imag);

	++pt ;
      }



}



/***  create a fake source for the sequential inversion ***/


void  fake_spin_wilson_src(spin_wilson_vector *in, int color ) 
{
  int pt = 0 ;
  int i , j , k=0 ;

  for(i = 0 ; i < 4 ; ++i)
    for(j = 0 ; j < 4 ; ++j)
      for(k = 0 ; k < 3 ; ++k)
      {
	in->d[i].d[j].c[k].imag = 0.0  ;
	in->d[i].d[j].c[k].real = 0.0  ;

      }




  for(i = 0 ; i < 4 ; ++i)
  {
    in->d[i].d[i].c[color].real = 1.0 ; 
  }


/*******
  if( this_node == 0 ) 
    printf("Here is the fake source for the HQET inverter\n"); 

  dump_spin_wilson_vector( in ) ; 
****************************************/

}



void  zero_spin_wilson(spin_wilson_vector *in ) 
{
  int pt = 0 ;
  int i , j , k=0 ;

  for(i = 0 ; i < 4 ; ++i)
    for(j = 0 ; j < 4 ; ++j)
      for(k = 0 ; k < 3 ; ++k)
      {
	in->d[i].d[j].c[k].imag = 0.0 ;
	in->d[i].d[j].c[k].real = 0.0 ;

      }


  if( this_node == 0 ) printf("DEBUG:: spin_wilson_vector set equal to zero  \n"); 

}


void fake_seq_src(int color)
{
  int i ;
  register site *s;


  FORALLSITES(i,s)
  {
    if( s->t == 0 && s->x == 0 && s->y == 0 && s->z == 0)
    {
      fake_spin_wilson_src((spin_wilson_vector *)&(s->quark_sequential), color )  ;
    }
    else 
      zero_spin_wilson((spin_wilson_vector *)&(s->quark_sequential))  ;



  } /*** end the loop over the lattice sites ***/


  tf  = 0.0  ;
  if( this_node == 0)
  {
    printf("DEBUG tf = %d\n\n",tf); 

  }


}





void gamma_five_sequential(int color)
{
  int i ;
  void mult_gamma5_sink(field_offset ans);
  void unit_quark_sequential( int color);


  unit_quark_sequential(color) ; 
  for( i = 0 ; i < 4 ; ++i)
    mult_gamma5_sink( F_OFFSET(quark_sequential.d[ i  ]) ) ; 
  


}


/**
    Multiply gamma_5 into the sinl of a wilson quark propagator
    at all timeslices. This is useful to help debug the sequential source
    routine.

**/


void mult_gamma5_sink(field_offset ans)
{
  register site *s; 
  int i ; 
  wilson_vector gammafive_quark ; 

  FORALLSITES(i,s) 
  {
      mult_by_gamma((wilson_vector *)F_PT(s,ans), 
		    &gammafive_quark , GAMMAFIVE );

      *((wilson_vector *)F_PT(s,ans)) = gammafive_quark  ; 

  } /** end the loop over the sites ***/




}




/*
 *  Set the zonked propagator equal to the unit
 *  matrix
 */

void unit_quark_zonked(int color)
{
  register site *s; 
  int i ; 
  int ic, i_spin, j_spin ; 

  FORALLSITES(i,s) 
  {
    for(ic= 0 ; ic < 3 ; ++ic)
      for(i_spin = 0 ; i_spin < 4 ; ++ i_spin ) 
	for(j_spin = 0 ; j_spin < 4 ; ++j_spin ) 
	{
	  s-> quark_zonked.d[i_spin].d[j_spin].c[ ic ].real = 0.0 ; 
	  s-> quark_zonked.d[i_spin].d[j_spin].c[ ic ].imag = 0.0 ; 
	}




      for(i_spin = 0 ; i_spin < 4 ; ++ i_spin ) 
      {
	s-> quark_zonked.d[i_spin].d[i_spin].c[color].real = 1.0 ; 
      }



  } /** end the loop over the sites ***/

  if( this_node == 0) printf("quark_zonked set to the UNIT matrix\n"); 


}









/*
 *  Set the sequential propagator equal to the unit
 *  matrix
 */

void unit_quark_sequential( int color)
{
  register site *s; 
  int i ; 
  int ic, i_spin, j_spin ; 

  FORALLSITES(i,s) 
  {

    for(ic= 0 ; ic < 3 ; ++ic)
      for(i_spin = 0 ; i_spin < 4 ; ++ i_spin ) 
	for(j_spin = 0 ; j_spin < 4 ; ++j_spin ) 
	{
	  s->quark_sequential.d[i_spin].d[j_spin].c[ic].real = 0.0 ; 
	  s->quark_sequential.d[i_spin].d[j_spin].c[ic].imag = 0.0 ; 
	}



      for(i_spin = 0 ; i_spin < 4 ; ++ i_spin ) 
      {
	s->quark_sequential.d[i_spin].d[i_spin].c[color].real = 1.0 ; 
      }



  } /** end the loop over the sites ***/

  if( this_node == 0) printf("quark_sequential set to the UNIT matrix\n"); 


}



void dump_quark_sequential(void)
{
  register site *s; 
  int i ; 

  if( this_node ==0 ) printf("Here is the sequential quark propagator \n"); 


  FORALLSITES(i,s) 
  {
    printf("xyzt = [%d %d %d %d]\n",s->x ,s->y ,s->z, s->t   ); 
    dump_spin_wilson_vector(&(s->quark_sequential)) ; 
  }


}



void dump_quark_zonked(void)
{
  register site *s; 
  int i ; 

  if( this_node ==0 ) printf("Here is the zonked quark propagator \n"); 


  FORALLSITES(i,s) 
  {
    printf("xyzt = [%d %d %d %d]\n",s->x ,s->y ,s->z, s->t   ); 
    dump_spin_wilson_vector(&(s->quark_zonked )) ; 
  }


}






/*
 *  For housekeeping purposes, zero the HQET propagator
 *
 */


void zero_hqet(void )
{
  register site *s;
  int i ; 
  int ic, jc ; 

  FORALLSITES(i,s)
  {
    for(ic = 0 ; ic < 3 ; ++ic)
      for(jc = 0 ; jc < 3 ; ++jc)
      {
	s->heavy_prop.e[ic][jc].real = 0.0 ; 
	s->heavy_prop.e[ic][jc].imag = 0.0 ; 
      }

  }


}





/*
 *  Set the HQET propagator equal to one 
 *
 */


void unit_hqet(void )
{
  register site *s;
  int i ; 
  int ic, jc ; 

  FORALLSITES(i,s)
  {
    for(ic = 0 ; ic < 3 ; ++ic)
      for(jc = 0 ; jc < 3 ; ++jc)
      {
	s->heavy_prop.e[ic][jc].real = 0.0 ; 
	s->heavy_prop.e[ic][jc].imag = 0.0 ; 
      }


      for(jc = 0 ; jc < 3 ; ++jc)
      {
	s->heavy_prop.e[jc][jc].real = 1.0 ; 
      }


  }


}





/**** end of debug routines *****/

#endif

/*** no debug options so no not include the code ***/


