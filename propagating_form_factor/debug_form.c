/**************************** debug_form.c ****************************/
/* MIMD version 6 */
/*
 *  This file contains a number of functions, thought
 *  to be useful to help debug the HQET code.
 *
 */


#include "prop_form_includes.h"
#ifdef DEBUGDEF
#include "debug_form.h"


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





void dump_heavy_smear_func(void)
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
    z = s->heavy_smear_func[0];
    printf("heavy_smear_func[%d %d %d %d] = %g %g\n",s->x,s->y,s->z,s->t,z.real,z.imag); 

  }



}






void dump_seq_smear_func(void)
{
  int i ;
  register site *s;
  int j,k ;
  double g_re ,g_im ;
  int t ;
  complex z ;
  int ip ;

  /*****------------------------------**********/

  if( this_node == 0 )
  {
    printf("Here are the \"seq_smear_func\" smearing functions\n"); 
  }

  for( ip = 0 ; ip < MAXPMOM ; ++ip )
  {
    
    if( this_node == 0 ) 
    {
      printf("---------------------------\n"); 
      printf("momentum pointer = %d\n",ip); 
      printf("---------------------------\n"); 
    }

    FORALLSITES(i,s)
    {
      z = s->seq_smear_func[ ip ]     ;
      printf("seq_smear_func[%d %d %d %d] = %g %g\n",s->x,s->y,s->z,s->t,z.real,z.imag); 
    }

  }  /*** end of the loops over the momentum values ****/



}






/**** end of debug routines *****/

#endif

/*** no debug options so no not include the code ***/



