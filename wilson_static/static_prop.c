/****************** static_prop.c ************************************/
/* MIMD version 7 */
/*
  $Header: /lqcdproj/detar/cvsroot/milc_qcd/wilson_static/static_prop.c,v 1.4 2007/10/09 21:02:32 detar Exp $
  Calculate the Wilson lines required for the static
  quark propgator -- this  is used for the 2-pt function.
 

  W(x,t)      =   U(x,0)_0 U(x,1)_0 .... U(x,t-1)_0        t > 0

  W(x,t = 0 ) = 1
                       dagger       dagger        dagger
  W(x,-t)     =   U(x,-1)_0 U(x,-2)_0 .... U(x,-t)_0        t > 0

 */

#include "w_static_includes.h"

static char rcs_id[] = "$Header: /lqcdproj/detar/cvsroot/milc_qcd/wilson_static/static_prop.c,v 1.4 2007/10/09 21:02:32 detar Exp $" ;


/*  
 *  Calculate the gauge part of the static propagator.
 *  This is used for the two point calculation.
 *
 * Note smear_w_line[0] is used as work space in this
 * routine.
 */


void static_prop() 
{
  register int i;
  register site *st;
  msg_tag *tag;
  int tloop ;
  int nthalf = nt/2 ;
  /*************---------**********-------------************/


  /* Initialise the gauge part of the  propagator  ***/
  setup_static_prop() ;

  /*
   *   Calculate the static propagator for positive time
   *
   *   W(t+1) = W(t) U_4(t)
   */

  for(tloop=1 ; tloop <= nthalf ; ++tloop)
  {

    /* The smear_w_line[0] object is used as work space ***/
    FORALLSITES(i,st)
    {
	mult_su3_nn(&(st->w_line),  &(st->link[TUP]), &(st->smear_w_line[0]));
    }

    /* Pull the w(t)*u(t)  from the previous time slice ***/
    tag=start_gather_site( F_OFFSET(smear_w_line[0]), sizeof(su3_matrix),
		     TDOWN, EVENANDODD, gen_pt[0] );
    wait_gather(tag);

    FORALLSITES(i,st)
    {
     if( st-> t == tloop )  
	su3mat_copy((su3_matrix *) gen_pt[0][i], &(st->w_line));
    }
    cleanup_gather(tag);


		     
  } /* end the loop over time slice ***/
  
  /*
   *   Calculate the static propagator for negative time
   *
   *   W(t) = W(t+1) U(t)^[ dagger ]
   *                    4
   */


  for(tloop=nt-2 ; tloop > nthalf ; --tloop)
  {
    /* Pull the Wilson lines from the previous time slice ***/

    tag=start_gather_site( F_OFFSET(w_line), sizeof(su3_matrix),
		     TUP, EVENANDODD, gen_pt[0] );
    wait_gather(tag);

    FORALLSITES(i,st)
    {
      if( st-> t == tloop )
	mult_su3_na((su3_matrix *) gen_pt[0][i] , &(st->link[TUP]),  &(st->w_line));
    }
    cleanup_gather(tag);

		     
  } /* end the loop over time slice ***/

}






/*  
 *  Calculate the gauge part of the static propagator.
 *
 *  This is used for the four fermion operator calculation
 *
 * Note smear_w_line[0] is used as work space in this
 * routine.
 */


void static_prop_bparam() {
  register int i;
  register site *st;
  msg_tag *tag;
  int tloop ;
  int nthalf = nt/2 ;
  /*************---------**********-------------************/


  /* Initialise the gauge part of the  propagator  ***/
  setup_static_prop() ;

  /*
   *   Calculate the static propagator for positive time
   *
   *   W(t+1) = W(t) U_4(t)
   */

  for(tloop=1 ; tloop <= nthalf ; ++tloop)
  {

    /* The smear_w_line[0] object is used as work space ***/
    FORALLSITES(i,st)
    {
	mult_su3_nn(&(st->w_line),  &(st->link[TUP]), &(st->smear_w_line[0]));
    }

    /* Pull the w(t)*u(t)  from the previous time slice ***/
    tag=start_gather_site( F_OFFSET(smear_w_line[0]), sizeof(su3_matrix),
		     TDOWN, EVENANDODD, gen_pt[0] );
    wait_gather(tag);

    FORALLSITES(i,st)
    {
     if( st-> t == tloop )  
	su3mat_copy((su3_matrix *) gen_pt[0][i], &(st->w_line));
    }
    cleanup_gather(tag);


		     
  } /* end the loop over time slice ***/
  
  /*
   *   Calculate the static propagator for negative time
   *
   *   W(t) = W(t+1) U(t)^[ dagger ]
   *                    4
   */


  for(tloop=nt-2 ; tloop > nthalf ; --tloop)
  {
    /* Pull the Wilson lines from the previous time slice ***/

    tag=start_gather_site( F_OFFSET(w_line), sizeof(su3_matrix),
		     TUP, EVENANDODD, gen_pt[0] );
    wait_gather(tag);

    FORALLSITES(i,st)
    {
      if( st-> t == tloop )
	mult_su3_na((su3_matrix *) gen_pt[0][i] , &(st->link[TUP]),  &(st->w_line));
    }
    cleanup_gather(tag);

		     
  } /* end the loop over time slice ***/

}



/*
 *  Initialize the static propagator -- set up
 *  the boundary conditions for the recursive equations used
 *  to generate the Wilson line.
 *
 *  This is useful for the two point calculation
 *
 *  W(t= 0 )  =  1
 *
 *  W(t = nt-1)  = U_4^[ dagger ]
 *
 *  W(t) = 0 every where else
 */

void setup_static_prop()
{
  int i,j,k ;
  register site *s;
   int last_slice = nt -1 ;


    FORALLSITES(i,s)
    {

      if ( s->t == 0)
      {
	/*** First time slice ***************/
	for(k=0; k < 3 ; ++k)
	{
	  for(j=k+1; j < 3 ; ++j)
	  {
	    s->w_line.e[k][j].real = 0.0 ;
	    s->w_line.e[k][j].imag = 0.0 ;
	    
	    s->w_line.e[j][k].real = 0.0 ;
	    s->w_line.e[j][k].imag = 0.0 ;
	  }
	  s->w_line.e[k][k].real = 1.0 ;
	  s->w_line.e[k][k].imag = 0.0 ;
	}

      } /**** end of the first time slice ****************/
      else if ( s->t == last_slice  )
      {
	su3_adjoint(&(s->link[TUP]), &(s->w_line) ) ;  
      } /**** end of the last time slice ****************/
      else
      {
	for(k=0; k < 3 ; ++k)
	  for(j=0; j < 3 ; ++j)
	  {
	    s->w_line.e[k][j].real = 0.0 ;
	    s->w_line.e[k][j].imag = 0.0 ;
	  }
      }  /**** end of middle of the time  ****************/



    } /** end of the loop over the sites ****/



}






/*
 *  Initialize the static propagator -- set up
 *  the boundary conditions for the recursive equations used
 *  to generate the Wilson line.
 *
 *  This is useful for the B parameter calculation
 *
 *  W(t= 0 )  =  1
 *
 *  W(t = nt-1)  = U_4(nt-1)
 *
 *  W(t) = 0 every where else
 */

void setup_static_prop_bparam()
{
  int i,j,k ;
  register site *s;
   int last_slice = nt -1 ;


    FORALLSITES(i,s)
    {

      if ( s->t == 0)
      {
	/*** First time slice ***************/
	for(k=0; k < 3 ; ++k)
	{
	  for(j=k+1; j < 3 ; ++j)
	  {
	    s->w_line.e[k][j].real = 0.0 ;
	    s->w_line.e[k][j].imag = 0.0 ;
	    
	    s->w_line.e[j][k].real = 0.0 ;
	    s->w_line.e[j][k].imag = 0.0 ;
	  }
	  s->w_line.e[k][k].real = 1.0 ;
	  s->w_line.e[k][k].imag = 0.0 ;
	}

      } /**** end of the first time slice ****************/
      else if ( s->t == last_slice  )
      {

	s->w_line = s->link[TUP] ;

      } /**** end of the last time slice ****************/
      else
      {
	for(k=0; k < 3 ; ++k)
	  for(j=0; j < 3 ; ++j)
	  {
	    s->w_line.e[k][j].real = 0.0 ;
	    s->w_line.e[k][j].imag = 0.0 ;
	  }
      }  /**** end of middle of the time  ****************/



    } /** end of the loop over the sites ****/



}



