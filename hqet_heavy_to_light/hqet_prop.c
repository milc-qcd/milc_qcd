/*************************** hqet_prop.c *****************************/ 
/* MIMD version 6 */
/* $Header: /lqcdproj/detar/cvsroot/milc_qcd/hqet_heavy_to_light/hqet_prop.c,v 1.2 2005/03/31 00:01:58 detar Exp $

  A collection of routines to generate hqet propagators

  If v1 == v2, then a two point function is 
   generated (and tcurrent is ignored).

  If v1 != v2 then a three point function is generated,
    for t < tcurrent then update hqet propagator with velocity v1


  v1 is the velocity of the hqet propagator for the
  first set of time slices

  Work space required in:   lattice.h
 
  work_space :: su3 matrix
                s->tempvec[dir])   dir = 1 -->3 XUP .. TUP

  This routine generates a possible hqet propagators from
  0 to tend, 


 This code has been modified from the 
 dslash for Kogut-Susskind fermions  


    Evolution equation for HQET propagators
    -----------------------------------------

      See Mandula and Ogilvie
      Phys.Rev D45 (1992) 2183
  
                   dagger
       G(x,t+1) = U(x,t)_[4]
               ( 1 +  \sum_{k=1}^{3} \frac{-i v_{\mu}{v_{0}}
               ( \Delta^{+}_{k} + \Delta^{-}_{k} )    )   G(x,t)

  Subroutine arguments

     hqet_prop :: (su3_matrix)  heavy quark propagator
     tend      :: final time slice to generate to
     tcurrent  :: timeslice to insert the current
     tstart    :: the timeslice the source is on
     v1        :: the velocity for t < tcurrent
     v2        :: the velocity for t > tcurrent

 */

#include "hqet_light_includes.h"

#include<assert.h>

#ifdef DEBUGDEF
#include DEBUGDEF
#endif

static char rcs_id[] = "$Header: /lqcdproj/detar/cvsroot/milc_qcd/hqet_heavy_to_light/hqet_prop.c,v 1.2 2005/03/31 00:01:58 detar Exp $" ;

void generate_hqet_prop(field_offset hqet_prop, int tstart, int tend , int tcurrent, 
int v1, int v2) 
{
  register int i,j,t;
  register site *s;
  msg_tag *tag[8];
  int tloop ;
  int nthalf = nt/2 ;
  int dir ;
  int jj ;
  register su3_vector *a,*b1,*b2,*b3,*b4;
  complex unorm[4] ;
  /*************---------**********-------------************/

  assert( tend > tstart ) ; 


  /** check to see whether to a two or three point function ***/
  if( v1 == v2 )
  {
    tcurrent = tend + 4 ;
  }


  /** set up the renormalized velocity ****/
  unorm[TUP].real =0.0 ; unorm[TUP].imag = 0.0 ;
  for(jj =XUP ; jj <= ZUP ; ++jj)
  {
    unorm[jj].real =0.0 ; 
    unorm[jj].imag = velocity[v1][jj] /( 2.0 * velocity[v1][TUP] );
  }

  ++tstart ; /*** tloop is the time slice of the required propagator ****/

  for(tloop=tstart ; tloop <= tend ; ++tloop)
  {

    if( tloop >  tcurrent ) 
    {
      for(jj =XUP ; jj <= ZUP ; ++jj)
      {
	unorm[jj].real =0.0 ; 
	unorm[jj].imag = velocity[v2][jj] /( 2.0 *  velocity[v2][TUP] ) ;
      }
    }


    /***-----  spatial parts of the evolution equation ------***/

    /* Start gathers of the propagator from positive directions */
    for(dir=XUP; dir<=ZUP; dir++)
    {
	tag[dir] = start_gather_site( hqet_prop, sizeof(su3_matrix), dir, EVENANDODD,
	    gen_pt[dir] );
    }

    /*
     *  tmpvec[i] = U_{\mu}^{\dagger}(x-i,t) s(x-i)
     */


    /* Multiply by adjoint matrix at other sites */
    FORALLSITES(i,s)
    {
      for(dir=XUP; dir<=ZUP; dir++)
      {
	mult_su3_an_z(unorm[dir] , &(s->link[dir]),
		    (su3_matrix *)F_PT(s,hqet_prop), &(s->tempvec[dir]) );
      }

    }


    /* Start gathers from negative directions */
    for( dir=XUP; dir <= ZUP; dir++)
    {
	tag[OPP_DIR(dir)] = start_gather_site( F_OFFSET(tempvec[dir]),
	    sizeof(su3_matrix), OPP_DIR( dir), EVENANDODD ,
	    gen_pt[OPP_DIR(dir)] );
    }



    /* Wait gathers from positive directions *******/
    for(dir=XUP; dir<=ZUP; dir++)
    {
      wait_gather(tag[dir]);
    }



    FORALLSITES(i,s)
    {
      su3mat_copy((su3_matrix *)F_PT(s,hqet_prop) , &(s->tempvec[ TUP ]) );

	for(dir=XUP; dir<=ZUP; dir++)
	{
	  mult_su3_nn_z_inc( unorm[dir] ,&(s->link[dir]),(su3_matrix *)(gen_pt[dir][i]),
			    &(s->tempvec[ TUP ]) );
	}

    }

    /* Wait gathers from negative directions */
    for(dir=XUP; dir<=ZUP; dir++)
    {
      wait_gather(tag[OPP_DIR(dir)]);
    }



    FORALLSITES(i,s)
    {
      sub3_su3_matrix((su3_matrix *)(gen_pt[XDOWN][i]) ,  
		      (su3_matrix *)(gen_pt[YDOWN][i]) ,
		      (su3_matrix *)(gen_pt[ZDOWN][i]) ,
		      &(s->tempvec[ TUP ])) ;

    }



    /*
     *   Calculate the hqet propagator for positive time
     *
     *   W(t+1) =  U_4(t)^{\dagger}   W(t)
     *
     *  Evolve the HQET propagator forward in time
     *
     *
     * Note: the structure tempvec[XUP] is used as new work space
     */


    FORALLSITES(i,s)
    {
      mult_su3_an(&(s->link[TUP]), (su3_matrix *) &(s->tempvec[ TUP ])  , &(s->tempvec[XUP]));
    }

    /* Pull the w(t)*u(t)  from the previous time slice ***/
    tag[TDOWN]=start_gather_site( F_OFFSET(tempvec[XUP]), sizeof(su3_matrix),
		     TDOWN, EVENANDODD, gen_pt[TDOWN] );
    wait_gather(tag[TDOWN]);

    FORALLSITES(i,s)
    {
     if( s-> t == tloop )  
	su3mat_copy((su3_matrix *) gen_pt[TDOWN][i], (su3_matrix *)   F_PT(s,hqet_prop)  );


    }


    /*** free up the message passing bufffers   ****/
    cleanup_gather(tag[TDOWN]);

    for(dir=XUP; dir<=ZUP; dir++)
    {
      cleanup_gather(tag[dir]);
      cleanup_gather(tag[OPP_DIR(dir)]);
    }


		     
  } /* end the loop over time slice ***/


}  /*** end of the function ****/

