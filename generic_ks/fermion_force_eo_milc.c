/****** fermion_force_eo_milc.c -- ******************/
/* MIMD version 7 */
/* Fermion force for a general quark action.  Not optimized.
* (Formerly fermion_force_general.c)
* D.T. 1/28/98, starting from gauge_stuff.c
* K.O. 3/99 Added optimized fattening for Asq actions
* D.T. 4/99 Combine force calculations for both mass quarks
* K.O. 4/99 Optimized force for Asq action
* S.G. 7/01, modified to use t_longlink and t_fatlink
* C.D. 10/02, consolidated quark_stuff.c and quark_stuff_tmp.c
*
* C.D. 3/05 Separated from quark_stuff.c

* In this directory, assume all paths connect even to odd sites, etc.
* Tabulate "backwards" paths (e.g. "XDOWN" is backward path to "XUP")
* as separate parity transforms of the fundamental paths.  They will
* generally need a negative sign in Dslash.  See bottom for a long
* comment on sign conventions.
*/

/*
 * 10/01/02, flopcount for ASQ_OPTIMIZED - C. DeTar
 * Fermion force: 253935 for eo_fermion_force_oneterm()
 * Fermion force: 433968 for eo_fermion_force_twoterms()
 */


#include "generic_ks_includes.h"	/* definitions files and prototypes */
#include "../include/fermion_links.h"
#include "../include/fermion_links_milc.h"
#include "../include/ks_action_paths.h"

#define GOES_FORWARDS(dir) (dir<=TUP)
#define GOES_BACKWARDS(dir) (dir>TUP)

/**********************************************************************/
/*   Version for a single set of degenerate flavors                   */
/**********************************************************************/

/* update the  momenta with the fermion force */
/* Assumes that the conjugate gradient has been run, with the answer in
   x_off, and dslash_site(x_off,x_off,ODD) has been run. (fills in x_off_odd) */
/* SEE LONG COMMENTS AT END */

/* update the  momenta with the fermion force */
/* Assumes that the conjugate gradient has been run, with the answer in
   x_off, and dslash_site(x_off,x_off,ODD) has been run. (fills in x_off_odd) */
/* SEE LONG COMMENTS AT END */

void eo_fermion_force_oneterm_site( Real eps, Real weight, field_offset x_off,
				    int prec, fermion_links_t *fl)
{
  /* Ignore prec for now */
  /* note CG_solution and Dslash * solution are combined in "x_off" */
  /* New version 1/21/99.  Use forward part of Dslash to get force */
  /* see long comment at end */
  /* For each link we need x_off transported from both ends of path. */
  /* For example weight = nflavors/4 */
  ks_action_paths *ap = get_action_paths(fl);
  register int i,dir,lastdir=-99,ipath,ilink;
  register site *s;
  int length;
  su3_matrix tmat,tmat2;
  Real ferm_epsilon, coeff;
  int num_q_paths = ap->p.num_q_paths;
  Q_path *q_paths = ap->p.q_paths;

#ifdef FFTIME
  int nflop = 0;
  double dtime;
#endif
  msg_tag *mtag0;
  half_wilson_vector *hw_tmp0,*hw_tmp1,*tmp_pt;

#ifdef FFTIME
dtime=-dclock();
#endif
  ferm_epsilon = 2.0*weight*eps;
  hw_tmp0 = (half_wilson_vector *)
    malloc(sites_on_node*sizeof(half_wilson_vector) );
  hw_tmp1 = (half_wilson_vector *)
    malloc(sites_on_node*sizeof(half_wilson_vector) );
  /* Use half_wilson_vectors to store x_off transported from ends of
     path.  0 component from forward end, 1 component from back end */

  /* loop over paths, and loop over links in path */
  for( ipath=0; ipath<num_q_paths; ipath++ ){
    if(q_paths[ipath].forwback== -1)continue;	/* skip backwards dslash */
    length = q_paths[ipath].length;

    /* path transport x_off and Dslash*x_off from far end.  Sometimes
	we need them at the start point of the path, and sometimes
	one link into the path - an optimization for later */
    path_transport_site( x_off, F_OFFSET(tempvec[0]),
      EVENANDODD, q_paths[ipath].dir, length );
    /* use tempvec[1] for transport from starting end */
    FORALLSITES(i,s){
      hw_tmp0[i].h[0]=s->tempvec[0];
      hw_tmp0[i].h[1]=*(su3_vector *)F_PT(s,x_off);
    }

    /* A path has (length+1) points, counting the ends.  At first
	 point, no "down" direction links have their momenta "at this
	 point". At last, no "up" ... */
    for( ilink=0; ilink<=length; ilink++ ){
      if(ilink<length)dir = q_paths[ipath].dir[ilink];
      else dir=NODIR;
      coeff = ferm_epsilon*q_paths[ipath].coeff;
      if( (ilink%2)==1 )coeff = -coeff;

      /* path transport x_off and Dslash*x_off from previous point */
      /* Use "half_wilson_vector" to handle pair of vectors -
	0 component is x_off from forward end, 1 component from back end */
      /* sometimes we don't need them */
      if( (ilink>0&&ilink<length) || 
	(ilink==length && GOES_BACKWARDS(lastdir)) ){

	if( GOES_FORWARDS(lastdir) ){
	  FORALLSITES(i,s){
            mult_adj_su3_mat_hwvec( &(s->link[lastdir]),
	      &(hw_tmp0[i]), &(hw_tmp1[i]) );
	  }
	  mtag0 = start_gather_field( hw_tmp1, 2*sizeof(su3_vector),
             OPP_DIR(lastdir), EVENANDODD, gen_pt[0] );
          wait_gather(mtag0);
          FORALLSITES(i,s){
	     hw_tmp0[i] = *(half_wilson_vector *)gen_pt[0][i];
	  }
          cleanup_gather(mtag0);
	}
	else{   /* GOES_BACKWARDS(lastdir) */
          mtag0 = start_gather_field( hw_tmp0, 2*sizeof(su3_vector),
               OPP_DIR(lastdir), EVENANDODD, gen_pt[0] );
          wait_gather(mtag0);
          FORALLSITES(i,s){
            mult_su3_mat_hwvec( &(s->link[OPP_DIR(lastdir)]),
                (half_wilson_vector *)(gen_pt[0][i]),
		&(hw_tmp1[i]) );
          }
	  tmp_pt = hw_tmp0; hw_tmp0 = hw_tmp1; hw_tmp1 = tmp_pt;
          cleanup_gather(mtag0);
	}
      }

      /* add in contribution to the force */
      /* Put antihermitian traceless part into momentum */
      FORALLSITES(i,s){
        if( ilink<length && GOES_FORWARDS(dir) ){
          uncompress_anti_hermitian( &(s->mom[dir]), &tmat2 );
	  su3_projector( &(hw_tmp0[i].h[0]), &(hw_tmp0[i].h[1]), &tmat );
	  if( s->parity==EVEN ){
	    scalar_mult_add_su3_matrix(&tmat2, &tmat,  coeff, &tmat2 );
	  }
	  else{
	    scalar_mult_add_su3_matrix(&tmat2, &tmat, -coeff, &tmat2 );
	  }
          make_anti_hermitian( &tmat2, &(s->mom[dir]) );
        }
	if( ilink>0 && GOES_BACKWARDS(lastdir) ){
          uncompress_anti_hermitian( &(s->mom[OPP_DIR(lastdir)]), &tmat2 );
	  su3_projector( &(hw_tmp0[i].h[0]), &(hw_tmp0[i].h[1]), &tmat );
	  if( s->parity==EVEN ){
	    scalar_mult_add_su3_matrix(&tmat2, &tmat, -coeff, &tmat2 );
	  }
	  else{
	    scalar_mult_add_su3_matrix(&tmat2, &tmat,  coeff, &tmat2 );
	  }
          make_anti_hermitian( &tmat2, &(s->mom[OPP_DIR(lastdir)]) );
        }
      }
      lastdir = dir;
    } /* end loop over links in path */
  } /* end loop over paths */

  free( hw_tmp0 ); free( hw_tmp1 );

#ifdef FFTIME
dtime += dclock();
node0_printf("FFTIME:  time = %e (1 mass) mflops = %e\n",dtime,
	     (Real)nflop*volume/(1e6*dtime*numnodes()) );
/**printf("TLENGTH: %d\n",tlength);**/
#endif
} /* eo_fermion_force_oneterm (version 7) */

/**********************************************************************/
/*   Version for two sets of flavors with distinct masses             */
/**********************************************************************/

void eo_fermion_force_twoterms_site( Real eps, Real weight1, Real weight2, 
				     field_offset x1_off, field_offset x2_off,
				     int prec, fermion_links_t *fl)
{
  /* Ignore prec for now */
  /* note CG_solution and Dslash * solution are combined in "x_off" */
  /* New version 1/21/99.  Use forward part of Dslash to get force */
  /* 4/14/99 combine force from two different mass quarks, (eg 2+1flavors) */
  /* see long comment at end */
  /* For each link we need x_off transported from both ends of path. */
  /* For example weight1 = nflavor1/4; weight2 = nflavor2/4 */
  ks_action_paths *ap = get_action_paths(fl);
  register int i,dir,lastdir=-99,ipath,ilink;
  register site *s;
  int length;
  su3_matrix tmat,tmat2;
  Real ferm_epsilon1, ferm_epsilon2, coeff1, coeff2;
#ifdef FFTIME
  int nflop = 0;
  double dtime;
#endif
  msg_tag *mtag0;
  wilson_vector *w_tmp0,*w_tmp1,*tmp_pt;
  int num_q_paths = ap->p.num_q_paths;
  Q_path *q_paths = ap->p.q_paths;

#ifdef FFTIME
dtime=-dclock();
#endif
  ferm_epsilon1 = 2.0*weight1*eps;
  ferm_epsilon2 = 2.0*weight2*eps;
  w_tmp0 = (wilson_vector *)
    malloc(sites_on_node*sizeof(wilson_vector) );
  w_tmp1 = (wilson_vector *)
    malloc(sites_on_node*sizeof(wilson_vector) );
  /* Use wilson_vectors to store x_off transported from ends of
     path.  0 and 1 components from forward end, 2 and 3 components
     from back end */

  /* loop over paths, and loop over links in path */
  for( ipath=0; ipath<num_q_paths; ipath++ ){
    if(q_paths[ipath].forwback== -1)continue;	/* skip backwards dslash */
    length = q_paths[ipath].length;

    /* path transport x_off and Dslash*x_off from far end.  Sometimes
	we need them at the start point of the path, and sometimes
	one link into the path - an optimization for later */
    /**
    path_transport_site( x1_off, F_OFFSET(tempvec[0]),
      EVENANDODD, q_paths[ipath].dir, length );
    path_transport_site( x2_off, F_OFFSET(tempvec[1]),
      EVENANDODD, q_paths[ipath].dir, length );
    **/
/** WARNING!! Assumes xxx1 and xxx2 contiguous **/
if( x2_off-x1_off != sizeof(su3_vector) ){node0_printf("BOTCH\n"); exit(0);}
    path_transport_hwv_site( x1_off, F_OFFSET(tempvec[0]),
      EVENANDODD, q_paths[ipath].dir, length );
    /* use tempvec[2] for transport from starting end */
    FORALLSITES(i,s){
      w_tmp0[i].d[0]=s->tempvec[0];
      w_tmp0[i].d[1]=s->tempvec[1];
      w_tmp0[i].d[2]=*(su3_vector *)F_PT(s,x1_off);
      w_tmp0[i].d[3]=*(su3_vector *)F_PT(s,x2_off);
    }

    /* A path has (length+1) points, counting the ends.  At first
	 point, no "down" direction links have their momenta "at this
	 point". At last, no "up" ... */
    for( ilink=0; ilink<=length; ilink++ ){
      if(ilink<length)dir = q_paths[ipath].dir[ilink];
      else dir=NODIR;
      coeff1 = ferm_epsilon1*q_paths[ipath].coeff;
      coeff2 = ferm_epsilon2*q_paths[ipath].coeff;
      if( (ilink%2)==1 ){ coeff1 = -coeff1; coeff2 = -coeff2;}

      /* path transport x_off and Dslash*x_off from previous point */
      /* Use "wilson_vector" to handle pair of vectors -
	0 component is x_off1 from forward end, 1 component is x_off2
	from forward end,  2 and 3  components are x_off1 and x_off2 
	from back end */
      /* sometimes we don't need them */
      if( (ilink>0&&ilink<length) || 
	(ilink==length && GOES_BACKWARDS(lastdir)) ){

	if( GOES_FORWARDS(lastdir) ){
	  FORALLSITES(i,s){
            mult_adj_mat_wilson_vec( &(s->link[lastdir]),
	      &(w_tmp0[i]), &(w_tmp1[i]) );
	  }
	  mtag0 = start_gather_field( w_tmp1, 4*sizeof(su3_vector),
             OPP_DIR(lastdir), EVENANDODD, gen_pt[0] );
          wait_gather(mtag0);
          FORALLSITES(i,s){
	     w_tmp0[i] = *(wilson_vector *)gen_pt[0][i];
	  }
          cleanup_gather(mtag0);
	}
	else{   /* GOES_BACKWARDS(lastdir) */
          mtag0 = start_gather_field( w_tmp0, 4*sizeof(su3_vector),
               OPP_DIR(lastdir), EVENANDODD, gen_pt[0] );
          wait_gather(mtag0);
          FORALLSITES(i,s){
            mult_mat_wilson_vec( &(s->link[OPP_DIR(lastdir)]),
                (wilson_vector *)(gen_pt[0][i]),
		&(w_tmp1[i]) );
          }
	  tmp_pt = w_tmp0; w_tmp0 = w_tmp1; w_tmp1 = tmp_pt;
          cleanup_gather(mtag0);
	}
      }

      /* add in contribution to the force */
      /* Put antihermitian traceless part into momentum */
      FORALLSITES(i,s){
        if( ilink<length && GOES_FORWARDS(dir) ){
	  if( s->parity==ODD ){coeff1 *= -1.0; coeff2 *= -1.0; }
          uncompress_anti_hermitian( &(s->mom[dir]), &tmat2 );
	  su3_projector( &(w_tmp0[i].d[0]), &(w_tmp0[i].d[2]), &tmat );
	  scalar_mult_add_su3_matrix( &tmat2, &tmat,  coeff1, &tmat2 );
	  su3_projector( &(w_tmp0[i].d[1]), &(w_tmp0[i].d[3]), &tmat );
	  scalar_mult_add_su3_matrix( &tmat2, &tmat,  coeff2, &tmat2 );
          make_anti_hermitian( &tmat2, &(s->mom[dir]) );
	  if( s->parity==ODD ){coeff1 *= -1.0; coeff2 *= -1.0; }
        }
	if( ilink>0 && GOES_BACKWARDS(lastdir) ){
	  if( s->parity==EVEN ){coeff1 *= -1.0; coeff2 *= -1.0; }
          uncompress_anti_hermitian( &(s->mom[OPP_DIR(lastdir)]), &tmat2 );
	  su3_projector( &(w_tmp0[i].d[0]), &(w_tmp0[i].d[2]), &tmat );
	  scalar_mult_add_su3_matrix( &tmat2, &tmat,  coeff1, &tmat2 );
	  su3_projector( &(w_tmp0[i].d[1]), &(w_tmp0[i].d[3]), &tmat );
	  scalar_mult_add_su3_matrix( &tmat2, &tmat,  coeff2, &tmat2 );
          make_anti_hermitian( &tmat2, &(s->mom[OPP_DIR(lastdir)]) );
	  if( s->parity==EVEN ){coeff1 *= -1.0; coeff2 *= -1.0; }
        }
      }
      lastdir = dir;
    } /* end loop over links in path */
  } /* end loop over paths */

  free( w_tmp0 ); free( w_tmp1 );

#ifdef FFTIME
dtime += dclock();
node0_printf("FFTIME:  time = %e (2 mass) mflops = %e\n",dtime,
	     (Real)nflop*volume/(1e6*dtime*numnodes()) );
/**printf("TLENGTH: %d\n",tlength);**/
#endif
} /* eo_fermion_force_twoterms (version 7) */

/* LONG COMMENTS
   Here we have combined "xxx", (offset "x_off")  which is
(M_adjoint M)^{-1} phi, with Dslash times this vector, which goes in the
odd sites of xxx.  Recall that phi is defined only on even sites.  In
computing the fermion force, we are looking at

< X |  d/dt ( Dslash_eo Dslash_oe ) | X >
=
< X | d/dt Dslash_eo | T > + < T | d/dt Dslash_oe | X >
where T = Dslash X.

The subsequent manipulations to get the coefficent of H, the momentum
matrix, in the simulation time derivative above look the same for
the two terms, except for a minus sign at the end, if we simply stick
T, which lives on odd sites, into the odd sites of X

 Each path in the action contributes terms when any link of the path
is the link for which we are computing the force.  We get a minus sign
for odd numbered links in the path, since they connect sites of the
opposite parity from what it would be for an even numbered link.
Minus signs from "going around" plaquette - ie KS phases, are supposed
to be already encoded in the path coefficients.
Minus signs from paths that go backwards are supposed to be already
encoded in the path coefficients.

Here, for example, are comments reproduced from the force routine for
the one-link plus Naik plus single-staple-fat-link action:

 The three link force has three contributions, where the link that
was differentiated is the first, second, or third link in the 3-link
path, respectively.  Diagramatically, where "O" represents the momentum,
the solid line the link corresponding to the momentum, and the dashed
lines the other links:
 

	O______________ x ............ x ...............
+
	x..............O______________x.................
+
	x..............x..............O________________
Think of this as
	< xxx | O | UUUxxx >		(  xxx, UUUX_p3 )
+
	< xxx U | O | UUxxx >		( X_m1U , UUX_p2 )
+
	< xxx U U | O | Uxxx >		( X_m2UU , UX_p1 )
where "U" indicates parallel transport, "X_p3" is xxx displaced
by +3, etc.
Note the second contribution has a relative minus sign
because it effectively contributes to the <odd|even>, or M_adjoint,
part of the force when we work on an even site. i.e., for M on
an even site, this three link path begins on an odd site.

The staple force has six contributions from each plane containing the
link direction:
Call these diagrams A-F:


	x...........x		O____________x
		    .			     .
		    .			     .
		    .			     .
		    .			     .
		    .			     .
	O___________x		x............x
	   (A)			    (B)



	x	    x		O____________x
	.	    .		.	     .
	.	    .		.	     .
	.	    .		.	     .
	.	    .		.	     .
	.	    .		.	     .
	O___________x		x	     x
	   (C)			    (D)



	x...........x		O____________x
	.			.
	.			.
	.			.
	.			.
	.			.
	O___________x		x............x
	   (E)			    (F)

As with the Naik term, diagrams C and D have a relative minus
sign because they connect sites of the other parity.

Also note an overall minus sign in the staple terms relative to the
one link term because, with the KS phase factors included, the fat
link is  "U - w3 * UUU", or the straight link MINUS w3 times the staples.

Finally, diagrams B and E get one more minus sign because the link
we are differentiating is in the opposite direction from the staple
as a whole.  You can think of this as this "U" being a correction to
a "U_adjoint", but the derivative of U is iHU and the derivative
of U_adjoint is -iHU_adjoint.

*/
/* LONG COMMENT on sign conventions
In most of the program, the KS phases and antiperiodic boundary
conditions are absorbed into the link matrices.  This greatly simplfies
multiplying by the fermion matrix.  However, it requires care in
specifying the path coefficients.  Remember that each time you
encircle a plaquette, you pick up a net minus sign from the KS phases.
Thus, when you have more than one path to the same point, you generally
have a relative minus sign for each plaquette in a surface bounded by
this path and the basic path for that displacement.

Examples:
  Fat Link:
    Positive:	X-------X

    Negative     --------
	 	|	|
		|	|
		X	X

  Naik connection, smeared
    Positive:	X-------x-------x-------X

    Negative:	---------
		|	|
		|	|
		X	x-------x-------X

    Positive:	--------x--------
		|		|
		|		|
		X		x-------X

    Negative:	--------x-------x-------x
		|			|
		|			|
		X			X
*/



/* Comment on acceptable actions.
   We construct the backwards part of dslash by reversing all the
   paths in the forwards part.  So, for example, in the p4 action
   the forwards part includes +X+Y+Y

		X
		|
		|
		X
		|
		|
	X---->--X

  so we put -X-Y-Y in the backwards part.  But this isn't the adjoint
  of U_x(0)U_y(+x)U_y(+x+y).  Since much of the code assumes that the
  backwards hop is the adjoint of the forwards (for example, in
  preventing going to 8 flavors), the code only works for actions
  where this is true.  Roughly, this means that the fat link must
  be symmetric about reflection around its midpoint.  Equivalently,
  the paths in the backwards part of Dslash are translations of the
  paths in the forwards part.  In the case of the "P4" or knight's move
  action, this means that we have to have both paths
   +X+Y+Y and +Y+Y+X to the same point, with the same coefficients.
  Alternatively, we could just use the symmetric path +Y+X+Y.
*/  /* fermion_forcee_eo_milc.c */
