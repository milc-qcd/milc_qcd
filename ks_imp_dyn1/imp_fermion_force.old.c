/***************** imp_ferion_force.old.c *****************************/
/* NOT MAINTAINED.  TEST BEFORE USE! */
/* MIMD version 6 */
/* update the  momenta with the fermion force */
/* Assumes that the conjugate gradient has been run, with the answer in
   xxx, and dslash(xxx,xxx,ODD) has been run. (fills in xxx_odd) */
/* SEE LONG COMMENT AT END */
/* OLD VERSION FOR OneFatNaik action!!! */
#ifdef FN	/* Actually, requirements much stricter!! */


void imp_fermion_force( Real eps ){
  /* note CG_solution and Dslash * solution are combined in "xxx" */
  /* see long comment at end */
  register int i, a,b;
  register site *s;
  register int dir,otherparity;
  msg_tag *tag[4];
  su3_vector tvec;
  su3_matrix tmat1, tmat2, tmat3;
  Real ferm_epsilon;
  double dtime, dclock();
  int dir2;

#ifdef FFTIME
dtime=-dclock();
#endif
  ferm_epsilon = (nflavors/2.0)*eps;

  for(dir=XUP; dir<=TUP; dir++) {

/* SINGLE LINK PARTS:  NEED UX_p1 and X_m1U */
    /* gather xxx and Dslash*xxx from + direction */
    tag[0] = start_gather( F_OFFSET(xxx), sizeof(su3_vector), dir, 
	EVENANDODD, gen_pt[0] );
    wait_gather(tag[0]);

    FORALLSITES(i,s) {
	/* finish parallel transorting from +1 */
	mult_su3_mat_vec(&(s->link[dir]),
	   (su3_vector *)gen_pt[0][i], &(s->UX_p1) );
	/* begin parallel transport from minus direction */
        mult_adj_su3_mat_vec(&(s->link[dir]), &(s->xxx), &(s->tempvec[0]) );
    }
    cleanup_gather(tag[0]);

/* FINISHED MAKING VECTORS FOR SINGLE LINK PART */
/* NOW ADD PIECES OF ACCELERATION FOR SINGLE LINK */
    FORALLSITES(i,s) {
      uncompress_anti_hermitian( &(s->mom[dir]), &tmat2);

      /* one link term in force */
      su3_projector( &(s->UX_p1), &(s->xxx), &tmat1 );
/**
if( dir==1 && i==node_index(1,0,0,0) ){
printf("IMP_FORCE TERM:  %d OneLink  0  %e\n",dir, c1*ferm_epsilon);
dumpvec( &(s->UX_p1) );
dumpvec( &(s->xxx) );
dumpmat( &tmat1 );
}
**/
      if(s->parity==EVEN) scalar_mult_add_su3_matrix(&tmat2,
	       &tmat1, c1*ferm_epsilon, &tmat2 );
      else 		  scalar_mult_add_su3_matrix(&tmat2,
	         &tmat1, -c1*ferm_epsilon, &tmat2 );

      make_anti_hermitian(&tmat2, &(s->mom[dir]) );
    }
/* FINISHED SINGLE LINK FORCE */

/* NAIK PARTS: UUX_p2, UUUX_p3, X_m1U, X_m2UU  */
    /* parallel transport xxx and dslash*xxx from minus */
    tag[1] = start_gather( F_OFFSET(tempvec[0]), sizeof(su3_vector),
	OPP_DIR(dir), EVENANDODD, gen_pt[1] );
    wait_gather(tag[1]);
    FORALLSITES(i,s) { 
	/* put result of p.t. from minus direction somewhere we can
	   gather it again */
        s->X_m1U = *(su3_vector *)gen_pt[1][i];
    }
    cleanup_gather(tag[1]);

    /* start parallel transport from +, for distance two */
    tag[0] = start_gather( F_OFFSET(UX_p1), sizeof(su3_vector), dir, 
			     EVENANDODD, gen_pt[0] );
    wait_gather(tag[0]);
    FORALLSITES(i,s) { 
	mult_su3_mat_vec(&(s->link[dir]),
	    (su3_vector *)gen_pt[0][i], &(s->UUX_p2) ); 
        mult_adj_su3_mat_vec(&(s->link[dir]), &(s->X_m1U), &(s->tempvec[0]) );
    }
    cleanup_gather(tag[0]);

    tag[2] = start_gather( F_OFFSET(tempvec[0]), sizeof(su3_vector),
	OPP_DIR(dir), EVENANDODD, gen_pt[2] );
    tag[0] = start_gather( F_OFFSET(UUX_p2), sizeof(su3_vector), dir, 
			     EVENANDODD, gen_pt[0] );
    wait_gather(tag[0]);
    wait_gather(tag[2]);
    FORALLSITES(i,s) {
	mult_su3_mat_vec(&(s->link[dir]),
	    (su3_vector *)gen_pt[0][i], &(s->UUUX_p3) );
        s->X_m2UU = *(su3_vector *)gen_pt[2][i]; 
    }
    cleanup_gather(tag[0]);
    cleanup_gather(tag[2]);
/* FINISHED MAKING VECTORS FOR NAIK PART */
/* NOW ADD PIECES OF ACCELERATION FOR NAIK PART */
    FORALLSITES(i,s) {
      uncompress_anti_hermitian( &(s->mom[dir]), &tmat2);

      /* three link (Naik) term, "our" link is first link */
      su3_projector( &(s->UUUX_p3), &(s->xxx), &tmat1 );
/**
if( dir==1 && i==node_index(1,0,0,0) ){
printf("IMP_FORCE TERM:  %d Naik  0  %e\n",dir, c3*ferm_epsilon);
dumpvec( &(s->UUUX_p3) );
dumpvec( &(s->xxx) );
dumpmat( &tmat1 );
}
**/
      /* three link (Naik) term, "our" link is second link */
      /* note different sign */
      su3_projector( &(s->UUX_p2), &(s->X_m1U), &tmat3 );
/**
if( dir==1 && i==node_index(1,0,0,0) ){
printf("IMP_FORCE TERM:  %d Naik  1  %e\n",dir, -c3*ferm_epsilon);
dumpvec( &(s->UUX_p2) );
dumpvec( &(s->X_m1U) );
dumpmat( &tmat3 );
}
**/
      sub_su3_matrix( &tmat1, &tmat3, &tmat1 );
      /* three link (Naik) term, "our" link is third link */
      su3_projector( &(s->UX_p1), &(s->X_m2UU), &tmat3 );
/**
if( dir==1 && i==node_index(1,0,0,0) ){
printf("IMP_FORCE TERM:  %d Naik  1  %e\n",dir, c3*ferm_epsilon);
dumpvec( &(s->UX_p1) );
dumpvec( &(s->X_m2UU) );
dumpmat( &tmat3 );
}
**/
      add_su3_matrix( &tmat1, &tmat3, &tmat1 );
      if(s->parity==EVEN)scalar_mult_add_su3_matrix(&tmat2,
	   &tmat1, c3*ferm_epsilon, &tmat2 );
      else		 scalar_mult_add_su3_matrix(&tmat2,
	   &tmat1, -c3*ferm_epsilon, &tmat2 );

      make_anti_hermitian(&tmat2, &(s->mom[dir]) );
    }
/* FINISHED NAIK PART OF FORCE */

/* NOW DO THE STAPLE PART */
    /* use tempmat1 to accumulate acceleration */
    FORALLSITES(i,s) {
      clear_su3mat( &(s->tempmat1) );
    }
    for( dir2=XUP; dir2<=TUP; dir2++ )if( dir2 != dir ){
	
	/* Diagram "A" (see long comments at end) */
        /* parallel transport from backwards "dir" direction */
        FORALLSITES(i,s) { /*  multiply by U_adj */
	    mult_adj_su3_mat_vec(&(s->link[dir]), &(s->xxx), &(s->tempvec[0]) );
        }
        tag[0] = start_gather( F_OFFSET(tempvec[0]), sizeof(su3_vector),
	    OPP_DIR(dir), EVENANDODD, gen_pt[0] );
        wait_gather(tag[0]);

	/* put result into place we can gather it, and parallel transport
	   from forwards "dir2" direction */
        FORALLSITES(i,s) {
	    s->tempvec[1] = *(su3_vector *)gen_pt[0][i];
	}
        cleanup_gather(tag[0]);

	tag[1] = start_gather( F_OFFSET(tempvec[1]), sizeof(su3_vector),
            dir2, EVENANDODD, gen_pt[1] );
	wait_gather(tag[1]);
        FORALLSITES(i,s) { /*  multiply by U */
	    mult_su3_mat_vec(&(s->link[dir2]), (su3_vector *)gen_pt[1][i],
		 &(s->tempvec[2]) );
        }
        cleanup_gather(tag[1]);
	/*now parallel transport from the +dir direction */
	tag[2] = start_gather( F_OFFSET(tempvec[2]), sizeof(su3_vector),
            dir, EVENANDODD, gen_pt[2] );
	wait_gather(tag[2]);
        FORALLSITES(i,s) { /*  multiply by U */
	    mult_su3_mat_vec(&(s->link[dir]), (su3_vector *)gen_pt[2][i],
		 &(s->tempvec[3]) );
        }
        cleanup_gather(tag[2]);
	/* and add to momentum */
	FORALLSITES(i,s){
            su3_projector( &(s->tempvec[3]), &(s->xxx), &tmat1 );
/**
if( dir==1 && i==node_index(1,0,0,0) ){
printf("IMP_FORCE TERM:  %d STAPLE  A  %e\n",dir, -w3*ferm_epsilon);
dumpvec( &(s->tempvec[3]) );
dumpvec( &(s->xxx) );
dumpmat( &tmat1 );
}
**/
            if(s->parity==EVEN)scalar_mult_add_su3_matrix(&(s->tempmat1),
	          &tmat1, -w3*ferm_epsilon, &(s->tempmat1) );
            else		 scalar_mult_add_su3_matrix(&(s->tempmat1),
	          &tmat1, +w3*ferm_epsilon, &(s->tempmat1) );
	}

	/* Diagram "B" */
        /* parallel transport from backwards "dir" direction */
        FORALLSITES(i,s) { /*  multiply by U_adj */
	    mult_adj_su3_mat_vec(&(s->link[dir]), &(s->xxx), &(s->tempvec[0]) );
        }
        tag[0] = start_gather( F_OFFSET(tempvec[0]), sizeof(su3_vector),
	    OPP_DIR(dir), EVENANDODD, gen_pt[0] );
        wait_gather(tag[0]);

	/*  parallel transport from backwards "dir2" direction */
        FORALLSITES(i,s) { /*  multiply by U_adj */
	    mult_adj_su3_mat_vec(&(s->link[dir2]), (su3_vector *)gen_pt[0][i],
		 &(s->tempvec[1]) );
        }
        cleanup_gather(tag[0]);
	tag[1] = start_gather( F_OFFSET(tempvec[1]), sizeof(su3_vector),
            OPP_DIR(dir2), EVENANDODD, gen_pt[1] );
	wait_gather(tag[1]);
        FORALLSITES(i,s) {
	    s->tempvec[2] = *(su3_vector *)gen_pt[1][i];
	}
        cleanup_gather(tag[1]);

	/*now parallel transport from the +dir direction */
	tag[2] = start_gather( F_OFFSET(tempvec[2]), sizeof(su3_vector),
            dir, EVENANDODD, gen_pt[2] );
	wait_gather(tag[2]);
        FORALLSITES(i,s) { /*  multiply by U */
	    mult_su3_mat_vec(&(s->link[dir]), (su3_vector *)gen_pt[2][i],
		 &(s->tempvec[3]) );
        }
        cleanup_gather(tag[2]);
	/* and add to momentum */
	FORALLSITES(i,s){
            su3_projector( &(s->tempvec[3]), &(s->xxx), &tmat1 );
/**
if( dir==1 && i==node_index(1,0,0,0) ){
printf("IMP_FORCE TERM:  %d STAPLE  B  %e\n",dir, w3*ferm_epsilon);
dumpvec( &(s->tempvec[3]) );
dumpvec( &(s->xxx) );
dumpmat( &tmat1 );
}
**/
            if(s->parity==EVEN)scalar_mult_add_su3_matrix(&(s->tempmat1),
	          &tmat1,  w3*ferm_epsilon, &(s->tempmat1) );
            else		 scalar_mult_add_su3_matrix(&(s->tempmat1),
	          &tmat1, -w3*ferm_epsilon, &(s->tempmat1) );
	}

	/* Diagram "C" */
	/* gather from +dir2 direction */
        tag[0] = start_gather( F_OFFSET(xxx), sizeof(su3_vector),
            dir2, EVENANDODD, gen_pt[0] );
        wait_gather(tag[0]);
        FORALLSITES(i,s) { /*  multiply by U */
            mult_su3_mat_vec(&(s->link[dir2]), (su3_vector *)gen_pt[0][i],
                 &(s->tempvec[0]) );
        }
	/* gather again from +dir direction */
        tag[1] = start_gather( F_OFFSET(tempvec[0]), sizeof(su3_vector),
            dir, EVENANDODD, gen_pt[1] );
        wait_gather(tag[1]);
        FORALLSITES(i,s) { /*  multiply by U */
            mult_su3_mat_vec(&(s->link[dir]), (su3_vector *)gen_pt[1][i],
                 &(s->tempvec[1]) );
        }
        /* and add to momentum */
        FORALLSITES(i,s){
            su3_projector( &(s->tempvec[1]), &(s->tempvec[0]), &tmat1 );
/**
if( dir==1 && i==node_index(1,0,0,0) ){
printf("IMP_FORCE TERM:  %d STAPLE  C  %e\n",dir, w3*ferm_epsilon);
dumpvec( &(s->tempvec[1]) );
dumpvec( &(s->tempvec[0]) );
dumpmat( &tmat1 );
}
**/
            if(s->parity==EVEN)scalar_mult_add_su3_matrix(&(s->tempmat1),
                  &tmat1, +w3*ferm_epsilon, &(s->tempmat1) );
            else                 scalar_mult_add_su3_matrix(&(s->tempmat1),
                  &tmat1, -w3*ferm_epsilon, &(s->tempmat1) );
        }
        cleanup_gather(tag[0]);
        cleanup_gather(tag[1]);

	/* Diagram "D" */
	/* gather from -dir2 direction */
        FORALLSITES(i,s) { /*  multiply by U_adj */
            mult_adj_su3_mat_vec(&(s->link[dir2]), &(s->xxx),
		&(s->tempvec[0]) );
        }
        tag[0] = start_gather( F_OFFSET(tempvec[0]), sizeof(su3_vector),
            OPP_DIR(dir2), EVENANDODD, gen_pt[0] );
        wait_gather(tag[0]);

        /* put result into place we can gather it, and parallel transport
           from forwards "dir" direction */
        FORALLSITES(i,s) {
            s->tempvec[1] = *(su3_vector *)gen_pt[0][i];
        }
        cleanup_gather(tag[0]);
	/* gather again from +dir direction */
        tag[1] = start_gather( F_OFFSET(tempvec[1]), sizeof(su3_vector),
            dir, EVENANDODD, gen_pt[1] );
        wait_gather(tag[1]);
        FORALLSITES(i,s) { /*  multiply by U */
            mult_su3_mat_vec(&(s->link[dir]), (su3_vector *)gen_pt[1][i],
                 &(s->tempvec[2]) );
        }
        /* and add to momentum */
        FORALLSITES(i,s){
            su3_projector( &(s->tempvec[2]), &(s->tempvec[1]), &tmat1 );
/**
if( dir==1 && i==node_index(1,0,0,0) ){
printf("IMP_FORCE TERM:  %d STAPLE  D  %e\n",dir, w3*ferm_epsilon);
dumpvec( &(s->tempvec[2]) );
dumpvec( &(s->tempvec[1]) );
dumpmat( &tmat1 );
}
**/
            if(s->parity==EVEN)scalar_mult_add_su3_matrix(&(s->tempmat1),
                  &tmat1,  w3*ferm_epsilon, &(s->tempmat1) );
            else                 scalar_mult_add_su3_matrix(&(s->tempmat1),
                  &tmat1, -w3*ferm_epsilon, &(s->tempmat1) );
        }
        cleanup_gather(tag[1]);

	/* Diagram "E" */
	/* gather from +dir direction */
        tag[0] = start_gather( F_OFFSET(xxx), sizeof(su3_vector),
            dir, EVENANDODD, gen_pt[0] );
        wait_gather(tag[0]);
        FORALLSITES(i,s) { /*  multiply by U */
            mult_su3_mat_vec(&(s->link[dir]), (su3_vector *)gen_pt[0][i],
                 &(s->tempvec[0]) );
        }
	/* gather again from +dir2 direction */
        tag[1] = start_gather( F_OFFSET(tempvec[0]), sizeof(su3_vector),
            dir2, EVENANDODD, gen_pt[1] );
        wait_gather(tag[1]);
        FORALLSITES(i,s) { /*  multiply by U */
            mult_su3_mat_vec(&(s->link[dir2]), (su3_vector *)gen_pt[1][i],
                 &(s->tempvec[1]) );
        }
        /* and add to momentum */
        FORALLSITES(i,s){
            su3_projector( &(s->tempvec[0]), &(s->tempvec[1]), &tmat1 );
/**
if( dir==1 && i==node_index(1,0,0,0) ){
printf("IMP_FORCE TERM:  %d STAPLE  E  %e\n",dir, w3*ferm_epsilon);
dumpvec( &(s->tempvec[0]) );
dumpvec( &(s->tempvec[1]) );
dumpmat( &tmat1 );
}
**/
            if(s->parity==EVEN)scalar_mult_add_su3_matrix(&(s->tempmat1),
                  &tmat1,  w3*ferm_epsilon, &(s->tempmat1) );
            else                 scalar_mult_add_su3_matrix(&(s->tempmat1),
                  &tmat1, -w3*ferm_epsilon, &(s->tempmat1) );
        }
        cleanup_gather(tag[0]);
        cleanup_gather(tag[1]);

	/* Diagram "F" */
	/* gather from +dir direction */
        tag[0] = start_gather( F_OFFSET(xxx), sizeof(su3_vector),
            dir, EVENANDODD, gen_pt[0] );
        wait_gather(tag[0]);
        FORALLSITES(i,s) { /*  multiply by U */
            mult_su3_mat_vec(&(s->link[dir]), (su3_vector *)gen_pt[0][i],
                 &(s->tempvec[0]) );
        }
	/* gather again from -dir2 direction */
        FORALLSITES(i,s) { /*  multiply by U_adj */
            mult_adj_su3_mat_vec(&(s->link[dir2]), &(s->tempvec[0]),
		&(s->tempvec[1]) );
        }
        tag[1] = start_gather( F_OFFSET(tempvec[1]), sizeof(su3_vector),
            OPP_DIR(dir2), EVENANDODD, gen_pt[1] );
        wait_gather(tag[1]);
        /* and add to momentum */
        FORALLSITES(i,s){
            su3_projector( &(s->tempvec[0]),
		(su3_vector *)gen_pt[1][i], &tmat1 );
/**
if( dir==1 && i==node_index(1,0,0,0) ){
printf("IMP_FORCE TERM:  %d STAPLE  F  %e\n",dir, -w3*ferm_epsilon);
dumpvec( &(s->tempvec[0]) );
dumpvec( (su3_vector *)gen_pt[1][i] );
dumpmat( &tmat1 );
}
**/
            if(s->parity==EVEN)scalar_mult_add_su3_matrix(&(s->tempmat1),
                  &tmat1, -w3*ferm_epsilon, &(s->tempmat1) );
            else                 scalar_mult_add_su3_matrix(&(s->tempmat1),
                  &tmat1,  w3*ferm_epsilon, &(s->tempmat1) );
        }
        cleanup_gather(tag[0]);
        cleanup_gather(tag[1]);

    }  /* end dir2 loop */
    FORALLSITES(i,s) {
      uncompress_anti_hermitian( &(s->mom[dir]), &tmat2);
      add_su3_matrix( &tmat2, &(s->tempmat1), &tmat2 );
      make_anti_hermitian(&tmat2, &(s->mom[dir]) );
/**
if( dir==1 && i==node_index(1,0,0,0) ){
anti_hermitmat tmat2;
printf("IMP_FORCE TOTAL:  %d\n",dir);
printf("%e  %e  %e\n",s->mom[dir].m00im,s->mom[dir].m11im,s->mom[dir].m22im);
printf("%e  %e  %e  %e  %e %e\n",s->mom[dir].m01.real,s->mom[dir].m01.imag,s->mom[dir].m02.real,s->mom[dir].m02.imag,s->mom[dir].m12.real,s->mom[dir].m12.imag);
}
**/

    }
/* DONE WITH STAPLE FORCE */

  } /* end dir loop */
#ifdef FFTIME
dtime += dclock();
node0_printf("FFTIME:  %e\n",dtime);
#endif
} /* imp_fermion_force */
#endif


/* LONG COMMENTS
   Here we have combined "xxx", which is  (M_adjoint M)^{-1} phi,
with Dslash times this vector, which goes in the odd sites of xxx.
Recall that phi is defined only on even sites.  In computing the fermion
force, we are looking at

< X |  d/dt ( Dslash_eo Dslash_oe ) | X >
=
< X | d/dt Dslash_eo | T > + < T | d/dt Dslash_oe | X >
where T = Dslash X.

The subsequent manipulations to get the coefficent of H, the momentum
matrix, in the simulation time derivative above look the same for
the two terms, except for a minus sign at the end, if we simply stick
T, which lives on odd sites, into the odd sites of X


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
