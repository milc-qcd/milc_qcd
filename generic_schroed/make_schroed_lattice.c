/******** make_schroed_lattice.c *********/

#include "generic_schroed_includes.h"

/********* set_boundary_fields() - set up boundary fields **********/
void set_boundary_fields()  {
/* Sets the link matrices at the boundaries t=0 and t=nt.
   The latter are stored in boundary at t=nt-1.
   Also computes the derivatives, to be stored in boundary at t=0,
   and for the other boundary at t=nt-2. */

register int i,j,k,dir;
register site *sit;
Real fact[3],b0r[3],b0i[3],btr[3],bti[3];
Real angle,c_t;

    c_t = 1.0 - (0.534+c_t11)/beta;

    if(bc_flag > 0){
	if(nx!=ny || nx!=nz){
	     printf("Boundary conditions require equal size spatial lattice!\n");
	     terminate(1);
	}
    }

    switch(bc_flag){
	case 0:
	    for(j=0; j<3; j++){
		fact[j] = 0.0;
		b0r[j] = c_t;
		b0i[j] = 0.0;
		btr[j] = c_t;
		bti[j] = 0.0;
	    }
	    break;

	case 1:
	    fact[2] = -1.0/(Real)(2*nx);
	    fact[1] = fact[2];
	    fact[0] = -2.0*fact[2];
	    angle = PI/(Real)(3*nx);
            b0r[2] = c_t*cos((double)angle);
            b0i[2] = c_t*sin((double)angle);
	    b0r[1] = c_t;
	    b0i[1] = 0.0;
	    b0r[0] = b0r[2];
            b0i[0] = -b0i[2];
            btr[2] = c_t*cos((double)(2.0*angle));
            bti[2] = c_t*sin((double)(2.0*angle));
            btr[1] = b0r[2];
            bti[1] = b0i[2];
            btr[0] = c_t*cos((double)(3.0*angle));
            bti[0] = -c_t*sin((double)(3.0*angle));
	    break;

	case 2:
	    fact[2] = -1.0/(Real)(2*nx);
	    fact[1] = fact[2];
	    fact[0] = -2.0*fact[2];
	    angle = PI/(Real)(6*nx);
            b0r[2] = c_t*cos((double)angle);
            b0i[2] = c_t*sin((double)angle);
	    b0r[1] = c_t;
	    b0i[1] = 0.0;
	    b0r[0] = b0r[2];
            b0i[0] = -b0i[2];
            btr[2] = c_t*cos((double)(3.0*angle));
            bti[2] = c_t*sin((double)(3.0*angle));
            btr[1] = c_t*cos((double)(2.0*angle));
            bti[1] = c_t*sin((double)(2.0*angle));
            btr[0] = c_t*cos((double)(5.0*angle));
            bti[0] = -c_t*sin((double)(5.0*angle));
	    break;

	default:
	    if(this_node==0)printf("Bad bc_flag %d\n",bc_flag);
	    terminate(1);
    }

    FORALLSITES(i,sit){
	if( sit->t == 0 ){
	    for(dir=XUP;dir<=ZUP;dir++){
		for(j=0; j<3; j++) for(k=0; k<3; k++)  {
		    if (j != k)  {
			sit->link[dir].e[j][k] = cmplx(0.0,0.0);
			sit->boundary[dir].e[j][k] = cmplx(0.0,0.0);
		    }
		    else  {
			sit->link[dir].e[j][k].real = b0r[j];
			sit->link[dir].e[j][k].imag = b0i[j];
			sit->boundary[dir].e[j][k].real = -fact[j]*b0i[j];
			sit->boundary[dir].e[j][k].imag = fact[j]*b0r[j];
		    }
		}
	    }
	}
	else if( sit->t == (nt-1) ){
	    for(dir=XUP;dir<=ZUP;dir++){
		for(j=0; j<3; j++) for(k=0; k<3; k++)  {
		    if (j != k)  {
			sit->boundary[dir].e[j][k] = cmplx(0.0,0.0);
		    }
		    else  {
			sit->boundary[dir].e[j][k].real = btr[j];
			sit->boundary[dir].e[j][k].imag = bti[j];
		    }
		}
	    }
	}
	else if( sit->t == (nt-2) ){
	    for(dir=XUP;dir<=ZUP;dir++){
		for(j=0; j<3; j++) for(k=0; k<3; k++)  {
		    if (j != k)  {
			sit->boundary[dir].e[j][k] = cmplx(0.0,0.0);
		    }
		    else  {
			sit->boundary[dir].e[j][k].real = fact[j]*bti[j];
			sit->boundary[dir].e[j][k].imag = -fact[j]*btr[j];
		    }
		}
	    }
	}
	else{
	    for(dir=XUP;dir<=ZUP;dir++){
		for(j=0; j<3; j++) for(k=0; k<3; k++)  {
		    sit->boundary[dir].e[j][k] = cmplx(0.0,0.0);
		}
	    }
	}
    }
    if(this_node==0)printf("Boundary values set\n");
}
