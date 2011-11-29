/******** wavefunc_t.c *************/
/* MIMD version5 */
/* NOT MAINTAINED.  TEST BEFORE USE! */

/* Compute wavefunctions with staggered quarks.  Fourier transform the
   2x2x2 cubes, do convolution, transform back, then construct hadrons.
   The "block Fourier transform" leaves the lowest order bit alone.
   This means the bit reverse is different than usual, and all the
   butterflies connect twice as far as usual (level 0 connects to second
   nearest neighbors, etc.)

   References:
	Fast Fourier Transform and Convolution Algorithms, H.J. Nussbaumer
	(Springer series in information sciences, 1982)  QA 403.5
*/
/* this is a  version which includes the baryon--we pin two
quarks on the same site and put the third quark everywhere. Thus we need
only one set of fft's. Other components of the baryon where two quarks
are offset with some other separation could be done with the
general_gather routine. This is more wasteful of time compared to
fft's if many relative separations were wanted, but if we just do a few,
it is easier to code.

Meson code is by D. T. */
#include "ks_dyn_includes.h"

/* Variables set by setup_block_fourier() and used by block_fourier() */
int dim[4];		/* dimensions */
int logdim[4];		/* log_2 of dimension */
int bitrev_dir;			/* index for bit reverse gather */
int bitrev_block_dir;		/* index for bit reverse gather */
int *butterfly_dir[4];	/* indices for butterflies.  First index is
			   direction, second is level.  The actual space
			   will be malloc'ed in setup_fourier(). */
/* Gathers for symmetry transforms */
int switch_dir,permute_dir,minus_dir[4],minus_block_dir[4];
Real sym_factor;	/* normalization factor for symmetries */



/* Hadron wave functions. */
void wavefunc_t() {
register int i,j,n;
register site *s;
register complex cc;
msg_tag *tag;
Real finalrsq,scale,x;
int tmin,tmax,cgn,color;
/* for baryon code */
int ca,ca1,ca2,cb,cb1,cb2;
void symmetry_combine(field_offset src,field_offset space,int size,int dir);
void block_fourier(
 field_offset src,	/* src is field to be transformed */
 field_offset space,	/* space is working space, same size as src */
 int size,		/* Size of field in bytes.  The field must
			   consist of size/sizeof(complex) consecutive
			   complex numbers.  For example, an su3_vector
			   is 3 complex numbers. */
 int isign);		/* 1 for x -> k, -1 for k -> x */
void fourier(
field_offset src,	/* src is field to be transformed */
field_offset space,	/* space is working space, same size as src */
int size,		/* Size of field in bytes.  The field must
			   consist of size/sizeof(complex) consecutive
			   complex numbers.  For example, an su3_vector
			   is 3 complex numbers. */
int isign);		/* 1 for x -> k, -1 for k -> x */
void write_wf(field_offset src,char *string,int tmin,int tmax);

    /* Fix TUP Coulomb gauge - gauge links only*/
    rephase( OFF );
    gaugefix(TUP,(Real)1.8,500,(Real)GAUGE_FIX_TOL);
    rephase( ON );

    for(color=0;color<3;color++){ /* Make wall source */
        FORALLSITES(i,s){
	    for(j=0;j<3;j++)s->phi.c[j]=cmplx(0.0,0.0);
	    if( s->x%2==0 && s->y%2==0 && s->z%2==0 && s->t==0 ){
                s->phi.c[color] = cmplx(-1.0,0.0);
	    }
        }
        /* do a C.G. (source in phi, result in xxx) */
	load_ferm_links(&fn_links);
        cgn = ks_congrad(F_OFFSET(phi),F_OFFSET(xxx),mass,
			 niter, rsqprop, PRECISION, EVEN, &finalrsq, 
			 &fn_links);
        /* Multiply by -Madjoint, result in propmat[color] */
        dslash_site( F_OFFSET(xxx), F_OFFSET(propmat[color]), ODD, &fn_links);
        scalar_mult_latvec( F_OFFSET(xxx), (Real)(-2.0*mass),
	    F_OFFSET(propmat[color]), EVEN);
    }


    /* construct the diquark propagator--uses tempmat1 and do this before
you fft the quark propagator */

    FORALLSITES(i,s){
	for(ca=0;ca<3;ca++)for(cb=0;cb<3;cb++){
		ca1= (ca+1)%3; ca2= (ca+2)%3;
		cb1= (cb+1)%3; cb2= (cb+2)%3;
		CMUL((s->propmat[ca1].c[cb1]),(s->propmat[ca2].c[cb2]),
			(s->tempmat1.e[ca][cb]));

	CMUL((s->propmat[ca1].c[cb2]),(s->propmat[ca2].c[cb1]),
			cc);


		CSUB((s->tempmat1.e[ca][cb]),cc,(s->tempmat1.e[ca][cb]));
	}
    }
/* complex conjugate the diquark prop */
    FORALLSITES(i,s){
	for(ca=0;ca<3;ca++)for(cb=0;cb<3;cb++){
		CONJG((s->tempmat1.e[ca][cb]),(s->tempmat1.e[ca][cb]));
	}
    }
    /* Transform the diquark propagator.  */
   block_fourier( F_OFFSET(tempmat1), F_OFFSET(tempvec[0]),
	3*sizeof(su3_vector), FORWARDS);
/* complex conjugate the diquark prop. Now we have D(-k) for convolution */
    FORALLSITES(i,s){
	for(ca=0;ca<3;ca++)for(cb=0;cb<3;cb++){
		CONJG((s->tempmat1.e[ca][cb]),(s->tempmat1.e[ca][cb]));
	}
    }

    /* Transform the propagator.  */
    block_fourier( F_OFFSET(propmat[0]), F_OFFSET(tempvec[0]),
	3*sizeof(su3_vector), FORWARDS);

/* CODE SPECIFIC TO PARTICULAR PARTICLES */

/* MESON CODE */


    /* Square the result, component by component, sum over source and
	sink colors, result in ttt.c[0] */
    FORALLSITES(i,s){
	s->ttt.c[0].real = s->ttt.c[0].imag = 0.0;
	for(color=0;color<3;color++){
	    s->ttt.c[0].real += magsq_su3vec( &(s->propmat[color]) );
	}
    }
	
    /* Transform back to coordinate space.  */
    block_fourier( F_OFFSET(ttt.c[0]), F_OFFSET(cg_p),
	sizeof(complex),BACKWARDS);

    /* make pion and rho wave functions by summing over block with
	signs */
    /* Site	sign		particle
	x,y,z
	0,0,0	+1		pi (PS)
	1,0,0	(-1)^x		rho-b1 (VT), gamma_x polarization
	0,1,0	(-1)^y		rho-b1 (VT), gamma_y polarization
	0,0,1	(-1)^z		rho-b1 (VT), gamma_z polarization
	1,1,0	(-1)^(x+y)	rho-a1 (PV), gamma_z polarization
	1,0,1	(-1)^(x+z)	rho-a1 (PV), gamma_y polarization
	0,1,1	(-1)^(y+z)	rho-a1 (PV), gamma_x polarization
	1,1,1	(-1)^(x+y+z)	pi_2, sigma (SC)
    */
    for(i=XUP;i<=ZUP;i++){
	tag = start_gather_site( F_OFFSET(ttt.c[0]), sizeof(complex),
	    butterfly_dir[i][0], EVENANDODD, gen_pt[0]);
	wait_gather(tag);
        FORALLSITES(j,s){
	    s->cg_p.c[0] = *((complex *)(gen_pt[0][j]));
	}
	FORALLSITES(j,s){
            /* Find coordinate - treat s->x,y,z,t as array */
            n = ((short *)&(s->x))[i];
	    if( (n&0x01)==0 ){  /* eta[i]=0 (coordinate within block) */
		CSUM( s->ttt.c[0], s->cg_p.c[0] );
	    }
	    else{
		CSUB( s->cg_p.c[0], s->ttt.c[0], s->ttt.c[0] );
	    }
	}
	cleanup_gather(tag);
    } /* end loop over directions */

    /* Combine displacements that are symmetry transforms */
    if(nx==ny && ny==nz){	/* average over permutations, one interchange */
	permute_combine( F_OFFSET(ttt.c[0]), F_OFFSET(cg_p.c[0]),
	    sizeof(complex), permute_dir);
	    /* average over permutations of x,y,z */
        symmetry_combine(F_OFFSET(ttt.c[0]),F_OFFSET(cg_p.c[0]),
	    sizeof(complex), switch_dir); /* combine with x <-> y site */
    }
    else{	/* maybe one pair equal */
        if(nx==ny || ny==nz || nx==nz )
	     symmetry_combine( F_OFFSET(ttt.c[0]),
	    F_OFFSET(cg_p.c[0]), sizeof(complex), switch_dir);
	    /* combine with x <-> y site, y <-> z or x <-> z */
    }
    for(i=XUP;i<=ZUP;i++){	/* average over coord <-> -coord */
        symmetry_combine( F_OFFSET(ttt.c[0]),
	    F_OFFSET(cg_p.c[0]), sizeof(complex), minus_block_dir[i]);
    }
    /* t <-> nt-t uses ordinary minus_map */
    symmetry_combine( F_OFFSET(ttt.c[0]),
        F_OFFSET(cg_p.c[0]), sizeof(complex), minus_dir[TUP]);

    /* Dump wavefunction */
    tmin = 0; tmax=nt/2;
    scale = 1.0/(sym_factor*(nx/2)*(ny/2)*(nz/2));
	/* symmetries and volume factor from FFT's */
    FORALLSITES(i,s){
	CMULREAL( s->ttt.c[0], scale, s->ttt.c[0] );
    }
    write_wf( F_OFFSET(ttt.c[0]), "WF:",tmin,tmax);
 /* END MESONS CODE */

/* BARYON CODE */
    /* Multiply quark by diquark, sum over source and
	sink colors, result in ttt.c[0] */
    FORALLSITES(i,s){
	s->ttt.c[0].real = s->ttt.c[0].imag = 0.0;
	for(ca=0;ca<3;ca++)for(cb=0;cb<3;cb++){

		CMUL((s->propmat[ca].c[cb]),(s->tempmat1.e[ca][cb]),cc);
		CSUM((s->ttt.c[0]),cc);
	}
    }
	
    /* Transform back to coordinate space.  */
    block_fourier( F_OFFSET(ttt.c[0]), F_OFFSET(cg_p),
	sizeof(complex),BACKWARDS);

    /*  nucleon wave function is (apparently) defined on EVEN sites */

    /* Combine displacements that are symmetry transforms */
    if(nx==ny && ny==nz){	/* average over permutations, one interchange */
	permute_combine( F_OFFSET(ttt.c[0]), F_OFFSET(cg_p.c[0]),
	    sizeof(complex), permute_dir);
	    /* average over permutations of x,y,t */
        symmetry_combine(F_OFFSET(ttt.c[0]),F_OFFSET(cg_p.c[0]),
	    sizeof(complex), switch_dir); /* combine with x <-> y site */
    }
    else{	/* maybe one pair equal */
        if(nx==ny || ny==nz || nx==nz )
	     symmetry_combine( F_OFFSET(ttt.c[0]),
	    F_OFFSET(cg_p.c[0]), sizeof(complex), switch_dir);
	    /* combine with x <-> y site, y <-> z or x <-> z */
    }
    for(i=XUP;i<=ZUP;i++){	/* average over coord <-> -coord */
        symmetry_combine( F_OFFSET(ttt.c[0]),
	    F_OFFSET(cg_p.c[0]), sizeof(complex), minus_block_dir[i]);
    }

    /* Dump wavefunction */
    tmin = 0; tmax=nt/2; 
    scale = 1.0/(sym_factor*(nx/2)*(ny/2)*(nz/2));
	/* symmetries and volume factor from FFT's */
    FORALLSITES(i,s){
	CMULREAL( s->ttt.c[0], scale, s->ttt.c[0] );
    }
    write_wf( F_OFFSET(ttt.c[0]), "BARW:",tmin,tmax);



} /* end wavefunct_t() */

/* BEGIN USUAL CODE FOR WAVE FUNCTION PACKAGE     */


/************************* setup_wavefunc_t() *************************/
/* Set up gathers for wave functions */
setup_wavefunc_t(){
int i,key[4],arg[2];
void butterfly_map(int x,int y,int z,int t,int *arg,
		   int fb,int *xp,int *yp,int *zp,int *tp);
void switch_map(int x,int y,int z,int t,int *arg,
		int fb,int *xp,int *yp,int *zp,int *tp);
void permute_map(int x,int y,int z,int t,int *arg,
		 int fb,int *xp,int *yp,int *zp,int *tp);
void minus_map(int x,int y,int z,int t,int *arg,
	       int fb,int *xp,int *yp,int *zp,int *tp);
void minus_block_map(int x,int y,int z,int t,int *arg,
		     int fb,int *xp,int *yp,int *zp,int *tp);
void setup_block_fourier(int *key );

    key[XUP]=1; key[YUP]=1; key[ZUP]=1; key[TUP]=0;
    setup_block_fourier(key);

    /* If two dimensions are equal and the third different, switching
	the two coordinates with equal dimensions is a symmetry.  If all
	three spatial dimensions are equal, we average the site
	x,y,z with y,z,x and z,x,y, then switch any pair of coordinates.  */
    /* Note we switch coordinates, not coordinates of doubled
       lattice.  Thus when we switch x and y we also switch
       gamma_x and gamma_y (displacements within 2-cube) */
    sym_factor=1.0;
    if( nx==ny && ny==nz ){	/* all three equal */
        /* interchange x and y coordinates */
	arg[0]=XUP; arg[1]=YUP;
        switch_dir = make_gather( switch_map, arg, OWN_INVERSE,
	    NO_EVEN_ODD, SAME_PARITY);
        /* mapping y<-x, z<-y, x<-z.  Its inverse will be permute_dir+1 */
        permute_dir = make_gather( permute_map, arg, WANT_INVERSE,
	    NO_EVEN_ODD, SAME_PARITY);
	sym_factor *= 6.0;
    }
    else {	/* not all three equal */
	if(nx==ny){	/* but one pair (nx,ny) is equal */
	    arg[0]=XUP; arg[1]=YUP;
            switch_dir = make_gather( switch_map, arg, OWN_INVERSE,
	        NO_EVEN_ODD, SAME_PARITY);
	}
	if(ny==nz){
	    arg[0]=YUP; arg[1]=ZUP;
            switch_dir = make_gather( switch_map, arg, OWN_INVERSE,
	        NO_EVEN_ODD, SAME_PARITY);
	}
	if(nx==nz){
	    arg[0]=XUP; arg[1]=ZUP;
            switch_dir = make_gather( switch_map, arg, OWN_INVERSE,
	        NO_EVEN_ODD, SAME_PARITY);
	}
	sym_factor *= 2.0;
    }

    /* mappings x<->(-x), etc. */
    for(i=XUP;i<=TUP;i++){
        key[XUP]=0; key[YUP]=0; key[ZUP]=0; key[TUP]=0;
	key[i]=1;
	if(i!=TUP) minus_block_dir[i] = make_gather( minus_block_map, key,
	    OWN_INVERSE, NO_EVEN_ODD, SAME_PARITY);
	/* For full FFT method, need ordinary minus_map */
        minus_dir[i] = make_gather( minus_map, key, OWN_INVERSE,
	    NO_EVEN_ODD, SAME_PARITY);
    }
    sym_factor *= 16.0;
} /* end setup_wavefunc_t() */




/**************************** setup_block_fourier() ****************/
void setup_block_fourier(int *key ){
    /* block_fourier Fourier transforms 2x2... blocks instead of single
	lattice sites. */
    /* "key" is a four component array.  If a component is 1, the Fourier
	transform is done in that direction, if it is 0 that direction is
	left alone. */
register int dir,i,j;
int arg[2];
void bitrev_map(int x,int y,int z,int t,int *key,
		int fb,int *xp,int *yp,int *zp,int *tp);
void bitrev_block_map(int x,int y,int z,int t,int *key,
		      int fb,int *xp,int *yp,int *zp,int *tp);
void butterfly_map(int x,int y,int z,int t,int *arg,
		   int fb,int *xp,int *yp,int *zp,int *tp);

    /* Check that relevant dimensions are power of 2, remember their logs */
    /* Allocate space for butterfly_dir */
    dim[0]=nx;dim[1]=ny;dim[2]=nz;dim[3]=nt;
    FORALLUPDIR(dir){
	if(key[dir]==1){
	    for( j=0,i=dim[dir] ; (i&0x01) == 0; i>>=1 )j++;
	    if( i != 1 ){	/* Not a power of 2 */
		if(this_node==0)printf("Can't Fourier transform dimension %d\n",
		    dim[dir]);
		terminate(1);
	    }
	    logdim[dir]=j;
	    butterfly_dir[dir] = (int *)malloc( j*sizeof(int) );
	}
	else{  /* ignore this dimension */
	    logdim[dir]= -1; 
	    butterfly_dir[dir]=NULL;/* This will remind us to ignore it */
	}
    }

    /* Set up bit-reverse for block FFT */
    bitrev_block_dir = make_gather( bitrev_block_map, key, OWN_INVERSE,
	NO_EVEN_ODD, SAME_PARITY );
    /* Set up bit-reverse for full FFT */
    /* If we don't want full FFT, don't need this one */
    bitrev_dir = make_gather( bitrev_map, key, OWN_INVERSE,
	NO_EVEN_ODD, SCRAMBLE_PARITY );


    /* Set up butterflies */
    FORALLUPDIR(dir)if(key[dir]==1){
	/* Level zero is not used by the block FFT, but is used for
	   combining within a block.  The full FFT uses all levels. */
	for(i=0;i<logdim[dir];i++){
	    arg[0] = dir;
	    arg[1] = i;
	    if(i!=0) butterfly_dir[dir][i] = make_gather( butterfly_map,
		arg, OWN_INVERSE, NO_EVEN_ODD, SAME_PARITY );
	    else     butterfly_dir[dir][i] = make_gather( butterfly_map,
		arg, OWN_INVERSE, NO_EVEN_ODD, SWITCH_PARITY );
	}
    }
} /* end setup_block_fourier() */


/********************* average fields over lattice symmetries ***********/
/* Average together field at sites related by gather in direction "dir" */
void symmetry_combine(field_offset src,field_offset space,int size,int dir){
msg_tag *tag;
register int i,j;
register site *s;
    tag = start_gather_site( src, size, dir, EVENANDODD, gen_pt[0]);
    wait_gather(tag);
    FORALLSITES(i,s){
	for(j=0;j<size/sizeof(complex);j++){
	    ((complex *)F_PT(s,space))[j] = ((complex *)(gen_pt[0][i]))[j];
	}
    }
    FORALLSITES(i,s){
	for(j=0;j<size/sizeof(complex);j++){
	    CSUM( ((complex *)F_PT(s,src))[j], ((complex *)F_PT(s,space))[j] );
	}
    }
    cleanup_gather(tag);
}

/* Average together src with sites related by forwards and backwards
   gather. Used for combining with permutations of (x,y,z) */
permute_combine(src,space,size,dir) field_offset src,space; int size,dir; {
msg_tag *tag0,*tag1;
register int i,j;
register site *s;
    tag0 = start_gather_site( src, size, dir, EVENANDODD, gen_pt[0]);
    tag1 = start_gather_site( src, size, dir+1, EVENANDODD, gen_pt[1]);
    wait_gather(tag0);
    wait_gather(tag1);
    FORALLSITES(i,s){
	for(j=0;j<size/sizeof(complex);j++){
	    CADD( ((complex *)(gen_pt[0][i]))[j],
	          ((complex *)(gen_pt[1][i]))[j],
		  ((complex *)F_PT(s,space))[j] );
	}
    }
    FORALLSITES(i,s){
	for(j=0;j<size/sizeof(complex);j++){
	    CSUM( ((complex *)F_PT(s,src))[j],
		  ((complex *)F_PT(s,space))[j] );
	}
    }
    cleanup_gather(tag0);
    cleanup_gather(tag1);
}


/****** mapping functions for gathers ***********************************/
/* Function that defines switching two coordinates (arg tells which two) */
void switch_map(int x,int y,int z,int t,int *arg,
		int fb,int *xp,int *yp,int *zp,int *tp)
{
int tt,d[4];
    d[XUP]=x; d[YUP]=y; d[ZUP]=z; d[TUP]=t;
    tt = d[arg[0]];
    d[arg[0]]=d[arg[1]];
    d[arg[1]]=tt;
    *xp=d[XUP]; *yp=d[YUP]; *zp=d[ZUP]; *tp=d[TUP];
}

/* Function that defines permuting spatial coordinates */
/* y <- x,  z <- y,  x <- z for forwards mapping.  "arg" is ignored. */
void permute_map(int x,int y,int z,int t,int *arg,
		 int fb,int *xp,int *yp,int *zp,int *tp)
{
    if(fb==FORWARDS){
	*yp=x; *zp=y, *xp=z; *tp=t;
    }
    else {
	*xp=y; *yp=z; *zp=x; *tp=t;
    }
}

/* Function that changes "sign" of selected coordinates. (arg tells which) */
/* if arg[i]==1, coord[i] -> dim[i]-coord[i], except 0 -> 0 */
/* This works on the block coordinates !!! *** !!!  */
void minus_block_map(int x,int y,int z,int t,int *arg,
		     int fb,int *xp,int *yp,int *zp,int *tp)
{
register int i,j;
    *xp=x; *yp=y; *zp=z; *tp=t;
    if(arg[XUP]==1){
	i=x/2; j=x%2;
	if(i != 0) i = nx/2-i;
	*xp = 2*i+j;
    }
    if(arg[YUP]==1){
	i=y/2; j=y%2;
	if(i != 0) i = ny/2-i;
	*yp = 2*i+j;
    }
    if(arg[ZUP]==1){
	i=z/2; j=z%2;
	if(i != 0) i = nz/2-i;
	*zp = 2*i+j;
    }
    if(arg[TUP]==1){
	i=t/2; j=t%2;
	if(i != 0) i = nt/2-i;
	*tp = 2*i+j;
    }
}

/* Function that changes "sign" of selected coordinates. (arg tells which) */
/* if arg[i]==1, coord[i] -> dim[i]-coord[i], except 0 -> 0 */
void minus_map(int x,int y,int z,int t,int *arg,
	       int fb,int *xp,int *yp,int *zp,int *tp)
{
    if(arg[XUP]==1){ if(x==0)*xp=0; else *xp = nx-x; } else *xp=x;
    if(arg[YUP]==1){ if(y==0)*yp=0; else *yp = ny-y; } else *yp=y;
    if(arg[ZUP]==1){ if(z==0)*zp=0; else *zp = nz-z; } else *zp=z;
    if(arg[TUP]==1){ if(t==0)*tp=0; else *tp = nt-t; } else *tp=t;
}

/* Function that defines the bit reverse mapping */
void bitrev_map(int x,int y,int z,int t,int *key,
		int fb,int *xp,int *yp,int *zp,int *tp)
{
   if(key[XUP]==1) *xp=bitrev_one_int(x,logdim[XUP]); else *xp=x;
   if(key[YUP]==1) *yp=bitrev_one_int(y,logdim[YUP]); else *yp=y;
   if(key[ZUP]==1) *zp=bitrev_one_int(z,logdim[ZUP]); else *zp=z;
   if(key[TUP]==1) *tp=bitrev_one_int(t,logdim[TUP]); else *tp=t;
}

/* Bit reverse a single integer, which ranges from 0 to 2^n-1 ( n bits ) */
int bitrev_one_int(int i,int n){
register int j,k;
    for(k=0,j=0;j<n;j++){ /* j counts the bits, k is result we build up */
	/* This is a dumb way to do this, but it is only the setup code */
	k |= (i&0x01);	/* pull off lowest order bit of i */
	i >>= 1;	/* throw away lowest order bit of i */
	k <<= 1;	/* open up lowest order bit of k for next bit */
    }
    k >>= 1;	/* we overran by one bit */
    return(k);
}

/* Function that defines the bit reverse mapping for block FFT */
void bitrev_block_map(int x,int y,int z,int t,int *key,
		      int fb,int *xp,int *yp,int *zp,int *tp)
{
   if(key[XUP]==1) *xp=bitrev_block_one_int(x,logdim[XUP]); else *xp=x;
   if(key[YUP]==1) *yp=bitrev_block_one_int(y,logdim[YUP]); else *yp=y;
   if(key[ZUP]==1) *zp=bitrev_block_one_int(z,logdim[ZUP]); else *zp=z;
   if(key[TUP]==1) *tp=bitrev_block_one_int(t,logdim[TUP]); else *tp=t;
}

/* Bit reverse a single integer, which ranges from 0 to 2^n-1 ( n bits ) */
/* For block FFT, lowest order bit isn't included in reversal */
int bitrev_block_one_int(int i,int n){
register int itemp,j,k;
    itemp=i;
    for(k=0,j=1;j<n;j++){ /* j counts the bits, k is result we build up */
	/* This is a dumb way to do this, but it is only the setup code */
	k |= (itemp&0x02);	/* pull off second lowest order bit of i */
	itemp >>= 1;	/* throw away lowest order bit of i */
	k <<= 1;	/* open up lowest order bit of j for next bit */
    }
    k >>= 1;	/* we overran by one bit */
    k |= i & 0x01;	/* just copy lowest order bit */
/**printf("The %d block bit reverse of %x is %x\n",n,itemp,k);**/
    return(k);
}
 
/* Function that defines the butterfly mappings */
void butterfly_map(int x,int y,int z,int t,int *arg,
		   int fb,int *xp,int *yp,int *zp,int *tp)
{
int mask,dir,level;
    /* arg contains the direction of the butterfly and the level */
    dir=arg[0];  level=arg[1];	/* just so I can remember them */
    mask = 0x01 << level;	/* one in the bit we switch */

    *xp=x; *yp=y; *zp=z; *tp=t;
    if(arg[0]==XUP) *xp ^= mask;
    if(arg[0]==YUP) *yp ^= mask;
    if(arg[0]==ZUP) *zp ^= mask;
    if(arg[0]==TUP) *tp ^= mask;
/**if(x==0 && y==0 && z==0)
printf("The level %d dir %d butterfly of %d %d %d %d is %d %d %d %d\n",
level,dir,x,y,z,t,*xp,*yp,*zp,*tp);**/
}





/******************************* block_fourier() ***********************/
/* The actual Fourier transform routine.
   The algorithm is the "decimation in frequency" scheme depicted in
   Nussbaumer figure 4.2.
*/
/* block_fourier Fourier transforms 2x2... blocks instead of single
    lattice sites. In other words, we do a Fourier transform as if it
    were an array of dimension Size/2, with each element a 2x3.. block
*/
void block_fourier(
 field_offset src,	/* src is field to be transformed */
 field_offset space,	/* space is working space, same size as src */
 int size,		/* Size of field in bytes.  The field must
			   consist of size/sizeof(complex) consecutive
			   complex numbers.  For example, an su3_vector
			   is 3 complex numbers. */
 int isign)		/* 1 for x -> k, -1 for k -> x */
{
/* local variables */
register complex *space_pt,*src_pt;
register complex cc2;
register int mask,level,i,j,n,power,dir;
register site *s;
register Real theta_0;
complex *phase;	/* array of phase factors */
msg_tag *tag;
int ncomp;	/* number of complex numbers in field */

    ncomp = size/sizeof(complex);
    /* Danielson-Lanczos section */
    /* loop over all directions, and if we are supposed to transform in
	that direction, do it */
    FORALLUPDIR(dir)if(logdim[dir] != -1){
        /* The fundamental angle, others are multiples of this */
        theta_0 = -isign*4*PI/dim[dir];
	/* Make an array of phase factors */
	phase = (complex *)malloc( (dim[dir]/4)*sizeof(complex) );
	for(i=0;i<dim[dir]/4;i++)phase[i]=ce_itheta( i*theta_0 );

	for(level=logdim[dir]-2,mask=dim[dir]>>2; level>=0; level--,mask>>=1 ){
	    /* "mask" picks out the bit that is flipped to find the
		coordinate you are combining with */
            /* In blocked version, actual bit is one to the left of mask*/

	    /* Get the site at other end of butterfly */
	    tag = start_gather_site( src, size,
		butterfly_dir[dir][level+1], EVENANDODD, gen_pt[0]);
	    wait_gather(tag);
	    FORALLSITES(i,s){
		memcpy( F_PT(s,space), gen_pt[0][i], size );
	    }
	    cleanup_gather(tag);

	    FORALLSITES(i,s){
		/* Find coordinate - treat s->x,y,z,t as array */
		n = ((short *)&(s->x))[dir];
		src_pt = (complex *)F_PT(s,src);	/* pointer to source */
		space_pt = (complex *)F_PT(s,space);	/* pointer to partner */

		/* If this is the "top" site in the butterfly, just
		   add in the partner.  If it is the "bottom" site,
		   subtract site from its partner, then multiply by
		   the phase factor. */
/**if(s->x==0 && s->y==0 && s->z==0){
printf("Local = ( %.04f , %.04f ) partner = ( %.04f , %.04f )\n",
src_pt[0].real, src_pt[0].imag, space_pt[0].real, space_pt[0].imag );
}**/
		if( (n>>1)&mask ){	/* Bottom site - the bit is one */
		    if(level==0){
		        /* Special case level 0 - all phases are 1 */
			for(j=0;j<ncomp;j++){/* loop over complex numbers */
		            CSUB( space_pt[j], src_pt[j], src_pt[j] );
			}
		    }
		    else {	/* General level */
		    	power = ((n>>1)&(mask-1))<<(logdim[dir]-level-2);
/**if(s->x==0 && s->y==0 && s->z==0)
printf("Level %d site %d power is %d\n", level,n,power);**/
			for(j=0;j<ncomp;j++){/* loop over complex numbers */
		            CSUB( space_pt[j], src_pt[j], cc2 );
			    /* cc2  <-  partner - local */
		            CMUL( phase[power], cc2, src_pt[j] );
			}
		    } /* end general level */
		}
		else {		/* Top site */
		    for(j=0;j<ncomp;j++){	/* loop over complex numbers */
		        CSUM( src_pt[j], space_pt[j] );
		    }
		}
/**if(s->x==0 && s->y==0 && s->z==0){
printf("Level %d site %d Result = ( %.04f , %.04f )\n",
level,n,src_pt[0].real, src_pt[0].imag );
}**/
	    }	/* end loop over sites */
	}  /* end loop over level */
	free(phase);
    } /* for loop on direction */

    /* Bit reverse */
    tag = start_gather_site( src, size, bitrev_block_dir,
	EVENANDODD, gen_pt[0]);
    wait_gather(tag);
    FORALLSITES(i,s){
	memcpy( F_PT(s,space), gen_pt[0][i], size );
    }
    FORALLSITES(i,s){
	memcpy( F_PT(s,src), F_PT(s,space), size );
    }
    cleanup_gather(tag);
} /* end block_fourier() */

 

/* The actual Fourier transform routine.
   The algorithm is the "decimation in frequency" scheme depicted in
   Nussbaumer figure 4.2.
*/
void fourier(
field_offset src,	/* src is field to be transformed */
field_offset space,	/* space is working space, same size as src */
int size,		/* Size of field in bytes.  The field must
			   consist of size/sizeof(complex) consecutive
			   complex numbers.  For example, an su3_vector
			   is 3 complex numbers. */
int isign)		/* 1 for x -> k, -1 for k -> x */
{
/* local variables */
register complex *space_pt,*src_pt;
register complex cc2;
register int mask,level,i,j,n,power,dir;
register site *s;
register Real theta_0;
complex *phase;	/* array of phase factors */
msg_tag *tag;
int ncomp;	/* number of complex numbers in field */

    ncomp = size/sizeof(complex);
    /* Danielson-Lanczos section */
    /* loop over all directions, and if we are supposed to transform in
	that direction, do it */
    FORALLUPDIR(dir)if(logdim[dir] != -1){
        /* The fundamental angle, others are multiples of this */
        theta_0 = -isign*2*PI/dim[dir];
	/* Make an array of phase factors */
	phase = (complex *)malloc( (dim[dir]/2)*sizeof(complex) );
	for(i=0;i<dim[dir]/2;i++)phase[i]=ce_itheta( i*theta_0 );

	for(level=logdim[dir]-1,mask=dim[dir]>>1; level>=0; level--,mask>>=1 ){
	    /* "mask" picks out the bit that is flipped to find the
		coordinate you are combining with */

	    /* Get the site at other end of butterfly */
	    tag = start_gather_site( src, size,
		butterfly_dir[dir][level], EVENANDODD, gen_pt[0]);
	    wait_gather(tag);
	    FORALLSITES(i,s){
		memcpy( F_PT(s,space), gen_pt[0][i], size );
	    }
	    cleanup_gather(tag);

	    FORALLSITES(i,s){
		/* Find coordinate - treat s->x,y,z,t as array */
		n = ((short *)&(s->x))[dir];
		src_pt = (complex *)F_PT(s,src);	/* pointer to source */
		space_pt = (complex *)F_PT(s,space);	/* pointer to partner */

		/* If this is the "top" site in the butterfly, just
		   add in the partner.  If it is the "bottom" site,
		   subtract site from its partner, then multiply by
		   the phase factor. */
/**if(s->x==0 && s->y==0 && s->z==0){
printf("Local = ( %.04f , %.04f ) partner = ( %.04f , %.04f )\n",
src_pt[0].real, src_pt[0].imag, space_pt[0].real, space_pt[0].imag );
}**/
		if( n&mask ){	/* Bottom site - the bit is one */
		    if(level==0){
		        /* Special case level 0 - all phases are 1 */
			for(j=0;j<ncomp;j++){/* loop over complex numbers */
		            CSUB( space_pt[j], src_pt[j], src_pt[j] );
			}
		    }
		    else {	/* General level */
		    	power = (n&(mask-1))<<(logdim[dir]-level-1);
/**if(s->x==0 && s->y==0 && s->z==0)
printf("Level %d site %d power is %d\n", level,n,power);**/
			for(j=0;j<ncomp;j++){/* loop over complex numbers */
		            CSUB( space_pt[j], src_pt[j], cc2 );
			    /* cc2  <-  partner - local */
		            CMUL( phase[power], cc2, src_pt[j] );
			}
		    } /* end general level */
		}
		else {		/* Top site */
		    for(j=0;j<ncomp;j++){	/* loop over complex numbers */
		        CSUM( src_pt[j], space_pt[j] );
		    }
		}
/**if(s->x==0 && s->y==0 && s->z==0){
printf("Result = ( %.04f , %.04f )\n",
src_pt[0].real, src_pt[0].imag );
}**/
	    }	/* end loop over sites */
	}  /* end loop over level */
	free(phase);
    } /* for loop on direction */

    /* Bit reverse */
    tag = start_gather_site( src, size, bitrev_dir,
	EVENANDODD, gen_pt[0]);
    wait_gather(tag);
    FORALLSITES(i,s){
	memcpy( F_PT(s,space), gen_pt[0][i], size );
    }
    FORALLSITES(i,s){
	memcpy( F_PT(s,src), F_PT(s,space), size );
    }
    cleanup_gather(tag);
}  /* end fourier() */


void write_wf(field_offset src,char *string,int tmin,int tmax)
/* this subroutine writes a wave function from time tmin to time tmax */
/* The string is a tag which identifies the output lines */
/* Output format:  STRING  t   x y z   real  imag  */
{
int currentnode,newnode;
int l,x,y,z,t;
complex lbuf;
int node0=0;

    g_sync();
    currentnode=0;

    /* write the elements */
    for(t=tmin;t<=tmax;t++)for(z=0;z<=nz/2+1;z++)
	for(y=0;y<=ny/2+1;y++)for(x=0;x<=nx/2+1;x++) {

	/* Skip elements that are redundant due to symmetries */
	if(nx==ny && y<x )continue;
	if(nx==nz && z<x )continue;
	if(ny==nz && z<y )continue;

	newnode=node_number(x,y,z,t);
	if(newnode != currentnode) {	/* switch to another node */
	    g_sync();
	    currentnode=newnode;
	}

	if(this_node==0) {
	    if(currentnode==0) {
		l=node_index(x,y,z,t);
                lbuf = *( (complex *)F_PT( &(lattice[l]), src ) );
	    }
	    else {
		get_field((char *)&lbuf,sizeof(complex),currentnode);
	    }
	    if( (printf("%s  %d   %d %d %d\t%.7e\t%.7e\n",string,t,x,y,z,
		(double)lbuf.real, (double)lbuf.imag))== EOF) {
		printf("Write error in write_wf\n"); 
		terminate(1);
	    } 
	}
	else {	/* for nodes other than 0 */
	    if(this_node==currentnode) {
		if(this_node==0) printf("this node is zero\n");
		l=node_index(x,y,z,t);
                lbuf = *( (complex *)F_PT( &(lattice[l]), src ) );
		send_field((char *)&lbuf,sizeof(complex),node0);
	    }
	}
    }
    g_sync();
    if(this_node==0) fflush(stdout);
}

/********************** TEST ROUTINE ***********************************/
test_block_fourier(){
int i,j,k;
complex cc2;
site *s;
Real phase;
int key[4];


    /* Set up Fourier transform */
	/* T direction only */
    key[XUP] = 1;
    key[YUP] = 1;
    key[ZUP] = 0;
    key[TUP] = 0;
    setup_block_fourier(key);

    /* Set source equal to e^{ikx} */
    FORALLSITES(i,s){
        /* k = 0,0,0,1 */
        if( s->x %2 == 0 && s->y%2==0 ){
	    phase =  2*PI*((Real)((s->x)>>1))/(nx>>1) ;
	    phase += 2*PI*((Real)((s->y)>>1))/(ny>>1) ;
	    s->cg_p.c[0] = ce_itheta( phase );
	}
	else
	s->cg_p.c[0] = cmplx(0.0,0.0);
	/**
	if(s->t==2)s->cg_p.c[0].real = 1.0;
	else if(s->t==3)s->cg_p.c[0].real = 1.0;
	else s->cg_p.c[0].real = 0.0;
	s->cg_p.c[0].imag = 0.0;
	**/
if(s->t==0 && s->z==0 )
printf("At site %d %d %d %d started with ( %.04f , %.04f )\n",
s->x,s->y,s->z,s->t,s->cg_p.c[0].real,s->cg_p.c[0].imag);
    }

    /* Transform it.  Use tempvec[0] as workspace */
    block_fourier( F_OFFSET(cg_p), F_OFFSET(tempvec[0]),
	sizeof(complex), FORWARDS);

    /* Dump the result */
    FORALLSITES(i,s){
if(s->t==0 && s->z==0 )
	printf("At site %d %d %d %d got ( %.04f , %.04f )\n",
	    s->x,s->y,s->z,s->t,s->cg_p.c[0].real,s->cg_p.c[0].imag);
    }
}
