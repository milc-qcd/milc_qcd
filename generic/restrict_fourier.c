/************ restrict_fourier.c *****************************************/
/* MIMD version 7 */

/* MIMD version 7.000 by T.D., D.T. */
/* Modified to permit F.T. with some coordinates fixed by C.D. */
/* Modified to allow dimension p * 2^k with any p (up to storage limit) C.D. */

/* These routines set up and perform Fourier transforms on fields in
   the lattice.  The field consists of "size" consecutive complex
   numbers.  For example, an su3_vector is three consecutive
   complex numbers, and a wilson_vector is 12.

   The setup_restrict_fourier() routine makes all the tables needed for the
   communications.

   References:
	Fast Fourier Transform and Convolution Algorithms, H.J. Nussbaumer
	(Springer series in information sciences, 1982)  QA 403.5
*/

#include "generic_includes.h"

#define TRUE 1
#define FALSE 0
#define restrict rstrict /* C-90 T3D cludge */

/* Variables set by setup_fourier() and used by fourier() */
int dim[4];		/* dimensions */
int logdim[4];		/* log_2 of dimension */
int bitrev_dir;		/* index for bit reverse gather */
int *butterfly_dir[4];	/* indices for butterflies.  First index is
			   direction, second is level.  The actual space
			   will be malloc'ed in setup_fourier(). */
			/* Level 0 is the lowest order butterfly, i.e.
			   reverse the righmost (ones) bit.  Levels
			   range from 0 to n-1 when the dimension is 2^n */
int pbaserev_dir;       /* index for base p analog of bit reverse map */
int pcyclic_dir[4];     /* index for cyclic mod p gather */
int dmin[4],dmax[4];    /* Restrictions on range of FT */
int pfactor[4];         /* Residual prime factor p */
int notbase2;           /* True if we need to do a base p != 2 transform */
  
void bitrev_map(int x,int y,int z,int t,int *key,int fb,
		int *xp,int *yp,int *zp,int *tp);
void butterfly_map(int x,int y,int z,int t,int *arg,int fb,
		   int *xp,int *yp,int *zp,int *tp);
void pbaserev_map( int x,int y,int z,int t,int *key,int fb,
		   int *xp,int *yp,int *zp,int *tp);
void pcyclic_map( int x,int y,int z,int t,int *arg,int fb,
		  int *xp,int *yp,int *zp,int *tp);

void setup_restrict_fourier( int *key, int *slice){
  /* "key" is a four component array.  If a component is 1, the Fourier
     transform is done in that direction, if it is 0 that direction is
     left alone. 
     If it is 2, then do the FT on a subset of the lattice with a
     fixed value slice[dir] of the coordinate in that direction. 
     "slice" is a four component array.  Not used unless key[dir]=2. */
  register int dir,i,j;
  int arg[2];
  char myname[] = "setup_restrict_fourier";

  /* Initialize limits on coordinate range of FT */
  dmin[XUP] = dmin[YUP] = dmin[ZUP] = dmin[TUP] = 0;
  dmax[XUP] = nx; dmax[YUP] = ny; dmax[ZUP] = nz; dmax[TUP] = nt;
  
  /* Check that relevant dimensions are power of 2, remember their logs */
  /* Allocate space for butterfly_dir */
  /* Set restrictions on coordinates */
  dim[0]=nx;dim[1]=ny;dim[2]=nz;dim[3]=nt;
  notbase2 = FALSE;
  FORALLUPDIR(dir){
    pfactor[dir] = 1;
    if(key[dir]==1){
      for( j=0,i=dim[dir] ; (i&0x01) == 0; i>>=1 )j++;
      pfactor[dir] = i;
      notbase2 |= (i != 1);  /* If dim[dir] != 2^k and key[dir] != 1 */
      logdim[dir]=j;
      butterfly_dir[dir] = (int *)malloc( j*sizeof(int) );
      if(butterfly_dir[dir]==NULL){
	printf("%s(%d): no room for butterfly_dir\n",myname,
	       this_node);
	terminate(1);
      }

    }
    else{  /* ignore this dimension */
      logdim[dir]= -1; 
      butterfly_dir[dir]=NULL;	/* This will remind us to ignore it */
    }

    if(key[dir]==2)
      {
	dmin[dir] = slice[dir];
	dmax[dir] = slice[dir] + 1;
      }
  }
  
  /* Set up bit-reverse */
  bitrev_dir = make_gather( bitrev_map, key, OWN_INVERSE,
			   NO_EVEN_ODD, SCRAMBLE_PARITY );
  /* Set up butterflies */
  FORALLUPDIR(dir)if(key[dir]==1){
    for(i=0;i<logdim[dir];i++){
      arg[0] = dir;
      arg[1] = i;
      butterfly_dir[dir][i] = make_gather( butterfly_map, arg,
			  OWN_INVERSE, NO_EVEN_ODD, SCRAMBLE_PARITY );
    }
  }

  /* Set up base p cyclic gather */
  if(notbase2)FORALLUPDIR(dir)if(key[dir]==1)if(pfactor[dir]!=1)
    {
      arg[0] = dir;
      pcyclic_dir[dir] = make_gather( pcyclic_map, arg,
			 NO_INVERSE, NO_EVEN_ODD, SCRAMBLE_PARITY );
    }
    /* Set up analog of bit-reverse map for base p */
    if(notbase2)pbaserev_dir = make_gather( pbaserev_map, key,
			 NO_INVERSE, NO_EVEN_ODD, SCRAMBLE_PARITY );
}

void cleanup_restrict_fourier(void){
  int dir;

  /* Unfortunately, we have no way to unmake gathers,
     so we can't clean up the most important allocation here.
     This means there will be a memory leak */

  FORALLUPDIR(dir){
    if(butterfly_dir[dir] != NULL)
      free(butterfly_dir[dir]);
  }
}
  
/* Bit reverse a single integer, which ranges from 0 to 2^n-1 ( n bits ) */
/* This version allows for a dimension p * 2^n.
   In this case the sites are arranged in groups mod p and the bit reverse
   map is done in p classes, identified by the coordinate mod p.
   e.g. for dimension 12, (p = 3) the sites are grouped in 3's with the
   classification abcabcabcabcabc.  Then all sites labeled a participate in
   one bit reverse map, sites labeled b in a second bit reverse map, and the
   same for c.  The result of the bit reverse map is then the ordering
   0 1 2  6 7 8  3 4 5  9 10 11 */
int bitrev_one_int( int i,int n,int p ) {
  register int j,k;
  int idivp,imodp;
  if(p!=1){idivp = i/p;  imodp = i % p;} else{idivp = i; imodp = 0;}
  for(k=0,j=0;j<n;j++){ /* j counts the bits, k is result we build up */
    /* This is a dumb way to do this, but it is only the setup code */
    k |= (idivp&0x01);	/* pull off lowest order bit of idivp */
    idivp >>= 1;	/* throw away lowest order bit of idivp */
    k <<= 1;	/* open up lowest order bit of k for next bit */
  }
  k >>= 1;	/* we overran by one bit */
  return(k*p + imodp);
}

/* Function that defines the bit reverse mapping */
void bitrev_map(int x,int y,int z,int t,int *key,int fb,
		int *xp,int *yp,int *zp,int *tp)
{
  *xp = x; *yp = y; *zp = z; *tp = t;
  /* If any coordinate is outside a restricted range, do not touch the point */
  if((x < dmin[XUP])||(x >= dmax[XUP]))return;
  if((y < dmin[YUP])||(y >= dmax[YUP]))return;
  if((z < dmin[ZUP])||(z >= dmax[ZUP]))return;
  if((t < dmin[TUP])||(t >= dmax[TUP]))return;

  /* Otherwise, if we are doing an FT in that direction, bit reverse */
  if(key[XUP]==1) *xp=bitrev_one_int(x,logdim[XUP],pfactor[XUP]);
  if(key[YUP]==1) *yp=bitrev_one_int(y,logdim[YUP],pfactor[YUP]);
  if(key[ZUP]==1) *zp=bitrev_one_int(z,logdim[ZUP],pfactor[ZUP]);
  if(key[TUP]==1) *tp=bitrev_one_int(t,logdim[TUP],pfactor[TUP]);
}

int butterfly_map_one_int( int i, int mask, int p)
{
  if(p!=1)return ((i/p ^ mask) * p + (i % p)); else return (i^mask);
}

/* Function that defines the butterfly mappings */
/* This version allows for a dimension p * 2^n
   In this case the butterfly map is done for the 2^n part,
   treating the sites in groups of p */
void butterfly_map(int x,int y,int z,int t,int *arg,int fb,
		   int *xp,int *yp,int *zp,int *tp)
{
  int mask,level;
  /* arg contains the direction of the butterfly and the level */
  level=arg[1];	/* just so I can remember them */
  mask = 0x01 << level;	/* one in the bit we switch */
  
  *xp=x; *yp=y; *zp=z; *tp=t;
  /* If any coordinate is outside a restricted range, do not touch the point */
  if((x < dmin[XUP])||(x >= dmax[XUP]))return;
  if((y < dmin[YUP])||(y >= dmax[YUP]))return;
  if((z < dmin[ZUP])||(z >= dmax[ZUP]))return;
  if((t < dmin[TUP])||(t >= dmax[TUP]))return;

  /* Otherwise, map it */
  if(arg[0]==XUP) *xp = butterfly_map_one_int(*xp,mask,pfactor[XUP]);
  if(arg[0]==YUP) *yp = butterfly_map_one_int(*yp,mask,pfactor[YUP]);
  if(arg[0]==ZUP) *zp = butterfly_map_one_int(*zp,mask,pfactor[ZUP]);
  if(arg[0]==TUP) *tp = butterfly_map_one_int(*tp,mask,pfactor[TUP]);
}

/* Do single level base p reverse map on one coordinate */
int pbaserev_one_int(int i, int p, int fb, int n)
{
  int pow2;
  /* The dimension for this site factors as n = p * pow2
     e.g., the ordering for values for n = 12 (p = 3) (pow2 = 4) 
     before gathering is
     0,4,8,1,5,9,2,6,10,3,7,11
     These must be gathered to produce serial ordering
     i is the serially ordered site index and the return value is the
     site from which the value must be gathered */
  
  pow2 = n/p;
  if(fb == FORWARDS) return p * (i % pow2) + i/pow2;
  else                return pow2 * (i % p) + i/p;
}

/* Function that defines the analog of bit-reverse for base p */
void pbaserev_map( int x,int y,int z,int t,int *key,int fb,
		  int *xp,int *yp,int *zp,int *tp)
{
  *xp = x; *yp = y; *zp = z; *tp = t;

  /* If any coordinate is outside a restricted range, do not touch the point */
  if((x < dmin[XUP])||(x >= dmax[XUP]))return;
  if((y < dmin[YUP])||(y >= dmax[YUP]))return;
  if((z < dmin[ZUP])||(z >= dmax[ZUP]))return;
  if((t < dmin[TUP])||(t >= dmax[TUP]))return;

  /* Otherwise, map it */
  if((pfactor[XUP]!=1)&&(key[XUP]==1))
    *xp=pbaserev_one_int(x,pfactor[XUP],fb,nx);
  if((pfactor[YUP]!=1)&&(key[YUP]==1))
    *yp=pbaserev_one_int(y,pfactor[YUP],fb,ny);
  if((pfactor[ZUP]!=1)&&(key[ZUP]==1))
    *zp=pbaserev_one_int(z,pfactor[ZUP],fb,nz);
  if((pfactor[TUP]!=1)&&(key[TUP]==1))
    *tp=pbaserev_one_int(t,pfactor[TUP],fb,nt);
}

/* Function that defines the base p cyclic map */
void pcyclic_map( int x,int y,int z,int t,int *arg,int fb,
		 int *xp,int *yp,int *zp,int *tp)
{
  /* arg contains the direction of the cyclic map */
  /* This map permutes the sites cyclically in groups mod p */
  *xp=x; *yp=y; *zp=z; *tp=t;
  /* If any coordinate is outside a restricted range, do not touch the point */
  if((x < dmin[XUP])||(x >= dmax[XUP]))return;
  if((y < dmin[YUP])||(y >= dmax[YUP]))return;
  if((z < dmin[ZUP])||(z >= dmax[ZUP]))return;
  if((t < dmin[TUP])||(t >= dmax[TUP]))return;

  /* Otherwise, if we are doing an FT in that direction, map it */
  if(fb == FORWARDS)
    {
      if(arg[0]==XUP){*xp += 1; if((*xp % pfactor[XUP])==0)*xp -=pfactor[XUP];}
      if(arg[0]==YUP){*yp += 1; if((*yp % pfactor[YUP])==0)*yp -=pfactor[YUP];}
      if(arg[0]==ZUP){*zp += 1; if((*zp % pfactor[ZUP])==0)*zp -=pfactor[ZUP];}
      if(arg[0]==TUP){*tp += 1; if((*tp % pfactor[TUP])==0)*tp -=pfactor[TUP];}
    }
  else
    {
      if(arg[0]==XUP){if((*xp % pfactor[XUP])==0)*xp+=pfactor[XUP];*xp -= 1; }
      if(arg[0]==YUP){if((*yp % pfactor[YUP])==0)*yp+=pfactor[YUP];*yp -= 1; }
      if(arg[0]==ZUP){if((*zp % pfactor[ZUP])==0)*zp+=pfactor[ZUP];*zp -= 1; }
      if(arg[0]==TUP){if((*tp % pfactor[TUP])==0)*tp+=pfactor[TUP];*tp -= 1; }
    }
}

/* The actual Fourier transform routine.
   The algorithm is the "decimation in frequency" scheme depicted in
   Nussbaumer figure 4.2.
   Modified to incoporate possible restriction on coordinate
   and to handle a residual factor p that is not a power of 2
   */
void restrict_fourier_field(
     complex *src,	 /* src is field to be transformed */
     int size,		 /* Size of field in bytes.  The field must
			    consist of size/sizeof(complex) consecutive
			    complex numbers.  For example, an su3_vector
			    is 3 complex numbers. */
     int isign)		 /* 1 for x -> k, -1 for k -> x */
{
  /* local variables */
  register complex *space_pt,*src_pt;
  register complex cc2;
  register int mask,level,i,j,n=0,power,dir;
  register site *s;
  register Real theta_0;
  complex *phase;	/* array of phase factors */
  msg_tag *tag;
  int ncomp;	/* number of complex numbers in field */
  int x,y,z,t;
  int ndivp,nmodp,ncycmodp,dimdirdivp;
  int jcycle;
  complex *space,*space2;
  char myname[] = "restrict_fourier_field";

  ncomp = size/sizeof(complex);

  space  = (complex *)malloc(sizeof(complex)*ncomp*sites_on_node);
  space2 = (complex *)malloc(sizeof(complex)*ncomp*sites_on_node);

  /* Danielson-Lanczos section for factor of 2 levels */
  /* loop over all directions, and if we are supposed to transform in
     that direction, do it */
  FORALLUPDIR(dir)if(logdim[dir] != -1){
    /* The fundamental angle, others are multiples of this */
    dimdirdivp = dim[dir]/pfactor[dir];
    theta_0 = -isign*2*PI/dimdirdivp;
    /* Make an array of phase factors */
    phase = (complex *)malloc( (dim[dir]/2)*sizeof(complex) );
    if(phase == NULL){
      printf("%s(%d): no room for 2's phase\n",myname,
		   this_node);
            terminate(1);
	  }      
    for(i=0;i<dimdirdivp/2;i++)phase[i]=ce_itheta( i*theta_0 );
    
    for(level=logdim[dir]-1,mask=dimdirdivp>>1; level>=0; level--,mask>>=1 ){
      /* "mask" picks out the bit that is flipped to find the
	 coordinate you are combining with */
      
      /* Get the site at other end of butterfly */
      tag = start_gather_field( src, size,
			 butterfly_dir[dir][level], EVENANDODD, gen_pt[0]);
      wait_gather(tag);
      
      /* Loop over all sites involved in FT with restriction on range */
      for(x=dmin[XUP]; x < dmax[XUP]; x++)for(y=dmin[YUP]; y < dmax[YUP]; y++)
	for(z=dmin[ZUP]; z < dmax[ZUP]; z++)for(t=dmin[TUP]; t < dmax[TUP]; t++)
	  {
	    if( node_number(x,y,z,t) != mynode() )continue;
	    i = node_index(x,y,z,t);
	    s = &lattice[i];
	    memcpy( space+i*ncomp, gen_pt[0][i], size );
	  }
      cleanup_gather(tag);


      
      /* For residual factor p=pfactor[dir], the coordinate index is treated 
	 in groups of p */
      for(x=dmin[XUP]; x < dmax[XUP]; x++)for(y=dmin[YUP]; y < dmax[YUP]; y++)
	for(z=dmin[ZUP]; z < dmax[ZUP]; z++)for(t=dmin[TUP]; t < dmax[TUP]; t++)
	  {
	    if( node_number(x,y,z,t) != mynode() )continue;
	    i = node_index(x,y,z,t);
	    s = &lattice[i];

	    /* Find coordinate - treat s->x,y,z,t as array */
/*	    n = ((short *)&(s->x))[dir];*/ /* Doesn't work with the T3D /mpp/bin/cc July, 1994 */

	    if(dir==XUP)n = s->x;
	    else if(dir==YUP)n = s->y;
	    else if(dir==ZUP)n = s->z;
	    else if(dir==TUP)n = s->t;
	    /* Reduction of coordinate index */
	    ndivp = n/pfactor[dir];
	    src_pt = src + i*ncomp;	/* pointer to source */
	    space_pt = space + i*ncomp;	/* pointer to partner */
	    
	    /* If this is the "top" site in the butterfly, just
	       add in the partner.  If it is the "bottom" site,
	       subtract site from its partner, then multiply by
	       the phase factor. */

	    if( ndivp & mask ){ /* Bottom site - the bit is one */
	      if(level==0){
		/* Special case level 0 - all phases are 1 */
		for(j=0;j<ncomp;j++){/* loop over complex numbers */
		  CSUB( space_pt[j], src_pt[j], src_pt[j] );
		}
	      }
	      else {	/* General level */
		power = (ndivp&(mask-1))<<(logdim[dir]-level-1);
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
	  }	/* end loop over sites */
    }  /* end loop over level */

    free(phase);

  } /* for loop on direction */

  /* Bit reverse */
  tag = start_gather_field( src, size, bitrev_dir,
		     EVENANDODD, gen_pt[0]);
  wait_gather(tag);
  
  
  for(x=dmin[XUP]; x < dmax[XUP]; x++)for(y=dmin[YUP]; y < dmax[YUP]; y++)
    for(z=dmin[ZUP]; z < dmax[ZUP]; z++)for(t=dmin[TUP]; t < dmax[TUP]; t++)
      {
	if( node_number(x,y,z,t) != mynode() )continue;
	i = node_index(x,y,z,t);
	s = &lattice[i];
	
	memcpy( space + i*ncomp, gen_pt[0][i], size );
      }

  cleanup_gather(tag);

  /* Copy back to "src" */
  for(x=dmin[XUP]; x < dmax[XUP]; x++)
    for(y=dmin[YUP]; y < dmax[YUP]; y++)
      for(z=dmin[ZUP]; z < dmax[ZUP]; z++)
	for(t=dmin[TUP]; t < dmax[TUP]; t++)
	  {
	    if( node_number(x,y,z,t) != mynode() )continue;
	    i = node_index(x,y,z,t);
	    s = &lattice[i];
	    memcpy( src + i*ncomp, space + i*ncomp, size );
	  }

  /* Debug */

  /* Dump current values */
/**  for(x=0;x<nx;x++)for(y=0;y<ny;y++)for(z=0;z<nz;z++)for(t=0;t<nt;t++)
    {
      if( node_number(x,y,z,t) != mynode() )continue;
      i = node_index(x,y,z,t);
      s = &lattice[i];
      src_pt = (complex *)F_PT(s,src);
      printf("BPC%d %d %d %d %d %f %f\n",this_node,x,y,z,t,
	     src_pt[0].real,src_pt[0].imag);
      fflush(stdout);
    }**/

  /* Previous result is in "src" AND "space" */

  /* Extension of Danielson-Lanczos for residual factor */
  FORALLUPDIR(dir)if((logdim[dir] != -1)&&(pfactor[dir] != 1)){

    /* Start asynchronous cyclic gather from "space"*/
    tag = start_gather_field( space, size, pcyclic_dir[dir],
		       EVENANDODD, gen_pt[0]);

    /* The fundamental angle, others are multiples of this */
    theta_0 = -isign*2*PI/dim[dir];
    /* Make an array of phase factors */
    phase = (complex *)malloc( dim[dir]*sizeof(complex) );
    if(phase == NULL){
      printf("%s(%d): no room for odd phase\n",myname,
	     this_node);
      terminate(1);
    }      
    for(i=0;i<dim[dir];i++)phase[i]=ce_itheta( i*theta_0 );

    /* Initialize by multiplying by appropriate phase factor */
    /* Values come from "space" and end up in "src" */
    for(x=dmin[XUP]; x < dmax[XUP]; x++)for(y=dmin[YUP]; y < dmax[YUP]; y++)
      for(z=dmin[ZUP]; z < dmax[ZUP]; z++)for(t=dmin[TUP]; t < dmax[TUP]; t++)
	{
	  if( node_number(x,y,z,t) != mynode() )continue;
	  i = node_index(x,y,z,t);
	  s = &lattice[i];
	  /* Find coordinate - treat s->x,y,z,t as array */
	  /* n = ((short *)&(s->x))[dir]; */ /* Doesn't work on T3D /mpp/bin/cc  July, 1994 */
	  if(dir==XUP)n = s->x;
	  else if(dir==YUP)n = s->y;
	  else if(dir==ZUP)n = s->z;
	  else if(dir==TUP)n = s->t;

	  ndivp = n/pfactor[dir]; nmodp = n % pfactor[dir];
	  src_pt = src + i*ncomp;	/* pointer to source */
	  space_pt = space + i*ncomp;	/* pointer to partner */
	  power = (nmodp * dim[dir]/pfactor[dir] + ndivp)*nmodp;
	  power %= dim[dir];
	  for(j=0;j<ncomp;j++){	/* loop over complex numbers */
	    CMUL( phase[power], space_pt[j], src_pt[j]);
	  }
	}

    /* Then for each of p - 1 cyclic permutations, accumulate next factor */
    for(jcycle=1;jcycle<pfactor[dir];jcycle++)
      {
	/* Wait for gather that was started before the initialization pass */
	wait_gather(tag);
	
	/* Gather goes to "space2" */
	for(x=dmin[XUP]; x < dmax[XUP]; x++)
	  for(y=dmin[YUP]; y < dmax[YUP]; y++)
	    for(z=dmin[ZUP]; z < dmax[ZUP]; z++)
	      for(t=dmin[TUP]; t < dmax[TUP]; t++)
		{
		  if( node_number(x,y,z,t) != mynode() )continue;
		  i = node_index(x,y,z,t);
		  s = &lattice[i];
		  memcpy( space2 + i*ncomp, gen_pt[0][i], size );
		}
	cleanup_gather(tag);
	
	/* Copy gathered values from "space2" to "space" */
	for(x=dmin[XUP]; x < dmax[XUP]; x++)
	  for(y=dmin[YUP]; y < dmax[YUP]; y++)
	    for(z=dmin[ZUP]; z < dmax[ZUP]; z++)
	      for(t=dmin[TUP]; t < dmax[TUP]; t++)
		{
		  if( node_number(x,y,z,t) != mynode() )continue;
		  i = node_index(x,y,z,t);
		  s = &lattice[i];
		  memcpy( space + i*ncomp, space2 + i*ncomp, size );
		}

	/* Start next asynchronous cyclic gather from "space" if needed */
	if(pfactor[dir]-jcycle>1)
	  tag = start_gather_field( space, size, pcyclic_dir[dir],
			     EVENANDODD, gen_pt[0]);

	/* Accumulate result from "space" to "src" */
	for(x=dmin[XUP]; x < dmax[XUP]; x++)
	  for(y=dmin[YUP]; y < dmax[YUP]; y++)
	    for(z=dmin[ZUP]; z < dmax[ZUP]; z++)
	      for(t=dmin[TUP]; t < dmax[TUP]; t++)
		{
		  if( node_number(x,y,z,t) != mynode() )continue;
		  i = node_index(x,y,z,t);
		  s = &lattice[i];
		  
		  /* Find coordinate - treat s->x,y,z,t as array */
		  /* n = ((short *)&(s->x))[dir]; */ /* Doesn't work on T3D /mpp/bin/cc  July, 1994 */
		  if(dir==XUP)n = s->x;
		  else if(dir==YUP)n = s->y;
		  else if(dir==ZUP)n = s->z;
		  else if(dir==TUP)n = s->t;

		  ndivp = n/pfactor[dir]; nmodp = n % pfactor[dir];
		  /* j-Cyclic transform on coordinate */
		  ncycmodp = (n + jcycle) % pfactor[dir];
		  src_pt = src + i*ncomp;	/* pointer to source */
		  space_pt = space + i*ncomp;	/* pointer to partner */
		  power = (nmodp * dim[dir]/pfactor[dir] + ndivp)*ncycmodp;
		  power %= dim[dir];
		  for(j=0;j<ncomp;j++){	/* loop over complex numbers */
		    CMUL( phase[power], space_pt[j], cc2);
		    CSUM( src_pt[j],cc2);
		  }
		}
      } /* loop over cyclic permutations */
    
    free(phase);

    /* Copy gathered values from "src" to "space" for next dir */
    for(x=dmin[XUP]; x < dmax[XUP]; x++)
      for(y=dmin[YUP]; y < dmax[YUP]; y++)
	for(z=dmin[ZUP]; z < dmax[ZUP]; z++)
	  for(t=dmin[TUP]; t < dmax[TUP]; t++)
	    {
	      if( node_number(x,y,z,t) != mynode() )continue;
	      i = node_index(x,y,z,t);
	      s = &lattice[i];
	      memcpy( space + i*ncomp, src + i*ncomp, size );
	    }

  } /* loop over directions for cyclic base p transform */

  /* Previous result is in "src" */

  /* Debug */

  /* Dump current values */
/**  for(x=0;x<nx;x++)for(y=0;y<ny;y++)for(z=0;z<nz;z++)for(t=0;t<nt;t++)
    {
      if( node_number(x,y,z,t) != mynode() )continue;
      i = node_index(x,y,z,t);
      s = &lattice[i];
      src_pt = src + i*ncomp;
      printf("BPB%d %d %d %d %d %f %f\n",this_node,x,y,z,t,
	     src_pt[0].real,src_pt[0].imag);
      fflush(stdout);
    }**/

  /* Previous result is in "src" AND "space" */

  /* Do p base analog of bit-reverse transform */
  if(notbase2)
    {
	tag = start_gather_field( src, size, pbaserev_dir,
			   EVENANDODD, gen_pt[0]);
	wait_gather(tag);

	/* Gather into space */

	for(x=dmin[XUP]; x < dmax[XUP]; x++)
	  for(y=dmin[YUP]; y < dmax[YUP]; y++)
	    for(z=dmin[ZUP]; z < dmax[ZUP]; z++)
	      for(t=dmin[TUP]; t < dmax[TUP]; t++)
		{
		  if( node_number(x,y,z,t) != mynode() )continue;
		  i = node_index(x,y,z,t);
		  s = &lattice[i];
		  memcpy( space + i*ncomp, gen_pt[0][i], size );
		}
	
	cleanup_gather(tag);
	
	/* Then copy back to src */
	
	for(x=dmin[XUP]; x < dmax[XUP]; x++)
	  for(y=dmin[YUP]; y < dmax[YUP]; y++)
	    for(z=dmin[ZUP]; z < dmax[ZUP]; z++)
	      for(t=dmin[TUP]; t < dmax[TUP]; t++)
		{
		  if( node_number(x,y,z,t) != mynode() )continue;
		  i = node_index(x,y,z,t);
		  s = &lattice[i];
		  
		  memcpy( src + i*ncomp, space + i*ncomp, size );
		}
      } /* If not base 2 */

  free(space); free(space2);
  
  /* Final result is in "src" */
}

void restrict_fourier_site(
     field_offset src,	 /* src is field to be transformed */
     int size,		 /* Size of field in bytes.  The field must
			    consist of size/sizeof(complex) consecutive
			    complex numbers.  For example, an su3_vector
			    is 3 complex numbers. */
     int isign)		 /* 1 for x -> k, -1 for k -> x */
{
  int i;
  site *s;
  complex *t_src;
  int ncomp = size/sizeof(complex);

  t_src  = (complex *) malloc(sites_on_node*size);

  if(t_src == NULL){
    printf("restrict_fourier_site(%d): Can't allocate src\n",this_node);
    terminate(1);
  }

  /* copy src to temporary */
  FORALLSITES(i,s) {
    memcpy( t_src + i*ncomp, F_PT(s,src), size );
  }

  restrict_fourier_field( t_src, size, isign);

  /* copy src back */
  FORALLSITES(i,s) {
    memcpy( F_PT(s,src), t_src + i*ncomp, size );
  }

  free(t_src);
}
