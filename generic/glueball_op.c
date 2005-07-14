/****** glueball_op.c  -- ******************/
/* MIMD version 7 */

/* C.D. 11/18/01 Patterned after gauge_stuff.c */

/* Measures glueball operators at zero momentum 

   Glueball operators are constructed from loop paths, defined in
   analogy with Tom DeGrand's general gauge action scheme.  Here, it
   is necessary to keep track of the structure of the cubic group so
   we can project onto irreps of O_h.

   We start with a set of basic loop shapes, defined in glueball_op.h.
   These can be changed to suit.  The make_glueball_ops procedure
   applies O_h transformations on the basic shapes, leading to a set
   of unique orientated loops.  These are measured in the
   measure_glueball_ops procedure.  

   A projection onto O_h irreps is done using the generating
   machine:

      Op(R) = \sum_g char(R,g) U(g) Op


   where Op is the basic shape, R labels the irrep, char(R,g) is
   character of the irrep R for group element g, and U(g) Op is the
   transformed basic shape.

   The set U(g) Op is organized into a list of unique orientations
   resulting from transfomations of the basic shape Op.  Usually there
   are fewer unique orientations than group elements g.  Only the
   unique shapes are tabulated for measuring.  It is convenient to sum
   the characters of all group elements leading to the same unique
   orientation.  This sum is called the "weight" of the orientation.
   
   The expectation values of Op(R) are written for each starting shape
   Op and for each time slice, so correlators can be computed offline.  */


#include "generic_includes.h"	     /* generally used prototypes for
                                        code in generic directory */
#include "../include/cubic_group.h"  /* defines cubic group */
#include <glueball_op.h>             /* defines glueball operator
                                        shapes */

#define LOOPEND  /* We make it standard */
#ifdef LOOPEND
#undef FORALLSITES
#define FORALLSITES(i,s) \
{ register int loopend; loopend=sites_on_node; \
for( i=0,  s=lattice ; i<loopend; i++,s++ )
#define END_LOOP }
#else
#define END_LOOP        /* define it to be nothing */
#endif

/*----------------------------------------------------------------------*/
void loop_print( const link_path *path ){
  const char dir_names[NDIRS][3] =
  {"+x", "+y", "+z", "+t", "-t", "-z", "-y", "-x"};
 
  register int i;
  node0_printf("\t(");
  for(i = 0; i < path->length; i++)
    node0_printf(" %s",dir_names[path->dir[i]]);
  node0_printf(",  L = %d )\n", path->length );
}

/*----------------------------------------------------------------------*/
void loop_reverse( link_path *dest, const link_path *src ){
  /* Reverses link path.  i.e. reflects dir and runs path backwards */
  register int j;
  register int length = src->length;

  for( j = 0; j < length; j++)
    dest->dir[j] = OPP_DIR(src->dir[length-j-1]);  
  dest->length = length;
}

/*----------------------------------------------------------------------*/
/* Convert loop path into a number with digits equal to
   the direction code in dirs.h, starting with position "start"
   and running cyclicly through the path */

int loop_digits( const link_path *path, int start){
  int base = 10;   /* Must be not less than the maximum direction code */
  int length = path->length;
  int digits;
  int j;

  digits = path->dir[(length-1+start) % length];
  for(j = length - 2; j >= 0; j--) 
    digits = digits*base + path->dir[(j + start) % length];

  return digits;
}

/*----------------------------------------------------------------------*/
/* Compute a number (signature) uniquely identifying the loop path up
   to cyclic permutations and path reversal.  The number identifies a
   canonical ordering of the path directions.  The signature returned
   is negative, if it was necessary to reverse the path to put it in
   canonical order 

   Code based on gauge_stuff:char_num()
*/

int loop_signature( const link_path *forwards){
  int j;
  register int length = forwards->length;
  link_path reverse;
  int digits,sign;
  int newv;

  /* forwards is array of directions.  reverse is the reversed path
     (i.e. reverse order of reflected directions). */
  
  sign = 1;

  /* digits = decimal digits of forward path, counting from 0*/
  digits = loop_digits( forwards, 0);

  /* forward */
  /* Make "digits" be the smallest number obtained by cyclic perms */

  for(j = 1; j < length; j++){
    newv = loop_digits ( forwards, j);
    if(newv < digits) digits = newv;
  }
  
  /* backward*/
  loop_reverse( &reverse, forwards );

  for(j = 0; j < length; j++){
    newv = loop_digits ( &reverse, j);
    if(newv < digits){
      digits = newv;
      /* Change sign, because the reverse direction is canonical */
      sign = -1;
    }
  }
  
  return sign*digits;
  
} /* loop_signature */

/*----------------------------------------------------------------------*/
/* Apply cubic group transformation on loop path: dest = U_g src */

void loop_cubic_group_transform(link_path *dest, const cubic_group *g, 
				const link_path *src){
  int pp[NDIRS];
  int length = src->length;
  int j;

  /* The O_h transformation defines a coordinate permutation
     and reflection */
  pp[TUP] = TUP;  pp[TDOWN] = TDOWN;  /* No transformations on t */
  for(j = XUP; j < TUP; j++){
    pp[j] = g->p[j];
    pp[OPP_DIR(j)] = OPP_DIR(pp[j]);
  }
  
  /* Apply the O_h transformation on the initial shape
     resulting in the new orientation */
  dest->length = length;
  for(j = 0; j <length; j++) 
    dest->dir[j] = pp[src->dir[j]];
}


/*----------------------------------------------------------------------*/
/* Make the table of oriented loops to be measured */
void make_glueball_ops() {
  
  char my_name[] = "make_glueball_ops";
  int length,iloop,i,j,signature;
  link_path xfm_shape;
  int orientation, flag;
  int kloop,irrep,save_orient,kclass,ielem;
  /* for each rotation/reflection, an integer distinct for each
     starting point, or each cyclic permutation of the links */
  int loop_signatures[MAX_ORIENTATION];

  /* Dump list of glueball operator shapes in use */

  node0_printf("\nGlueball operator shapes:\n");
  for(iloop=0;iloop<NGLUEBALL_BASIC_SHAPES;iloop++){
    node0_printf("%2d %6s ",iloop,glueball_shape[iloop].name);     
    loop_print(&glueball_shape[iloop]);
  }

  /* Run through all shapes */
  for(iloop=0;iloop<NGLUEBALL_BASIC_SHAPES;iloop++){
    length=glueball_shape[iloop].length;
    /* Count distinct orientations for this loop */
    orientation=0;

    /* Apply all O_h transformations on the current operator shape */
    for(ielem=0; ielem<N_CUBIC_GROUP; ielem++)
      {

	/* Do cubic group transform on loop resulting in xfm_shape */
	loop_cubic_group_transform(&xfm_shape,
				   &g[ielem],
				   &glueball_shape[iloop]);

	/* See if this xfm_shape is already in our list of unique
           orientations or a reflection of one of them */
	signature = loop_signature(&xfm_shape);

	flag = 0;
	for(save_orient=0;save_orient<orientation;save_orient++) 
	  {
	    if(signature == loop_signatures[save_orient])
	      { flag =  1; break; }  /* a match without reflection */
	    else if(signature == -loop_signatures[save_orient])
	      { flag = -1; break; }  /* a match with reflection */
	  }
	
	/* If not found, add new loop path to tables */
	if(flag == 0){
	  save_orient = orientation;
	  loop_signatures[orientation] = signature;
	  glueball_op[iloop].loop[orientation].path.length = length;
	  for(j=0;j<length;j++)	
	    glueball_op[iloop].loop[orientation].path.dir[j] 
	      = xfm_shape.dir[j];
	  for(irrep = 0; irrep < N_CUBIC_IRREP; irrep++){
	    glueball_op[iloop].loop[orientation].weight_r[irrep] = 0;
	    glueball_op[iloop].loop[orientation].weight_i[irrep] = 0;
	  }
	  /* Allocate space for trace of operator */
	  glueball_op[iloop].loop[orientation].trace = 
	    (double_complex *)malloc(nt*sizeof(double_complex));
	  if(glueball_op[iloop].loop[orientation].trace == NULL){
	    printf("%s: No room for operator trace\n");
	    terminate(1);
	  }

	  /** node0_printf("ADD LOOP: %d %d",iloop,signature); 
	      loop_print( &xfm_shape );**/
	  
	  orientation++;
	  if(orientation>MAX_ORIENTATION){
	    printf("%s: MAX_ORIENTATION too small\n",my_name);
	    terminate(1);
	  }
	}

	/* Update irrep weights = sum of characters */
	/* kclass is the class of the current O_h element */
	kclass = g[ielem].class;   
	/* Add character for each irrep to weight table */
	for(irrep = 0; irrep < N_CUBIC_IRREP; irrep++){
	  /* weight for real part is sum of charcters */
	  glueball_op[iloop].loop[save_orient].weight_r[irrep] += 
	    cubic_char[irrep][kclass];
	  /* weight for imaginary part is modified by sign for path reversal */
	  if(flag == -1)
	    glueball_op[iloop].loop[save_orient].weight_i[irrep] -= 
	      cubic_char[irrep][kclass];
	  else
	    glueball_op[iloop].loop[save_orient].weight_i[irrep] += 
	      cubic_char[irrep][kclass];
	}

	/* Update orientations */
	glueball_op[iloop].norientations = orientation;

      } /* group element */
    
  } /* end iloop */
  
  /* print out the glueball ops constructed here */

  node0_printf("\nGlueball operators: \n index orientation path\n");
  for(iloop = 0; iloop < NGLUEBALL_BASIC_SHAPES; iloop++)
    for(orientation = 0; 
	orientation < glueball_op[iloop].norientations; 
	orientation++){
      node0_printf(" %d         %d     ",
		   iloop,orientation);
      loop_print(&glueball_op[iloop].loop[orientation].path);
      node0_printf("                 ");
      for(j = 0; j < N_CUBIC_IRREP; j++){
	if(glueball_op[iloop].loop[orientation].weight_r[j] != 0)
	  node0_printf(" %s+ %d, ",cubic_irrep_name[j],
             glueball_op[iloop].loop[orientation].weight_r[j]);
	if(glueball_op[iloop].loop[orientation].weight_i[j] != 0)
	  node0_printf(" %s- %d, ",cubic_irrep_name[j],
             glueball_op[iloop].loop[orientation].weight_i[j]);
      }
      node0_printf("\n");
    }
  
  /* print table of cubic group irreps - needed for interpreting
     indices in glueball measurements below */

  node0_printf("\nCubic group irreps\n");
  for(irrep = 0; irrep < N_CUBIC_IRREP; irrep++){
    node0_printf("%3d %s\n",irrep,cubic_irrep_name[irrep]);
  }
  
} /* make_glueball_ops */


/*----------------------------------------------------------------------*/
/* Measure glueball operators and write results projected onto O_h
   irreps */
void measure_glueball_ops() {
    register int i;
    int irrep;
    register site *s;
    int length;
    int iloop,orientation,t,weight_r,weight_i,nonnull_r,nonnull_i;
    double_complex *glueball_irrep;
    complex trace;
    Real norm;

    glueball_irrep = (double_complex *)malloc(nt*sizeof(double_complex));
    if(glueball_irrep == NULL){
      printf("measure_glueball_op: No room for glueball_irrep\n");
      terminate(1);
    }

    /* traces of operators */
    for(iloop=0; iloop < NGLUEBALL_BASIC_SHAPES; iloop++)
      for(orientation = 0; 
	  orientation < glueball_op[iloop].norientations; 
	  orientation++){
	path_product( glueball_op[iloop].loop[orientation].path.dir,
		      glueball_op[iloop].loop[orientation].path.length);
	
	for(t = 0; t < nt; t++)
	  glueball_op[iloop].loop[orientation].trace[t] = dcmplx(0.,0.);
	
	FORALLSITES(i,s){
	  trace=trace_su3( &s->tempmat1 );
	  CSUM(glueball_op[iloop].loop[orientation].trace[s->t],trace);
	} END_LOOP /* sites */
    } /* iloop */
    
    /* Project all irrep contributions and write them out */
    for(irrep = 0; irrep < N_CUBIC_IRREP; irrep++){
      /* Run through all operators, projecting onto
	 the specified irrep */
      for(iloop=0; iloop < NGLUEBALL_BASIC_SHAPES; iloop++){
	for(t = 0; t < nt; t++)glueball_irrep[t] = dcmplx(0.,0.);
	nonnull_r = nonnull_i = 0;
	trace.real = trace.imag = 0.;
	/* Run through all orientations of the basic operator */
	for(orientation = 0; 
	    orientation < glueball_op[iloop].norientations; 
	    orientation++){
	  /* Operator weights for real and imag parts for this irrep */
	  weight_r = glueball_op[iloop].loop[orientation].weight_r[irrep];  
	  weight_i = glueball_op[iloop].loop[orientation].weight_i[irrep];  
	  if(weight_r != 0){
	    nonnull_r = 1;
	    for(t = 0; t < nt; t++){
	      trace.real = weight_r *
		glueball_op[iloop].loop[orientation].trace[t].real;
	      glueball_irrep[t].real += trace.real;
	    }
	  }
	  if(weight_i != 0){
	    nonnull_i = 1;
	    for(t = 0; t < nt; t++){
	      trace.imag = weight_i *
		glueball_op[iloop].loop[orientation].trace[t].imag;
	      glueball_irrep[t].imag += trace.imag;
	    }
	  }
	} /* orientation */

	/* Write result if nonnull */
	norm = 1./((Real)N_CUBIC_GROUP*nx*ny*nz);
	if(nonnull_r == 1 || nonnull_i == 1)
	  g_vecdcomplexsum(glueball_irrep,nt);
	
	if(nonnull_r == 1)
	  for(t = 0; t < nt; t++)
	    printf("GBALL %s+ %d %d %e\n",cubic_irrep_name[irrep],iloop,t,
		   norm*glueball_irrep[t].real);
	if(nonnull_i == 1)
	  for(t = 0; t < nt; t++)
	    printf("GBALL %s- %d %d %e\n",cubic_irrep_name[irrep],iloop,t,
		   norm*glueball_irrep[t].imag);
	
	for(t = 0; t < nt; t++)glueball_irrep[t] = dcmplx(0.,0.);
	nonnull_r = nonnull_i = 0;
	trace.real = trace.imag = 0.;
      } /* iloop */
    } /* irrep */
    
    free(glueball_irrep);
    
} /* measure_glueball_ops */



