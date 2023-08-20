/************************* fmunu.c *******************************/
/* Calculates the temporal and spatial field strength components */
/* and the topological charge.                                   */

/* Includes */
#include "wilson_flow_includes.h"
#include "../include/field_strength.h"

/* Computes the real trace of the su3 matrix product: ReTr(A*B) */
Real
real_trace_nn( su3_matrix *a, su3_matrix *b )
{
  register int i,j;
  register complex x;
  register Real sum;

  sum = 0.0;
  for( i=0; i<3; i++ ) for( j=0; j<3; j++ ) {
    CMUL( a->e[i][j], b->e[j][i], x );
    sum += x.real;
  }

  return sum;
}

/* Computes the field strength components and topological charge */
void
fmunu_fmunu( double *time, double *space, double *charge )
{
  /* Site variables */
  register int i;
  register site *s;

  /* Temporary component storage */
  su3_matrix *ft, *fs;

  /* Initialize sums */
  *time = *space = *charge = 0;

  /* Compute 8*F_mu,nu at each site */
  make_field_strength( F_OFFSET(link), F_OFFSET(fieldstrength) );

  /* Loop over each site to sum F_mu,nu components */
  FORALLSITES(i, s) {

    fs = &(s->fieldstrength[FS_XY]);
    ft = &(s->fieldstrength[FS_ZT]);
    *time -= real_trace_nn(ft, ft);
    *space -= real_trace_nn(fs, fs);
    *charge -= real_trace_nn(fs, ft);

    fs = &(s->fieldstrength[FS_XZ]);
    ft = &(s->fieldstrength[FS_YT]);
    *time -= real_trace_nn(ft, ft);
    *space -= real_trace_nn(fs, fs);
    *charge -= realtrace_su3(fs, ft); /* ReTr{ fs.dag * ft } */

    fs = &(s->fieldstrength[FS_YZ]);
    ft = &(s->fieldstrength[FS_XT]);
    *time -= real_trace_nn(ft, ft);
    *space -= real_trace_nn(fs, fs);
    *charge -= real_trace_nn(fs, ft);
  }

  /* Sum over all nodes */
  g_doublesum(time);
  g_doublesum(space);
  g_doublesum(charge);

  /* Normalizations */
  *time /= (volume*64.0);
  *space /= (volume*64.0);
  *charge *= 0.0003957858736028819197; /* normalization of 1/(8^2 * 4 * PI^2) */
}


/* Compute loops: 1x1 -- plaquette, 1x2 + 2x1 -- rectangle
   for Wilson (one-plaquette) and
   Symanzik tree-level (plaquette and rectangle) action,
   temporal and spatial part separately */
void
gauge_action_w_s( double *wl1x1t, double *wl1x1s,
                  double *wl1x2t, double *wl1x2s ) {

#define NTEMP_STORAGE 6
  register int i, dir1, dir2;
  register site *s;
  msg_tag *tag0, *tag1, *tag2, *tag3, *tag4, *tag5, *tag6, *tag7, *tag8;
  su3_matrix *su3mat[NTEMP_STORAGE];
  su3_matrix tempmat;
  double tt;

  for( i=0; i<NTEMP_STORAGE; i++ ) {
    su3mat[i] = (su3_matrix *)malloc( sizeof(su3_matrix)*sites_on_node );
    if(su3mat[i] == NULL) {
        printf( "gauge_action_w_s: can't malloc su3mat[%d]\n", i );
        fflush(stdout); terminate(1);
    }
  }

  // prepare accumulators
  *wl1x1s = 0;
  *wl1x1t = 0;
  *wl1x2s = 0;
  *wl1x2t = 0;

  for( dir1=YUP; dir1<=TUP; dir1++ ) {
    for( dir2=XUP; dir2<dir1; dir2++ ) {
      // request link[dir2] from direction dir1
      tag0 = start_gather_site( F_OFFSET(link[dir2]) , sizeof(su3_matrix),
                                dir1, EVENANDODD, gen_pt[0]);

      // request link[dir1] from direction dir2
      tag1 = start_gather_site( F_OFFSET(link[dir1]), sizeof(su3_matrix),
                                dir2, EVENANDODD, gen_pt[1]);

      // while waiting for gathers multiply two two links on the site
      FORALLSITES(i, s)
        mult_su3_an( &(s->link[dir2]), &(s->link[dir1]), &(su3mat[0][i]) );

      wait_gather(tag0);

      // form a staple for "right" link in the plaquette
      FORALLSITES(i, s)
        mult_su3_nn( &(su3mat[0][i]), (su3_matrix *)(gen_pt[0][i]), &(su3mat[3][i]) );

      wait_gather(tag1);

//fflush(stdout);printf("dir1=%d dir2=%d\n",dir1,dir2);fflush(stdout);
      // form a staple for "left" link in the plaquette
      FORALLSITES(i, s) {
        mult_su3_nn( &(s->link[dir2]), (su3_matrix *)(gen_pt[1][i]), &tempmat );
        mult_su3_na( &tempmat, (su3_matrix *)(gen_pt[0][i]), &(su3mat[5][i]) );
      }

      // request staple for "left" = su3mat[5] from dir2
      tag2 = start_gather_field( su3mat[5], sizeof(su3_matrix),
                                 dir2, EVENANDODD, gen_pt[2]);

//fflush(stdout);printf("dir1=%d dir2=%d\n",dir1,dir2);fflush(stdout);
      // form a staple for "bottom" link in the plaquette
      FORALLSITES(i, s) {
        mult_su3_na( (su3_matrix *)(gen_pt[1][i]), (su3_matrix *)(gen_pt[0][i]), &(su3mat[1][i]) );
        mult_su3_na( &(su3mat[1][i]), &(s->link[dir1]), &(su3mat[4][i]) );
      }

      // request staple for "bottom" = su3mat[4] from dir1
      tag3 = start_gather_field( su3mat[4], sizeof(su3_matrix),
                                 dir1, EVENANDODD, gen_pt[3]);

      wait_gather(tag2);

      FORALLSITES(i, s) {
        // form a staple for "top" link in the plaquette
        mult_su3_an( (su3_matrix *)(gen_pt[1][i]), &(su3mat[0][i]), &(su3mat[2][i]) );
        // get the contribution of 1x2 rectangle extended in dir2
        // to the accumulator
        if( dir1==TUP ) *wl1x2t +=
            realtrace_su3( (su3_matrix *)(gen_pt[2][i]), &(su3mat[3][i]) );
        else            *wl1x2s +=
            realtrace_su3( (su3_matrix *)(gen_pt[2][i]), &(su3mat[3][i]) );
      }

      wait_gather(tag3);

      FORALLSITES(i, s) {
        // get the contribution of 1x2 rectangle extended in dir1
        // and of 1x1 plaquette to the accumulators
        if( dir1==TUP ) {
          *wl1x2t +=
            realtrace_su3( (su3_matrix *)(gen_pt[3][i]), &(su3mat[2][i]) );
          *wl1x1t +=
            realtrace_su3( &(su3mat[1][i]), &(su3mat[0][i]) );
        }
        else {
          *wl1x2s +=
            realtrace_su3( (su3_matrix *)(gen_pt[3][i]), &(su3mat[2][i]) );
          *wl1x1s +=
            realtrace_su3( &(su3mat[1][i]), &(su3mat[0][i]) );
        }
      }

      // clean up all gathers
      cleanup_gather(tag0);
      cleanup_gather(tag1);
      cleanup_gather(tag2);
      cleanup_gather(tag3);

    } // dir2
  } // dir1

  // global sum
  g_doublesum( wl1x1s );
  g_doublesum( wl1x1t );
  g_doublesum( wl1x2s );
  g_doublesum( wl1x2t );

  // get densities
  *wl1x1s /= volume;
  *wl1x1s /= 3;
  *wl1x1t /= volume;
  *wl1x1t /= 3;
  *wl1x2s /= volume;
  *wl1x2s /= 6;
  *wl1x2t /= volume;
  *wl1x2t /= 6;

  // deallocate temporary storage
  for( i=0; i<NTEMP_STORAGE; i++ ) {
    free( su3mat[i] );
  }

#undef NTEMP_STORAGE
}
