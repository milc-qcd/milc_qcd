/************************** d_plaq4.c *******************************/
/* MIMD version 7 */
/* This version mallocs the temporary su3_matrix */

/* AB 11/29/7 - dump all plaquettes measured on fat links on disk */

#include "generic_includes.h"

void d_plaquette_field_dump(su3_matrix **U_field, char *file_name_prefix) {
/* su3mat is scratch space of size su3_matrix */
su3_matrix *su3mat;
register int i,dir1,dir2;
register site *s;
register su3_matrix *m1,*m4;
su3_matrix mtmp;
double ss_sum,st_sum;
double rtrace, rtrace3;
msg_tag *mtag0,*mtag1;
    ss_sum = st_sum = 0.0;
FILE *fp;
char plaq_file_name[300];

    su3mat = (su3_matrix *)malloc(sizeof(su3_matrix)*sites_on_node);
    if(su3mat == NULL)
      {
	printf("plaquette: can't malloc su3mat\n");
	fflush(stdout); terminate(1);
      }

    sprintf( plaq_file_name, "%s_node%04d.dat", file_name_prefix, this_node );
    fp = fopen( plaq_file_name, "wt" );


    for(dir1=YUP;dir1<=TUP;dir1++){
	for(dir2=XUP;dir2<dir1;dir2++){

	    mtag0 = start_gather_field( U_field[dir2], sizeof(su3_matrix),
		dir1, EVENANDODD, gen_pt[0] );
	    mtag1 = start_gather_field( U_field[dir1], sizeof(su3_matrix),
		dir2, EVENANDODD, gen_pt[1] );

	    FORALLSITES(i,s){
		m1 = &(U_field[dir1][i]);
		m4 = &(U_field[dir2][i]);
		mult_su3_an(m4,m1,&su3mat[i]);
	    }

	    wait_gather(mtag0);
	    wait_gather(mtag1);

	    FORALLSITES(i,s){
		mult_su3_nn( &su3mat[i], (su3_matrix *)(gen_pt[0][i]),
		    &mtmp);

		if(dir1==TUP ) {
                  rtrace = (double)
		    realtrace_su3((su3_matrix *)(gen_pt[1][i]),&mtmp);
                  st_sum += rtrace;
                }
		else {
                  rtrace = (double)
		    realtrace_su3((su3_matrix *)(gen_pt[1][i]),&mtmp);
                  ss_sum += rtrace;
                }

                /* output plaquette into file */
                fprintf( fp, "%18.12g\n", rtrace );
            }

            cleanup_gather(mtag0);
	    cleanup_gather(mtag1);
	}
    }


    fclose( fp );

    free(su3mat);
} /* d_plaquette4 */

