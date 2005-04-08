/************************** monte_block_ape_b_many.c *******************************/
/* MIMD version 6 */
/* Determinisitc ape blocking  WITH printout*/
/* MIMD version 3 */
/* T. DeGrand May 97--
This version ultimately replaces link by an N-level APE-blocked link,
with blocking parameter alpha. It holds the final
link in blocked_link[4+dir] and uses blocked_link[8+dir] as temporary
storage


 */



#include "arb_dirac_eig_includes.h"


void monte_block_ape_b(int NumStp1)
{
int NumTrj,Nhit, index1, ina, inb,ii,cb;
int parity;
Real b3;
register int dir,i;
register site *st;
void dsdu_ape(register int dir1, int parity);
su3_matrix tmat1;
Real a0,a1,a2,a3,asq;
int index,ind1,ind2,step;
su2_matrix h;
Real alpha;

int NumStp;
NumStp=10 ;

Nhit=5;

alpha=0.45;
b3=alpha/(1.0-alpha)/6.0;
b3=1.0/b3;

if(this_node==0)printf("pure APE blocking with alpha  %e N %d\n",alpha,NumStp);

/*  set bb_link=link */
 
                FORALLSITES(i,st)for(dir=XUP;dir<=TUP;dir++){
                st->blocked_link[8+dir]= st->link[dir];
               }
/* ape blocking steps_rg levels*/
 
        for(step=1;step<=NumStp;step++){


        for(parity=ODD;parity<=EVEN;parity++)
	for(dir=XUP;dir<=TUP;dir++){

              /* compute the gauge force */
		dsdu_ape(dir,parity); 

        FORSOMEPARITY(i,st,parity){
	{
/* set blocked_link=bb_link temporarily*/
                 st->blocked_link[4+dir]= st->blocked_link[8+dir];
/* add the staple to the blocked link. ``staple'' will become the new
blocked_link after normalization */
                scalar_mult_add_su3_matrix(&(st->tempmat2),
	 &(st->blocked_link[8+dir]),b3,&(st->tempmat2));
  
/* if(i==0&&step==1){
	printf("\n\n step=%d i=%d dir=%d\n",step,i,dir);
	dumpmat(&(st->blocked_link[8+dir]));
	dumpmat(&(st->tempmat2));
}*/

        /* Now do hits in the SU(2) subgroup to "normalize" staple */

        for(index=0;index<3*Nhit;index++){

             /*  pick out an SU(2) subgroup */
                        ind1=(index) % 3;
                        ind2=(index+1) % 3;

                        if(ind1 > ind2){ ii=ind1; ind1=ind2; ind2=ii;}



                mult_su3_na( &(st->blocked_link[4+dir]), &(st->tempmat2), &tmat1 );

                /* Extract SU(2) subgroup in Pauli matrix representation,
                   a0 + i * sum_j a_j sigma_j, from the SU(3) matrix tmat1 */
                a0 = tmat1.e[ind1][ind1].real + tmat1.e[ind2][ind2].real;
                a1 = tmat1.e[ind1][ind2].imag + tmat1.e[ind2][ind1].imag;
                a2 = tmat1.e[ind1][ind2].real - tmat1.e[ind2][ind1].real;
                a3 = tmat1.e[ind1][ind1].imag - tmat1.e[ind2][ind2].imag;

                /* Normalize and put complex conjugate into u */
                asq = a0*a0 + a1*a1 + a2*a2 + a3*a3;
                asq = sqrt((double)asq);
                a0 = a0/asq; a1 = a1/asq; a2 = a2/asq; a3 = a3/asq;
                h.e[0][0] = cmplx( a0,-a3);
                h.e[0][1] = cmplx(-a2,-a1);
                h.e[1][0] = cmplx( a2,-a1);
                h.e[1][1] = cmplx( a0, a3);

                /* Do the SU(2) hit */
                left_su2_hit_n( &h, ind1, ind2, &(st->blocked_link[4+dir]));
		} /* indices */
            } /* end loop over sites */

	}} /*  direction and parity */

/* if you're not at the end copy b_link into bb_link */
        if(step != NumStp)
                FORALLSITES(i,st)for(dir=XUP;dir<=TUP;dir++){
                  st->blocked_link[8+dir]= st->blocked_link[4+dir];
                        }


        }/* step  */




                FORALLSITES(i,st)for(dir=XUP;dir<=TUP;dir++){
                st->link[dir] = st->blocked_link[4+dir];
/* printf("\n %d %d %d %d     %d\n",st->x,st->y,st->z,st->t,dir);
	dumpmat(&(st->link[dir])); */
                }


} /* monte */




/* dsdu_ape_m.c  -- compute the staple, using bb_link for the staple **/

/* MIMD version 3 */

void dsdu_ape(register int dir1, int parity) 
{
register int i,dir2;
register site *st;
msg_tag *tag0,*tag1,*tag2,*tag3,*tag4;
int start;
su3_matrix tmat1,tmat2;
int disp[4];	/* displacement vector for general gather */
	/* Loop over other directions, computing force from plaquettes in
	   the dir1,dir2 plane */
	start=1; /* indicates staple sum not initialized */
	for(dir2=XUP;dir2<=TUP;dir2++)if(dir2 != dir1)
	{
	    /* displacement vector for bb_link 2 sites away */
	    for(i=XUP;i<=TUP;i++)disp[i]=0;
	    disp[dir1] = 1;
	    disp[dir2] = -1;

	    /* get blocked_link[8+dir2] from direction dir1 */
	    tag0 = start_gather_site( F_OFFSET(blocked_link[8+dir2]), sizeof(su3_matrix),
		dir1, parity, gen_pt[0] );

	    /* get blocked_link[8+dir1] from direction dir2 */
	    tag1 = start_gather_site( F_OFFSET(blocked_link[8+dir1]), sizeof(su3_matrix),
		dir2, parity, gen_pt[1] );

	    /* get blocked_link[8+dir2] from direction -dir2 */
	    tag2 = start_gather_site( F_OFFSET(blocked_link[8+dir2]), sizeof(su3_matrix),
		OPP_DIR(dir2), parity, gen_pt[2] );

	    /* get blocked_link[8+dir1] from direction -dir2 */
	    tag3 = start_gather_site( F_OFFSET(blocked_link[8+dir1]), sizeof(su3_matrix),
		OPP_DIR(dir2), parity, gen_pt[3] );

	    /* get blocked_link[8+dir2] from displacement +dir1-dir2 */
	    tag4 = start_general_gather_site( F_OFFSET(blocked_link[8+dir2]),
		sizeof(su3_matrix), disp, parity, gen_pt[4] );

	    /* Upper staple */
	    wait_gather(tag0);
	    wait_gather(tag1);
          if(start){  /* this is the first contribution to staple */
	FORSOMEPARITY(i,st,parity){
	        mult_su3_nn( &(st->blocked_link[8+dir2]), (su3_matrix *)gen_pt[1][i], &tmat1 );
		mult_su3_na( &tmat1, (su3_matrix *)gen_pt[0][i], &(st->tempmat2) );

		}
		start=0; 
		}
		else{
	FORSOMEPARITY(i,st,parity){
		mult_su3_nn( &(st->blocked_link[8+dir2]), (su3_matrix *)gen_pt[1][i], &tmat1 );
		mult_su3_na( &tmat1, (su3_matrix *)gen_pt[0][i], &tmat2 );
		add_su3_matrix( &(st->tempmat2), &tmat2, &(st->tempmat2));
		}
	    } /* upper tempmat2 */
	    cleanup_gather(tag0);
	    cleanup_gather(tag1);

	    /* Lower tempmat2 */
	    wait_gather(tag2);
	    wait_gather(tag3);
	    wait_general_gather(tag4);
	FORSOMEPARITY(i,st,parity){
	        mult_su3_an( (su3_matrix *)gen_pt[2][i], (su3_matrix *)gen_pt[3][i], &tmat1 );
	        mult_su3_nn( &tmat1, (su3_matrix *)gen_pt[4][i], &tmat2 );
		add_su3_matrix( &(st->tempmat2), &tmat2, &(st->tempmat2));
	    }  /* lower tempmat2 */
	    cleanup_gather(tag2);
	    cleanup_gather(tag3);
	    cleanup_general_gather(tag4);
	}
}
