/* setup_offset.c  **/

/*  label[k] labels the offset of the kth connection between psibar and psi
by an integer 1 to 4, the separation of the paths in ``new york metric.''
It corresponds to our labeling of rho and lambda, lambda_1 = lambda(1,0,0,0).
offset[k] is a four component array with the actual offset of the kth
connection.  This program only generates the ``positive offsets'', a set of
40 of the 80 offsets. The paths for the others can be found by using
adjoints of these paths (i.e., -1-1-1-1 is the adjoint of 1111)

*/

#include "arb_dirac_eig_includes.h"

void setup_offset()
{
int j,n[4];
int k,ii,flag,mu;

/* construct the first list */
if(this_node==0)printf("only enough offsets for planar action\n");
off_max=0;
for(n[0]=1; n[0]>= -1; n[0]--)for(n[1]=1; n[1]>= -1; n[1]--)
for(n[2]=1; n[2]>= -1; n[2]--)for(n[3]=1; n[3]>= -1; n[3]--){

k=abs(n[0])+abs(n[1])+abs(n[2])+abs(n[3]);
	if(k==1 || k==2){
/* run the list of previous entries and see if the new one is here, reversed */
		flag=0;
		for(ii=0;ii<off_max;ii++){
		if(
		n[0] == -offset[ii][0]  &&
		n[1] == -offset[ii][1]  &&
		n[2] == -offset[ii][2]  &&
		n[3] == -offset[ii][3]  ) flag=1;
		}
			if(flag==0){
			for(mu=0;mu<4;mu++) offset[off_max][mu]=n[mu];
			label[off_max]=k;
			off_max++;

			}

	}
}
/*
printf(" there are %d distinct paths\n",off_max);
for(ii=0;ii<off_max;ii++){
        printf("     %d: ",label[ii]);
	for(mu=0;mu<4;mu++) printf("%d ",offset[ii][mu]);
	printf("\n");
}
*/
}
