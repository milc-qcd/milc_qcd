/*************** setup_offset.c ****************/
/* MIMD version 7 */

/*  label[k] labels the offset of the kth connection between psibar and psi
by an integer 1 to 4, the separation of the paths in ``new york metric.''
It corresponds to our labeling of rho and lambda, lambda_1 = lambda(1,0,0,0).
offset[k] is a four component array with the actual offset of the kth
connection.  This program only generates the ``positive offsets'', a set of
40 of the 80 offsets. The paths for the others can be found by using
adjoints of these paths (i.e., -1-1-1-1 is the adjoint of 1111)

*/

#include "arb_ov_includes.h"

void setup_offset()
{
int n[4];
int k,ii,flag,mu;
int jj, dir1, dir2;

/* for gather tables */
void cubic_neighbor(int x,int y,int z,int t,int *arg, int forw_back,
int *xpt,int *ypt,int *zpt,int *tpt);
int ipath;

/* construct the first list */
    if(this_node==0)printf("only enough offsets for planar action\n");
    off_max=0;
    for(n[0]=1; n[0]>= -1; n[0]--)for(n[1]=1; n[1]>= -1; n[1]--)
    for(n[2]=1; n[2]>= -1; n[2]--)for(n[3]=1; n[3]>= -1; n[3]--)
    {

        k=abs(n[0])+abs(n[1])+abs(n[2])+abs(n[3]);
	if(k==1 || k==2){
/* run the list of previous entries and see if the new one is here, reversed */
	    flag=0;
		
	    for(ii=0;ii<off_max;ii++){
		if(n[0] == -offset[ii][0]  &&
		   n[1] == -offset[ii][1]  &&
		   n[2] == -offset[ii][2]  &&
		   n[3] == -offset[ii][3]  ) flag=1;
	    }
			
	    if(flag==0){
		for(mu=0;mu<4;mu++) offset[off_max][mu]=n[mu];
		label[off_max]=k;
		
/* fill translation table for paths with length 2,
 * which contains elementary steps -> offsets 
 */ 
		if (k==1){
		    for (ii=0;ii<4;ii++)
			if (n[ii]!=0) break;
		    if (n[ii]==1) dir1=ii;
		    else dir1=OPP_DIR(ii);
		    
		    diroff1[dir1]=off_max;
		    diroff1[OPP_DIR(dir1)]=off_max;
		    signoff1[dir1]=+1;
		    signoff1[OPP_DIR(dir1)]=-1;
		}
		
		if (k==2){
		    for (ii=0;ii<4;ii++)
			if (n[ii]!=0) break;
		    if (n[ii]==1) dir1=ii;
		    else dir1=OPP_DIR(ii);
		    
		    for (jj=ii+1;jj<4;jj++)
			if (n[jj]!=0) break;
		    if (n[jj]==1) dir2=jj;
		    else dir2=OPP_DIR(jj);

		    diroff[dir1][dir2]=off_max;
		    diroff[OPP_DIR(dir1)][OPP_DIR(dir2)]=off_max;
		    signoff[dir1][dir2]=+1;
		    signoff[OPP_DIR(dir1)][OPP_DIR(dir2)]=-1;

		    diroff[dir2][dir1]=off_max;
		    diroff[OPP_DIR(dir2)][OPP_DIR(dir1)]=off_max;
		    signoff[dir2][dir1]=+1;
		    signoff[OPP_DIR(dir2)][OPP_DIR(dir1)]=-1;
		    /*
		    printf("so %i %i %i %i %i %i %i %i \n",dir1,dir2,n[0],n[1],n[2],n[3],diroff[dir1][dir2],diroff[OPP_DIR(dir1)][OPP_DIR(dir2)]);
		    */
		}





		
		off_max++;
	    }

	
	}/* k == 1 or 2 */

    }
/*
printf(" there are %d distinct paths\n",off_max);
for(ii=0;ii<off_max;ii++){
        printf("     %d: ",label[ii]);
	for(mu=0;mu<4;mu++) printf("%d ",offset[ii][mu]);
	printf("\n");
}
*/
/* now make the gather tables */

for(ipath=0;ipath<off_max;ipath++){
  goffset[ipath]=make_gather(cubic_neighbor,offset[ipath],
			 WANT_INVERSE,NO_EVEN_ODD,SCRAMBLE_PARITY);


/*
        node0_printf(" ahead    %d: ",ipath);
	for(mu=0;mu<4;mu++) node0_printf("%d ",offset[ipath][mu]);
        node0_printf("     %d: ",goffset[ipath]);
	node0_printf("\n");
*/




}
}
void cubic_neighbor(int x,int y,int z,int t,int *arg, int forw_back,
int *xpt,int *ypt,int *zpt,int *tpt)
{
  if(forw_back==FORWARDS){
*xpt = (x+nx+arg[0])%nx;
*ypt = (y+ny+arg[1])%ny;
*zpt = (z+nz+arg[2])%nz;
*tpt = (t+nt+arg[3])%nt;
  }
else
{
*xpt = (x+nx-arg[0])%nx;
*ypt = (y+ny-arg[1])%ny;
*zpt = (z+nz-arg[2])%nz;
*tpt = (t+nt-arg[3])%nt;
  }
}
