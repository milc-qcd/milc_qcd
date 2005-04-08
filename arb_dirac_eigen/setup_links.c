/************************** setup_links.c *******************************/
/*
 working version: same paths for rho and lambda terms,
paths hardwired in here. The (unused) 'iflag' is present in case you
want several options for paths.
*/

#include "arb_dirac_eig_includes.h"

void setup_links(int iflag)
{
register int i;
register site *s;
int n[4],k,l1,l2,index;
int ioffset,ipath;

/* data for paths: [k][npath[k]][along each path] */
int dir[5][24][4],sign[5][24][4],length[5][24],npath[5];
Real cc[5][24];

msg_tag *tag;


/* we organize the calculation by the ``label'' variable,equal to
an integer 1 to 4, which identifies the separation of psi and psibar */
        for(ioffset=0;ioffset<off_max;ioffset++){
                k=label[ioffset];


	path_stuff(iflag,ioffset,k,npath,dir,sign,length,cc);


/* this fragment loops over all the paths. It constructs the product of
links along the ith path, weights it by the factor cc[k][ipath],
and adds it to the total, which is then put in blocked_link[ioffset]. */	

		FORALLSITES(i,s){
		  for(l1=0;l1<3;l1++)for(l2=0;l2<3;l2++){
		  s->staple.e[l1][l2].real=0.0;
		  s->staple.e[l1][l2].imag=0.0;}
		  }

		  for(ipath=0;ipath<npath[k];ipath++){
		  path(dir[k][ipath],sign[k][ipath],length[k][ipath]);


			FORALLSITES(i,s){
			scalar_mult_add_su3_matrix(&(s->staple),
				&(s->tempmat1),cc[k][ipath],&(s->staple));
			}
		}

/* path.c puts the end of the path at location ``zero.''
We must gather so that the beginning of the path is at zero */
		tag=start_general_gather_site(F_OFFSET(staple),sizeof(su3_matrix),
                offset[ioffset], EVENANDODD, gen_pt[0] );

		wait_general_gather(tag);

		FORALLSITES(i,s){
			su3mat_copy((su3_matrix *)(gen_pt[0][i]),
				&(s->blocked_link[ioffset]));
/* printf("%d %d %d %d  : %d\n",s->x,s->y,s->z,s->t,ioffset); 
dumpmat(&(s->blocked_link[ioffset])); */
		}
		cleanup_general_gather(tag);
	}

}

/* here we put all the grubby path-generation stuff */

void path_stuff(int iflag, int ioffset, int k, int npath[5],
 int dir[5][24][4], int sign[5][24][4], int length[5][24],
Real cc[5][24])
{
int l1,l2,l3;
int perm[4];
Real sum;


/*  a simple set of paths for debugging: npath=1 (for all paths), 
and the length of the path is equal to k (i.e. one minimum length path)


	npath[k]=1;
	length[k][0]=k;
	cc[k][0]=1.0;
	l1=0;
	for(l2=0;l2<4;l2++)if(offset[ioffset][l2]!=0){
	dir[k][0][l1]= abs(l2);
	sign[k][0][l1]= offset[ioffset][l2];
	l1++;
	}

*/

/* a complicated set of paths (average over all paths in an offset
over the hypercube) */

if(k==1){


	npath[k]=1;
        length[k][0]=1;
	cc[k][0]=1.0;
	l1=0;
        for(l2=0;l2<4;l2++)if(offset[ioffset][l2]!=0){
        dir[k][0][l1]= abs(l2);
        sign[k][0][l1]= offset[ioffset][l2];
        l1++;
        }

}

if(k==2){
	npath[k]=2;
	for(l1=0;l1<2;l1++) {length[k][l1]=2;cc[k][l1]=0.5;}
	l1=0;
	for(l2=0;l2<4;l2++)if(offset[ioffset][l2]!=0){
	dir[k][0][l1]= abs(l2);
	sign[k][0][l1]= offset[ioffset][l2];
	l1++;
	}
	dir[k][1][0]=dir[k][0][1];
	dir[k][1][1]=dir[k][0][0];
	sign[k][1][0]=sign[k][0][1];
	sign[k][1][1]=sign[k][0][0];

}


if(k==3){
	for(l1=0;l1<6;l1++) {length[k][l1]=3;cc[k][l1]=1.0/6.0;}
	l1=0;
	for(l2=0;l2<4;l2++)if(offset[ioffset][l2]!=0){
	dir[k][0][l1]= abs(l2);
	sign[k][0][l1]= offset[ioffset][l2];
	l1++;
	}


	l1=0;
  for(perm[0]=0;perm[0]<3;perm[0]++)
  for(perm[1]=0;perm[1]<3;perm[1]++)
  for(perm[2]=0;perm[2]<3;perm[2]++)
     if(perm[0] != perm[1] && perm[0] != perm[2] 
         && perm[1] != perm[2]  ) {
	for(l2=0;l2<3;l2++){
	  dir[k][l1][l2]=dir[k][0][perm[l2]] ;
	  sign[k][l1][l2]=sign[k][0][perm[l2]] ;
	}
  l1++;
  }
	npath[k]=l1;

}


if(k==4){
	npath[k]=24;
	for(l1=0;l1<24;l1++) {length[k][l1]=4;cc[k][l1]=1.0/24.0;}
	l1=0;
	for(l2=0;l2<4;l2++)if(offset[ioffset][l2]!=0){
	dir[k][0][l1]= abs(l2);
	sign[k][0][l1]= offset[ioffset][l2];
	l1++;
	}

	l1=0;
  for(perm[0]=0;perm[0]<4;perm[0]++)
  for(perm[1]=0;perm[1]<4;perm[1]++)
  for(perm[2]=0;perm[2]<4;perm[2]++)
  for(perm[3]=0;perm[3]<4;perm[3]++)
     if(perm[0] != perm[1] && perm[0] != perm[2] && perm[0] != perm[3]
         && perm[1] != perm[2] && perm[1] != perm[3] && 
            perm[2] != perm[3] ) {
	for(l2=0;l2<4;l2++){
	  dir[k][l1][l2]=dir[k][0][perm[l2]] ;
	  sign[k][l1][l2]=sign[k][0][perm[l2]] ;
	}
  l1++;
  }
	npath[k]=l1;
}

/*
	printf("ioffset: %d: ",ioffset);
	for(l2=0;l2<4;l2++)printf("%d ",offset[ioffset][l2]);
	printf(": ");
	printf("%d paths:",npath[k]);
	sum=0.0;
	for(l2=0;l2<npath[k];l2++) sum += cc[k][l2];
	printf(" norm= %e\n",sum);
*/

/*
	for(l2=0;l2<npath[k];l2++){
	printf("		    ");
	  for(l1=0;l1<length[k][l2];l1++)
		printf(" %d",sign[k][l2][l1]*dir[k][l2][l1]);
	  printf(":   %e\n",cc[k][l2]);
	}
*/

}
