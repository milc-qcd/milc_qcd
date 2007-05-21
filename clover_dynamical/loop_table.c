/************************** loop_table.c *******************************/
/* MIMD version 7 */
/* WE NEED TO NORMALIZE THIS PROCEDURE WITH generic/gauge_stuff.c -CD! */

#include "cl_dyn_includes.h"

void make_loop_table2()
{
/*
static int loop_ind[nloop][10] = {
{0,1,7,6,0,0,0,0,0,0},
{0,1,1,7,6,6,0,0,0,0},
{0,1,2,6,7,5,0,0,0,0},
{0,1,2,7,6,5,0,0,0,0},
{  0,  1,  0,  1,  7,  6,  7,  6,0 ,0},
{  0,  1,  0,  1,  7,  7,  6,  6,0 ,0},
{  0,  1,  0,  2,  6,  7,  7,  5,0 ,0},
{  0,  1,  0,  2,  7,  5,  7,  6,0 ,0},
{  0,  1,  0,  2,  7,  6,  7,  5,0 ,0},
{  0,  1,  0,  2,  7,  7,  5,  6,0 ,0},
{  0,  1,  0,  6,  7,  1,  7,  6,0 ,0},
{  0,  1,  1,  0,  6,  7,  7,  6,0 ,0},
{  0,  1,  1,  2,  6,  5,  7,  6,0 ,0},
{  0,  1,  1,  2,  6,  6,  7,  5,0 ,0},
{  0,  1,  1,  2,  6,  7,  5,  6,0 ,0},
{  0,  1,  1,  2,  7,  6,  6,  5,0 ,0},
{  0,  1,  1,  7,  2,  6,  6,  5,0 ,0},
{  0,  1,  1,  7,  7,  6,  6,  0,0 ,0},
{  0,  1,  2,  3,  5,  4,  7,  6,0 ,0},
{  0,  1,  2,  3,  5,  6,  7,  4,0 ,0},
{  0,  1,  2,  3,  5,  7,  4,  6,0 ,0},
{  0,  1,  2,  3,  5,  7,  6,  4,0 ,0},
{  0,  1,  2,  3,  6,  7,  4,  5,0 ,0},
{  0,  1,  2,  3,  6,  7,  5,  4,0 ,0},
{  0,  1,  2,  3,  7,  6,  5,  4,0 ,0},
{  0,  1,  2,  6,  5,  1,  7,  6,0 ,0},
{  0,  1,  2,  6,  7,  1,  5,  6,0 ,0},
{  0,  1,  7,  6,  0,  1,  7,  6,0 ,0}
};
static int loop_length_in[nloop] = {4,6,6,6,
                       8,8,8,8,8,8,8,8,8,8,
                       8,8,8,8,8,8,8,8,8,8,
                       8,8,8,8};
*/

int perm[8],pp[8],ir[4];
int length,iloop,i,j,chr;
int vec[10];
int count,flag;
void char_num(int dig[],int *chr,int *ch,int length);
void make_loop_term();



/*
The Plaquette 
static int loop_ind[nloop][10] = {
{0,1,7,6,0,0,0,0,0,0}
};
static int loop_length_in[nloop] = {4};
*/

/* Plaquette plus O_3 promoted to O_1
static int loop_ind[nloop][10] = {
{0,1,7,6,0,0,0,0,0,0},
{0,1,2,7,6,5,0,0,0,0},
};
static int loop_length_in[nloop] = {4,6};
if(this_node == 0) printf("\n\nNote O_1 is really O_3 in the notes! \n\n\n");
*/


/* Plaquette plus O_1 plus O_3 promoted to O_2 */
static int loop_ind[nloop][10] = {
{0,1,7,6,0,0,0,0,0,0},
{0,1,1,7,6,6,0,0,0,0},
{0,1,2,7,6,5,0,0,0,0},
};
static int loop_length_in[nloop] = {4,6,6};
if(this_node == 0) printf("\n\nNote O_2 is really O_3 in the notes! \n\n\n");


for(j=0;j<nloop;j++){loop_num[j]=0;loop_length[j]=loop_length_in[j];}

/* set up the loop coefficients */

for(i=0;i<nreps;i++)for(j=0;j<nloop;j++)loop_coeff[j][i]=0.0;


/* the plaquette 
loop_coeff[0][0]= 1.0;
*/


/* loop action 9 Sept  

loop_coeff[0][0]= 5.230000e-01;
loop_coeff[1][0]= 5.970000e-02;
loop_coeff[0][1]= 2.109197e-03;
loop_coeff[1][1]= 5.408617e-03;
loop_coeff[0][2]= 5.341197e-03;
loop_coeff[1][2]= 5.091361e-03;
loop_coeff[0][3]= 1.668338e-02;
loop_coeff[1][3]= -6.360714e-04;
*/

/* a**2 improved Symanzik  with  funny couplings */
loop_coeff[0][0]= 1.0000000;
loop_coeff[1][0]= -1.00/(20.0*u0*u0) * (1.00 - 0.6264*log(u0) );
loop_coeff[2][0]= 1.00/(u0*u0) * 0.04335 * log(u0);


for(iloop=0;iloop<nloop;iloop++){
  length=loop_length[iloop];

  count=0;
/* permutations */
  for(perm[0]=0;perm[0]<4;perm[0]++)
  for(perm[1]=0;perm[1]<4;perm[1]++)
  for(perm[2]=0;perm[2]<4;perm[2]++)
  for(perm[3]=0;perm[3]<4;perm[3]++)
     {if(perm[0] != perm[1] && perm[0] != perm[2] && perm[0] != perm[3]
         && perm[1] != perm[2] && perm[1] != perm[3] && 
            perm[2] != perm[3] ) {
/* reflections*/
      for(ir[0]=0;ir[0]<2;ir[0]++)
      for(ir[1]=0;ir[1]<2;ir[1]++)
      for(ir[2]=0;ir[2]<2;ir[2]++)
      for(ir[3]=0;ir[3]<2;ir[3]++){
         for(j=0;j<4;j++){ pp[j]=perm[j];
                         {if(ir[j] == 1) pp[j]=7-pp[j];}
                           pp[7-j]=7-pp[j];}
/* create new vector*/
         for(j=0;j<length;j++) vec[j]=pp[loop_ind[iloop][j]];
                
         char_num(vec,&chr,&ch,length);
         flag=0;
/* check if it's a new set: */
         for(j=0;j<count;j++) if(chr == loop_char[j])flag=1;
         if(flag == 0 ){loop_char[count]=chr;
                        loop_ch[iloop][count]=ch;
                        for(j=0;j<length;j++)
                          loop_table[iloop][count][j]=vec[j];
                        count++; }
          loop_num[iloop]=count;

                                   } /* end reflection*/
                     }} /* end permutation */
/*
if(this_node == 0) printf("iloop,loop_num %d %d \n",iloop,loop_num[iloop]);
*/
      } /* end iloop */

/*
for(iloop=0;iloop<nloop;iloop++)
for(count=0;count<loop_num[iloop];count++){
if(this_node == 0) printf(" %d %d   %d  %d ",iloop,count,loop_char[count]
,loop_ch[iloop][count]);
for(j=0;j<length;j++)if(this_node == 0) printf(" %d",loop_table[iloop][count][j]);
if(this_node == 0) printf("\n");}
*/

 /* print out the loop coefficients */
if(this_node == 0) printf("loop coefficients: nloop rep loop_coeff   multiplicity\n");
for(i=0;i<nreps;i++)
for(j=0;j<nloop;j++) {
if(this_node == 0) printf("                    %d %d      %e     %d\n",
j,i,loop_coeff[j][i],loop_num[j]);}

make_loop_term();
}

void make_loop_term()
{
int dirs[10];
int dir,ln,iloop,l_num,length,k,irep;


l_num= 0;dir=0;

/*  loop over all the different contributions to plaquettes */
 
   for(iloop=0;iloop<nloop;iloop++){
/*
if(this_node == 0) printf("iloop=%d\n",iloop);
*/
      length=loop_length[iloop];
   /* loop over rotations and reflections */
      for(ln=0;ln<loop_num[iloop];ln++)if(loop_ch[iloop][ln]==1){
 /* set up dirs and sign  */
       for(k=0;k<length;k++){
                if(loop_table[iloop][ln][k] < 4 ){
                       dirs[k]=(dir+loop_table[iloop][ln][k] )% 4;
                }
                else {
                       dirs[k]=(7+dir-loop_table[iloop][ln][k] )% 4;
                }
}
        for(k=0;k<length;k++)if(dirs[k]==dir) {
		for(irep=0;irep<nreps;irep++) {
		loop_term[l_num][irep]=loop_coeff[iloop][irep];
/*
	  	if(this_node == 0) printf("%d %d %e\n",l_num,irep,loop_term[l_num][irep]);
*/
		}
 
        l_num++;
        } /* k (location in path) */
        }} /* ln and iloop */
 

}

void char_num(int dig[],int *chr,int *ch,int length)
{
int j;
int bdig[10],tenl,newv,old;

  *ch=0;
  tenl=1;
  for(j=0;j<length-1;j++) tenl=tenl*10;

  *chr=dig[length-1];
  for(j=length-2;j>=0;j--) *chr= *chr*10+dig[j];

/* forward*/
  old=*chr;
  for(j=length-1;j>=1;j--){
       newv=old-tenl*dig[j];
       newv=newv*10+dig[j];
       if(newv < *chr) *chr=newv;
       old=newv;           }

/* backward*/
   for(j=0;j<length;j++)bdig[j]=7-dig[length-j-1];
   old=bdig[length-1];
   for(j=length-2;j>=0;j--) old=old*10+bdig[j];
   if(old < *chr ) *chr=old;
   for(j=length-1;j>=1;j--){
       newv=old-tenl*bdig[j];
       newv=newv*10+bdig[j];
       if(newv < *chr) *chr=newv;
       old=newv;           }

  if(*chr < tenl) *ch=1;
}


