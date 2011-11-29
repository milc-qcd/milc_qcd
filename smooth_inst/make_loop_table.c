/************************** make_loop_table.c *******************************/
/* MIMD version 7 */
/* WE NEED TO NORMALIZE THIS PROCEDURE WITH generic/gauge_stuff.o -CD! */

#include "smooth_inst_includes.h"

void make_loop_table2(void)
{
   void char_num(int dig[13], int chr[2], int *ch, const int length);

   /* Paths and coefficients after DeGrand, Hasenfratz, Kovacs,
      hep-lat/9705009 */

   static int inst_ind[nist][max_inst_length] =
   {
      {0,1,2,6,7,3,0,4,7,5,0,0,0,0,0,0},
      {0,1,2,7,3,5,0,4,7,6,0,0,0,0,0,0}
/*
      {0,1,2,3,7,4,0,6,7,5,0,0,0,0,0,0},
      {0,1,7,6,2,3,0,1,7,6,5,4,0,0,0,0}
*/
   };

   static int inst_length_in[nist] = {10,10};

   int perm[8],pp[8],ir[4];
   int length,iloop,i,j,k,chr[2];
   int vec[max_inst_length];
   int count,flag;
   int jsave = 0;
   int pnum,ppnum;

   for(j=0;j<nist;j++){inst_num[j]=0;inst_length[j]=inst_length_in[j];}

   for(i=0;i<inreps;i++)for(j=0;j<nist;j++)loop_coeff_inst[j][i]=0.0;
   loop_coeff_inst[0][0]=1.0;

   loop_coeff_inst[0][0]= 7.872507e-02;
   loop_coeff_inst[0][1]= 3.173630e-01;
   loop_coeff_inst[1][0]=-1.888383e-01;
   loop_coeff_inst[1][1]= 2.854577e-01;
/*
   loop_coeff_inst[0][0]= 9.029782e-02;
   loop_coeff_inst[0][1]= 4.884618e-01;
   loop_coeff_inst[1][0]=-1.786298e-01;
   loop_coeff_inst[1][1]= 4.126987e-01;
*/


   for(iloop=0;iloop<nist;iloop++)
   {
      length=inst_length[iloop];

      count=0;
      /* permutations */
      for(perm[0]=0;perm[0]<4;perm[0]++)
        for(perm[1]=0;perm[1]<4;perm[1]++)
          for(perm[2]=0;perm[2]<4;perm[2]++)
            for(perm[3]=0;perm[3]<4;perm[3]++)
            {
               if(   perm[0] != perm[1] && perm[0] != perm[2]
                  && perm[0] != perm[3]
                  && perm[1] != perm[2] && perm[1] != perm[3]
                  && perm[2] != perm[3] )
               {
                  /* permutation number */
                  pnum=0;
                  for(k=0;k<4;k++)
                  {
                     pp[k]=perm[k];
                  }

                  for(k=0;k<4;k++)
                  {
                     if ( pp[k] != k )
                     {
                        for(j=k+1;j<4;j++)
                        {
                           if( pp[j] == k )
                           {
                              pp[j]=pp[k];
                              pp[k]=k;
                              pnum +=1;
                           }
                        }
                     }
                  }


                  /* reflections*/
                  for(ir[0]=0;ir[0]<2;ir[0]++)
                    for(ir[1]=0;ir[1]<2;ir[1]++)
                      for(ir[2]=0;ir[2]<2;ir[2]++)
                        for(ir[3]=0;ir[3]<2;ir[3]++)
                        {
                           for(j=0;j<4;j++)
                           {
                              pp[j]=perm[j];
                              if(ir[j] == 1)
                              {
                                 pp[j]=7-pp[j];
                              }
                              pp[7-j]=7-pp[j];
                           }
                           /* create new vector*/
                           for(j=0;j<length;j++) vec[j]=pp[inst_ind[iloop][j]];

                           char_num(vec,chr,&ch,length);

                           flag=0;
                           /* check if it's a new set: */
                           for(j=0;j<count;j++)
                             if(chr[0] == loop_char[j][0] &&
                                chr[1] == loop_char[j][1])
                             {
                                flag=1; jsave=j; break;
                             }

                           if( flag == 0 )
                           {
                              loop_char[count][0]=chr[0];
                              loop_char[count][1]=chr[1];
                              for(j=0;j<length;j++)
                                inst_table[iloop][count][j]=vec[j];
                           }
                           /* permutation number */
                           ppnum=0;
                           for(j=0;j<4;j++)
                           {
                              ppnum +=ir[j];
                           }
                           ppnum= (ppnum+pnum) % 2;
                           if (flag == 0)
                           {
                              eps[iloop][count] = 1 - 2*ppnum;
                              count++;
                           }

                           if (flag == 1) eps[iloop][jsave] += 1 - 2*ppnum;

                           inst_num[iloop]=count;

                        } /* end reflection*/
               }
            } /* end permutation */
      if(this_node==0)
        printf("iloop,loop_num %d %d %d \n",iloop,inst_num[iloop],
               inst_length[iloop]);

/* for(j=0;j<count;j++)   printf("epsilon= %d \n", eps[iloop][j]); */

   } /* end iloop */

   if(this_node==0)
   {
      printf("Instanton coefficients \n");
      for(iloop=0;iloop<nist;iloop++)
      {
         printf(" inst %d ",iloop);
         for(i=0;i<inreps;i++)printf(" %e ",loop_coeff_inst[iloop][i]);
         printf("\n");
      }
   }

}
