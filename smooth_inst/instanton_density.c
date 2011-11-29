/************************** instanton_density.c ****************************/
/* MIMD version 7 */
/* Formerly named set_of_paths.c */

#include "smooth_inst_includes.h"
#include <string.h>

#define NORM (1./4./32/PI/PI)

void instanton_density(void)
{
   void path(int *dir, int *sign, const int length);

   register int i,k;
   register site *s;
   complex trace;
   Real plaq1,plaq2,plaq3,plaq4,plaq5,tttr;
   Real sum;
   int dirs[max_inst_length],sign[max_inst_length],length,epsloop;

   Real loop_pwr,sump,summ,ch_d;

   /* these are for inst_table  */
   int ln,iloop;

   FORALLSITES(i,s)
   {
      s->ch_dens=0.0;
   }
   sum=0.0;

   for(iloop=0;iloop<nist;iloop++)
   {
      plaq1=0.0;
      plaq2=0.0;
      plaq3=0.0;
      plaq4=0.0;
      plaq5=0.0;

      length=inst_length[iloop];


      /* loop over rotations and reflections */
      for(ln=0;ln<inst_num[iloop];ln++)
      {
         /* set up dirs and sign  */
         for(k=0;k<length;k++)
         {
            if( inst_table[iloop][ln][k] < 4 )
            {
               sign[k]=1;
               dirs[k]=( inst_table[iloop][ln][k] )% 4;
            }
            else
            {
               sign[k]=-1;
               dirs[k]=( 7-inst_table[iloop][ln][k] )% 4;
            }
         }

         path(dirs, sign, length);

         epsloop = eps[iloop][ln];

         FORALLSITES(i,s)
         {
            trace=trace_su3( &(s->tempmat1) );
            tttr = trace.real - 3.0;

            loop_pwr=tttr;
            plaq1 += epsloop*loop_pwr;
            s->ch_dens += loop_coeff_inst[iloop][0]*loop_pwr*epsloop;
            loop_pwr *= tttr;
            plaq2 += epsloop*loop_pwr;
            s->ch_dens += loop_coeff_inst[iloop][1]*loop_pwr*epsloop;
            loop_pwr *= tttr;
            plaq3 += epsloop*loop_pwr; loop_pwr *= tttr;
            plaq4 += epsloop*loop_pwr; loop_pwr *= tttr;
            plaq5 += epsloop*loop_pwr;

         } /* for sites */
      } /* ln */

      plaq1 *= 1./4./32/PI/PI;
      plaq2 *= 1./4./32/PI/PI;
      plaq3 *= 1./4./32/PI/PI;
      plaq4 *= 1./4./32/PI/PI;
      plaq5 *= 1./4./32/PI/PI;

      if(this_node == 0) printf("loop %d completed \n",iloop);
      fflush(stdout);

   } /* iloop */

   FORALLSITES(i,s)
   {
      s->ch_dens *= NORM;
      sum += s->ch_dens;
   }

   g_floatsum( &sum );

   sump=0.0;
   summ=0.0;

   FORALLSITES(i,s)
   {
      ch_d=s->ch_dens;
      if (ch_d<0.0) summ += ch_d;
      if (ch_d>0.0) sump += ch_d;
   }

   g_floatsum( &summ );
   g_floatsum( &sump );

   if(this_node == 0) printf("Topological charge  %e \n",sum);
   if(this_node == 0) printf("Sum_plus  %e    Sum_minus  %e \n", sump,summ);

} /* instanton_density */
