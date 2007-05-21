/************************** char_num.c *******************************/
/* MIMD version 7 */

#include "smooth_inst_includes.h"

void char_num(int dig[13], int chr[2], int *ch, const int length)
{
   int j;
   int tenl;
   int base;
   int save;

   base=8;
   *ch=0;
   tenl=1;
   for (j=0;j<length/2;j++) tenl=tenl*base;

   chr[0]=0;
   chr[1]=0;

   save=1;
   for (j=0;j<length/2;j++)
   {
      chr[0] += dig[j]*save;
      save *= base;
   }

   save=1;
   for (j=length/2;j<length;j++)
   {
      chr[1] += dig[j]*save;
      save *= base;
   }

   return;
}
