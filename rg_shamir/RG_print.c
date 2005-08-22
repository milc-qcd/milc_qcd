#include <stdio.h>
#include <qdp.h>
#include <math.h>
#include "RG_Shamir_includes.h"
#include "RG_include.h"
#define TOL 10e-6


void check_cv(QLA_ColorVector *df, int coords[])
{
  QLA_Complex z;
  int i,j;

  for(i=0; i<QLA_Nc; i++) 
   {
      z = QLA_elem_V(*df,i);
      if ( (fabs(QLA_real(z)) > TOL) || fabs((QLA_imag(z)) > TOL) )
       {
       fprintf(stderr,"node=%d \n",this_node);
       fprintf(stderr,"x=something wrong for color=%d \n",i);
       fprintf(stderr,"%g\t%g\n", QLA_real(z), QLA_imag(z));
       }
   }
     

return;
}

void print_cv(QLA_ColorVector *df, int coords[])
{
  QLA_Complex z;
  int i,j;

 // fprintf(stderr,"x=(%d %d %d %d)\n",coords[0],coords[1],coords[2],coords[3]);
// if((coords[0]==0)&&(coords[1]==0)&&(coords[2]==0)&&(coords[3]==0))
    for(j=0; j<QLA_Nc; j++) 
    {
      z = QLA_elem_V(*df,j);
     if ( (fabs(QLA_real(z)) > TOL) || fabs((QLA_imag(z)) > TOL) )
     {
      printf("%d %d %d %d %g\t%g\n",coords[0],coords[1],coords[2],coords[3],QLA_real(z), QLA_imag(z));fflush(stdout);
      //printf("x=(%d %d %d %d) this node %d\n",coords[0],coords[1],coords[2],coords[3],this_node);fflush(stdout);
      //printf("%g\t%g this node %d\n", QLA_real(z), QLA_imag(z),this_node);fflush(stdout);
      }
    }

return;
}


void print_cv_node(QLA_ColorVector *df, int coords[])
{
  QLA_Complex z;
  char file[10] = "RG_file";
  char filename[20];
  FILE *fw;
  int i,j;

  sprintf(file,"%s.%d",file,this_node);
  fw = fopen(file,"a");

  z = QLA_elem_V(*df,0);
//  if ( (fabs(QLA_real(z)) > TOL) || fabs((QLA_imag(z)) > TOL) )
    fprintf(fw,"%d %d %d %d %d %d %g\t%g\n",S,T,coords[0],coords[1],coords[2],coords[3],QLA_real(z), QLA_imag(z));fflush(stdout);
  

  fclose(fw);

return;
}

void print_cv_node_tot(QLA_ColorVector *df, int coords[])
{
  QLA_Complex z;
  char file[10] = "RG_file";
  char filename[20];
  FILE *fw;
  int i,j;

  sprintf(file,"%s.%d",file,this_node);
  fw = fopen(file,"a");

    for(j=0; j<QLA_Nc; j++) 
    {
     z = QLA_elem_V(*df,j);
    fprintf(fw,"%d %d %d %d %d %d %d %g\t%g\n",S,T,j,coords[0],coords[1],coords[2],coords[3],QLA_real(z), QLA_imag(z));fflush(stdout);
    }

  fclose(fw);

return;
}


void print_gl(QLA_ColorMatrix *gl, int coords[])
{
  QLA_Complex z;
  int i,j;

//  fprintf(stderr,"x=(%d %d %d %d)\n",coords[0],coords[1],coords[2],coords[3]);
//if((coords[0]==0)&&(coords[1]==0)&&(coords[2]==0)&&(coords[3]==0))
//{
  printf("x=(%d %d %d %d) this node %d\n",coords[0],coords[1],coords[2],coords[3],this_node);
  fflush(stdout);
  for(i=0; i<QLA_Nc; i++) 
  for(j=0; j<QLA_Nc; j++) 
   {
      z = QLA_elem_M(*gl,i,j);
      printf("%g\t%g this node %d\n ", QLA_real(z), QLA_imag(z),this_node);
      fflush(stdout);
   }  
//}
return;
}

void check_gl(QLA_ColorMatrix *gl, int coords[])
{
  QLA_Complex z;
  int i,j;

  for(i=0; i<QLA_Nc; i++) 
  for(j=0; j<QLA_Nc; j++) 
   {
      z = QLA_elem_M(*gl,i,j);
      if ( (fabs(QLA_real(z)) > TOL) || fabs((QLA_imag(z)) > TOL) )
      {
      printf("node=%d \n",this_node);
      printf("x=something wrong for color=%d,%d \n",i,j);
      printf("%g\t%g this node %d\n ", QLA_real(z), QLA_imag(z),this_node);
      fflush(stdout);
     }
   }  

return;
}


