/*********************** save_topo.c *************************/
/* MIMD version 7 */

/* routine for FFdual output. */
/* This works for both Intel and Ncube, but other machines may need
   special treatment */

#include "smooth_inst_includes.h"

void save_topo(char *filenam)
{
   FILE *fp = NULL;
   int currentnode,newnode;
   int i,x,y,z,t;
   Real lbuf;
  /* Hack to distinguish single and double precision files */
   int32type topo_magic_number = TOPO_VERSION_NUMBER ;
   int32type tmp;

   /* node 0 does all the writing */
   if(this_node==0)
   {
      fp = fopen(filenam,"wb");
      if( fwrite(&topo_magic_number, sizeof(int32type), 1, fp) != 1 )
      {
         printf("Write error in save_topo\n"); terminate(1);
      }
      tmp = (int32type)nx;
      if( fwrite(&tmp, sizeof(int32type), 1, fp) != 1 )
      {
         printf("Write error in save_topo\n"); terminate(1);
      }
      tmp = (int32type)ny;
      if( fwrite(&tmp, sizeof(int32type), 1, fp) != 1 )
      {
         printf("Write error in save_topo\n"); terminate(1);
      }
      tmp = (int32type)nz;
      if( fwrite(&tmp, sizeof(int32type), 1, fp) != 1 )
      {
         printf("Write error in save_topo\n"); terminate(1);
      }
      tmp = (int32type)nt;
      if( fwrite(&tmp, sizeof(int32type), 1, fp) != 1 )
      {
         printf("Write error in save_topo\n"); terminate(1);
      }
      if( fwrite((char *)startfile, sizeof(char), 80, fp) != 80 )
      {
         printf("Write error in save_topo\n"); terminate(1);
      }
      tmp = (int32type)total_sweeps;
      if( fwrite(&tmp, sizeof(int32type), 1, fp) != 1 )
      {
         printf("Write error in save_topo\n"); terminate(1);
      }
      if( fwrite((Real *)&ape_weight, sizeof(Real), 1, fp) != 1 )
      {
         printf("Write error in save_topo\n"); terminate(1);
      }
   }
   g_sync();
   currentnode=0;

   for(t=0;t<nt;t++)for(z=0;z<nz;z++)for(y=0;y<ny;y++)for(x=0;x<nx;x++)
   {
      newnode=node_number(x,y,z,t);
      if(newnode != currentnode)
      { /* switch to another node */
         g_sync();
         currentnode=newnode;
      }

      if(this_node==0)
      {
         if(currentnode==0)
         {
            i=node_index(x,y,z,t);
            lbuf = (Real)lattice[i].ch_dens;
         }
         else
         {
            get_field((char *)&lbuf, sizeof(Real), currentnode);
         }
         if( fwrite(&lbuf, sizeof(Real), 1, fp) != 1 )
         {
            printf("Write error in save_topo\n"); terminate(1);
         }
      }
      else       /* for nodes other than 0 */
      {
         if(this_node==currentnode)
         {
            i=node_index(x,y,z,t);
            lbuf=lattice[i].ch_dens;
            send_field((char*)&lbuf, sizeof(Real), 0);
         }
      }
   }
   g_sync();
   if(this_node==0)
   {
      fflush(fp);
      printf("Saved FFD in file  %s\n",filenam);
      fclose(fp);
      fflush(stdout);
   }
}
