/*  File: inst_size.c
    Code to identify and measure instantons on a lattice
    Authors: A. Hasenfratz and T. Kovacs
    */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include "../../include/config.h"
#include "../../include/int32type.h"
#define MAXFILENAME 256   /* ASCII string length for all file names */
#define MAXINSTCOUNT 200

int NX;                /* lttice size */
int NY;
int NZ;
int NT;

/* hex 00010203 for Real 00010223 for double */
#define TOPO_VERSION_NUMBER (66051+8*(sizeof(Real)-4)  

void read_in(char *fname, int32type dim[4]);
void dump_lat(char *filename);
void clust(int x,int y,int z,int t,Real *radi,Real *dradi,int inst);
void fluct();
void aver_fluct();               /* averages the lattice */
void distance(int x,int y,int z,int t,Real rho,Real drho,
	      int *fl,Real *dist_like,Real *dist_opp);
void find_max_small(int *xmax,int *ymax,int *zmax,int *tmax,
		    int x,int y,int z,int t,int *inst);
void find_max(int *x,int *y,int *z,int *t,int *inst);
void subtr(int x,int y,int z,int t,Real r,int inst);
void root(Real q,Real *rh,int ir);
void byterevn(int32type w[], int n);

Real ****array4d_Real(int nx, int ny, int nz, int nt);
int ****array4d_int(int nx, int ny, int nz, int nt);

Real ****latt;
int ****label;
Real stop_min,stop_max;

int loc_i[MAXINSTCOUNT][4],loc_a[MAXINSTCOUNT][4];
Real rad_i[MAXINSTCOUNT],drad_i[MAXINSTCOUNT],
  height_i[MAXINSTCOUNT],charge_i[MAXINSTCOUNT];
Real rad_a[MAXINSTCOUNT],drad_a[MAXINSTCOUNT],
  height_a[MAXINSTCOUNT],charge_a[MAXINSTCOUNT];
int count_i=0,count_a=0;
int inst;
int step[3]={0,1,-1};

/* Instanton charge density profile: To be precise, the total
   instanton charge between spherical shells of radius r and r + dr.
   Here dr = 0.5.  The first shell is exceptional and covers 0 <= r < 1.
   The second shell covers 1 <= r < 1.5.  Successive shells start at the
   previous radius plus dr.

   First index counts spherical shells.
   Second index is the instanton size parameter

   */
/* charge[radius=1,1.5,...4.0...][rho=1,2,..8]  */
Real inst_charge[12][8] = {
  {0.054,  0.0183, 0.0052, 0.00196, 0.000877, 0.00045, 0.00025, 0.00014 },
  {0.55,   0.255 , 0.0993, 0.04451, 0.022171, 0.01206, 0.00699, 0.00410 },
  {0.20,   0.136 , 0.0685, 0.03504, 0.018795, 0.01069, 0.00637, 0.00382 },
  {0.14,   0.253 , 0.1943, 0.12156, 0.073544, 0.04527, 0.02822, 0.01760 },
  {0.0267, 0.063 , 0.0646, 0.04681, 0.030923, 0.02014, 0.01301, 0.00837 },
  {0.0090, 0.128 , 0.1716, 0.14876, 0.110894, 0.07868, 0.05297, 0.0357  },
  {0.0088, 0.055 , 0.0964, 0.10074, 0.085374, 0.06634, 0.04715, 0.0335  },
  {0.0045, 0.039 , 0.0812, 0.09864, 0.094020, 0.08014, 0.05829, 0.0436  },
  {0.0022, 0.021 , 0.0510, 0.06942, 0.072360, 0.06614, 0.04988, 0.0389  },
  {0.0022, 0.023 , 0.0626, 0.09678, 0.112827, 0.11356, 0.08534, 0.0707  },
  {0.0022, 0.009 , 0.0273, 0.04684, 0.059777, 0.06497, 0.04994, 0.0437  },
  {0.0022, 0.010 , 0.0328, 0.06129, 0.085341, 0.10082, 0.07609, 0.0673  }
};
/* Size corresponding to the second index in inst_charge */
Real inst_rho[8]={1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0};

int main(int argc, char *argv[])
{
   int ix,iy,iz,it;
   int32type dim[4];
   int xmax,ymax,zmax,tmax;     /* local max to start search for instanton */
   int x,y,z,t;
   int c_i,c_a;
   Real rad,drad,rad_test,charge,height;
   Real d_like,d_opp;
   int shift,i,j,k,flag=1;
   Real dist;
   /* of the instanton to be taken       */
   Real max_height[7]=
   {0.0182985, 0.0052485, 0.00195799, 0.000876759, 0.000446378, 0.00024912,
      1.412184e-04 };


   char fname[FILENAME_MAX];

   if(argc != 3)
   {
      printf("Usage: %s infile shift\n", argv[0]);
      exit(1);
   }
   else
   {
      strcpy(fname, argv[1]);
      shift=atoi(argv[2]);
   }


   /* read data into array */
   read_in(fname, dim);            

   /** dump_lat("foo"); **/

   label = array4d_int(NX, NY, NZ, NT);

   /*aver_fluct();*/

   for(it=0; it<NT; it++)   /* all pts are initialized to be */
                            /* unassociated with any instanton */
     for(iz=0; iz<NZ; iz++)
       for(iy=0; iy<NY; iy++)
         for(ix=0; ix<NX; ix++)
         {
	   label[ix][iy][iz][it] = -1; /* -1 means unassociated */
         }


   /* Mark low-fluctuating points for exclusion */
   fluct();

   /* find global max on lattice among unexcluded points */
   find_max(&x,&y,&z,&t,&inst);    
   /** printf("global_max %d %d %d %d\n",x,y,z,t); **/



   while(x != -1)
   {
     /* find the largest charge in the vicinity of x,y,z,t.  xmax = -1
        if none was found above the threshold */
      find_max_small(&xmax,&ymax,&zmax,&tmax,x,y,z,t,&inst);

      /** printf("find_max %d %d %d %d\n",xmax,ymax,zmax,tmax); **/

      while(xmax != -1)
      {
         label[xmax][ymax][zmax][tmax] =1;

	 /* Match the instanton charge density profile with the charge
	   distribution around xmax, ymax, zmax, tmax.  If OK, remove
	   the instanton from the distribution.  We we get its size
	   rad < NX/2 and add it to the list. If not OK, we go walk
	   around among the nearest neighbors, trying to match the
	   profile.  Eventually we either succeed or give up on the
	   site x,y,z,t. */
         if(inst == -1){
            height=latt[xmax][ymax][zmax][tmax];
            clust(xmax,ymax,zmax,tmax,&rad,&drad,inst);
            if(rad<NX/2){
               rad_test=rad;
               rad_a[count_a]=rad;
               drad_a[count_a]=drad;
               height_a[count_a]=height;
               loc_a[count_a][0]=xmax;
               loc_a[count_a][1]=ymax;
               loc_a[count_a][2]=zmax;
               loc_a[count_a][3]=tmax;
               count_a++;
            }

         }else{
            height=latt[xmax][ymax][zmax][tmax];
            clust(xmax,ymax,zmax,tmax,&rad,&drad,inst);
            if(rad<NX/2){
               rad_test=rad;
               rad_i[count_i]=rad;
               drad_i[count_i]=drad;
               height_i[count_i]=height;
               loc_i[count_i][0]=xmax;
               loc_i[count_i][1]=ymax;
               loc_i[count_i][2]=zmax;
               loc_i[count_i][3]=tmax;
               count_i++;
            } }

	 /* Find next maximum in vicinity of x,y,z,t */
         find_max_small(&xmax,&ymax,&zmax,&tmax,x,y,z,t,&inst);

	 /** printf("find_max %d %d %d %d\n",xmax,ymax,zmax,tmax); **/

      } /* end of while xmax */


      /*fluct();*/
      if(rad==NX/2+1){flag=0;x=-1;}else{flag=1;}

      /* Find the global max among all remaining sites not already checked
	 for an instanton profile */
      if(flag == 1 && x !=-1){
         find_max(&x,&y,&z,&t,&inst);
      }

   }/* end of while */
   /* printout */
   c_i=count_i;
   for(i=0;i<c_i;i++){
      count_i=i+1;
      inst=1;flag=0;
      for(k=0;k<i;k++){
         dist=0.;
         for(j=0;j<4;j++)dist+=(loc_i[i][j]-loc_i[k][j])*(loc_i[i][j]-loc_i[k][j]);
         if(dist==0)flag=1; }
      if(flag==0){
         distance(
                  loc_i[i][0],loc_i[i][1],loc_i[i][2],loc_i[i][3],
                  rad_i[i],drad_i[i],
                  &flag,&d_like,&d_opp);
         if(flag == 0)
           printf("I %d  %.2f  %d %d %d %d\n"
                  ,i+1, rad_i[i],
                  (loc_i[i][3]+shift+NT)%NT, (loc_i[i][2]+shift+NZ)%NZ,
                  (loc_i[i][1]+shift+NY)%NY, (loc_i[i][0]+shift+NX)%NX);
/*
           printf("I %d  %e  %e  %e  %e   %d %d %d %d     %d %d %d %d \n"
                  ,i+1, rad_i[i],drad_i[i], d_opp,d_like,
                  (loc_i[i][3]+shift+NT)%NT, (loc_i[i][2]+shift+NZ)%NZ,
                  (loc_i[i][1]+shift+NY)%NY, (loc_i[i][0]+shift+NX)%NX,
                  loc_i[i][3],loc_i[i][2],loc_i[i][1],loc_i[i][0]);
*/
      }}

   c_a=count_a;
   for(i=0;i<c_a;i++){
      count_a=i+1;
      inst=-1;flag=0;
      for(k=0;k<i;k++){
         dist=0.;
         for(j=0;j<4;j++)dist+=(loc_a[i][j]-loc_a[k][j])*(loc_a[i][j]-loc_a[k][j]);
        if(dist==0)flag=1;
      }
      if(flag==0){
         distance(
                  loc_a[i][0],loc_a[i][1],loc_a[i][2],loc_a[i][3],
                  rad_a[i],drad_a[i],
                  &flag,&d_like,&d_opp);
         if(flag == 0)
           printf("A %d  %.2f  %d %d %d %d\n"
                  ,i+1, rad_a[i],
                  (loc_a[i][3]+shift+NT)%NT, (loc_a[i][2]+shift+NZ)%NZ,
                  (loc_a[i][1]+shift+NY)%NY, (loc_a[i][0]+shift+NX)%NX);
/*
           printf("AI %d  %e  %e  %e %e   %d %d %d %d     %d %d %d %d \n"
                  ,-i-1, rad_a[i],drad_a[i], d_opp,d_like,
                  (loc_a[i][3]+shift+NT)%NT, (loc_a[i][2]+shift+NZ)%NZ,
                  (loc_a[i][1]+shift+NY)%NY, (loc_a[i][0]+shift+NX)%NX,
                  loc_a[i][3],loc_a[i][2],loc_a[i][1],loc_a[i][0]);
*/
      }}
}

/* Sets label = 1 for points that have densities fluctuating below 2
   sigma */

void fluct()
{               /* calculates the average flutuation */
   int all;
   int ix,iy,iz,it;
   Real val;
   Real sum2=0.;
   Real sum=0.;
   all=0;

   /* Calculate variance for sites with label < 0 */
   for(it=0; it<NT; it++)
     for(iz=0; iz<NZ; iz++)
       for(iy=0; iy<NY; iy++)
         for(ix=0; ix<NX; ix++)
           if(label[ix][iy][iz][it] <0 )
           {
              val=latt[ix][iy][iz][it];
              sum+=val;
              sum2+=val*val;
              all++;
           }
   if(all!=0)
   {
     /** printf(" fluct %d %e %e  \n",all,sum,sqrt(sum2/all)); **/
      stop_min = sqrt(sum2/all);
      if(stop_min<5.e-5)stop_min=5.0e-5;
      stop_max = stop_min*2.;
      /* looks for instantons that are at least 2 times above the fluctuations */
   }
   /* set the label to 1 for those that are below stop_max */
   for(it=0; it<NT; it++)
     for(iz=0; iz<NZ; iz++)
       for(iy=0; iy<NY; iy++)
         for(ix=0; ix<NX; ix++)
           if(sqrt(latt[ix][iy][iz][it]*latt[ix][iy][iz][it]) <stop_max )
/*
           if(latt[ix][iy][iz][it] >-stop_max )
*/
             label[ix][iy][iz][it]=1;
}

void aver_fluct()               /* averages the lattice */
{
   Real ****latt_av;
   int jx,jy,jz,jt;
   int kx,ky,kz,kt;
   int ix,iy,iz,it;

   latt_av = array4d_Real(NX, NY, NZ, NT);

   for(it=0; it<NT; it++)
     for(iz=0; iz<NZ; iz++)
       for(iy=0; iy<NY; iy++)
         for(ix=0; ix<NX; ix++)
         {
            latt_av[ix][iy][iz][it]=0;
            for(jt=-1; jt<=1; jt++)
              for(jz=-1; jz<=1; jz++)
                for(jy=-1; jy<=1; jy++)
                  for(jx=-1; jx<=1; jx++)
                  {
                     kx = (ix+jx+NX)%NX;
                     ky = (iy+jy+NY)%NY;
                     kz = (iz+jz+NZ)%NZ;
                     kt = (it+jt+NT)%NT;
                     latt_av[ix][iy][iz][it]+=latt[kx][ky][kz][kt];
                  }
         }

   for(it=0; it<NT; it++)
     for(iz=0; iz<NZ; iz++)
       for(iy=0; iy<NY; iy++)
         for(ix=0; ix<NX; ix++)
         {
            latt[ix][iy][iz][it]=latt_av[ix][iy][iz][it]/81;
         }

}


void read_in(char *fname, int32type dim[4])       /* reads data into array */
{
   FILE *fp;
   int32type int32typebuf[1];
   int32type topo_magic_number;
   char startfile[MAXFILENAME];
   int32type total_sweeps[1];
   Real ape_weight[1];
   int byterevflag;
   int ix,iy,iz,it;
   Real val[1];

   if ( (fp = fopen(fname, "rb")) == NULL )
   {
      fprintf(stderr, "Unable to open input file %s, ", fname);
      fprintf(stderr, "error %d.\n", errno);
      fprintf(stderr, "Exiting program.\n");
      exit(1);
   }
   fread(&int32typebuf[0], sizeof(int32type), 1, fp);
   if(int32typebuf[0] == TOPO_VERSION_NUMBER)
   {
      topo_magic_number = int32typebuf[0];
      byterevflag = 0;
   }
   else
   {
      byterevn(int32typebuf, 1);
      if(int32typebuf[0] == TOPO_VERSION_NUMBER)
      {
         topo_magic_number = int32typebuf[0];
         byterevflag = 1;
      }
      else
      {
         printf("Wrong version number. Terminating program.\n");
         exit(1);
      }
   }
   fread(dim, sizeof(int32type), 4, fp);
   if (byterevflag == 1)
   {
      byterevn(dim, 4);
   }
   NX = dim[0];
   NY = dim[1];
   NZ = dim[2];
   NT = dim[3];
   latt = array4d_Real(NX, NY, NZ, NT);
   fread(startfile, sizeof(char), 80, fp);
   fread(total_sweeps, sizeof(int32type), 1, fp);
   if (byterevflag == 1)
   {
      byterevn(total_sweeps, 1);
   }
   fread(ape_weight, sizeof(Real), 1, fp);
   if (byterevflag == 1)
   {
      byterevn((int32type *)ape_weight, 1);
   }

   for(it=0; it<NT; it++)
     for(iz=0; iz<NZ; iz++)
       for(iy=0; iy<NY; iy++)
         for(ix=0; ix<NX; ix++)
         {
            fread(&val[0], sizeof(Real), 1, fp);
            if (byterevflag == 1)
            {
               byterevn((int32type *)val, 1);
            }
            latt[ix][iy][iz][it] = val[0];
         }
   fclose(fp);
}

void dump_lat(char *filename)
{
   FILE *fp;
   int ix,iy,iz,it;

   if ( (fp = fopen(filename, "w")) == NULL )
   {
      fprintf(stderr, "Unable to open output file %s, ", filename);
      fprintf(stderr, "Exiting program.\n");
      exit(1);
   }

   for(it=0; it<NT; it++)
     for(iz=0; iz<NZ; iz++)
       for(iy=0; iy<NY; iy++)
         for(ix=0; ix<NX; ix++)
	   fprintf(fp,"%d %d %d %d %f\n",ix,iy,iz,it,latt[ix][iy][iz][it]);

   fclose(fp);
}

/* finds largest charge among the sites x,y,z,t and its nearest
   neighbors that has not yet been placed in an instanton
   Return its coordinates as xmax, ymax, zmax, tmax */
void find_max_small(int *xmax,int *ymax,int *zmax,int *tmax,
		    int x,int y,int z,int t,int *inst)
{
   int jx,jy,jz,jt;
   int kx,ky,kz,kt;
   int ix,iy,iz,it;
   int x0,y0,z0,t0;
   int flag = 0;
   Real curr_max;

   *xmax=*ymax=*zmax=*tmax=-1;

   curr_max = stop_max;

   /* Iterates over nearest neighbors */
   for(it=-1; it<=1; it++)
     for(iz=-1; iz<=1; iz++)
       for(iy=-1; iy<=1; iy++)
         for(ix=-1; ix<=1; ix++){
            x0 = (ix+x+NX)%NX;
            y0 = (iy+y+NY)%NY;
            z0 = (iz+z+NZ)%NZ;
            t0 = (it+t+NT)%NT;
            if( label[x0][y0][z0][t0] < 0 ){
               if(curr_max < latt[x0][y0][z0][t0] )
               {
                  *xmax=x0; *ymax=y0; *zmax=z0; *tmax=t0;
                  curr_max = latt[x0][y0][z0][t0];
                  *inst=1;
               }

               if(curr_max < -latt[x0][y0][z0][t0] )
               {
                  *xmax=x0; *ymax=y0; *zmax=z0; *tmax=t0;
                  curr_max = -latt[x0][y0][z0][t0];
                  *inst=-1;
               }
            }}

}


/* finds largest value on lattice that */
/* is not yet in an instanton and not a low-lying fluctuation */
void find_max(int *x,int *y,int *z,int *t,int *inst)
{
   int jx,jy,jz,jt;
   int kx,ky,kz,kt;
   int ix,iy,iz,it;
   int x0,y0,z0,t0;
   int flag = 0;
   Real curr_max;

   x0 = y0 = z0 = t0 = -1;

   curr_max = stop_max;

   for(it=0; it<NT; it++)
     for(iz=0; iz<NZ; iz++)
       for(iy=0; iy<NY; iy++)
         for(ix=0; ix<NX; ix++)
           if( label[ix][iy][iz][it] < 0 ){
              if(curr_max < latt[ix][iy][iz][it] )
              {
                 x0=ix; y0=iy; z0=iz; t0=it;
                 curr_max = latt[x0][y0][z0][t0];
                 *inst=1;
              }

              if(curr_max < -latt[ix][iy][iz][it] )
              {
                 x0=ix; y0=iy; z0=iz; t0=it;
                 curr_max = -latt[x0][y0][z0][t0];
                 *inst=-1;
              }
           }

   *x = x0; *y = y0; *z = z0; *t = t0;
}

/* Starting from x,y,z,t add up the charge inside a series of
   increasingly larger spherical shells.  Compare the charge against a
   standard profile to determine the instanton size parameter.  If the
   comparison succeeds, remove the instanton.  Return the size "radi"
   or NX/2 if the comparison fails.  dradi is the largest shell radius
   tried */
void clust(int x,int y,int z,int t,Real *radi,Real *dradi,int inst)
{
   int ix,iy,iz,it;
   int jx,jy,jz,jt;
   int kx,ky,kz,kt;
   int flag=0,NN;

   Real q_charge = 0.;
   Real q_charge_old=0;
   Real partial_sum,partial_sum_old;
   Real radius,dr,rad,RR[20];
   Real rms;
   int radnum;
   int ir,all,keep,count;
   int dist_flag;
   Real dist_like,dist_opp;
   Real rho,rho_s,rho_s2,rho_old;


   for(ix=0;ix<20;ix++)RR[ix]=0;
   q_charge=inst*latt[x][y][z][t];
   partial_sum=latt[x][y][z][t];


   radius=1.;
   radnum=0;
   dr=.5;
   flag=0;
   root(q_charge,&rho,radnum);
   RR[0]=rho;
   count=1;rho_s=rho;rho_s2=rho*rho;
   rho_old=rho;




   all=1;
/*
        printf("inst %e  %e  %e   %e %d \n",radius,rho,q_charge,partial_sum*all,all);
*/
   while(flag == 0 ){
      partial_sum_old=partial_sum;
      partial_sum=0.;
      q_charge_old=q_charge;
      rms=0.;
      all=0; keep=0;

      /* Compute the total charge "partial_sum" inside a spherical
	 shell with radius rad in the interval [radius,radius+dr) */
      NN=radius+dr+1;
      if(NN>NX/2)NN=NX/2;
      for(it=-NN; it<=NN; it++)
        for(iz=-NN; iz<=NN; iz++)
          for(iy=-NN; iy<=NN; iy++)
            for(ix=-NN; ix<=NN; ix++) {
               rad=sqrt(1.0*(it*it+iz*iz+iy*iy+ix*ix));
               if(rad>=radius && rad < radius+dr){
                  kx = (ix+x+NX)%NX;
                  ky = (iy+y+NY)%NY;
                  kz = (iz+z+NZ)%NZ;
                  kt = (it+t+NT)%NT;
                  all++;
		  /* Looks like "keep" and "all" are the same now */
                {
                   keep += 1;
                   partial_sum+=latt[kx][ky][kz][kt];}
               }} /* end of it,iy,iz,ix */

      radius +=dr;
      radnum++;
      if(radnum >=11)flag=1;

      /* Continue? */
      if(keep != 0 ){partial_sum=inst*partial_sum/keep;
                     q_charge+=partial_sum*(all);}
      if(q_charge_old >= q_charge || q_charge > 1. ) flag=1;
      if(radius >NX/2 || radius > rho*2)flag=1;

      root(partial_sum*all,&rho,radnum);
      RR[count]=rho;
      if(rho >NX/2)flag=1;
      if(count==1)rho_old=(rho_old+rho*32)/33.;
      /*if(count==2)rho_old=(rho_old+rho)/2.;*/
      if( rho>rho_old*1.05 || rho<rho_old/1.05 )flag=1;
      /*if(radius > rho)flag=1;*/



      if(flag==0  ){rho_s +=rho;rho_s2 +=rho*rho;count++;}
      else
      {
         if(count != 0 )
         {
            rho_s=rho_s/count;
            rho_s2=rho_s2/count-rho_s*rho_s;
            rho_s2=sqrt(rho_s2); 
            *radi=rho_s;
            *dradi=radius;
         }
         else
         {
            *radi=NX/2;
         }
      }

/*
        printf("inst %e  %e  %e   %e %d \n",radius,rho,q_charge,partial_sum*all,all);
*/
   } /* end of while */

   dist_flag=1;
   radius=radius-dr/2;
   if( rho_s <NX/2  && radius >= 1.5  &&
      (rho<2.0 || rho_s<radius*1.5 || radius>4.0) )
   {
      distance(x,y,z,t,rho_s,rho_s2,&dist_flag,&dist_like,&dist_opp);
   }

   /* Which one to keep? If rho>NX, clear choice;
      If dist_flag=1, one object is closer than its radius -
      keep the smaller radius,i.e. the one found before;
      If rho<4.5, a relatively small instanton - keep it if it's
      not killed by dist_flag;
      If rho>4.5, things can get murky. Keep only those that are fairly
      stable, i.e. they were explored to radius > rho/1.5
*/
   if( dist_flag==0 )
   {
      subtr(x,y,z,t,rho_s,inst);
      /** printf("dist %d %e %e \n",dist_flag,dist_like,dist_opp);
      printf(" %d %d %d %d %e %e\n   ",t,z,y,x,*radi,*dradi);
      for(ix=0;ix<count;ix++)printf("%e %e \n",dr*ix+1,RR[ix]); **/
   }
   else
   {
      *radi=NX/2;
   }
}

/* Given an instanton position x,y,z,t and charge "inst" and 
   size "rho", find the nearest like charge and unlike charge
   instanton and return the distances dist_like and dist_opp.
   If another instanton is closer than the size rho of either
   instanton, return fl = 1 */
void distance(int x,int y,int z,int t,Real rho,Real drho,
	      int *fl,Real *dist_like,Real *dist_opp)
{
   int i,j,k;
   int r[4],flag;
   int xx;
   Real dist;
   Real di=NT,da=NT;

   flag=0;
   r[0]=x;r[1]=y;r[2]=z;r[3]=t;
   for(i=0;i<count_i;i++){
      dist=0;
      xx = (loc_i[i][0]-r[0]+NX)%NX;
      if(xx>NX/2)xx=xx-NX;
      dist+=xx*xx*1.0;
      for(k=1;k<4;k++){
         xx = (loc_i[i][k]-r[k]+NY)%NY;
         if(xx>NY/2)xx=xx-NY;
         dist+=xx*xx*1.0;}
      dist=sqrt(dist);
      if(di>dist && dist >1.){di=dist;}
      if((dist<rho || dist <rad_i[i] || dist<(rad_i[i]+rho)*0.7)
         && inst==1 && dist != 0)
      {
         flag=1;
      }
   }
   for(i=0;i<count_a;i++){
      dist=0.;
      xx = (loc_a[i][0]-r[0]+NX)%NX;
      if(xx>NX/2)xx=xx-NX;
      dist+=xx*xx*1.0;
      for(k=1;k<4;k++){
         xx = (loc_a[i][k]-r[k]+NY)%NY;
         if(xx>NY/2)xx=xx-NY;
         dist+=xx*xx*1.0;}
      dist=sqrt(dist);
      if(da>dist && dist>1.){da=dist;}
      if((dist<rho || dist <rad_a[i] || dist<(rad_a[i]+rho)*0.7)
         &&  inst==- 1 && dist!=0){
         flag=1;
      }
                }
   *fl=flag;
   if(inst==1){*dist_like=di;*dist_opp=da;}
   if(inst==-1){*dist_like=da;*dist_opp=di;}
}


void subtr(int x,int y,int z,int t,Real r,int inst)
{
   int i,j,k,ir;
   int ix,iy,iz,it;
   int kx,ky,kz,kt;
   int all;
   Real radius,dr,rad;
   int radnum;
   Real q[12];

   i=r-1;
   if(i>=6)i=5;
   if(i<0)i=0;
   for(ir=0;ir<12;ir++)
     q[ir]=inst_charge[ir][i]+(inst_charge[ir][i]-inst_charge[ir][i+1])*(r-inst_rho[i])/(inst_rho[i]-inst_rho[i+1]);
   latt[x][y][z][t] -= q[0];
   label[x][y][z][t]++;

   radius=1.;
   radnum=1;
   dr=.5;
   while( radnum <= 11)
   {
      /* first count the number of points in the shell */
      all=0;
      for(it=-NT/2; it<NT/2; it++)
        for(iz=-NZ/2; iz<NZ/2; iz++)
          for(iy=-NY/2; iy<NY/2; iy++)
            for(ix=-NX/2; ix<NX/2; ix++) {
               rad=sqrt(1.0*(it*it+iz*iz+iy*iy+ix*ix));
               if(rad>=radius && rad < radius+dr)all++;}

      for(it=-NT/2; it<NT/2; it++)
        for(iz=-NZ/2; iz<NZ/2; iz++)
          for(iy=-NY/2; iy<NY/2; iy++)
            for(ix=-NX/2; ix<NX/2; ix++) {
               rad=sqrt(1.0*(it*it+iz*iz+iy*iy+ix*ix));
               if(rad>=radius && rad < radius+dr){
                  kx = (ix+x+NX)%NX;
                  ky = (iy+y+NY)%NY;
                  kz = (iz+z+NZ)%NZ;
                  kt = (it+t+NT)%NT;
                  latt[kx][ky][kz][kt] = latt[kx][ky][kz][kt]-inst*q[radnum]/all;
                  if((radius+dr)/r <1.0 || radius == 1.) label[kx][ky][kz][kt]++;
               }}
      radius+=dr;radnum++;
   }


}



/* Given the instanton charge q inside a shell and a distance index
   ir, find the instanton size parameter by interpolating in a table
   of known lattice instanton profiles for various size parameters.
   Return rh = the radius or return rh = NX/2 if not determined from
   the table. */
void root(Real q,Real *rh,int ir)
{
   int i,j;
   Real y1,y2,al;
   Real r1,r2,r;
   Real q1,q2;


   if(q<1)
   {
      r1=0;r2=0;
      q1=0;q2=0;
      for(i=1;i<8;i++)
        if((q>= inst_charge[ir][i] || i==7) &&
           inst_charge[ir][i] < inst_charge[ir][i-1] &&  r1 == 0)
        {
           r1=inst_rho[i-1];r2=inst_rho[i];
           q1=inst_charge[ir][i-1];q2=inst_charge[ir][i];
        }

      /* Interpolate to get instanton size */
      if(q2!= q1)
      {
         *rh=r2+(r2-r1)*(q-q2)/(q2-q1);
      }
      else
      {
         *rh=NX/2;
      }
   }
   else
   {
      *rh=NX/2;
   }

}

/*----------------------------------------------------------------------*/

/* For doing byte reversal on 32-bit words */
/* From io_lat4.c */

void byterevn(int32type w[], int n)
{
  register int32type old,new;
  int j;

  for(j=0; j<n; j++)
    {
      old = w[j];
      new = old >> 24 & 0x000000ff;
      new |= old >> 8 & 0x0000ff00;
      new |= old << 8 & 0x00ff0000;
      new |= old << 24 & 0xff000000;
      w[j] = new;
    }
} /* byterevn */

/*----------------------------------------------------------------------*/

/* allocates a 4-dimensional Real array */
Real ****array4d_Real(int nx, int ny, int nz, int nt)
{
   Real ****q;                  /* 4-d array */
   long x;
   long y;
   long z;
   long t;

   q = (Real ****)malloc((size_t)(nx * sizeof(Real ***)));
   if (q == NULL)
   {
      printf("Lattice nx too big.\n");
      exit(1);
   }
   for (x = 0; x < nx; x++)
   {
      q[x] = (Real ***)malloc((size_t)(ny * sizeof(Real **)));
      if (q[x] == NULL)
      {
         printf("Lattice ny too big.\n");
         exit(1);
      }
      for (y = 0; y < ny; y++)
      {
         q[x][y] = (Real **)malloc((size_t)(nz * sizeof(Real *)));
         if (q[x][y] == NULL)
         {
            printf("Lattice nz too big.\n");
            exit(1);
         }
         for (z = 0; z < nz; z++)
         {
            q[x][y][z] = (Real *)malloc((size_t)(nt * sizeof(Real)));
            if (q[x][y][z] == NULL)
            {
               printf("Lattice nt too big.\n");
               exit(1);
            }
         }
      }
   }
   return q;
} /* array4d_Real */

/*----------------------------------------------------------------------*/

/* allocates a 4-dimensional int array */
int ****array4d_int(int nx, int ny, int nz, int nt)
{
   int ****q;                  /* 4-d array */
   long x;
   long y;
   long z;
   long t;

   q = (int ****)malloc((size_t)(nx * sizeof(int ***)));
   if (q == NULL)
   {
      printf("Lattice nx too big.\n");
      exit(1);
   }
   for (x = 0; x < nx; x++)
   {
      q[x] = (int ***)malloc((size_t)(ny * sizeof(int **)));
      if (q[x] == NULL)
      {
         printf("Lattice ny too big.\n");
         exit(1);
      }
      for (y = 0; y < ny; y++)
      {
         q[x][y] = (int **)malloc((size_t)(nz * sizeof(int *)));
         if (q[x][y] == NULL)
         {
            printf("Lattice nz too big.\n");
            exit(1);
         }
         for (z = 0; z < nz; z++)
         {
            q[x][y][z] = (int *)malloc((size_t)(nt * sizeof(int)));
            if (q[x][y][z] == NULL)
            {
               printf("Lattice nt too big.\n");
               exit(1);
            }
         }
      }
   }
   return q;
} /* array4d_int */
