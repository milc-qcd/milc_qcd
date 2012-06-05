/* ************************************************************ */
/*								*/
/*			    U1POT.C				*/
/*								*/
/* 1. u1ploop(): Polyakov loop with U(1) links			*/
/* 3. hvy_u1pot(): Static q-q_bar potential in U(1) field	*/
/*								*/
/* Last updated on 06.21.07					*/
/*								*/
/* ************************************************************ */

#include "generic_u1_includes.h"

#define MAX_R nx-2
#define MAX_T (nt/4)
#define MAX_X (nx/2-2)

void shift_cfield(field_offset src,field_offset dest,int dir);

complex u1ploop(void)
{

  int d[4];
  Real sum;
  complex plp;
  Real *tempmat1;
  int i,t;
  site *s;
  msg_tag *mtag;

  tempmat1 = create_r_field();

  sum = 0.0;
  d[XUP]=d[YUP]=d[ZUP]=0;

  mtag = declare_strided_gather(u1_A+TUP, 4*sizeof(Real), sizeof(Real), 
				TUP, EVEN, gen_pt[0]);
  prepare_gather(mtag);
  do_gather(mtag);

  //  mtag=start_gather_site(F_OFFSET(u1link[TUP]),sizeof(complex),
  //		 TUP,EVEN,gen_pt[0]);
  wait_gather(mtag);
  FOREVENSITES(i,s){
    //      l1=&(s->u1link[TUP]);
    //    l2 = (complex *)gen_pt[0][i];
    tempmat1[i] = u1_A[4*i+TUP] + *((Real *)gen_pt[0][i]);
  }
  cleanup_gather(mtag);

  for(t=2;t<nt;t+=2){
      d[TUP]=t;
      mtag=start_general_gather_field(tempmat1,sizeof(Real),
				      d,EVEN,gen_pt[0]);
      wait_general_gather(mtag);
      FOREVENSITES(i,s){
	  if(s->t>1) continue;	/* only compute on first two slices */
	  //l1=(complex *)gen_pt[0][i];
	  //CMUL(tempmat1[i],*l1,tmat);
	  //          tempmat1[i].real=tmat.real;
	  //          tempmat1[i].imag=tmat.imag;
	  tempmat1[i] += *((Real *)gen_pt[0][i]);
         }
      cleanup_general_gather(mtag);
     }
  FOREVENSITES(i,s){
      if(s->t>1) continue;
      //      plp.real=tempmat1[i].real;
      //      plp.imag=tempmat1[i].imag;
      plp.real = cos(tempmat1[i]);
      sum += plp.real;
      s->loop.real=plp.real;
      s->loop.imag=0.0;
     }

  g_floatsum(&sum);
  plp.real=sum/((Real)(nx*ny*nz));
  plp.imag=0.0;
  
  destroy_r_field(tempmat1);
  return(plp);

} /* end of u1ploop() */

void hvy_u1plpot(void)
{

  complex plp;
  Real norml,plpcorr[MAX_R];
  register int i,r;
  register site *s;

  for(r=0;r<MAX_R;r++) plpcorr[r]=0.0;
  norml=(Real)(nx*ny*nz/2);

  /* polyakov loop correlation function */
  for(r=1;r<=MAX_R;r++){
      FOREVENSITES(i,s){
	  if(s->t>0) continue;
	  CMUL_J(lattice[i].loop,lattice[i+r].loop,plp);
	  plpcorr[r-1]+=plp.real;
	 }
      plpcorr[r-1]*=norml;
     }

  node0_printf("Polyakov loop correlation:\n");
  for(r=0;r<MAX_R;r++)
      node0_printf("%2d \t %e\n",r+1,plpcorr[r]);

} /* end of hvy_u1pot() */

void hvy_u1wlpot(field_offset links)
{

  int t_dist,x_dist,y_dist,z_dist;
  complex *l1,*l2,*l3;
  double wloop;
  msg_tag *mtag0;
  register int i;
  register site *s;
  field_offset oldmat,newmat,tt;

  node0_printf("hvy_u1wlpot(): MAX_T = %d, MAX_X = %d\n",MAX_T,MAX_X);

  for(t_dist=1;t_dist<=MAX_T;t_dist ++){
      if(t_dist==1)
	{
        FORALLSITES(i,s){
	    s->u1tmp.real=((complex *)F_PT(s,links))[TUP].real;
	    s->u1tmp.imag=((complex *)F_PT(s,links))[TUP].imag;
	   }
	}
      else
	{
	mtag0=start_gather_site(F_OFFSET(u1tmp),sizeof(complex),
				TUP,EVENANDODD,gen_pt[0]);
	wait_gather(mtag0);
	FORALLSITES(i,s){
	    l1=&(((complex *)F_PT(s,links))[TUP]);
	    l2=(complex *)gen_pt[0][i];
	    l3=&(s->loop);
	    CMUL(*l1,*l2,*l3);
	   }
	cleanup_gather(mtag0);
	FORALLSITES(i,s){
	    s->u1tmp.real=s->loop.real;
	    s->u1tmp.imag=s->loop.imag;
	   }
	}
      oldmat=F_OFFSET(gftmp2);
      newmat=F_OFFSET(loop);

      for(x_dist=0;x_dist<=MAX_X;x_dist++){
	  for(y_dist=0;y_dist<=MAX_X;y_dist++){
	      FORALLSITES(i,s){
		  ((complex *)F_PT(s,oldmat))->real=s->u1tmp.real;
		  ((complex *)F_PT(s,oldmat))->imag=s->u1tmp.imag;
		 }
	      for(i=0;i<x_dist;i++){
		  shift_cfield(oldmat,newmat,XUP);
		  tt=oldmat;oldmat=newmat;newmat=tt;
		 }
	      for(i=0;i<y_dist;i++){
		  shift_cfield(oldmat,newmat,YUP);
		  tt=oldmat;oldmat=newmat;newmat=tt;
		 }
	      for(z_dist=0;z_dist<=MAX_X;z_dist++){
		  wloop=0.0;
		  FORALLSITES(i,s){
		      wloop+=(double)(s->u1tmp.real*
				((complex *)F_PT(s,oldmat))->real
				     +s->u1tmp.imag*
				((complex *)F_PT(s,oldmat))->imag);
		     }
		  g_doublesum(&wloop);
		  node0_printf("POT_LOOP: %d %d %d %d \t %e\n",
			x_dist,y_dist,z_dist,t_dist,wloop/volume);
		  shift_cfield(oldmat,newmat,ZUP);
		  tt=oldmat;oldmat=newmat;newmat=tt;
		 } /* z_dist-loop ends */
	     } /* y_dist-loop ends */
	 } /* x_dist-loop ends */
     } /* t_dist-loop ends */

} /* end of hvy_u1plpot() */

void shift_cfield(field_offset src,field_offset dest,int dir)
{

  register int i;
  register site *s;
  msg_tag *mtag;

  mtag=start_gather_site(src,sizeof(complex),dir,EVENANDODD,gen_pt[0]);
  wait_gather(mtag);
  FORALLSITES(i,s){
      ((complex *)F_PT(s,dest))->real=((complex *)gen_pt[0][i])->real;
      ((complex *)F_PT(s,dest))->imag=((complex *)gen_pt[0][i])->imag;
     }
  cleanup_gather(mtag);

} /* end of shift_field() */

/* ************************************************************ */

