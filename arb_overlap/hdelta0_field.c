/*************** hdelta0_field.c ****************************/
/* MIMD version 7 */

#include "arb_ov_includes.h"

#ifdef SSE
#define SSE_SUBS
#include "../sse/include/inline_sse.h"
#endif




void hdelta0_field(wilson_vector* src,wilson_vector* dest, int ikind)
{
        register int i;
        register site *s;
	Real r02;
	int chb,che,chbo,cheo,ii,jj;
	wilson_vector wtmp,wtmp2;
	wilson_vector * tmpvec;


	//	double source_norm;

	tmpvec=(wilson_vector*)malloc(sites_on_node*sizeof(wilson_vector));


        /* D^\dagger D */
	if(ikind == HZERO){
	  delta0_field(src,tmpvec,PLUS);
	  delta0_field(tmpvec,dest,MINUS);
	}

	else if(ikind == HTEST){
	  delta0_field(src,tmpvec,PLUS);
	  delta0_field(tmpvec,dest,MINUS);
	  
	  if(chirality_flag==1){chb=0;che=2;chbo=2;cheo=4;}
	  else{chb=2,che=4;chbo=0;cheo=2;}
	    FORALLSITES(i,s){
	      for(ii=chbo;ii<cheo;ii++)for(jj=0;jj<3;jj++){
		dest[i].d[ii].c[jj]=cmplx(0.0,0.0);
	      }
	    }
	  
	  /* check norm of dest 

	    source_norm=0.0;
            FORALLSITES(i,s){
              source_norm += (double)magsq_wvec(&(dest[i])  );
            }
            g_doublesum( &source_norm );
	    node0_printf("XXX %e\n",source_norm);
	  */
	}
	else if(ikind == HOVERLAP){
	  /* this does zero mass h_overlap^2  */
	  step_field(src,tmpvec);
	  /* for chiral src, dest = 2*R0*R0*(1+ chirality*eps)*src */
	  if(chirality_flag==1){chb=0;che=2;chbo=2;cheo=4;}
	  else{chb=2,che=4;chbo=0;cheo=2;}
	  r02=2.0*R0*R0;
	  if(chirality_flag== 1){
	    FORALLSITES(i,s){
	      copy_wvec(&src[i],&wtmp2);
	      for(ii=chb;ii<che;ii++)for(jj=0;jj<3;jj++){
		CADD(wtmp2.d[ii].c[jj],tmpvec[i].d[ii].c[jj],wtmp.d[ii].c[jj]);
		CMULREAL(wtmp.d[ii].c[jj],r02,wtmp.d[ii].c[jj]); 
	      }
	      for(ii=chbo;ii<cheo;ii++)for(jj=0;jj<3;jj++){
		wtmp.d[ii].c[jj]=cmplx(0.0,0.0);
	      }
	    copy_wvec(&wtmp,&(dest[i]));
	    }
	  }
	  else if(chirality_flag== -1){
	    FORALLSITES(i,s){
	      copy_wvec(&src[i],&wtmp2);
	      for(ii=chb;ii<che;ii++)for(jj=0;jj<3;jj++){
		CSUB(wtmp2.d[ii].c[jj],tmpvec[i].d[ii].c[jj],wtmp.d[ii].c[jj]);
		CMULREAL(wtmp.d[ii].c[jj],r02,wtmp.d[ii].c[jj]);
	    }
	    for(ii=chbo;ii<cheo;ii++)for(jj=0;jj<3;jj++){
	      wtmp.d[ii].c[jj]=cmplx(0.0,0.0);
	    }
	    copy_wvec(&wtmp,&(dest[i]));
	    }
	  }
	  else{node0_printf("chirality_flag undefined\n"); exit(1);}
	} /* kind_of_h0 flag */
	  else{
	    node0_printf("error in input: wrong kind_of_h0\n"); exit(0);
	  }

	  free(tmpvec);
}

#ifdef FIELDx
void hdelta0_field(wilson_vector *src,wilson_vector *dest,  int ikind)
{
        register int i;
        register site *s;
	Real r02;
	wilson_vector wtmp,wtmp2;
	complex ctmp,cn;
	Real chirality,cd;

	wilson_vector *vr; 

        /* D^\dagger D */

	if(ikind == HZERO){

	  vr=(wilson_vector*)malloc(sites_on_node*sizeof(wilson_vector));

	  r02=1.0/scalez/scalez;
	  delta0_field(src,vr,PLUS);
	  delta0_field(vr,dest,MINUS);

	  FORALLSITES(i,s)
                        scalar_mult_wvec(&(dest[i]),r02,&(dest[i]));

	  free(vr) ;

	}

	else{
        node0_printf("error in input: wrong kind_of_h0 in hdelta0_field\n"); exit(0);
      }
}

#endif
void g5delta0_field(wilson_vector *src,wilson_vector *dest,int sign)
{
    register int i;
    register site* s;
    
    delta0_field(src,dest,sign);
    FORALLSITES(i,s) mult_by_gamma(&(dest[i]),&(dest[i]),GAMMAFIVE);
    
}
