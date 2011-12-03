/****** gauge_force_symzk1_qdp.c   ************/
/* MIMD version 7 */
/* Gauge action Symanzik improved 1x1 + 1x2 + 1x1x1 */
/* 9/2006 started by L.L. */
/* 10/06  CD standalone QDP */

#include "generic_includes.h"
#include "../include/generic_qdp.h"
#include <qdp.h>

static QDP_ColorMatrix *fblink[8];
extern Real beta;

void
imp_gauge_force_qdp(QDP_ColorMatrix *flinks[], QDP_ColorMatrix *force[], 
		    QLA_Real coeff[], Real eps){
  
  Real Plaq, Rect, Pgm ;
  QDP_ColorMatrix *tempmom_qdp[4];
  QDP_ColorMatrix *Amu[6]; // products of 2 links Unu(x)*Umu(x+nu)
  QDP_ColorMatrix *tmpmat;
  QDP_ColorMatrix *tmpmat1;
  QDP_ColorMatrix *tmpmat2;
  QDP_ColorMatrix *staples;
  QDP_ColorMatrix *tmpmat3;
  QDP_ColorMatrix *tmpmat4;

  int i, k;
  int mu, nu, sig;
  double dtime;
  Real eb3 = -eps*beta/3.0;
#ifdef GFTIME
  int nflop = 96720;
#endif

  int j[3][2] = {{1,2},
                 {0,2},
                 {0,1}};
  
  dtime = -dclock();

  for(mu=0; mu<4; mu++) {
    tempmom_qdp[mu] = QDP_create_M();
    QDP_M_eq_zero(tempmom_qdp[mu], QDP_all);
  }

  tmpmat = QDP_create_M();
  for(i=0; i< 4; i++) {
    fblink[i] = flinks[i];
    fblink[OPP_DIR(i)] = QDP_create_M();
    QDP_M_eq_sM(tmpmat, fblink[i], QDP_neighbor[i], QDP_backward, QDP_all);
    QDP_M_eq_Ma(fblink[OPP_DIR(i)], tmpmat, QDP_all);
  }
  

  for(i=0; i<6; i++) {
    Amu[i] = QDP_create_M();
    //Smu[i] = QDP_create_M();
  }

  staples = QDP_create_M();
  tmpmat1 = QDP_create_M();
  tmpmat2 = QDP_create_M();
  tmpmat3 = QDP_create_M();
  tmpmat4 = QDP_create_M();

  Plaq = coeff[0];
  Rect = coeff[1];
  Pgm  = coeff[2];

  for(mu=0; mu<4; mu++) {
    i=0;
    for(nu=0; nu<4; nu++) {
      if(nu!=mu){
	// tmpmat1 = Umu(x+nu)
	QDP_M_eq_sM(tmpmat1, fblink[mu], QDP_neighbor[nu], QDP_forward, QDP_all); 
        QDP_M_eq_M_times_M(Amu[i], fblink[nu], tmpmat1, QDP_all);

        //tmpmat2 = Umu(x-nu)
	QDP_M_eq_sM(tmpmat2, fblink[mu], QDP_neighbor[nu], QDP_backward, QDP_all);
        QDP_M_eq_M_times_M(Amu[i+3], fblink[OPP_DIR(nu)], tmpmat2, QDP_all);
       

 
	//tmpmat = U_{nu}(x+mu)
        QDP_M_eq_sM(tmpmat, fblink[nu], QDP_neighbor[mu], QDP_forward, QDP_all);
        QDP_M_eq_M_times_Ma(staples, Amu[i], tmpmat, QDP_all);        
        QDP_M_peq_r_times_M(tempmom_qdp[mu], &Plaq, staples, QDP_all);
 
        //tmpmat = U_{-nu}(x+mu)
        QDP_M_eq_sM(tmpmat, fblink[OPP_DIR(nu)], QDP_neighbor[mu], QDP_forward, QDP_all);
        QDP_M_eq_Ma_times_M(tmpmat3, fblink[OPP_DIR(nu)], staples, QDP_all);
        QDP_M_eq_M_times_M(tmpmat4, tmpmat3, tmpmat, QDP_all);
        QDP_M_eq_sM(tmpmat, tmpmat4, QDP_neighbor[nu], QDP_forward, QDP_all);
        QDP_M_peq_r_times_M(tempmom_qdp[mu], &Rect, tmpmat, QDP_all);

        QDP_M_eq_Ma_times_M(tmpmat4, tmpmat2, tmpmat3, QDP_all);
        QDP_M_eq_sM(tmpmat, tmpmat4, QDP_neighbor[nu], QDP_forward, QDP_all);
        QDP_M_eq_sM(tmpmat3, tmpmat, QDP_neighbor[mu], QDP_backward, QDP_all);
        QDP_M_peq_r_times_M(tempmom_qdp[nu], &Rect, tmpmat3, QDP_all);

        //tmpmat = U_{-nu}(x+mu)
        QDP_M_eq_sM(tmpmat, fblink[OPP_DIR(nu)], QDP_neighbor[mu], QDP_forward, QDP_all);
        QDP_M_eq_M_times_Ma(tmpmat3, tmpmat2, tmpmat, QDP_all);
        QDP_M_eq_M_times_Ma(tmpmat, tmpmat3, staples, QDP_all);        
        QDP_M_eq_sM(tmpmat3, tmpmat, QDP_neighbor[nu], QDP_forward, QDP_all);
        QDP_M_peq_r_times_M(tempmom_qdp[nu], &Rect, tmpmat3, QDP_all);




        //tmpmat = U_{-nu}(x+mu) 
        QDP_M_eq_sM(tmpmat, fblink[OPP_DIR(nu)], QDP_neighbor[mu], QDP_forward, QDP_all);
        QDP_M_eq_M_times_Ma(staples, Amu[i+3], tmpmat, QDP_all);        
        QDP_M_peq_r_times_M(tempmom_qdp[mu], &Plaq, staples, QDP_all);

        QDP_M_eq_Ma_times_M(tmpmat3, fblink[nu], staples, QDP_all);
        QDP_M_eq_sM(tmpmat, fblink[nu], QDP_neighbor[mu], QDP_forward, QDP_all);
        QDP_M_eq_M_times_M(tmpmat4, tmpmat3, tmpmat, QDP_all);
        QDP_M_eq_sM(tmpmat, tmpmat4, QDP_neighbor[nu], QDP_backward, QDP_all);
        QDP_M_peq_r_times_M(tempmom_qdp[mu], &Rect, tmpmat, QDP_all);

        QDP_M_eq_Ma_times_M(tmpmat, tmpmat3, tmpmat1, QDP_all);
        QDP_M_eq_sM(tmpmat4, tmpmat, QDP_neighbor[mu], QDP_backward, QDP_all);
        QDP_M_peq_r_times_M(tempmom_qdp[nu], &Rect, tmpmat4, QDP_all);

        QDP_M_eq_sM(tmpmat, fblink[nu], QDP_neighbor[mu], QDP_forward, QDP_all);
        QDP_M_eq_M_times_M(tmpmat3, staples, tmpmat, QDP_all);
        QDP_M_eq_M_times_Ma(tmpmat4, tmpmat3, tmpmat1, QDP_all);
        QDP_M_peq_r_times_M(tempmom_qdp[nu], &Rect, tmpmat4, QDP_all);
        i++;
      }
      
    }

    // Construct the 5 pgm staples and add them to force
    QDP_M_eq_zero(staples, QDP_all);
    i=0;
    for(nu=0; nu<4; nu++){
      if(nu!=mu){
        k=0;
	for(sig=0; sig<4;sig ++){
	  if(sig!=mu && nu!=sig){
	    
	    // the nu_sig_mu ... staple and 3 reflections
            //tmpmat = Amu["sig"](x+nu)
	    QDP_M_eq_sM(tmpmat, Amu[j[i][k]], QDP_neighbor[nu], QDP_forward, QDP_all);
            //tmpmat1 = Unu(x)*Amu["sig"](x+nu)
            QDP_M_eq_M_times_M(tmpmat1, fblink[nu], tmpmat, QDP_all);   
            //tmpmat3 = Unu(x+mu+sig)
            QDP_M_eq_sM(tmpmat, fblink[nu], QDP_neighbor[mu], QDP_forward, QDP_all);
	    QDP_M_eq_sM(tmpmat3, tmpmat, QDP_neighbor[sig], QDP_forward, QDP_all); // HERE?
            //tmpmat2 = Unu(x)*Amu["sig"](x+nu)*adj(Unu(x+mu+sig))
	    QDP_M_eq_M_times_Ma(tmpmat2, tmpmat1, tmpmat3, QDP_all);
            //tmpmat = Usig(x+mu)
	    QDP_M_eq_sM(tmpmat, fblink[sig], QDP_neighbor[mu], QDP_forward, QDP_all);
            //tmpmat1 = Unu(x)*Amu["sig"](x+nu)*adj(Unu(x+mu+sig))*adj(Usig(x+mu))
	    QDP_M_eq_M_times_Ma(tmpmat1, tmpmat2, tmpmat, QDP_all);

	    QDP_M_peq_M(staples, tmpmat1, QDP_all);


            //tmpmat = Amu["sig"](x-nu)
	    QDP_M_eq_sM(tmpmat, Amu[j[i][k]], QDP_neighbor[nu], QDP_backward, QDP_all);
            //tmpmat1 = U_{-nu}(x)*Amu["sig"](x-nu)
            QDP_M_eq_M_times_M(tmpmat1, fblink[OPP_DIR(nu)], tmpmat, QDP_all);   
            //tmpmat3 = U_{-nu}(x+mu+sig)
            QDP_M_eq_sM(tmpmat, fblink[OPP_DIR(nu)], QDP_neighbor[mu], QDP_forward, QDP_all);
	    QDP_M_eq_sM(tmpmat3, tmpmat, QDP_neighbor[sig], QDP_forward, QDP_all); // HERE?
            //tmpmat2 = U_{-nu}nu(x)*Amu["sig"](x-nu)*adj(Unu(x+mu+sig))
	    QDP_M_eq_M_times_Ma(tmpmat2, tmpmat1, tmpmat3, QDP_all);
            //tmpmat = Usig(x+mu)
	    QDP_M_eq_sM(tmpmat, fblink[sig], QDP_neighbor[mu], QDP_forward, QDP_all);
            //tmpmat1 = U_{-nu}(x)*Amu["sig"](x-nu)*adj(Unu(x+mu+sig))*adj(Usig(x+mu))
	    QDP_M_eq_M_times_Ma(tmpmat1, tmpmat2, tmpmat, QDP_all);

	    QDP_M_peq_M(staples, tmpmat1, QDP_all);


            //tmpmat = Amu["-sig"](x-nu)
	    QDP_M_eq_sM(tmpmat, Amu[j[i][k]+3], QDP_neighbor[nu], QDP_backward, QDP_all);
            //tmpmat1 = U_{-nu}(x)*Amu["-sig"](x-nu)
            QDP_M_eq_M_times_M(tmpmat1, fblink[OPP_DIR(nu)], tmpmat, QDP_all);   
            //tmpmat = U_{-nu}(x+mu-sig)
            QDP_M_eq_sM(tmpmat, fblink[OPP_DIR(nu)], QDP_neighbor[mu], QDP_forward, QDP_all);
	    QDP_M_eq_sM(tmpmat3, tmpmat, QDP_neighbor[sig], QDP_backward, QDP_all); // HERE?
            //tmpmat2 = U_{-nu}nu(x)*Amu["-sig"](x-nu)*adj(Unu(x+mu-sig))
	    QDP_M_eq_M_times_Ma(tmpmat2, tmpmat1, tmpmat3, QDP_all);
            //tmpmat = U_{-sig}(x+mu)
	    QDP_M_eq_sM(tmpmat, fblink[OPP_DIR(sig)], QDP_neighbor[mu], QDP_forward, QDP_all);
            //tmpmat1 = U_{-nu}(x)*Amu["-sig"](x-nu)*adj(Unu(x+mu-sig))*adj(U_{-sig}(x+mu))
	    QDP_M_eq_M_times_Ma(tmpmat1, tmpmat2, tmpmat, QDP_all);

	    QDP_M_peq_M(staples, tmpmat1, QDP_all);

            


            //tmpmat = Amu["-sig"](x+nu)
	    QDP_M_eq_sM(tmpmat, Amu[j[i][k]+3], QDP_neighbor[nu], QDP_forward, QDP_all);
            //tmpmat1 = Unu(x)*Amu["-sig"](x+nu)
            QDP_M_eq_M_times_M(tmpmat1, fblink[nu], tmpmat, QDP_all);   
            //tmpmat3 = Unu(x+mu-sig)
            QDP_M_eq_sM(tmpmat, fblink[nu], QDP_neighbor[mu], QDP_forward, QDP_all);
	    QDP_M_eq_sM(tmpmat3, tmpmat, QDP_neighbor[sig], QDP_backward, QDP_all); // HERE?
            //tmpmat2 = Unu(x)*Amu["-sig"](x+nu)*adj(Unu(x+mu-sig))
	    QDP_M_eq_M_times_Ma(tmpmat2, tmpmat1, tmpmat3, QDP_all);
            //tmpmat = U_{-sig}(x+mu)
	    QDP_M_eq_sM(tmpmat, fblink[OPP_DIR(sig)], QDP_neighbor[mu], QDP_forward, QDP_all);
            //tmpmat1 = Unu(x)*Amu["sig"](x+nu)*adj(Unu(x+mu+sig))*adj(Usig(x+mu))
	    QDP_M_eq_M_times_Ma(tmpmat1, tmpmat2, tmpmat, QDP_all);

	    QDP_M_peq_M(staples, tmpmat1, QDP_all);

	    k++;
	  }//close if sig!=nu ...
	}//close sig loop
	i++;
      }// close if nu!=mu
    }//close the pgm nu loop

    QDP_M_peq_r_times_M(tempmom_qdp[mu], &Pgm, staples, QDP_all);
   

    
  }// closes the mu loop

  for(mu=0; mu<4; mu++){
    QDP_M_eq_M_times_Ma(tmpmat, fblink[mu], tempmom_qdp[mu], QDP_all);
    QDP_M_eq_r_times_M_plus_M( tempmom_qdp[mu], &eb3, tmpmat, force[mu], QDP_all);
    QDP_M_eq_antiherm_M(force[mu], tempmom_qdp[mu], QDP_all);
   }


  //DESTROY various fields

  QDP_destroy_M(tmpmat);
  QDP_destroy_M(tmpmat1);
  QDP_destroy_M(tmpmat2);
  QDP_destroy_M(tmpmat3);
  QDP_destroy_M(staples);
  QDP_destroy_M(tmpmat4);

  for(mu=0; mu<4; mu++){
    QDP_destroy_M(tempmom_qdp[mu]);
  }
  for(i=0; i<6; i++) {
    QDP_destroy_M(Amu[i]);
  }

  for(i=4; i<8; i++) {
    QDP_destroy_M(fblink[i]);
  }

  dtime += dclock();

#ifdef GFTIME
  node0_printf("GFTIME:   time = %e (Symanzik1) mflops = %e\n",dtime,
	       nflop*(double)volume/(1e6*dtime*numnodes()) );
#endif
} 


void 
imp_gauge_force( Real eps, field_offset mom_off ){

  int i,dir;
  site *s;
  Real **loop_coeff = get_loop_coeff();
  su3_matrix tmat;
  anti_hermitmat* momentum;
  QLA_Real coeff[3];
  QLA_ColorMatrix *temp;
  QDP_ColorMatrix *gf[4];
  QDP_ColorMatrix *force[4];
  double remaptime = -dclock();

  /* Allocate space for gauge field and force */
  FORALLUPDIR(dir){
    gf[dir] = QDP_create_M();
    force[dir]  = QDP_create_M();
  }

  /* Map gauge links to QDP */
  set4_M_from_site(gf, F_OFFSET(link), EVENANDODD);

  /* The force requires a special conversion from the antihermit type */
  FORALLUPDIR(dir){
    temp = QDP_expose_M (force[dir]);
    FORALLSITES(i,s) {
      momentum = (anti_hermitmat *)F_PT(s,mom_off);
      uncompress_anti_hermitian( momentum+dir, &tmat );
      memcpy((void *)&temp[i], (void *)&tmat, sizeof(QLA_ColorMatrix));
    }
    QDP_reset_M (force[dir]);
  }

  /* Map gauge coefficients to QLA */
  coeff[0] = loop_coeff[0][0]; /* plaq */
  coeff[1] = loop_coeff[1][0]; /* rect */
  coeff[2] = loop_coeff[2][0]; /* pgm */

  remaptime += dclock();

  /* Evaluate the gauge force */
  imp_gauge_force_qdp(gf, force, coeff, eps);

  remaptime -= dclock();

  /* Map the force back to MILC */
  /* It requires a special conversion to the antihermitian type */
  FORALLUPDIR(dir){
    temp = QDP_expose_M (force[dir]);
    FORALLSITES(i,s) {
      memcpy((void *)&tmat, (void *)&temp[i], sizeof(su3_matrix));
      momentum = (anti_hermitmat *)F_PT(s,mom_off);
      make_anti_hermitian( &tmat, momentum + dir); 
    }
    QDP_reset_M (force[dir]);
  }
  
  /* Clean up */
  FORALLUPDIR(dir){
    QDP_destroy_M(gf[dir]);
    QDP_destroy_M(force[dir]);
  }

  remaptime += dclock();
#ifdef GFTIME
  node0_printf("GFREMAP:  time = %e\n",remaptime);
#endif
}


