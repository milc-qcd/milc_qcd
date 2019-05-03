/***************** symmetrize.c *****************************************/

/* Symmetrize the current-current correlator over hypercubic group transformations */

/* MIMD version 7 */

/* 04/05/15 C. DeTar */

#include "rcorr_includes.h"

/*---------------------------------------------------------------*/
static void 
map_perm_xy(int x, int y, int z, int t, int *args, int fb,
	    int *xp, int *yp, int *zp, int *tp){
  *xp = y;  *yp = x;  *zp = z;  *tp = t; 
}

static int 
make_perm_xy(void){
  int dir;
  dir = make_gather(map_perm_xy, NULL, OWN_INVERSE, ALLOW_EVEN_ODD, SAME_PARITY);
  return dir;
}

/*---------------------------------------------------------------*/
static void 
map_perm_xz(int x, int y, int z, int t, int *args, int fb,
	    int *xp, int *yp, int *zp, int *tp){
  *xp = z;  *yp = y;  *zp = x;  *tp = t; 
}

static int 
make_perm_xz(void){
  int dir;
  dir = make_gather(map_perm_xz, NULL, OWN_INVERSE, ALLOW_EVEN_ODD, SAME_PARITY);
  return dir;
}

/*---------------------------------------------------------------*/
static void 
map_perm_yz(int x, int y, int z, int t, int *args, int fb,
	    int *xp, int *yp, int *zp, int *tp){
  *xp = x;  *yp = z;  *zp = y;  *tp = t; 
}

static int 
make_perm_yz(void){
  int dir;
  dir = make_gather(map_perm_yz, NULL, OWN_INVERSE, ALLOW_EVEN_ODD, SAME_PARITY);
  return dir;
}
/*---------------------------------------------------------------*/
void map_refl_x(int x, int y, int z, int t, int *args, int fb,
		int *xp, int *yp, int *zp, int *tp){
  *xp = (nx-x)%nx;  *yp = y;  *zp = z;  *tp = t; 
}

static int 
make_refl_x(void){
  int dir;
  dir = make_gather(map_refl_x, NULL, OWN_INVERSE, ALLOW_EVEN_ODD, SAME_PARITY);
  return dir;
}
/*---------------------------------------------------------------*/
static void 
map_refl_y(int x, int y, int z, int t, int *args, int fb,
	   int *xp, int *yp, int *zp, int *tp){
  *xp = x;  *yp = (ny-y)%ny;  *zp = z;  *tp = t; 
}

static int 
make_refl_y(void){
  int dir;
  dir = make_gather(map_refl_y, NULL, OWN_INVERSE, ALLOW_EVEN_ODD, SAME_PARITY);
  return dir;
}
/*---------------------------------------------------------------*/
static void 
map_refl_z(int x, int y, int z, int t, int *args, int fb,
	   int *xp, int *yp, int *zp, int *tp){
  *xp = x;  *yp = y;  *zp = (nz-z)%nz;  *tp = t; 
}

static int 
make_refl_z(void){
  int dir;
  dir = make_gather(map_refl_z, NULL, OWN_INVERSE, ALLOW_EVEN_ODD, SAME_PARITY);
  return dir;
}
/*---------------------------------------------------------------*/
static void 
map_refl_t(int x, int y, int z, int t, int *args, int fb,
	   int *xp, int *yp, int *zp, int *tp){
  *xp = x;  *yp = y;  *zp = z;  *tp = (nt-t)%nt; 
}

static int 
make_refl_t(void){
  int dir;
  dir = make_gather(map_refl_t, NULL, OWN_INVERSE, ALLOW_EVEN_ODD, SAME_PARITY);
  return dir;
}

/* Symmetrize the correlator */
/* This procedure assumes nx = ny = nz */
void
symmetrize(Real *q, Real *q2){
  Real *qxy, *qxz, *qyz;
  Real *qxy2, *qxz2, *qyz2;
  int dirxy, dirxz, diryz;
  msg_tag *mtagxy, *mtagxz, *mtagyz;
  msg_tag *mtagxy2, *mtagxz2, *mtagyz2;
  int i;

  /*******************************************************/
  /* make xy and xz permutation mappings                 */
  /*******************************************************/
  dirxy = make_perm_xy();
  dirxz = make_perm_xz();

  qxy = create_r_field();
  qxz = create_r_field();

  qxy2 = create_r_field();
  qxz2 = create_r_field();
  
  /* Permutation x-y */
  /* Defines qxy to be the x-y permutation of q */

  mtagxy = start_gather_field(q, sizeof(Real), dirxy, EVENANDODD, gen_pt[0]);
  mtagxz = start_gather_field(q, sizeof(Real), dirxz, EVENANDODD, gen_pt[1]);

  mtagxy2 = start_gather_field(q2, sizeof(Real), dirxy, EVENANDODD, gen_pt[2]);
  mtagxz2 = start_gather_field(q2, sizeof(Real), dirxz, EVENANDODD, gen_pt[3]);

  wait_gather(mtagxy);
  FORALLFIELDSITES(i){
    qxy[i] = *((Real *)gen_pt[0][i]);
  }

  wait_gather(mtagxz);
  FORALLFIELDSITES(i){
    qxz[i] = *((Real *)gen_pt[1][i]);
  }
  wait_gather(mtagxy2);
  FORALLFIELDSITES(i){
    qxy2[i] = *((Real *)gen_pt[2][i]);
  }

  wait_gather(mtagxz2);
  FORALLFIELDSITES(i){
    qxz2[i] = *((Real *)gen_pt[3][i]);
  }

  FORALLFIELDSITES(i){
    q[i] = (q[i] + qxy[i] + qxz[i])/3.;
    q2[i] = (q2[i] + qxy2[i] + qxz2[i])/3.;
  }

  cleanup_gather(mtagxy);
  cleanup_gather(mtagxz);

  cleanup_gather(mtagxy2);
  cleanup_gather(mtagxz2);

  destroy_r_field(qxy);
  destroy_r_field(qxz);

  destroy_r_field(qxy2);
  destroy_r_field(qxz2);

  /*******************************************************/
  /* yz permutations                                     */
  /*******************************************************/

  diryz = make_perm_yz();
  qyz = create_r_field();
  qyz2 = create_r_field();
  mtagyz = start_gather_field(q, sizeof(Real), diryz, EVENANDODD, gen_pt[0]);
  mtagyz2 = start_gather_field(q2, sizeof(Real), diryz, EVENANDODD, gen_pt[1]);

  wait_gather(mtagyz);
  FORALLFIELDSITES(i){
    qyz[i] = *((Real *)gen_pt[0][i]);
  }

  wait_gather(mtagyz2);
  FORALLFIELDSITES(i){
    qyz2[i] = *((Real *)gen_pt[1][i]);
  }

  FORALLFIELDSITES(i){
    q[i] = (q[i] + qyz[i])/2.;
    q2[i] = (q2[i] + qyz2[i])/2.;
  }

  cleanup_gather(mtagyz);
  cleanup_gather(mtagyz2);

  destroy_r_field(qyz);
  destroy_r_field(qyz2);

  /*******************************************************/
  /* x reflections                                       */
  /*******************************************************/

  Real *qxp, *qyp, *qzp, *qtp;
  Real *qxp2, *qyp2, *qzp2, *qtp2;
  int dirxp, diryp, dirzp, dirtp;
  msg_tag *mtagxp, *mtagyp, *mtagzp, *mtagtp;
  msg_tag *mtagxp2, *mtagyp2, *mtagzp2, *mtagtp2;

  /* Reflection x */
  dirxp = make_refl_x();
  qxp = create_r_field();
  qxp2 = create_r_field();
  mtagxp = start_gather_field(q, sizeof(Real), dirxp, EVENANDODD, gen_pt[0]);
  mtagxp2 = start_gather_field(q2, sizeof(Real), dirxp, EVENANDODD, gen_pt[1]);
  wait_gather(mtagxp);
  FORALLFIELDSITES(i){
    qxp[i] = *((Real *)gen_pt[0][i]);
  }
  wait_gather(mtagxp2);
  FORALLFIELDSITES(i){
    qxp2[i] = *((Real *)gen_pt[1][i]);
  }
  FORALLFIELDSITES(i){
    q[i] = (q[i] + qxp[i])/2.;
    q2[i] = (q2[i] + qxp2[i])/2.;
  }
  cleanup_gather(mtagxp);
  cleanup_gather(mtagxp2);
  destroy_r_field(qxp);
  destroy_r_field(qxp2);

  /* Reflection y */
  diryp = make_refl_y();
  qyp = create_r_field();
  qyp2 = create_r_field();
  mtagyp = start_gather_field(q, sizeof(Real), diryp, EVENANDODD, gen_pt[0]);
  mtagyp2 = start_gather_field(q2, sizeof(Real), diryp, EVENANDODD, gen_pt[1]);
  wait_gather(mtagyp);
  FORALLFIELDSITES(i){
    qyp[i] = *((Real *)gen_pt[0][i]);
  }
  wait_gather(mtagyp2);
  FORALLFIELDSITES(i){
    qyp2[i] = *((Real *)gen_pt[1][i]);
  }
  FORALLFIELDSITES(i){
    q[i] = (q[i] + qyp[i])/2.;
    q2[i] = (q2[i] + qyp2[i])/2.;
  }
  cleanup_gather(mtagyp);
  cleanup_gather(mtagyp2);
  destroy_r_field(qyp);
  destroy_r_field(qyp2);

  /* Reflection z */
  dirzp = make_refl_z();
  qzp = create_r_field();
  qzp2 = create_r_field();
  mtagzp = start_gather_field(q, sizeof(Real), dirzp, EVENANDODD, gen_pt[0]);
  mtagzp2 = start_gather_field(q2, sizeof(Real), dirzp, EVENANDODD, gen_pt[1]);
  wait_gather(mtagzp);
  FORALLFIELDSITES(i){
    qzp[i] = *((Real *)gen_pt[0][i]);
  }
  wait_gather(mtagzp2);
  FORALLFIELDSITES(i){
    qzp2[i] = *((Real *)gen_pt[1][i]);
  }
  FORALLFIELDSITES(i){
    q[i] = (q[i] + qzp[i])/2.;
    q2[i] = (q2[i] + qzp2[i])/2.;
  }
  cleanup_gather(mtagzp);
  cleanup_gather(mtagzp2);
  destroy_r_field(qzp);
  destroy_r_field(qzp2);

  /* Reflection t */
  dirtp = make_refl_t();
  qtp = create_r_field();
  qtp2 = create_r_field();
  mtagtp = start_gather_field(q, sizeof(Real), dirtp, EVENANDODD, gen_pt[0]);
  mtagtp2 = start_gather_field(q2, sizeof(Real), dirtp, EVENANDODD, gen_pt[1]);
  wait_gather(mtagtp);
  FORALLFIELDSITES(i){
    qtp[i] = *((Real *)gen_pt[0][i]);
  }
  wait_gather(mtagtp2);
  FORALLFIELDSITES(i){
    qtp2[i] = *((Real *)gen_pt[1][i]);
  }
  FORALLFIELDSITES(i){
    q[i] = (q[i] + qtp[i])/2.;
    q2[i] = (q2[i] + qtp2[i])/2.;
  }
  cleanup_gather(mtagtp);
  cleanup_gather(mtagtp2);
  destroy_r_field(qtp);
  destroy_r_field(qtp2);
} /* symmetrize.c */


