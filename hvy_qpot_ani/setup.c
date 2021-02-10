/******** setup.c *********/
/* MIMD version 7 */
/* original code by UMH */
/* 2/19/98 Version 5 port CD */

#define IF_OK if(status==0)

#include "hvy_qpot_includes.h"
#include "../include/dirs.h"
#include <string.h>
#include <ctype.h>

/* Each node has a params structure for passing simulation parameters */
#include "params.h"
params param;


/* Forward declarations */
static int initial_set(void);
static void next_local_lattice(int x, int y, int z, int t, int *dirpt, int FB,
                               int *xp, int *yp, int *zp, int *tp);
static void make_nll_gathers(void);

int  setup()   {
  /* print banner, get volume */
  int prompt=initial_set();
  if(prompt == 2)return prompt;
  /* Initialize the layout functions, which decide where sites live, and 
     initialize the node random number generator */
  setup_layout();
  this_node = mynode();
  /* allocate space for lattice, set up coordinate fields */
  make_lattice();
  /* set up nearest neighbor pointers and comlink structures */
  make_nn_gathers();

  /* set up pointers for passing full local lattices and 
     comlink structures code for this routine is below */
  make_nll_gathers();
  return(prompt);
}

/* SETUP ROUTINES */
static int initial_set(){
  int prompt=0,status;

  /* On node zero, read lattice size, seed, nflavors and send to others */
  if(mynode()==0){
    /* print banner */
    printf("Pure gauge SU3 Wilson loops\n");
    printf("MIMD version 7\n");
    printf("Machine = %s, with %d nodes\n",machine_type(),numnodes());
    time_stamp("start");

    status=get_prompt(stdin, &prompt);
    IF_OK status += get_i(stdin, prompt,"nx",&param.nx);
    IF_OK status += get_i(stdin, prompt,"ny",&param.ny);
    IF_OK status += get_i(stdin, prompt,"nz",&param.nz);
    IF_OK status += get_i(stdin, prompt,"nt",&param.nt);

    if(status>0) param.stopflag=1; else param.stopflag=0;
  } /* end if(mynode()==0) */

  /* Node 0 broadcasts parameter buffer to all other nodes */
  broadcast_bytes((char *)&param,sizeof(param));

  if( param.stopflag != 0 )
    normal_exit(0);

  nx=param.nx;
  ny=param.ny;
  nz=param.nz;
  nt=param.nt;
    
  number_of_nodes = numnodes();
  volume=(size_t)(nx*ny*nz*nt);
  return(prompt);
}

/* read in parameters and coupling constants    */
int readin(int prompt){
/* read in parameters for su3 monte carlo       */
/* argument "prompt" is 1 if prompts are to be given for input  */

  int status;
  int mu;
  char savebuf[128];

    /* On node zero, read parameters and send to all other nodes */
  if(this_node==0){

    printf("\n\n");
    status = 0;
    
    /* find out what kind of starting lattice to use */
    IF_OK status += ask_starting_lattice(stdin,  prompt, &(param.startflag),
        param.startfile );
    IF_OK status += get_f(stdin, prompt,"u0", &param.u0 );
#ifdef COULOMB
    IF_OK if (prompt==1)
      printf("enter 'no_gauge_fix', or 'coulomb_gauge_fix'\n");
    IF_OK scanf("%s",savebuf);
    IF_OK printf("%s\n",savebuf);
    IF_OK {
      if (strcmp("coulomb_gauge_fix",savebuf) == 0 ){
        param.fixflag = COULOMB_GAUGE_FIX;
      } else {
        if(strcmp("no_gauge_fix",savebuf) == 0 ) {
          param.fixflag = NO_GAUGE_FIX;
        } else{
          printf("error in input: fixing_command %s is invalid\n",savebuf); status++;
        }
      }
    }
#else
    IF_OK param.fixflag = NO_GAUGE_FIX;
#endif
    /* find out what to do with lattice at end */
    IF_OK status += ask_ending_lattice(stdin,  prompt, &(param.saveflag),
        param.savefile );
    IF_OK status += ask_ildg_LFN(stdin,  prompt, param.saveflag,
				      param.stringLFN );

#ifdef ANISOTROPY
    /* Direction of anisotropy */
    IF_OK status += get_s(stdin, prompt,"ani_dir",savebuf);
    IF_OK param.ani_dir = dirchar2index( savebuf[0], &status);
    /* Heavy quark anisotropy */
    IF_OK status += get_f(stdin, prompt, "ani_xiq", &param.ani_xiq);
#ifdef ONEDIM_ANISO_TEST
    /* Heavy quark isotropic link factor for debugging */
    IF_OK status += get_f(stdin, prompt, "iso_xiq", &param.iso_xiq);
#endif
#endif

    /* cor-time direction for the heavy-quark correlation function */
    IF_OK status += get_s(stdin, prompt,"cor_dir",savebuf);
    IF_OK param.cor_dir = dirchar2index( savebuf[0], &status);

#ifndef GFIXONLY
#if ( ( ( defined COULOMB || defined WLOOP_MEAS) \
     && ( ! ( defined PLANE_AVERAGED_PLC ) ) ) )
    /* minimum cor-time value for loops */
    IF_OK status += get_i(stdin, prompt,"min_ct",&param.min_ct);
#endif // END #if (defined COULOMB || defined WLOOP_MEAS)
#if (defined COULOMB || defined WLOOP_MEAS || defined PLOOPCOR_EXCLUSIVE)
    /* maximum cor-spatial distance */
    IF_OK status += get_i(stdin, prompt,"max_x",&param.max_x);
    IF_OK status += get_i(stdin, prompt,"max_y",&param.max_y);
    IF_OK status += get_i(stdin, prompt,"max_z",&param.max_z);
    IF_OK status += get_i(stdin, prompt,"max_t",&param.max_t);
#endif // END #if (defined COULOMB || defined WLOOP_MEAS || defined PLOOPCOR_EXCLUSIVE)

#if (defined COULOMB || defined WLOOP_MEAS || defined PLOOPCOR_MEAS)
    IF_OK status += get_f(stdin, prompt,"max_r",&param.max_r);
#endif // END #if (defined COULOMB || defined WLOOP_MEAS || defined PLOOPCOR_EXCLUSIVE)

#ifdef WLOOP_MEAS
    /* off_axis_flag : do off-axis Wilson loops? */
    IF_OK status += get_i(stdin, prompt,"off_axis_flag",&param.off_axis_flag);
#endif // END #ifdef WLOOP_MEAS

#ifndef SMEARING        
    IF_OK param.no_smear_level=1;
    IF_OK param.smear_num[0] = 0;
    IF_OK param.smear_num[1] = 0;
    IF_OK param.smear_num[2] = 0;
    IF_OK param.smear_num[3] = 0;
    IF_OK param.smear_num[4] = 0;
#else 
    IF_OK status += get_i(stdin, prompt,"no_smear_level",
                              &param.no_smear_level);
    if (param.no_smear_level > 0) {
      IF_OK status += get_i(stdin, prompt,"smear_num[0]",&param.smear_num[0]);
      IF_OK status += get_i(stdin, prompt,"smear_num[1]",&param.smear_num[1]);
      IF_OK status += get_i(stdin, prompt,"smear_num[2]",&param.smear_num[2]);
      IF_OK status += get_i(stdin, prompt,"smear_num[3]",&param.smear_num[3]);
      IF_OK status += get_i(stdin, prompt,"smear_num[4]",&param.smear_num[4]);

#ifdef HYP_4D_SMEARING /* ONLY in 4D HYP smearing */
      IF_OK status += get_f(stdin, prompt,"hyp_alpha1",&param.hyp_alpha1); 
#elif defined HYP_3D_SMEARING /* SET TO ZERO for 3D HYP smearing */
      IF_OK param.hyp_alpha1=0.;
#endif // END #ifdef HYP_4D_SMEARING
#ifdef HYP_SMEARING /* BOTH in 3D and 4D HYP smearing */
      IF_OK status += get_f(stdin, prompt,"hyp_alpha2",&param.hyp_alpha2);
      IF_OK status += get_f(stdin, prompt,"hyp_alpha3",&param.hyp_alpha3);
#else /* spatial APE smearing */
      IF_OK status += get_f(stdin, prompt,"smear_fac",&param.smear_fac);
#if (defined APE_1D_SMEARING || defined APE_1D2_SMEARING) 
      IF_OK status += get_s(stdin, prompt,"stap_dir",savebuf);
      IF_OK param.stap_dir = dirchar2index( savebuf[0],&status);
#endif
#endif // END #ifdef HYP_SMEARING 
    }
#endif // END #ifndef SMEARING
#ifdef NEW_HVY_POT 
    IF_OK status += get_i(stdin, prompt,"hqp_alg",&param.hqp_alg);
#endif
#endif // END #ifndef GFIXONLY

    if( status > 0)param.stopflag=1; else param.stopflag=0;
  } /* end if(this_node==0) */

  /* Node 0 broadcasts parameter buffer to all other nodes */
  broadcast_bytes((char *)&param,sizeof(param));

  if( param.stopflag != 0 ) normal_exit(0);

  if(prompt==2)return prompt;

  fixflag = param.fixflag;
  u0 = param.u0;

#if (defined COULOMB || defined WLOOP_MEAS)
  min_ct = param.min_ct;
#endif // END #if (defined COULOMB || defined WLOOP_MEAS)
#if (defined COULOMB || defined WLOOP_MEAS || defined PLOOPCOR_EXCLUSIVE)
  max_x = param.max_x;
  max_y = param.max_y;
  max_z = param.max_z;
  max_t = param.max_t;
#endif // END #if (defined COULOMB || PLOOPCOR_MEAS)

  /* Set initial values in xc[], nc[] and maxc[] arrays assuming standard geometry */
  xc[XUP]=XUP; xc[YUP]=YUP; xc[ZUP]=ZUP; xc[TUP]=TUP;
  nc[XUP]=nx; nc[YUP]=ny; nc[ZUP]=nz; nc[TUP]=nt;
  maxc[XUP]=max_x; maxc[YUP]=max_y; maxc[ZUP]=max_z; maxc[TUP]=max_t;

  /* Initialize cor-time direction */
  cor_dir = param.cor_dir;
  if ( cor_dir != TUP ) { 
    xc  [TUP] = cor_dir; 
    nc  [TUP] = nc[cor_dir]; 
    maxc[TUP] = maxc[cor_dir]; 
    xc  [cor_dir] = TUP; 
    nc  [cor_dir] = nt; 
    maxc[cor_dir] = max_t;
  }
#ifdef ANISOTROPY
  ani_dir = param.ani_dir;
  if ( ani_dir != cor_dir && ani_dir != TUP ) { 
    int tmp;
    tmp = xc[XUP];
    xc  [XUP] = ani_dir; 
    xc[ani_dir] = tmp; 
    tmp = nc[XUP];
    nc  [XUP] = nc[ani_dir]; 
    nc  [ani_dir] = tmp; 
    tmp =maxc[XUP];
    maxc[XUP] = maxc[ani_dir]; 
    maxc[ani_dir] = tmp;
  }
  ani_xiq = param.ani_xiq;
#ifdef ONEDIM_ANISO_TEST
  iso_xiq = param.iso_xiq;
#endif
#endif
  for ( mu=XUP; mu <=TUP; mu++) xc[OPP_DIR(mu)] = OPP_DIR(xc[mu]);
  
#if (defined COULOMB || defined WLOOP_MEAS || defined PLOOPCOR_EXCLUSIVE)
  max_r = param.max_r;
#endif // END #if (defined COULOMB || PLOOPCOR_MEAS)

#ifdef WLOOP_MEAS
  off_axis_flag = param.off_axis_flag;
#endif

  no_smear_level = param.no_smear_level;
  smear_num[0] = param.smear_num[0];
  smear_num[1] = param.smear_num[1];
  smear_num[2] = param.smear_num[2];
  smear_num[3] = param.smear_num[3];
  smear_num[4] = param.smear_num[4];
#ifdef HYP_SMEARING
  hyp_alpha1 = param.hyp_alpha1;
  hyp_alpha2 = param.hyp_alpha2;
  hyp_alpha3 = param.hyp_alpha3;
#else
  staple_weight = param.staple_weight;
  ape_iter = param.ape_iter;
#endif
  smear_fac = param.smear_fac;
#if (defined APE_1D_SMEARING || defined APE_1D2_SMEARING) 
  stap_dir = param.stap_dir;
  if ( maxc[stap_dir] > 0 ) { 
    node0_printf("readin: inconsistent choice of stap_dir=%d and maxc[%d]=%d\n",stap_dir,stap_dir,maxc[stap_dir]);
    terminate(1);
  }
#endif
#ifdef NEW_HVY_POT 
  hqp_alg = param.hqp_alg;
#endif

  startflag = param.startflag;
  saveflag = param.saveflag;
  strcpy(startfile,param.startfile);
  strcpy(savefile,param.savefile);
  strcpy(stringLFN, param.stringLFN);

  /* Do whatever is needed to get lattice */
  if( startflag != CONTINUE )
    startlat_p = reload_lattice( startflag, startfile );

  return(0);
}

/* Set up comlink structures for nearest local lattice gather pattern; 
   make_lattice() and  make_nn_gathers() must be called first, 
   preferably just before calling make_3n_gathers().
*/
static void
make_nll_gathers(void)
{
  int i;
  int x,y,z,t;
  int of=0;
  int ox,oy,oz,ot;
  int ix[4]={0,0,0,0};

  /* get the displacement of local lattices, assumes 
     that all local lattices have the same size */
  for(t=0;t<nt;t++)for(z=0;z<nz;z++)for(y=0;y<ny;y++)for(x=0;x<nx;x++){
    if(node_number(x,y,z,t)==mynode()){
      if (!of) { of=1; ox=x; oy=y; oz=z; ot=t; }
      if (ox==x&&oy==y&&oz==z) { ix[TUP]++;}
      if (ot==t&&ox==x&&oy==y) { ix[ZUP]++;}
      if (oz==z&&ot==t&&ox==x) { ix[YUP]++;}
      if (oy==y&&oz==z&&ot==t) { ix[XUP]++;}
    }
  }
  for(i=XUP; i<=TUP; i++) { locx[i]=ix[i]; }
  
  /* Wouldn't need the inverse if the shifts traveled around
     the lattice always in positive directions only */
  /* However, that is pretty inefficient for max_x << nxyz/2,
     for which I should introduce the inverse */
  /* keep in mind that parity in SAME for even displaments 
     and SWITCH for odd displacements */
  for(i=XUP; i<=TUP; i++) {
    if (locx[i]&1) {
      make_gather(next_local_lattice, &i, WANT_INVERSE,
                ALLOW_EVEN_ODD, SWITCH_PARITY);
    }
    else {
      make_gather(next_local_lattice, &i, WANT_INVERSE,
                ALLOW_EVEN_ODD, SAME_PARITY);
    }
  }
  
  /* need that only due to use of inverse */
  sort_eight_gathers(NLLXUP);
}

/* this routine uses only fundamental directions (XUP..TDOWN) as directions */
/* returning the displaccement of the origin of the next sublattice in that direction */

static void
next_local_lattice(int x, int y, int z, int t, int *dirpt, int FB,
               int *xp, int *yp, int *zp, int *tp)
     /* int x,y,z,t,*dirpt,FB;  coordinates of site, direction (eg XUP), and
        "forwards/backwards"  */
     /* int *xp,*yp,*zp,*tp;    pointers to coordinates of neighbor */
{
  int dir;
  dir = (FB==FORWARDS) ? *dirpt : OPP_DIR(*dirpt);
  *xp = x; *yp = y; *zp = z; *tp = t;
  switch(dir){
  case XUP: *xp = (x+locx[XUP])%nx; break;
  case XDOWN: *xp = (x+(locx[XUP]+1)*nx-locx[XUP])%nx; break;
  case YUP: *yp = (y+locx[YUP])%ny; break;
  case YDOWN: *yp = (y+(locx[YUP]+1)*ny-locx[YUP])%ny; break;
  case ZUP: *zp = (z+locx[ZUP])%nz; break;
  case ZDOWN: *zp = (z+(locx[ZUP]+1)*nz-locx[ZUP])%nz; break;
  case TUP: *tp = (t+locx[TUP])%nt; break;
  case TDOWN: *tp = (t+(locx[TUP]+1)*nt-locx[TUP])%nt; break;
  default: printf("next_local_lattice: bad direction\n"); exit(1);
  }
}
