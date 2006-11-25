/************************* control.c *******************************/
/* MIMD version 7 */
/* Main procedure for RG blocking a la Shamir        */
/* general quark action, general gauge action */
/* F. Maresca Jul 2005 */

#define CONTROL
#include "RG_Shamir_includes.h"	/* definitions files and prototypes */
#include "RG_include.h"
EXTERN gauge_header start_lat_hdr;	/* Input gauge field header */

#define NTOL 1e-06


//EXTERN int rvs[16][4]={{0,0,0,0},{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1},{1,1,0,0},{1,1,1,0},{1,1,1,1},{0,1,1,0},{0,1,1,1},{0,0,1,1},{0,1,0,1},{1,0,0,1},{1,1,0,1},{1,0,1,1},{1,0,1,0}};

EXTERN int rvs[16][4]={{0,0,0,0},{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1},{1,1,0,0},{1,0,1,0},{1,0,0,1},{0,1,1,0},{0,1,0,1},{0,0,1,1},{1,1,1,0},{1,1,0,1},{1,0,1,1},{0,1,1,1},{1,1,1,1}};


void
initial_li(QLA_Int *li, int coords[])
{
  int i,t;

  t = coords[0];
  for(i=1; i<4; i++) {
    t = t*8 + coords[i];
  }
  *li = t;
}

int
main( int argc, char **argv )
{
  int i,j,etap,eta,cmp[3];
  QDP_ColorVector *phi_ST,*phi,*phi_c,*tmp;
  QDP_ColorVector *phi_1,*phi_ST_f,*chi[RG_Ncn];
  QDP_ColorMatrix *wlink[NRG][RG_Ncn];
  QDP_Sub_Block QDP_block[NRG+1];
  QDP_Subset *test;
  QDP_Shift offset;
  int prompt,status;
  double dtimed, dclock();
  QDP_RandomState *rs;
  QDP_Int *li;
  
  
  initialize_machine(&argc,&argv);
#ifdef HAVE_QDP
  QDP_initialize(&argc, &argv);
#endif
  /* Remap standard I/O */
  if(remap_stdio_from_args(argc, argv) == 1)terminate(1);

  g_sync();
  /* set up */
  prompt = setup();
  
  /* loop over input sets */
  node0_printf("Start \n"); fflush(stdout);
  
  dtimed = -dclock();

  /* Read in the remaining parameters  */
  status = readin(prompt); 
  
  for (i=0;i<nrg+1;i++)
    RG_create_block(&QDP_block[i],i);
  
#ifdef CHECK
  node0_printf("Some checks........ \n");fflush(stdout);
#ifdef CHECK_PATH
  if ( RG_check_path() != 0 )
    printf("FAILED!!!\n");
#endif
  if ( RG_check(QDP_block) != 0 ) exit(-1);
  node0_printf(".......... Some checks end \n");fflush(stdout);
#else 
  
#ifndef NOTRANS
  //  for (i = 0; i < nrg; i++)
  //    for (j = 0; j < RG_Ncn; j++)
  //      wlink[i][j] = QDP_create_M();
#endif
  
  /* Define the parallel transporters for all RG blockings */
  RG_setup(QDP_block,wlink);
  
  phi = QDP_create_V();
  tmp = QDP_create_V();
  phi_c = QDP_create_V();
  phi_ST = QDP_create_V();
  phi_ST_f = QDP_create_V();

  for (i=0; i< RG_Ncn; i++)  
  {
   chi[i] = QDP_create_V();
   QDP_V_eq_zero(chi[i],QDP_all);
  }
  node0_printf("Color source set to one at the origin\n");
  fflush(stdout);

  SQDP_V_eq_func(phi, point_V, QDP_block[0]);
  SQDP_V_eq_func(phi, print_cv, QDP_block[0]);

  QDP_V_eq_zero(phi_c,QDP_all);
  /* Do inverse RG blocking on phi coarse to fine to get phi_c */
  RG_coarse_to_fine(phi_c,QDP_block,phi,wlink);
  node0_printf("wils-wils completed\n");
  fflush(stdout);
 
/* Check with a random source */
/*
  li = QDP_create_I();
  rs = QDP_create_S();
  QDP_I_eq_func(li, initial_li, QDP_all);
  QDP_S_eq_seed_i_I(rs, 987654321, li, QDP_all);
  QDP_V_eq_gaussian_S(phi_c,rs,QDP_all);
*/

  /* Compute chi[i] = M^{-1} phi_c|_i the field
     phi_c restricted to hypercube offsets i */
  RG_bulk(chi,phi_c);

  node0_printf("16 inversions completed\n");
  fflush(stdout);

#ifdef CHECK_INV
  RG_check_inversion(phi_c,chi);
#endif

  /* Iterate over all spin and taste labels S and T */
  for (S=0 ;S < RG_Ncn; S++)
    for (T=0 ;T < RG_Ncn; T++)
      {
	
	QDP_V_eq_zero(phi_ST,QDP_all);
	QDP_V_eq_zero(phi_ST_f,QDP_all);
	
	cmp[0]=S; cmp[1]=T; 
	
	/* Iterate over all hypercube sites */
	for (etap=0; etap < RG_Ncn; etap++)
	  {
	    cmp[2]=etap;
	    
	    QDP_V_eq_zero(tmp,QDP_all);
	    /* Project out the spin and taste component for S, T 
	       result in phi_ST at fine hypercube origins */
	    RG_decompose(tmp,chi,cmp,wlink[nrg-1],QDP_block[nrg-1]);
	    SQDP_V_peq_V(phi_ST,tmp,QDP_block[nrg-1]);
	  }
	
	//printf("S = %d T = %d:  this node %d\n",S,T,this_node);
	//SQDP_V_eq_func(phi_ST, print_cv, QDP_block[nrg-1]);
	/* RG block phi_ST from fine to coarse - result in phi_ST_f */
	RG_fine_to_coarse(phi_ST_f, QDP_block, phi_ST, wlink);
	
	/* Print result */
	node0_printf("S = %d T = %d  nrg = %d\n",S,T,nrg); fflush(stdout);
	//SQDP_V_eq_func(phi_ST_f, print_cv_node, QDP_block[0]);
	SQDP_V_eq_func(phi_ST_f, print_cv_node_tot, QDP_block[0]);
      }
  
  
  QDP_destroy_V(phi);
  QDP_destroy_V(tmp);
  QDP_destroy_V(phi_c);
  QDP_destroy_V(phi_ST);
  QDP_destroy_V(phi_ST_f);

  for (i=0; i< RG_Ncn; i++)  
    QDP_destroy_V(chi[i]);

#ifndef NOTRANS
  for (i = 0; i < nrg; i++)
  for (j = 0; j < RG_Ncn; j++)
   QDP_destroy_M(wlink[i][j]);
#endif

  dtimed += dclock();
  node0_printf("Timing: %lf\n",dtimed); fflush(stdout);
#endif

  printf("QDP finalize: this node:%d\n",this_node); fflush(stdout);
  QDP_finalize();
  printf("QDP finalize completed: this node %d\n",this_node); fflush(stdout);
  normal_exit(0);

  return 0;
}
