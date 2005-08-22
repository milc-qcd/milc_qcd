/************************* RG_read.c *******************************/
/* MIMD version 7 */
/* F. Maresca Jul 2005 */
/************************* control.c *******************************/
/* MIMD version 7 */
/* Main procedure for SU3 with dynamical staggered fermions        */
/* general quark action, general gauge action */

/* This version combines code for the PHI algorithm (approriate for 4
   flavors) and the R algorithm for "epsilon squared" updating of 
   1 to 4 flavors.  Compilation should occur with PHI_ALGORITHM defined
   for the former and not defined for the latter.  It also contains code
   for the hybrid Monte Carlo algorithm, for which HMC_ALGORITHM and
   PHI_ALGORITHM should be defined.  (Actually, the
   changes to control.c are minimal and the real differences will appear
   in update.c */

#define CONTROL
#include "RG_Shamir_includes.h"	/* definitions files and prototypes */
#include "RG_include.h"
EXTERN gauge_header start_lat_hdr;	/* Input gauge field header */

#define no_smear_level 6
#define NTOL 1e-06



EXTERN int rvs[16][4]={{0,0,0,0},{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1},{1,1,0,0},{1,1,1,0},{1,1,1,1},{0,1,1,0},{0,1,1,1},{0,0,1,1},{0,1,0,1},{1,0,0,1},{1,1,0,1},{1,0,1,1},{1,0,1,0}};

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
  int i,j,S,T,etap,eta,cmp[3];
  QDP_ColorVector *phi_ST,*phi,*phi_c,*tmp;
  QDP_ColorVector *phi_1,*phi_ST_f,*chi[RG_Ncn];
  QDP_ColorMatrix *wlink[NRG][RG_Ncn];
  QDP_Sub_Block QDP_block[NRG+1];
  QDP_Subset *test;
  QDP_Shift offset;
  int meascount,traj_donex,status;
  int prompt,todo,sm_lev;
  int s_iters, avs_iters, avspect_iters, avbcorr_iters;
  double dtime, dclock();
  QDP_String *xml_file_out, *xml_record_out;
  QDP_Writer *outfile;
  QDP_Reader *infile;
  QDP_RandomState *rs;
  QDP_Int *li;
  char file[80] = "file";


  /*  XML strings */
  char xml_write_file[] = "Shamir user file XML";
  char xml_write_record[] = "Value of the field";
 
  initialize_machine(argc,argv);
  QDP_initialize(&argc, &argv);
  g_sync();
  /* set up */
  prompt = setup();

  /* loop over input sets */
  node0_printf("Start \n"); fflush(stdout);

  status = readin(prompt); 

  for (i=0;i<NRG+1;i++)
  RG_create_block(&QDP_block[i],i);


#ifdef CHECK
   node0_printf("Some checks........ \n");fflush(stdout);
   if ( RG_check(QDP_block) != 0 ) exit(-1);
   node0_printf(".......... Some checks end \n");fflush(stdout);
#else 
#ifdef PRINT
  /* Create the file XML and the records XML*/
  xml_file_out = QDP_string_create();
  QDP_string_set(xml_file_out,xml_write_file);
  xml_record_out = QDP_string_create();
  QDP_string_set(xml_record_out,xml_write_record);

  node0_printf("Opened file for writing\n");
  outfile = QDP_open_write(xml_file_out, file, QDP_SINGLEFILE);
  node0_printf("Opened file for writing\n");
#endif

  for (i = 0; i < NRG; i++)
  for (j = 0; j < RG_Ncn; j++)
   wlink[i][j] = QDP_create_M();

  RG_setup(QDP_block,wlink);

  phi = QDP_create_V();
  tmp = QDP_create_V();
  phi_c = QDP_create_V();
  phi_ST = QDP_create_V();
  phi_ST_f = QDP_create_V();

  for (i=0; i< RG_Ncn; i++)  
   chi[i] = QDP_create_V();

  node0_printf("Color source set to one at the origin\n");

  SQDP_V_eq_func(phi, point_V, QDP_block[0]);
//  SQDP_V_eq_func(phi, point_2V, QDP_block[0]);
  SQDP_V_eq_func(phi, print_cv, QDP_block[0]);

  RG_coarse_to_fine(phi_c,QDP_block,phi,wlink);
  
  //SQDP_V_eq_func(phi_c, print_cv, QDP_block[NRG]);
  
  node0_printf("wils-wils completed\n");
  fflush(stdout);
 
 // valid_fatlinks = valid_longlinks = 0;

/* Check with a random source */
/*
  li = QDP_create_I();
  rs = QDP_create_S();
  QDP_I_eq_func(li, initial_li, QDP_all);
  QDP_S_eq_seed_i_I(rs, 987654321, li, QDP_all);
  QDP_V_eq_gaussian_S(phi_c,rs,QDP_all);
*/

  RG_bulk(chi,phi_c);

  node0_printf("16 inversions completed\n");
  fflush(stdout);
/*  
  for (eta=0; eta < RG_Ncn; eta++)
  {
  fprintf(stderr,"eta %d\n",eta);
  QDP_V_eq_func(chi[eta], print_cv, QDP_all);
  }
*/


#ifdef CHECK_INV
RG_check_inversion(phi_c,chi);
#endif

  for (S=0 ;S < RG_Ncn; S++)
  for (T=0 ;T < RG_Ncn; T++)
  {

  SQDP_V_eq_zero(phi_ST,QDP_block[NRG-1]);
  SQDP_V_eq_zero(phi_ST_f,QDP_block[NRG]);

  cmp[0]=S; cmp[1]=T; 

  for (etap=0; etap < RG_Ncn; etap++)
  {
   cmp[2]=etap;

   SQDP_V_eq_zero(tmp,QDP_block[NRG-1]);
   RG_decompose(tmp,chi,cmp,wlink[NRG-1],QDP_block[NRG-1]);

   SQDP_V_peq_V(phi_ST,tmp,QDP_block[NRG-1]);
  }

  //printf("S = %d T = %d:  this node %d\n",S,T,this_node);
  //SQDP_V_eq_func(phi_ST, print_cv, QDP_block[NRG-1]);
  RG_fine_to_coarse(phi_ST_f, QDP_block, phi_ST, wlink);
 
#ifdef PRINT
   status = QDP_write_V(outfile, xml_record_out, phi_ST_f);
   if (status != 0) exit(-1);
#endif

  //node0_printf("S = %d T = %d\n",S,T);
  //SQDP_V_eq_func(phi_ST_f, print_cv, QDP_block[0]);

  }


  QDP_destroy_V(phi);
  QDP_destroy_V(tmp);
  QDP_destroy_V(phi_c);
  QDP_destroy_V(phi_ST);
  QDP_destroy_V(phi_ST_f);

  for (i=0; i< RG_Ncn; i++)  
    QDP_destroy_V(chi[i]);

  for (i = 0; i < NRG; i++)
  for (j = 0; j < RG_Ncn; j++)
   QDP_destroy_M(wlink[i][j]);

#ifdef PRINT
  QDP_close_write(outfile);

  phi = QDP_create_V();

  infile = QDP_open_read(xml_file_out, file);
  if(infile == NULL){
    printf("Quitting... error\n");
    return 1;
  }

  node0_printf("Opened file for reading\n");

  for ( S = 0; S < RG_Ncn; S++)
  for ( T = 0; T < RG_Ncn; T++)
  {
   //printf("Record %d %d: this node\n",S,T,this_node);
   //printf("Start reading: this node %d\n",this_node);
   //QIO_verbose(QIO_VERB_DEBUG);
   status = QDP_read_V(infile, xml_record_out, phi);
   //printf("End reading , status %d: this node %d\n",status,this_node);
   if (status != 0)
   {
   printf("Error on record %d %d : status %d: this node %d\n",S,T,status,this_node);
   fflush(stdout);
   terminate(-1);
   }
   
   node0_printf("Record %d %d\n",S,T);
   fflush(stdout);
   SQDP_V_eq_func(phi, print_cv, QDP_block[0]);
  }

  QDP_close_read(infile);



  QDP_string_destroy(xml_file_out);
  QDP_string_destroy(xml_record_out);
  QDP_destroy_V(phi);
#endif

#endif 


  node0_printf("QDP finalize\n"); fflush(stdout);
  QDP_finalize();
  node0_printf("QDP finalize completed\n"); fflush(stdout);
  normal_exit(0);

  return 0;
}
