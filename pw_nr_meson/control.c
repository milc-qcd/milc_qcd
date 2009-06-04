/***************** control_pw_mes.c ****************************************/

/* Main procedure for P-wave nonrelativistic meson */
/* MIMD version 7 */

/* Comment out if you want to suppress detailed timing */
#define PRTIME

#define CONTROL
#include "pw_nr_meson_includes.h"
#ifdef HAVE_QDP
#include <qdp.h>
#endif

int main(int argc,char *argv[])
{
  int prompt,dir1,dir2;
  int ksnk, ksrc;
  double starttime,endtime,dtime;
  Real space_vol;
  int key[4];
  int slice[4];

  initialize_machine(&argc,&argv);
#ifdef HAVE_QDP
  QDP_initialize(&argc, &argv);
#endif
  /* Remap standard I/O */
  if(remap_stdio_from_args(argc, argv) == 1)terminate(1);

  g_sync();
  /* set up */
  prompt = setup(); 

  space_vol = (Real)(nx*ny*nz);

  while( readin(prompt) == 0){
    starttime=dclock();
    total_iters = 0;
    
    /**************************************************************/
    /* Set up gauge field */

    if( fixflag == COULOMB_GAUGE_FIX)
      {
	if(this_node == 0) 
	  printf("Fixing to Coulomb gauge\n");
	
	STARTTIME;
	gaugefix(TUP,(Real)1.5,500,GAUGE_FIX_TOL);
	ENDTIME("gauge fix");
	invalidate_this_clov(gen_clov);
      }
    else
      if(this_node == 0)printf("COULOMB GAUGE FIXING SKIPPED.\n");
    
    /* save lattice if requested */
    if( saveflag != FORGET ){
      savelat_p = save_lattice( saveflag, savefile, stringLFN );
    }

    /**************************************************************/
    /* Allocate space for quark propagators */

    quark_prop = 
      (block_pauli_propagator *)malloc(sites_on_node*
				       sizeof(block_pauli_propagator));
    if(quark_prop == NULL){
      printf("main(%d): No room for prop\n",this_node);
      terminate(1);
    }

    quark_prop_smear = 
      (block_pauli_propagator *)malloc(sites_on_node*
				       sizeof(block_pauli_propagator));
    if(quark_prop_smear == NULL){
      printf("main(%d): No room for prop\n",this_node);
      terminate(1);
    }

    antiquark_prop = 
      (block_pauli_propagator *)malloc(sites_on_node*
				       sizeof(block_pauli_propagator));
    if(antiquark_prop == NULL){
      printf("main(%d): No room for prop\n",this_node);
      terminate(1);
    }

    /**************************************************************/
    /* Read and/or generate antiquark propagator in Pauli basis */

    total_iters += get_wprop_to_field(a_startflag_w, a_startfile_w, 
				      a_saveflag_w, a_savefile_w,
				      antiquark_prop, &a_wqs, &qic, &dcp);
    
    /* Fourier transform antiquark */

    key[XUP] = 1;
    key[YUP] = 1;
    key[ZUP] = 1;
    key[TUP] = 0;
    
    setup_restrict_fourier(key, slice);

    STARTTIME;
    restrict_fourier_field((complex *)antiquark_prop,
			   sizeof(block_pauli_propagator), FORWARDS);
    ENDTIME("Fourier transform antiquark");

    /**************************************************************/
    /* Set up smearing functions for source and sink */

    STARTTIME;
    load_smearing(source_wqs, sink_wqs, num_smear);
    ENDTIME("load smearing");

    //_______________________________________________________________________

    /* Generate the meson propagators for each pairing of source and sink
       smearing and each pairing of source and sink P-wave components */

    for(ksnk=0; ksnk<num_smear;ksnk++){
      for(ksrc=0; ksrc<num_smear; ksrc++){

	/* Clear the meson propagators*/
	init_pw_prop();

	for(dir1=0;dir1<3;dir1++){
	  
	  /* Construct P-wave source for quark */
	  STARTTIME;
	  make_pwave_source(source_wqs, ksrc, dir1);
	  ENDTIME("make P-wave source");
	  
	  /* Read and/or generate quark prop in Pauli basis */
	  total_iters += 
	    get_wprop_to_field(startflag_w[ksrc][dir1],startfile_w[ksrc][dir1], 
			       saveflag_w[ksrc][dir1], savefile_w[ksrc][dir1],
			       quark_prop, &source_wqs[ksrc], &qic, &dcp);
	  
	  /* Fourier transform */
	  STARTTIME;
	  restrict_fourier_field((complex *)quark_prop,
				 sizeof(block_pauli_propagator), FORWARDS);
	  ENDTIME("Fourier transform quark");
	      
	  STARTTIME;
	  for(dir2=0;dir2<3;dir2++){
	    smear_quark(ksnk,dir1,dir2);
	  
	    /**************************************************************/
	    /*Calculate and print out the pw spectrum */
	    
	    all_pw_prop(dir1, dir2);
	  }
	  ENDTIME("compute raw meson props");

	}//dir1 loop
	    
	STARTTIME;
	assemble_pw_prop();
	ENDTIME("assemble meson props");

	print_pw_prop(ksrc, ksnk);
	
	av_pw_prop_out();
	
      }/*loop over ksrc*/
    }/*loop ksnk */
    
    free_smearing(num_smear);
    free(quark_prop);
    free(quark_prop_smear);
    free(antiquark_prop);

    endtime = dclock();
    node0_printf("RUNNING COMPLETED\n");
    node0_printf("Time = %e seconds\n",(double)(endtime-starttime));
    node0_printf("total_iters= %d\n", total_iters);
    node0_printf("\n\n\n");
    

  } /* while(readin) */

  return 0;
}

