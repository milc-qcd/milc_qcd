#ifndef _PARAMS_H
#define _PARAMS_H

#include "../include/macros.h"  /* For MAXFILENAME */
#include "../include/generic_quark_types.h"

/* structure for passing simulation parameters to each node */
typedef struct {
  int stopflag;   /* 1 if it is time to stop */
  /* INITIALIZATION PARAMETERS */
  int nx,ny,nz,nt;  /* lattice dimensions */
  int tf ; /** Fixed timeslice in the sequentiaol source ***/
  /*  REPEATING BLOCK */
  Real clov_c,u0;
  int verbose_flag  ; 
  int startflag ;      /* beginning lattice: CONTINUE, RELOAD, FRESH */
  char startfile[MAXFILENAME] ; /*** file containing the gauge configurations ***/
  
  int no_spectator ;  /** number of spectator kappa values ***/
  int niter_spectator ;   /* number of cg iterations for spectator
			     quark */
  int nrestart_spectator ; /* number of restarts for spectator
			      inversion */
  Real resid_spectator;   /* resid error for cg inversion */
  quark_source wqs_spectator[MAX_KAPPA]; /* Spectator source parameters */
  quark_source wqs_zonked_light[MAX_KAPPA]; /* Zonked_Light source parameters */
  
  Real kappa_spectator[MAX_KAPPA] ; /** Kappa values for the spectator
					 quarks inversion **/
  int startflag_spectator[MAX_KAPPA] ;  /** type of IO read for the
					    spectator quark ***/
  char qfile_spectator[MAX_KAPPA][MAXFILENAME];  

  int no_zonked_light ; /*** the number of light kappa values of
			  zonked quarks ***/
  int niter_zonked ;   /* number of cg iterations for zonked
			     quark */
  int nrestart_zonked ; /* number of restarts for zonked
			      inversion */
  Real resid_zonked;   /* resid error for cg inversion */

  Real kappa_zonked_light[MAX_KAPPA] ; /** Kappa values for the light
					    zonked quark **/
  int startflag_zonked[MAX_KAPPA] ; /** type of IO read for the zonked quark **/
  
  char qfile_zonked[MAX_KAPPA][MAXFILENAME];
  
  char hqet_smear_file[MAXVEL ][MAXFILENAME]; /** File containing the sources
					 for the sequential inversion ****/
  char heavy_light_out[MAXFILENAME] ; /** File to write the heavy --> light
				 form factors to ***/
  char twopt_out[MAXFILENAME] ;
  char seq_out[MAXFILENAME] ;       /** File to write the sequential two point
			       functions ***/
  
  int no_q_values ;
  int quark_type  ;
  
  int q_momstore[MAXMOM][3] ;
  
  Real velocity[MAXVEL][4] ;  /*** The store of velocities ******/
  int novel ;   /** The number of velocities   ****/
  
}  params;



#endif /* _PARAMS_H */
