#ifndef _PARAMS_H
#define _PARAMS_H

#include "../include/macros.h"  /* For MAXFILENAME */
#include "../include/generic_quark_types.h" /* For quark_source, etc. */


/* structure for passing simulation parameters to each node */
typedef struct {
  int stopflag;   /* 1 if it is time to stop */
  /* INITIALIZATION PARAMETERS */
  int nx,ny,nz,nt;  /* lattice dimensions */
  /*  REPEATING BLOCK */
  int verbose_flag  ; 
  Real clov_c,u0;
  int startflag;  /* what to do for beginning lattice */
  char startfile[MAXFILENAME] ; /*** file containing the gauge configurations ***/
  
  int no_spectator ;  /** number of spectator kappa values ***/
  int niter_spectator ;   /* number of cg iterations for spectator
			     quark */
  int nrestart_spectator ; /* number of restarts for spectator
			      inversion */
  Real resid_spectator;   /* resid error for cg inversion */
  
  Real kappa_spectator[MAX_KAPPA] ; /** Kappa values for the
  				       spectator quarks inversion **/
  quark_source wqs_spectator[MAX_KAPPA]; /* Spectator source parameters */
  int startflag_spectator[MAX_KAPPA]   ;  /** type of IO read for the
                                            spectator quark **/
  char qfile_spectator[MAX_KAPPA][MAXFILENAME];  

  int no_zonked_light ; /*** the number of light kappa values of
                          zonked quarks ***/
  int niter_zonked_light; 	/* maximum number of MR or cg iterations */
  int nrestart_zonked_light;	/* maximum number of MR or cg restarts */
  Real resid_zonked_light;   /* resid error for cg inversion */
  Real kappa_zonked_light[MAX_KAPPA] ; /** Kappa values for the light
                                          zonked quark **/
  quark_source wqs_zonked_light[MAX_KAPPA]; /* Zonked_Light source parameters */
  int startflag_zonked_light[MAX_KAPPA]   ;  /** type of IO read for the zonked light quark **/
  int saveflag_zonked_heavy[MAX_KAPPA]   ;  /** type of IO write for the zonked heavy quark **/

  int saveflag_zonked_light_ssink;  /* For zonked light smeared props */
  int saveflag_zonked_heavy_ssink;  /* For zonked heavy smeared props */
  char qfile_zonked_light[MAX_KAPPA][MAXFILENAME];
  char qfile_suffix_zonked_light[MAXFILENAME];
  char qfile_zonked_heavy[MAX_KAPPA][MAXFILENAME];
  char qfile_suffix_zonked_heavy[MAXFILENAME];
  
  int no_zonked_heavy ; /*** the number of light kappa values of
  			  zonked quarks ***/
  int niter_zonked_heavy; 	/* maximum number of MR or cg iterations */
  int nrestart_zonked_heavy;	/* maximum number of MR or cg restarts */
  Real resid_zonked_heavy;   /* resid error for cg inversion */
  Real kappa_zonked_heavy[MAX_KAPPA] ; /** Kappa values for the heavy
                                          zonked quark **/
  int inverter_type_zonked_heavy[MAX_KAPPA];
  int inverter_type_sequential[MAX_KAPPA];
  quark_source wqs_zonked_heavy[MAX_KAPPA]; /* Zonked_heavy source parameters */
  
  int  no_sequential  ;
  Real kappa_sequential[MAX_KAPPA] ; /** Kappa values for the
                                        sequential inversion **/
  int tf ; /** Fixed timeslice in the sequentiaol source ***/
  char filename_HL3[MAXFILENAME] ;   
  char filename_HH3[MAXFILENAME] ;   
  
  char seq_smear_file[MAXPMOM][MAXFILENAME];
  
  int no_p_values ;
  int no_q_values ;
  int no_k_values ;
  
  int p_momstore[MAXPMOM][3] ;
  int q_momstore[MAXMOM][3] ;
  int k_momstore[MAXMOM][3] ;
  
  char filename_HH2_GL[MAXFILENAME] ; /** File to write the heavy-heavy two **/
  char filename_LL2_GG[MAXFILENAME] ; /** File to write the light-light two
                               point functions to ***/
  char filename_HL2_GG[MAXFILENAME] ; /** File to write the heavy-light two
                                   point functions to ***/
  char filename_HL2_GE[MAXFILENAME] ; /** File to write the heavy-light two
  				   point functions to ***/

  char filename_HL2_GL[MAXFILENAME] ; /** File to write the heavy-light two
  				   point functions to ***/

  int saveflag_HH3;
  int saveflag_HL3;
  int saveflag_HH2_GL;
  int saveflag_LL2_GG;           /* Type of output file */
  int saveflag_HL2_GG;
  int saveflag_HL2_GE;
  int saveflag_HL2_GL;

}  params;


#endif /* _PARAMS_H */
