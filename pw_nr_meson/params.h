#ifndef _PARAMS_H
#define _PARAMS_H

/* structure for passing simulation parameters to each node */
typedef struct {
  int stopflag;	/* 1 if it is time to stop */
  /* INITIALIZATION PARAMETERS */
  int nx,ny,nz,nt;  /* lattice dimensions */
  int nflavors;	/* the number of flavors */
  
  /*  REPEATING BLOCK */
  int startflag;  /* what to do for beginning lattice */
  int saveflag;   /* what to do with lattice at end */
  char startfile[MAXFILENAME],savefile[MAXFILENAME];
  char stringLFN[MAXFILENAME];  /** ILDG LFN if applicable ***/
  int fixflag;   /* Gauge fixing */

  quark_invert_control qic;   /* Inversion parameters */
  dirac_clover_param dcp;     /* Clover parameters */
  
  char ensemble_id[MAXFILENAME];
  int sequence_number;
  
  int a_startflag_w;
  char a_startfile_w[MAXFILENAME];
  int a_saveflag_w;
  char a_savefile_w[MAXFILENAME];
  wilson_quark_source a_wqs;
  
  int startflag_w[NSM][MAXDIR];
  char startfile_w[NSM][MAXDIR][MAXFILENAME];
  
  wilson_quark_source source_wqs[NSM];  /* source parameters */
  int saveflag_w[NSM][MAXDIR];
  char savefile_w[NSM][MAXDIR][MAXFILENAME];
  wilson_quark_source sink_wqs[NSM];  /* sink parameters */
  char source_wf_label[NSM][MAXFILENAME];
  char sink_wf_label[NSM][MAXFILENAME];
  int  num_smear; /* number of smearings */
  int sequence;
  char a0_file[MAXFILENAME];
  char b1_file[MAXFILENAME];
  char a1_file[MAXFILENAME];
  char a2_t2_file[MAXFILENAME];
  char a2_e_file[MAXFILENAME];
  char a2_file[MAXFILENAME];
}  params;

#endif /* _PARAMS_H */
