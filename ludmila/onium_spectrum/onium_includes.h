#include "onium_generic.h"
#include "../include/io_wb.h"
#include "../include/io_prop_ks.h"

complex  two_pt_trace(wilson_propagator *antiquark, wilson_propagator *quark, 
		      int *g_snk, int n_snk, int *g_src, int n_src,
		      int *p,site *s);
  
void two_pt_func(field_offset snk, field_offset src, int *g_snk, int n_snk,
		 int *g_src, int n_src, int *p, complex *prop);

void All_KS_hl_prop(field_offset snk, field_offset src, complex **propagator);

int calculate_stag_prop();

void rotate_w_quark(field_offset, field_offset, float);

void get_smearings(char *filename);

void r_prop_w_fm(char *filename,field_offset);


void r_prop_w_fm1(char *filename);
