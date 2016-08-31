#ifndef _OPENMP_DEFS_H
#define _OPENMP_DEFS_H
#include "../include/loopend.h"
// Definitions for putting OpenMP into the MILC code
//
// Note including loopend.h before this is required
//
// To use this
// #include "../include/openmp_defs.h"
//
//  Then replace loops over sites (FORALLSITES, FORSOMEPARITY ... )
//  by FORALLSITES_OMP etc, and end each loop with END_LOOP_OMP
//
//  These ..._OMP macros take one more argument, which is the OMP clauses
//  that specify how variables are treated: private, firstprivate, lastprivate,
//  reduction, shared, ...
//  ONE ANNOYING FEATURE at the moment is that this can't be empty.  So if you
//  don't need to specify any such clauses, simply say "default(shared)",
//  which is the default anyway.
//
//	FORSOMEPARITY_OMP(i,s,l_parity,private(j) reduction(+:source_norm) ){
//	    source_norm += (double) magsq_su3vec( src+i );
//	    su3vec_copy( src+i, resid+i);
//	    su3vec_copy(resid+i, cg_p+i);
//	    clearvec(psim[j_low]+i);
//	    for(j=0;j<num_offsets;j++) if(j!=j_low){
//		clearvec(psim[j]+i);
//		su3vec_copy(resid+i, pm[j]+i);
//	    }
//       } END_LOOP_OMP
// where only the first and last lines have been modified from the current code
//
// Then, if "-DOMP" is defined at compile time, you get OMP, otherwise the
// extra argument is ignored and you get the code without OMP pragmas
//
// For the moment (01/30/13) do not use -DC_GLOBAL_INLINE in your compilation.
// It sometimes works, but the hidden inlining means that private variables
// won't really be properly identified.
//
#ifdef OMP
#define STR(s) #s  // copied from p. 222 of "C in a Nutshell"

// redefine FORALLSITES, FORSOMEPARITY, but must include specifying private and
// reduction variables

#define FORALLSITES_OMP(i,s,omp_args) \
  { \
_Pragma( STR(omp parallel for private(i,s) omp_args) )\
  for( i=0; i<sites_on_node; i++){ s= &(lattice[i]);
  
  //note extra bracket so all these loops can use the same END_LOOP_OMP

#define FORSOMESUBLATTICE_OMP(i,s,subl,omp_args) \
  { \
_Pragma( STR(omp parallel for private(i,s,subl,last_in_loop) omp_args) ) \
    for( i=(subl*subl_sites_on_node), s= &(lattice[i]), \
	 last_in_loop=((subl+1)*subl_sites_on_node); \
	 i< last_in_loop; i++,s++){
  //note extra bracket so all these loops can use the same END_LOOP_OMP

#define FORALLFIELDSITES_OMP(i,omp_args) \
  { \
_Pragma( STR(omp parallel for private(i) omp_args) )\
  for( i=0; i<sites_on_node; i++){
  //note extra bracket so all these loops can use the same END_LOOP_OMP

#define FOREVENFIELDSITES_OMP(i,omp_args) \
  { \
_Pragma( STR(omp parallel for private(i) omp_args) )\
  for( i=0; i<even_sites_on_node; i++){
  //note extra bracket so all these loops can use the same END_LOOP_OMP

#define FORODDFIELDSITES_OMP(i,omp_args) \
  { \
_Pragma( STR(omp parallel for private(i) omp_args) )\
  for( i=even_sites_on_node; i<sites_on_node; i++){
  //note extra bracket so all these loops can use the same END_LOOP_OMP

#define FORSOMEPARITY_OMP(i,s,parity,omp_args) \
  {register int loopend,loopstart;\
  loopend= (parity)==EVEN ? even_sites_on_node : sites_on_node ;\
  loopstart=((parity)==ODD ? even_sites_on_node : 0 );\
_Pragma( STR(omp parallel for private(i,s) omp_args) )\
  for( i=loopstart; i<loopend; i++){ s= &(lattice[i]);

#define FORSOMEFIELDPARITY_OMP(i,parity,omp_args) \
  {register int loopend,loopstart;\
  loopend= (parity)==EVEN ? even_sites_on_node : sites_on_node ;\
  loopstart=((parity)==ODD ? even_sites_on_node : 0 );\
_Pragma( STR(omp parallel for private(i) omp_args) )\
  for( i=loopstart; i<loopend; i++){


#ifdef SCHROED_FUN
#define FOREVENSITESDOMAIN_OMP(i,s,omp_args)	\
  FOREVENSITES_OMP(i,s,omp_args) if(s->t > 0)
#define FORODDSITESDOMAIN_OMP(i,s,omp_args)		\
  FORODDSITES_OMP(i,s,omp_args) if(s->t > 0)
#define FORALLSITESDOMAIN_OMP(i,s,omp_args)		\
  FORALLSITES_OMP(i,s,omp_args) if(s->t > 0)
#define FORSOMEPARITYDOMAIN_OMP(i,s,parity,omp_args)	\
  FORSOMEPARITY_OMP(i,s,parity,omp_args) if(s->t > 0)
#else
#define FOREVENSITESDOMAIN_OMP FOREVENSITES_OMP
#define FORODDSITESDOMAIN_OMP FORODDSITES_OMP
#define FORALLSITESDOMAIN_OMP FORALLSITES_OMP
#define FORSOMEPARITYDOMAIN_OMP FORSOMEPARITY_OMP
#endif

#define END_LOOP_OMP }} // need to close the bracket at beginning of definition above
		    // and bracket  in "{ s= ..."

#else // Not OMP


#define FORSOMESUBLATTICE_OMP(i,s,subl,omp_args) FORSOMESUBLATTICE(i,s,subl){
#define FORALLSITES_OMP(i,s,omp_args) FORALLSITES(i,s){
#define FOREVENSITES_OMP(i,s,omp_args) FOREVENSITES(i,s){
#define FORODDSITES_OMP(i,s,omp_args) FORODDSITES(i,s){
#define FORALLFIELDSITES_OMP(i,omp_args) FORALLFIELDSITES(i){
#define FOREVENFIELDSITES_OMP(i,omp_args) FOREVENFIELDSITES(i){
#define FORODDFIELDSITES_OMP(i,omp_args) FORODDFIELDSITES(i){
#define FORSOMEPARITY_OMP(i,s,parity,omp_args) FORSOMEPARITY(i,s,parity)
#define FORSOMEFIELDPARITY_OMP(i,parity,omp_args) FORSOMEFIELDPARITY(i,parity)
#define FOREVENSITESDOMAIN_OMP(i,s,omp_args) FOREVENSITESDOMAIN(i,s)
#define FORODDSITESDOMAIN_OMP(i,s,omp_args) FORODDSITESDOMAIN(i,s)
#define FORALLSITESDOMAIN_OMP(i,s,omp_args) FORALLSITESDOMAIN(i,s)
#define FORSOMEPARITYDOMAIN_OMP(i,s,parity,omp_args) FORSOMEPARITYDOMAIN(i,s,parity)
#define END_LOOP_OMP END_LOOP


#endif //OMP
#endif // _OPENMP_DEFS_H
