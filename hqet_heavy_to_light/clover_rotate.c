/************* clover_rotate.c **************************/
/* MIMD version 4  ***/

/* Modifications
   C. DeTar 5/24/97 Original version
 */

/* Right multiplies spin Wilson vector by Dslash to prepare for 
   building SW rotation correction for current */

#include "hqet_light_includes.h"

/* NOTE! Requires two scratch wilson_vector members of the site structure */

#define TMP1 fft_one 
#define TMP2 fft_two

void clover_rotate(field_offset src, field_offset dest)
{
  /* src and dest must be offsets for spin_wilson_vectors in the site structure */
  register int i;
  register site *s; 
  int si;
  field_offset src_wv,dest_wv;  /* Wilson vector offsets  */

  for(si = 0; si< 4; si++)
    {
      /* Compute field offset of wilson_vector si in the source spin Wilson vector */
      /* (We didn't want to assume that the compiler uses 
	 F_OFFSET(src.d[si]) = F_OFFSET(src) + sizeof(wilson_vector)*si),
	 since it may be adding extra bytes between structure members 
	 for alignment or other reasons -- so we do it this way...) */

      src_wv  = (char *)(&((spin_wilson_vector *)F_PT(&lattice[0],src ))->d[si])
	-((char *)&(lattice[0]));
      dest_wv = (char *)(&((spin_wilson_vector *)F_PT(&lattice[0],dest))->d[si])
	-((char *)&(lattice[0]));
      
      /* Do Wilson "dslash" on the src_wv field with (1 +/- gamma) factor */
      dslash_w_site(src_wv,F_OFFSET(TMP1), PLUS,EVENANDODD);
      dslash_w_site(src_wv,F_OFFSET(TMP2),MINUS,EVENANDODD);

      /* Then subtraction leaves 2*Dslash (i.e. 2*gamma factor) */
      FORALLSITES(i,s){
	sub_wilson_vector(&s->TMP1,&s->TMP2,
			  (wilson_vector *)F_PT(s, dest_wv));
      }
    }

}
