/* make_channel_wfs.c */
/* MIMD version 6 */

/* Collection of routines for defining three baryon-baryon channels
   coupling to the H dibaryon */

/* C. DeTar 5 April 1997 */

/* Defines wavefunction for 6 quark state constructed from the spin singlet
   lambda-lambda channel with all quarks in the same spatial orbital */

#include "h_dibaryon_includes.h"
#include <string.h>

void lam_lam(diqkwf *wf, int *nterm, Real *norm)
{

  *nterm = 90;
  *norm = 12.;

  /* Weight was divided by 3 */

  wf[ 0].dq[0] =  5;  wf[ 0].dq[1] = 14;  wf[ 0].dq[2] =  2;  wf[ 0].wt =  -2;
  wf[ 1].dq[0] =  5;  wf[ 1].dq[1] = 12;  wf[ 1].dq[2] =  9;  wf[ 1].wt =  -1;
  wf[ 2].dq[0] =  5;  wf[ 2].dq[1] = 12;  wf[ 2].dq[2] =  4;  wf[ 2].wt =  -1;
  wf[ 3].dq[0] =  5;  wf[ 3].dq[1] = 13;  wf[ 3].dq[2] =  6;  wf[ 3].wt =   1;
  wf[ 4].dq[0] =  5;  wf[ 4].dq[1] = 13;  wf[ 4].dq[2] =  3;  wf[ 4].wt =   1;
  wf[ 5].dq[0] =  7;  wf[ 5].dq[1] = 11;  wf[ 5].dq[2] =  2;  wf[ 5].wt =  -2;
  wf[ 6].dq[0] =  7;  wf[ 6].dq[1] =  9;  wf[ 6].dq[2] =  9;  wf[ 6].wt =   1;
  wf[ 7].dq[0] =  7;  wf[ 7].dq[1] =  9;  wf[ 7].dq[2] =  4;  wf[ 7].wt =   1;
  wf[ 8].dq[0] =  7;  wf[ 8].dq[1] =  4;  wf[ 8].dq[2] =  9;  wf[ 8].wt =   1;
  wf[ 9].dq[0] =  7;  wf[ 9].dq[1] =  4;  wf[ 9].dq[2] =  4;  wf[ 9].wt =   1;
  wf[10].dq[0] =  7;  wf[10].dq[1] =  2;  wf[10].dq[2] = 11;  wf[10].wt =  -2;
  wf[11].dq[0] =  8;  wf[11].dq[1] =  8;  wf[11].dq[2] =  2;  wf[11].wt =   2;
  wf[12].dq[0] =  8;  wf[12].dq[1] =  6;  wf[12].dq[2] =  9;  wf[12].wt =  -1;
  wf[13].dq[0] =  8;  wf[13].dq[1] =  6;  wf[13].dq[2] =  4;  wf[13].wt =  -1;
  wf[14].dq[0] =  8;  wf[14].dq[1] =  4;  wf[14].dq[2] =  6;  wf[14].wt =  -1;
  wf[15].dq[0] =  8;  wf[15].dq[1] =  4;  wf[15].dq[2] =  3;  wf[15].wt =  -1;
  wf[16].dq[0] =  8;  wf[16].dq[1] =  2;  wf[16].dq[2] =  8;  wf[16].wt =   1;
  wf[17].dq[0] =  8;  wf[17].dq[1] =  2;  wf[17].dq[2] = 10;  wf[17].wt =   1;
  wf[18].dq[0] =  6;  wf[18].dq[1] =  8;  wf[18].dq[2] =  9;  wf[18].wt =  -1;
  wf[19].dq[0] =  6;  wf[19].dq[1] =  8;  wf[19].dq[2] =  4;  wf[19].wt =  -1;
  wf[20].dq[0] =  6;  wf[20].dq[1] =  6;  wf[20].dq[2] = 11;  wf[20].wt =   2;
  wf[21].dq[0] =  6;  wf[21].dq[1] = 11;  wf[21].dq[2] =  6;  wf[21].wt =   1;
  wf[22].dq[0] =  6;  wf[22].dq[1] = 11;  wf[22].dq[2] =  3;  wf[22].wt =   1;
  wf[23].dq[0] =  6;  wf[23].dq[1] =  9;  wf[23].dq[2] =  8;  wf[23].wt =  -1;
  wf[24].dq[0] =  6;  wf[24].dq[1] =  9;  wf[24].dq[2] = 10;  wf[24].wt =  -1;
  wf[25].dq[0] = 10;  wf[25].dq[1] = 10;  wf[25].dq[2] =  2;  wf[25].wt =   2;
  wf[26].dq[0] = 10;  wf[26].dq[1] =  9;  wf[26].dq[2] =  6;  wf[26].wt =  -1;
  wf[27].dq[0] = 10;  wf[27].dq[1] =  9;  wf[27].dq[2] =  3;  wf[27].wt =  -1;
  wf[28].dq[0] = 10;  wf[28].dq[1] =  3;  wf[28].dq[2] =  9;  wf[28].wt =  -1;
  wf[29].dq[0] = 10;  wf[29].dq[1] =  3;  wf[29].dq[2] =  4;  wf[29].wt =  -1;
  wf[30].dq[0] = 10;  wf[30].dq[1] =  2;  wf[30].dq[2] =  8;  wf[30].wt =   1;
  wf[31].dq[0] = 10;  wf[31].dq[1] =  2;  wf[31].dq[2] = 10;  wf[31].wt =   1;
  wf[32].dq[0] = 11;  wf[32].dq[1] =  7;  wf[32].dq[2] =  2;  wf[32].wt =  -2;
  wf[33].dq[0] = 11;  wf[33].dq[1] =  6;  wf[33].dq[2] =  6;  wf[33].wt =   1;
  wf[34].dq[0] = 11;  wf[34].dq[1] =  6;  wf[34].dq[2] =  3;  wf[34].wt =   1;
  wf[35].dq[0] = 11;  wf[35].dq[1] =  3;  wf[35].dq[2] =  6;  wf[35].wt =   1;
  wf[36].dq[0] = 11;  wf[36].dq[1] =  3;  wf[36].dq[2] =  3;  wf[36].wt =   1;
  wf[37].dq[0] = 11;  wf[37].dq[1] =  2;  wf[37].dq[2] =  7;  wf[37].wt =  -2;
  wf[38].dq[0] =  9;  wf[38].dq[1] =  7;  wf[38].dq[2] =  9;  wf[38].wt =   1;
  wf[39].dq[0] =  9;  wf[39].dq[1] =  7;  wf[39].dq[2] =  4;  wf[39].wt =   1;
  wf[40].dq[0] =  9;  wf[40].dq[1] =  6;  wf[40].dq[2] =  8;  wf[40].wt =  -1;
  wf[41].dq[0] =  9;  wf[41].dq[1] =  6;  wf[41].dq[2] = 10;  wf[41].wt =  -1;
  wf[42].dq[0] =  9;  wf[42].dq[1] = 10;  wf[42].dq[2] =  6;  wf[42].wt =  -1;
  wf[43].dq[0] =  9;  wf[43].dq[1] = 10;  wf[43].dq[2] =  3;  wf[43].wt =  -1;
  wf[44].dq[0] =  9;  wf[44].dq[1] =  9;  wf[44].dq[2] =  7;  wf[44].wt =   2;
  wf[45].dq[0] =  0;  wf[45].dq[1] = 14;  wf[45].dq[2] =  9;  wf[45].wt =  -1;
  wf[46].dq[0] =  0;  wf[46].dq[1] = 14;  wf[46].dq[2] =  4;  wf[46].wt =  -1;
  wf[47].dq[0] =  0;  wf[47].dq[1] = 12;  wf[47].dq[2] = 11;  wf[47].wt =  -2;
  wf[48].dq[0] =  0;  wf[48].dq[1] = 13;  wf[48].dq[2] =  8;  wf[48].wt =   1;
  wf[49].dq[0] =  0;  wf[49].dq[1] = 13;  wf[49].dq[2] = 10;  wf[49].wt =   1;
  wf[50].dq[0] =  1;  wf[50].dq[1] = 14;  wf[50].dq[2] =  6;  wf[50].wt =   1;
  wf[51].dq[0] =  1;  wf[51].dq[1] = 14;  wf[51].dq[2] =  3;  wf[51].wt =   1;
  wf[52].dq[0] =  1;  wf[52].dq[1] = 12;  wf[52].dq[2] =  8;  wf[52].wt =   1;
  wf[53].dq[0] =  1;  wf[53].dq[1] = 12;  wf[53].dq[2] = 10;  wf[53].wt =   1;
  wf[54].dq[0] =  1;  wf[54].dq[1] = 13;  wf[54].dq[2] =  7;  wf[54].wt =  -2;
  wf[55].dq[0] =  3;  wf[55].dq[1] = 10;  wf[55].dq[2] =  9;  wf[55].wt =  -1;
  wf[56].dq[0] =  3;  wf[56].dq[1] = 10;  wf[56].dq[2] =  4;  wf[56].wt =  -1;
  wf[57].dq[0] =  3;  wf[57].dq[1] = 11;  wf[57].dq[2] =  6;  wf[57].wt =   1;
  wf[58].dq[0] =  3;  wf[58].dq[1] = 11;  wf[58].dq[2] =  3;  wf[58].wt =   1;
  wf[59].dq[0] =  3;  wf[59].dq[1] =  3;  wf[59].dq[2] = 11;  wf[59].wt =   2;
  wf[60].dq[0] =  3;  wf[60].dq[1] =  4;  wf[60].dq[2] =  8;  wf[60].wt =  -1;
  wf[61].dq[0] =  3;  wf[61].dq[1] =  4;  wf[61].dq[2] = 10;  wf[61].wt =  -1;
  wf[62].dq[0] =  4;  wf[62].dq[1] =  7;  wf[62].dq[2] =  9;  wf[62].wt =   1;
  wf[63].dq[0] =  4;  wf[63].dq[1] =  7;  wf[63].dq[2] =  4;  wf[63].wt =   1;
  wf[64].dq[0] =  4;  wf[64].dq[1] =  8;  wf[64].dq[2] =  6;  wf[64].wt =  -1;
  wf[65].dq[0] =  4;  wf[65].dq[1] =  8;  wf[65].dq[2] =  3;  wf[65].wt =  -1;
  wf[66].dq[0] =  4;  wf[66].dq[1] =  3;  wf[66].dq[2] =  8;  wf[66].wt =  -1;
  wf[67].dq[0] =  4;  wf[67].dq[1] =  3;  wf[67].dq[2] = 10;  wf[67].wt =  -1;
  wf[68].dq[0] =  4;  wf[68].dq[1] =  4;  wf[68].dq[2] =  7;  wf[68].wt =   2;
  wf[69].dq[0] =  2;  wf[69].dq[1] =  7;  wf[69].dq[2] = 11;  wf[69].wt =  -2;
  wf[70].dq[0] =  2;  wf[70].dq[1] =  8;  wf[70].dq[2] =  8;  wf[70].wt =   1;
  wf[71].dq[0] =  2;  wf[71].dq[1] =  8;  wf[71].dq[2] = 10;  wf[71].wt =   1;
  wf[72].dq[0] =  2;  wf[72].dq[1] = 10;  wf[72].dq[2] =  8;  wf[72].wt =   1;
  wf[73].dq[0] =  2;  wf[73].dq[1] = 10;  wf[73].dq[2] = 10;  wf[73].wt =   1;
  wf[74].dq[0] =  2;  wf[74].dq[1] = 11;  wf[74].dq[2] =  7;  wf[74].wt =  -2;
  wf[75].dq[0] = 14;  wf[75].dq[1] =  5;  wf[75].dq[2] =  2;  wf[75].wt =  -2;
  wf[76].dq[0] = 14;  wf[76].dq[1] =  0;  wf[76].dq[2] =  9;  wf[76].wt =  -1;
  wf[77].dq[0] = 14;  wf[77].dq[1] =  0;  wf[77].dq[2] =  4;  wf[77].wt =  -1;
  wf[78].dq[0] = 14;  wf[78].dq[1] =  1;  wf[78].dq[2] =  6;  wf[78].wt =   1;
  wf[79].dq[0] = 14;  wf[79].dq[1] =  1;  wf[79].dq[2] =  3;  wf[79].wt =   1;
  wf[80].dq[0] = 12;  wf[80].dq[1] =  5;  wf[80].dq[2] =  9;  wf[80].wt =  -1;
  wf[81].dq[0] = 12;  wf[81].dq[1] =  5;  wf[81].dq[2] =  4;  wf[81].wt =  -1;
  wf[82].dq[0] = 12;  wf[82].dq[1] =  0;  wf[82].dq[2] = 11;  wf[82].wt =  -2;
  wf[83].dq[0] = 12;  wf[83].dq[1] =  1;  wf[83].dq[2] =  8;  wf[83].wt =   1;
  wf[84].dq[0] = 12;  wf[84].dq[1] =  1;  wf[84].dq[2] = 10;  wf[84].wt =   1;
  wf[85].dq[0] = 13;  wf[85].dq[1] =  5;  wf[85].dq[2] =  6;  wf[85].wt =   1;
  wf[86].dq[0] = 13;  wf[86].dq[1] =  5;  wf[86].dq[2] =  3;  wf[86].wt =   1;
  wf[87].dq[0] = 13;  wf[87].dq[1] =  0;  wf[87].dq[2] =  8;  wf[87].wt =   1;
  wf[88].dq[0] = 13;  wf[88].dq[1] =  0;  wf[88].dq[2] = 10;  wf[88].wt =   1;
  wf[89].dq[0] = 13;  wf[89].dq[1] =  1;  wf[89].dq[2] =  7;  wf[89].wt =  -2;
}

/*********************************************************************/
/* Defines wavefunction for 6 quark state constructed from the 
   spin singlet, isospin singlet nucleon-xi channel with all 
   quarks in the same spatial orbital */

void nuc_xi(diqkwf *wf, int *nterm, Real *norm)
{

  *nterm = 132;
  *norm = 12*sqrt(3.);

  /* Weight was divided by 4 */

  wf[ 0].dq[0] =  5;  wf[ 0].dq[1] =  6;  wf[ 0].dq[2] = 13;  wf[ 0].wt =   3;
  wf[ 1].dq[0] =  5;  wf[ 1].dq[1] =  9;  wf[ 1].dq[2] = 12;  wf[ 1].wt =  -3;
  wf[ 2].dq[0] =  5;  wf[ 2].dq[1] =  3;  wf[ 2].dq[2] = 13;  wf[ 2].wt =  -1;
  wf[ 3].dq[0] =  5;  wf[ 3].dq[1] =  4;  wf[ 3].dq[2] = 12;  wf[ 3].wt =   1;
  wf[ 4].dq[0] =  5;  wf[ 4].dq[1] =  2;  wf[ 4].dq[2] = 14;  wf[ 4].wt =  -2;
  wf[ 5].dq[0] =  5;  wf[ 5].dq[1] = 12;  wf[ 5].dq[2] =  9;  wf[ 5].wt =  -2;
  wf[ 6].dq[0] =  5;  wf[ 6].dq[1] = 12;  wf[ 6].dq[2] =  4;  wf[ 6].wt =   2;
  wf[ 7].dq[0] =  5;  wf[ 7].dq[1] = 13;  wf[ 7].dq[2] =  6;  wf[ 7].wt =   2;
  wf[ 8].dq[0] =  5;  wf[ 8].dq[1] = 13;  wf[ 8].dq[2] =  3;  wf[ 8].wt =  -2;
  wf[ 9].dq[0] =  7;  wf[ 9].dq[1] =  9;  wf[ 9].dq[2] =  9;  wf[ 9].wt =   1;
  wf[10].dq[0] =  7;  wf[10].dq[1] =  9;  wf[10].dq[2] =  4;  wf[10].wt =  -1;
  wf[11].dq[0] =  7;  wf[11].dq[1] =  1;  wf[11].dq[2] = 13;  wf[11].wt =  -2;
  wf[12].dq[0] =  7;  wf[12].dq[1] =  4;  wf[12].dq[2] =  9;  wf[12].wt =  -1;
  wf[13].dq[0] =  7;  wf[13].dq[1] =  4;  wf[13].dq[2] =  4;  wf[13].wt =   1;
  wf[14].dq[0] =  7;  wf[14].dq[1] = 13;  wf[14].dq[2] =  1;  wf[14].wt =  -2;
  wf[15].dq[0] =  8;  wf[15].dq[1] =  9;  wf[15].dq[2] =  6;  wf[15].wt =  -1;
  wf[16].dq[0] =  8;  wf[16].dq[1] =  9;  wf[16].dq[2] =  3;  wf[16].wt =   1;
  wf[17].dq[0] =  8;  wf[17].dq[1] =  0;  wf[17].dq[2] = 13;  wf[17].wt =   3;
  wf[18].dq[0] =  8;  wf[18].dq[1] =  1;  wf[18].dq[2] = 12;  wf[18].wt =  -1;
  wf[19].dq[0] =  8;  wf[19].dq[1] =  3;  wf[19].dq[2] =  9;  wf[19].wt =   1;
  wf[20].dq[0] =  8;  wf[20].dq[1] =  3;  wf[20].dq[2] =  4;  wf[20].wt =  -1;
  wf[21].dq[0] =  8;  wf[21].dq[1] =  2;  wf[21].dq[2] =  8;  wf[21].wt =   1;
  wf[22].dq[0] =  8;  wf[22].dq[1] =  2;  wf[22].dq[2] = 10;  wf[22].wt =  -1;
  wf[23].dq[0] =  8;  wf[23].dq[1] = 12;  wf[23].dq[2] =  1;  wf[23].wt =  -1;
  wf[24].dq[0] =  8;  wf[24].dq[1] = 13;  wf[24].dq[2] =  0;  wf[24].wt =   3;
  wf[25].dq[0] =  6;  wf[25].dq[1] =  5;  wf[25].dq[2] = 13;  wf[25].wt =   3;
  wf[26].dq[0] =  6;  wf[26].dq[1] = 10;  wf[26].dq[2] =  9;  wf[26].wt =  -1;
  wf[27].dq[0] =  6;  wf[27].dq[1] = 10;  wf[27].dq[2] =  4;  wf[27].wt =   1;
  wf[28].dq[0] =  6;  wf[28].dq[1] = 11;  wf[28].dq[2] =  6;  wf[28].wt =   1;
  wf[29].dq[0] =  6;  wf[29].dq[1] = 11;  wf[29].dq[2] =  3;  wf[29].wt =  -1;
  wf[30].dq[0] =  6;  wf[30].dq[1] =  1;  wf[30].dq[2] = 14;  wf[30].wt =  -1;
  wf[31].dq[0] =  6;  wf[31].dq[1] =  4;  wf[31].dq[2] =  8;  wf[31].wt =  -1;
  wf[32].dq[0] =  6;  wf[32].dq[1] =  4;  wf[32].dq[2] = 10;  wf[32].wt =   1;
  wf[33].dq[0] =  6;  wf[33].dq[1] = 14;  wf[33].dq[2] =  1;  wf[33].wt =  -1;
  wf[34].dq[0] =  6;  wf[34].dq[1] = 13;  wf[34].dq[2] =  5;  wf[34].wt =   3;
  wf[35].dq[0] = 10;  wf[35].dq[1] =  6;  wf[35].dq[2] =  9;  wf[35].wt =  -1;
  wf[36].dq[0] = 10;  wf[36].dq[1] =  6;  wf[36].dq[2] =  4;  wf[36].wt =   1;
  wf[37].dq[0] = 10;  wf[37].dq[1] =  0;  wf[37].dq[2] = 13;  wf[37].wt =  -1;
  wf[38].dq[0] = 10;  wf[38].dq[1] =  1;  wf[38].dq[2] = 12;  wf[38].wt =   3;
  wf[39].dq[0] = 10;  wf[39].dq[1] =  4;  wf[39].dq[2] =  6;  wf[39].wt =   1;
  wf[40].dq[0] = 10;  wf[40].dq[1] =  4;  wf[40].dq[2] =  3;  wf[40].wt =  -1;
  wf[41].dq[0] = 10;  wf[41].dq[1] =  2;  wf[41].dq[2] =  8;  wf[41].wt =  -1;
  wf[42].dq[0] = 10;  wf[42].dq[1] =  2;  wf[42].dq[2] = 10;  wf[42].wt =   1;
  wf[43].dq[0] = 10;  wf[43].dq[1] = 12;  wf[43].dq[2] =  1;  wf[43].wt =   3;
  wf[44].dq[0] = 10;  wf[44].dq[1] = 13;  wf[44].dq[2] =  0;  wf[44].wt =  -1;
  wf[45].dq[0] = 11;  wf[45].dq[1] =  6;  wf[45].dq[2] =  6;  wf[45].wt =   1;
  wf[46].dq[0] = 11;  wf[46].dq[1] =  6;  wf[46].dq[2] =  3;  wf[46].wt =  -1;
  wf[47].dq[0] = 11;  wf[47].dq[1] =  0;  wf[47].dq[2] = 12;  wf[47].wt =  -2;
  wf[48].dq[0] = 11;  wf[48].dq[1] =  3;  wf[48].dq[2] =  6;  wf[48].wt =  -1;
  wf[49].dq[0] = 11;  wf[49].dq[1] =  3;  wf[49].dq[2] =  3;  wf[49].wt =   1;
  wf[50].dq[0] = 11;  wf[50].dq[1] = 12;  wf[50].dq[2] =  0;  wf[50].wt =  -2;
  wf[51].dq[0] =  9;  wf[51].dq[1] =  5;  wf[51].dq[2] = 12;  wf[51].wt =  -3;
  wf[52].dq[0] =  9;  wf[52].dq[1] =  7;  wf[52].dq[2] =  9;  wf[52].wt =   1;
  wf[53].dq[0] =  9;  wf[53].dq[1] =  7;  wf[53].dq[2] =  4;  wf[53].wt =  -1;
  wf[54].dq[0] =  9;  wf[54].dq[1] =  8;  wf[54].dq[2] =  6;  wf[54].wt =  -1;
  wf[55].dq[0] =  9;  wf[55].dq[1] =  8;  wf[55].dq[2] =  3;  wf[55].wt =   1;
  wf[56].dq[0] =  9;  wf[56].dq[1] =  0;  wf[56].dq[2] = 14;  wf[56].wt =   1;
  wf[57].dq[0] =  9;  wf[57].dq[1] =  3;  wf[57].dq[2] =  8;  wf[57].wt =   1;
  wf[58].dq[0] =  9;  wf[58].dq[1] =  3;  wf[58].dq[2] = 10;  wf[58].wt =  -1;
  wf[59].dq[0] =  9;  wf[59].dq[1] = 14;  wf[59].dq[2] =  0;  wf[59].wt =   1;
  wf[60].dq[0] =  9;  wf[60].dq[1] = 12;  wf[60].dq[2] =  5;  wf[60].wt =  -3;
  wf[61].dq[0] =  0;  wf[61].dq[1] =  8;  wf[61].dq[2] = 13;  wf[61].wt =   3;
  wf[62].dq[0] =  0;  wf[62].dq[1] = 10;  wf[62].dq[2] = 13;  wf[62].wt =  -1;
  wf[63].dq[0] =  0;  wf[63].dq[1] = 11;  wf[63].dq[2] = 12;  wf[63].wt =  -2;
  wf[64].dq[0] =  0;  wf[64].dq[1] =  9;  wf[64].dq[2] = 14;  wf[64].wt =   1;
  wf[65].dq[0] =  0;  wf[65].dq[1] =  4;  wf[65].dq[2] = 14;  wf[65].wt =  -3;
  wf[66].dq[0] =  0;  wf[66].dq[1] = 14;  wf[66].dq[2] =  9;  wf[66].wt =   2;
  wf[67].dq[0] =  0;  wf[67].dq[1] = 14;  wf[67].dq[2] =  4;  wf[67].wt =  -2;
  wf[68].dq[0] =  0;  wf[68].dq[1] = 13;  wf[68].dq[2] =  8;  wf[68].wt =   2;
  wf[69].dq[0] =  0;  wf[69].dq[1] = 13;  wf[69].dq[2] = 10;  wf[69].wt =  -2;
  wf[70].dq[0] =  1;  wf[70].dq[1] =  7;  wf[70].dq[2] = 13;  wf[70].wt =  -2;
  wf[71].dq[0] =  1;  wf[71].dq[1] =  8;  wf[71].dq[2] = 12;  wf[71].wt =  -1;
  wf[72].dq[0] =  1;  wf[72].dq[1] =  6;  wf[72].dq[2] = 14;  wf[72].wt =  -1;
  wf[73].dq[0] =  1;  wf[73].dq[1] = 10;  wf[73].dq[2] = 12;  wf[73].wt =   3;
  wf[74].dq[0] =  1;  wf[74].dq[1] =  3;  wf[74].dq[2] = 14;  wf[74].wt =   3;
  wf[75].dq[0] =  1;  wf[75].dq[1] = 14;  wf[75].dq[2] =  6;  wf[75].wt =  -2;
  wf[76].dq[0] =  1;  wf[76].dq[1] = 14;  wf[76].dq[2] =  3;  wf[76].wt =   2;
  wf[77].dq[0] =  1;  wf[77].dq[1] = 12;  wf[77].dq[2] =  8;  wf[77].wt =  -2;
  wf[78].dq[0] =  1;  wf[78].dq[1] = 12;  wf[78].dq[2] = 10;  wf[78].wt =   2;
  wf[79].dq[0] =  3;  wf[79].dq[1] =  5;  wf[79].dq[2] = 13;  wf[79].wt =  -1;
  wf[80].dq[0] =  3;  wf[80].dq[1] =  8;  wf[80].dq[2] =  9;  wf[80].wt =   1;
  wf[81].dq[0] =  3;  wf[81].dq[1] =  8;  wf[81].dq[2] =  4;  wf[81].wt =  -1;
  wf[82].dq[0] =  3;  wf[82].dq[1] = 11;  wf[82].dq[2] =  6;  wf[82].wt =  -1;
  wf[83].dq[0] =  3;  wf[83].dq[1] = 11;  wf[83].dq[2] =  3;  wf[83].wt =   1;
  wf[84].dq[0] =  3;  wf[84].dq[1] =  9;  wf[84].dq[2] =  8;  wf[84].wt =   1;
  wf[85].dq[0] =  3;  wf[85].dq[1] =  9;  wf[85].dq[2] = 10;  wf[85].wt =  -1;
  wf[86].dq[0] =  3;  wf[86].dq[1] =  1;  wf[86].dq[2] = 14;  wf[86].wt =   3;
  wf[87].dq[0] =  3;  wf[87].dq[1] = 14;  wf[87].dq[2] =  1;  wf[87].wt =   3;
  wf[88].dq[0] =  3;  wf[88].dq[1] = 13;  wf[88].dq[2] =  5;  wf[88].wt =  -1;
  wf[89].dq[0] =  4;  wf[89].dq[1] =  5;  wf[89].dq[2] = 12;  wf[89].wt =   1;
  wf[90].dq[0] =  4;  wf[90].dq[1] =  7;  wf[90].dq[2] =  9;  wf[90].wt =  -1;
  wf[91].dq[0] =  4;  wf[91].dq[1] =  7;  wf[91].dq[2] =  4;  wf[91].wt =   1;
  wf[92].dq[0] =  4;  wf[92].dq[1] =  6;  wf[92].dq[2] =  8;  wf[92].wt =  -1;
  wf[93].dq[0] =  4;  wf[93].dq[1] =  6;  wf[93].dq[2] = 10;  wf[93].wt =   1;
  wf[94].dq[0] =  4;  wf[94].dq[1] = 10;  wf[94].dq[2] =  6;  wf[94].wt =   1;
  wf[95].dq[0] =  4;  wf[95].dq[1] = 10;  wf[95].dq[2] =  3;  wf[95].wt =  -1;
  wf[96].dq[0] =  4;  wf[96].dq[1] =  0;  wf[96].dq[2] = 14;  wf[96].wt =  -3;
  wf[97].dq[0] =  4;  wf[97].dq[1] = 14;  wf[97].dq[2] =  0;  wf[97].wt =  -3;
  wf[98].dq[0] =  4;  wf[98].dq[1] = 12;  wf[98].dq[2] =  5;  wf[98].wt =   1;
  wf[99].dq[0] =  2;  wf[99].dq[1] =  5;  wf[99].dq[2] = 14;  wf[99].wt =  -2;
  wf[100].dq[0] =  2;  wf[100].dq[1] =  8;  wf[100].dq[2] =  8;  wf[100].wt =   1;
  wf[101].dq[0] =  2;  wf[101].dq[1] =  8;  wf[101].dq[2] = 10;  wf[101].wt =  -1;
  wf[102].dq[0] =  2;  wf[102].dq[1] = 10;  wf[102].dq[2] =  8;  wf[102].wt =  -1;
  wf[103].dq[0] =  2;  wf[103].dq[1] = 10;  wf[103].dq[2] = 10;  wf[103].wt =   1;
  wf[104].dq[0] =  2;  wf[104].dq[1] = 14;  wf[104].dq[2] =  5;  wf[104].wt =  -2;
  wf[105].dq[0] = 14;  wf[105].dq[1] =  6;  wf[105].dq[2] =  1;  wf[105].wt =  -1;
  wf[106].dq[0] = 14;  wf[106].dq[1] =  9;  wf[106].dq[2] =  0;  wf[106].wt =   1;
  wf[107].dq[0] = 14;  wf[107].dq[1] =  0;  wf[107].dq[2] =  9;  wf[107].wt =   2;
  wf[108].dq[0] = 14;  wf[108].dq[1] =  0;  wf[108].dq[2] =  4;  wf[108].wt =  -2;
  wf[109].dq[0] = 14;  wf[109].dq[1] =  1;  wf[109].dq[2] =  6;  wf[109].wt =  -2;
  wf[110].dq[0] = 14;  wf[110].dq[1] =  1;  wf[110].dq[2] =  3;  wf[110].wt =   2;
  wf[111].dq[0] = 14;  wf[111].dq[1] =  3;  wf[111].dq[2] =  1;  wf[111].wt =   3;
  wf[112].dq[0] = 14;  wf[112].dq[1] =  4;  wf[112].dq[2] =  0;  wf[112].wt =  -3;
  wf[113].dq[0] = 14;  wf[113].dq[1] =  2;  wf[113].dq[2] =  5;  wf[113].wt =  -2;
  wf[114].dq[0] = 12;  wf[114].dq[1] =  5;  wf[114].dq[2] =  9;  wf[114].wt =  -2;
  wf[115].dq[0] = 12;  wf[115].dq[1] =  5;  wf[115].dq[2] =  4;  wf[115].wt =   2;
  wf[116].dq[0] = 12;  wf[116].dq[1] =  8;  wf[116].dq[2] =  1;  wf[116].wt =  -1;
  wf[117].dq[0] = 12;  wf[117].dq[1] = 10;  wf[117].dq[2] =  1;  wf[117].wt =   3;
  wf[118].dq[0] = 12;  wf[118].dq[1] = 11;  wf[118].dq[2] =  0;  wf[118].wt =  -2;
  wf[119].dq[0] = 12;  wf[119].dq[1] =  9;  wf[119].dq[2] =  5;  wf[119].wt =  -3;
  wf[120].dq[0] = 12;  wf[120].dq[1] =  1;  wf[120].dq[2] =  8;  wf[120].wt =  -2;
  wf[121].dq[0] = 12;  wf[121].dq[1] =  1;  wf[121].dq[2] = 10;  wf[121].wt =   2;
  wf[122].dq[0] = 12;  wf[122].dq[1] =  4;  wf[122].dq[2] =  5;  wf[122].wt =   1;
  wf[123].dq[0] = 13;  wf[123].dq[1] =  5;  wf[123].dq[2] =  6;  wf[123].wt =   2;
  wf[124].dq[0] = 13;  wf[124].dq[1] =  5;  wf[124].dq[2] =  3;  wf[124].wt =  -2;
  wf[125].dq[0] = 13;  wf[125].dq[1] =  7;  wf[125].dq[2] =  1;  wf[125].wt =  -2;
  wf[126].dq[0] = 13;  wf[126].dq[1] =  8;  wf[126].dq[2] =  0;  wf[126].wt =   3;
  wf[127].dq[0] = 13;  wf[127].dq[1] =  6;  wf[127].dq[2] =  5;  wf[127].wt =   3;
  wf[128].dq[0] = 13;  wf[128].dq[1] = 10;  wf[128].dq[2] =  0;  wf[128].wt =  -1;
  wf[129].dq[0] = 13;  wf[129].dq[1] =  0;  wf[129].dq[2] =  8;  wf[129].wt =   2;
  wf[130].dq[0] = 13;  wf[130].dq[1] =  0;  wf[130].dq[2] = 10;  wf[130].wt =  -2;
  wf[131].dq[0] = 13;  wf[131].dq[1] =  3;  wf[131].dq[2] =  5;  wf[131].wt =  -1;
}

void make_channel_wfs(diqkwf channel_wf[MAXCHANNEL][MAXWF],
		      int channel_terms[MAXCHANNEL],
		      char channel_label[MAXCHANNEL][MAXLABEL], 
		      Real channel_norm[MAXCHANNEL],
		      int *nchannel)
{
  lam_lam(channel_wf[0],&channel_terms[0],&channel_norm[0]);
  strcpy(channel_label[0],"LL");
  nuc_xi (channel_wf[1],&channel_terms[1],&channel_norm[1]);
  strcpy(channel_label[1],"NX");
  *nchannel = 2;
} /* make_channel_wfs */
