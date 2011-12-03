/************** extract_gauge_v5.c ***********************/
/* MIMD version 7 */

/* Extract the binary payload from a version 5 gauge field file */
/* No checksums are verified */

#include <stdio.h>

int main(int argc, char *argv[]){

  unsigned int magicno;
  int nx, ny, nz, nt;
  char time_stamp[64];
  int order;
  unsigned int cksuma, cksumb;
  unsigned long long site,vol;
  float su3mat[72];

  fread(&magicno, 4, 1, stdin);
  fprintf(stderr, "magicno = %x\n",magicno);

  fread(&nx, 4, 1, stdin);
  fprintf(stderr, "nx      = %d\n",nx);

  fread(&ny, 4, 1, stdin);
  fprintf(stderr, "ny      = %d\n",ny);

  fread(&nz, 4, 1, stdin);
  fprintf(stderr, "nz      = %d\n",nz);

  fread(&nt, 4, 1, stdin);
  fprintf(stderr, "nt      = %d\n",nt);

  vol = nx*ny*nz*nt;
  fprintf(stderr, "vol     = %lld\n",vol);

  fread(time_stamp, 64, 1, stdin);
  time_stamp[64] = '\0';
  fprintf(stderr, "time    = %s\n",time_stamp);

  fread(&order, 4, 1, stdin);
  fprintf(stderr, "order   = %d\n",order);

  fread(&cksuma, 4, 1, stdin);
  fprintf(stderr, "cksuma  = %x\n",cksuma);

  fread(&cksumb, 4, 1, stdin);
  fprintf(stderr, "cksumb  = %x\n",cksumb);

  for(site = 0; site < 4*vol; site++){
    if(fread(su3mat,72,1,stdin) != 1){
      fprintf(stderr,"Unexpected EOF\n");
      break;
    }
    fwrite(su3mat,72,1,stdout);
  }

  return 0;
}
