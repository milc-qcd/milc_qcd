/*************************** stuff_gauge_v5.c ********************/
/* MIMD version 7 */

/* Reconstruct a version 5 gauge field file from the binary payload */
/* No checksums are computed */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(int argc, char *argv[]){

  unsigned int magicno = 0x4e87;
  int nx, ny, nz, nt;
  time_t time_stamp;
  int order = 0;
  unsigned int cksuma = 0, cksumb = 0;
  unsigned long long site,vol;
  float su3mat[72];

  if(argc < 5){
    fprintf(stderr,"Usage %s nx ny nz nt\n",argv[0]);
    return 1;
  }

  nx = atoi(argv[1]);
  ny = atoi(argv[2]);
  nz = atoi(argv[3]);
  nt = atoi(argv[4]);

  time(&time_stamp);

  fwrite(&magicno, 4, 1, stdout);
  fprintf(stderr, "magicno = %x\n",magicno);

  fwrite(&nx, 4, 1, stdout);
  fprintf(stderr, "nx      = %d\n",nx);

  fwrite(&ny, 4, 1, stdout);
  fprintf(stderr, "ny      = %d\n",ny);

  fwrite(&nz, 4, 1, stdout);
  fprintf(stderr, "nz      = %d\n",nz);

  fwrite(&nt, 4, 1, stdout);
  fprintf(stderr, "nt      = %d\n",nt);

  vol = nx*ny*nz*nt;
  fprintf(stderr, "vol     = %lld\n",vol);

  fwrite(ctime(&time_stamp), 64, 1, stdout);
  fprintf(stderr, "time    = %s",ctime(&time_stamp));

  fwrite(&order, 4, 1, stdout);
  fprintf(stderr, "order   = %d\n",order);

  fwrite(&cksuma, 4, 1, stdout);
  fprintf(stderr, "cksuma  = %x\n",cksuma);

  fwrite(&cksumb, 4, 1, stdout);
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
