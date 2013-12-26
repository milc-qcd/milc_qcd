/*************************** lattice_to_scidac.c ************************/
/* MIMD version 7 */
/* Make modulation file for modulating a quark source */
/* MPP or single-processor */
/* C. DeTar 10/23/15 */

/* Usage ...
   
   Modify the functions f and f_type below to suit
   make make_modulation

   make_modulation nx ny nz nt t0 modulation_filename

*/

#define CONTROL

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <math.h>
#include "../include/complex.h"
#include "../include/su3.h"
#include <lattice.h>
#include "../include/macros.h"
#include "../include/comdefs.h"
#include "../include/generic.h"
#include "../include/io_scidac.h"

#include <qio.h>

/*----------------------------------------------------------------------*/

complex f(int x, int y, int z, int t, int t0){
  if(t == t0 || t0 == -1){
    Real chi = 1. + cos(2.*PI*x/nx) + cos(2.*PI*y/ny) + cos(2.*PI*z/nz);
    return cmplx(chi,0.);
  } else {
    return cmplx(0.,0.);
  }
}

char *f_type(void){
  return "000 + 100 + 010 + 001 + -100 + -010 + -001";
}

/*----------------------------------------------------------------------*/
void make_lattice(){
  register int i;               /* scratch */
  int x,y,z,t;            /* coordinates */
  /* allocate space for lattice, fill in parity, coordinates and index.  */
  lattice = (site *)malloc( sites_on_node * sizeof(site) );
  if(lattice==NULL){
    printf("NODE %d: no room for lattice\n",this_node);
    terminate(1);
  }
  for(t=0;t<nt;t++)for(z=0;z<nz;z++)for(y=0;y<ny;y++)for(x=0;x<nx;x++){
	  if(node_number(x,y,z,t)==mynode()){
            i=node_index(x,y,z,t);
            lattice[i].x=x;     lattice[i].y=y; lattice[i].z=z; lattice[i].t=t;
            lattice[i].index = x+nx*(y+ny*(z+nz*t));
            if( (x+y+z+t)%2 == 0)lattice[i].parity=EVEN;
            else                 lattice[i].parity=ODD;
	  }
	}
}
/*----------------------------------------------------------------------*/

void setup() {
  setup_layout();
  make_lattice();
}

/*----------------------------------------------------------------------*/

int main(int argc, char *argv[])
{
  
  char *mod_file;
  char *mod_type;
  complex *chi;
  char fileinfo[] = "Modulation file";
  char *recinfo;
  int i, t0;
  site *s;
  
  if(argc < 6)
    {
      fprintf(stderr,"Usage %s nx ny nz nt t0 modulation_filename\n",argv[0]);
      exit(1);
    }
  
  /* Unpack args */
  nx = atoi(argv[1]);
  ny = atoi(argv[2]);
  nz = atoi(argv[3]);
  nt = atoi(argv[4]);
  volume = nx*ny*nz*nt;
  t0 = atoi(argv[5]);
  mod_file = argv[6];
  
  mod_type = f_type();
  
  if(this_node == 0)printf("Generating SciDAC-formatted modulation file %s of type %s\n",
			   mod_file, mod_type);
  
  initialize_machine(&argc,&argv);
  
  this_node = mynode();
  number_of_nodes = numnodes();
  
  setup();

  /* Create the modulation field */
  chi = create_c_field();
  FORALLSITES(i,s){
    chi[i] = f(s->x, s->y, s->z, s->t, t0);
  }

  /* Write file in SciDAC format */
  recinfo = mod_type;
  save_complex_scidac_from_field(mod_file, fileinfo, recinfo, 
				 QIO_SINGLEFILE, QIO_SERIAL, chi, 1);

  normal_exit(0);

  return 0;
}
