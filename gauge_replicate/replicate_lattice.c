/************************  replicate_lattice.c  -- ******************/
/* MIMD version 7 */

/* Enlarge the gauge field copying it in the specified directions */

#include "gauge_replicate_includes.h"
#include <qmp.h>
#include <qio.h>
#include "../include/io_lat.h"
#include "../include/io_scidac.h"
#include <string.h>
#define LATDIM 4

/* layout for replicated lattice */

static int squaresize2[4];
static size_t sites_on_node2;
static size_t volume2;

/* layout for replicated lattice */
/* Return the node number (MPI rank) that has the data from the
   original lattice that will be written for the replica-lattice site
   coordinates coords2 */

int node_number2(const int coords2[]) {

  int x = coords2[XUP] % nx ;
  int y = coords2[YUP] % ny ;
  int z = coords2[ZUP] % nz ;
  int t = coords2[TUP] % nt ;

  return node_number(x, y, z, t);
}

/*------------------------------------------------------------------*/
/* Convert PE rank to coordinates */
static void lex_coords(int coords[], const int dim, const int size[], 
	   const size_t rank)
{
  int d;
  size_t r = rank;

  for(d = 0; d < dim; d++){
    coords[d] = r % size[d];
    r /= size[d];
  }
}

/*------------------------------------------------------------------*/
/* Parity of the coordinate */
static int coord_parity(const int r[]){
  return (r[0] + r[1] + r[2] + r[3]) % 2;
}

/*------------------------------------------------------------------*/
/* Convert coordinate to linear lexicographic PE rank (inverse of
   lex_coords) */

static size_t lex_rank(const int coords[], int dim, int size[])
{
  int d;
  size_t rank = coords[dim-1];

  for(d = dim-2; d >= 0; d--){
    rank = rank * size[d] + coords[d];
  }
  return rank;
}

/*------------------------------------------------------------------*/
/* Return the index of the site on the replicated lattice that has
   data belonging to the replicated lattice at coords2 */

static int 
node_index2(const int coords2[]) {
  register int i,xr,yr,zr,tr;
  
  xr = coords2[XUP] % squaresize2[XUP]; 
  yr = coords2[YUP] % squaresize2[YUP];
  zr = coords2[ZUP] % squaresize2[ZUP]; 
  tr = coords2[TUP] % squaresize2[TUP];
  
  i = xr + squaresize2[XUP]*( yr + squaresize2[YUP]*( zr + squaresize2[ZUP]*tr));
  if(coord_parity(coords2) == 0 ){	/* even site */
    return( i/2 );
  }
  else {
    return( (i + sites_on_node2)/2 );
  }
}

/*------------------------------------------------------------------*/
/* Map PE-rank number and serialize site index of the replicated
   lattice to coordinates of the original lattice  */
/* Assumes even sites come first */
static void 
get_coords2(int coords2[], int node, size_t index2){
  int mc[4];
  int ir;
  int meo, neven2, xeo;
  int k = node;

  /* mc = the machine coordinates for node k */
#ifdef HAVE_QMP
  QMP_comm_get_logical_coordinates_from2(QMP_comm_get_default(), mc, k);
#else
  lex_coords(mc, 4, nsquares, k);
#endif

  /* meo = the parity of the machine coordinate */
  meo = coord_parity(mc);

  /* neven2 = the number of even sites on node k */
  neven2 = (sites_on_node2 + 1 - meo)/2;
  
  /* ir = the even part of the lexicographic index within the
     sublattice on node k */
  if(index2 < neven2){
    ir = 2*index2;
    xeo = 0;
  } else {
    ir = 2*(index2 - neven2);
    xeo = 1;
  }

  /* coords2 = the sublattice coordinate */
  lex_coords(coords2, 4, squaresize2, ir);

  /* Add offset to get full lattice coordinate (still a 2-fold ambiguity) */
  coords2[XUP] += mc[XUP]*squaresize2[XUP];
  coords2[YUP] += mc[YUP]*squaresize2[YUP];
  coords2[ZUP] += mc[ZUP]*squaresize2[ZUP];
  coords2[TUP] += mc[TUP]*squaresize2[TUP];

  /* Adjust coordinate according to parity */
  if( coord_parity(coords2) != xeo ){
    coords2[XUP]++;
    if(coords2[XUP] >= squaresize2[XUP]*(mc[XUP]+1)){
      coords2[XUP] -= squaresize2[XUP];
      coords2[YUP]++;
      if(coords2[YUP] >= squaresize2[YUP]*(mc[YUP]+1)){
	coords2[YUP] -= squaresize2[YUP];
	coords2[ZUP]++;
	if(coords2[ZUP] >= squaresize2[ZUP]*(mc[ZUP]+1)){
	  coords2[ZUP] -= squaresize2[ZUP];
	  coords2[TUP]++;
	}
      }
    }
  }
}

/*------------------------------------------------------------------*/
int num_sites2(int node) {
    return sites_on_node2;
}

/*------------------------------------------------------------------*/
static void 
p2f_mat(fsu3_matrix *dest, su3_matrix *src){
  int i,j;
  
  for(i = 0; i < 3; i++)for(j = 0; j < 3; j++){
    dest->e[i][j].real = src->e[i][j].real;
    dest->e[i][j].imag = src->e[i][j].imag;
  }
}

/*------------------------------------------------------------------*/
/* Factory function for replicating data from a MILC field to
   the output buffer.  The index presented is for the replicated
   field.  We convert it to the index of the original field */

void vget2_F3_M_from_site(char *buf, size_t index2, int count,
			  void *arg)
{
  int i;
  int coords[4];
  fsu3_matrix *dest = (fsu3_matrix *)buf;
  /* arg contains pointer to the data array */
  field_offset src_off = *((field_offset *)arg); \

  /* Translate index of replicated lattice to index of original lattice */
  get_coords2(coords, this_node, index2);

  /* Convert to coordinates of the original lattice */
  coords[XUP] = coords[XUP] % nx;
  coords[YUP] = coords[YUP] % ny;
  coords[ZUP] = coords[ZUP] % nz;
  coords[TUP] = coords[TUP] % nt;

  size_t index = node_index(coords[XUP], coords[YUP], 
			    coords[ZUP], coords[TUP]);

  site *s = &lattice[index];
  su3_matrix *src = (su3_matrix *)F_PT(s,src_off);
  for(i = 0; i < count; i++)
    p2f_mat(dest+i, src+i);
}

int
setup_layout2(int reps[]){
  int rep = reps[0]*reps[1]*reps[2]*reps[3];
  const int *nsquares = get_logical_dimensions();

  volume2 = volume*rep;
  sites_on_node2 = sites_on_node*rep;

  squaresize2[0] = reps[0]*nx/nsquares[0];
  squaresize2[1] = reps[1]*ny/nsquares[1];
  squaresize2[2] = reps[2]*nz/nsquares[2];
  squaresize2[3] = reps[3]*nt/nsquares[3];
}

static void 
build_qio_layout2(QIO_Layout *layout, int reps[]){
  static int lattice_size2[LATDIM];

  lattice_size2[0] = nx*reps[0];
  lattice_size2[1] = ny*reps[1];
  lattice_size2[2] = nz*reps[2];
  lattice_size2[3] = nt*reps[3];

  layout->node_number     = node_number2;
  layout->node_index      = node_index2;
  layout->get_coords      = get_coords2;
  layout->num_sites       = num_sites2;
  layout->latsize         = lattice_size2;
  layout->latdim          = LATDIM;
  layout->volume          = volume2;
  layout->sites_on_node   = sites_on_node2;
  layout->this_node       = this_node;
  layout->number_of_nodes = number_of_nodes;
}
 
int write_F3_M_from_site2(QIO_Writer *outfile,
         QIO_String *xml_record_out, field_offset src, int count){
  int status;
  QIO_RecordInfo *rec_info;
  /* We assume output precision is single */
  char qdptype[] = "M";
  char prec[] = "F";
  int datum_size = sizeof(fsu3_matrix);
  int word_size = sizeof(float);

  /* Create the record info for the field */
  rec_info = QIO_create_record_info(QIO_FIELD, 0, 0, 0, qdptype, prec, 3,
  				    0, datum_size, count);

  /* Write the record for the field */
  status = QIO_write(outfile, rec_info, xml_record_out,
		     vget2_F3_M_from_site,
		     count*datum_size, word_size, (void *)&src);
 if(status != QIO_SUCCESS)return 1;

  QIO_destroy_record_info(rec_info);

  return 0;
}

void save_replicated_lattice(int reps[], int saveflag, 
			     char filename[], char stringLFN[],
			     int volfmt, int serpar, int ildgstyle,
			     char ildgLFN[])
{
  int ndim[4] = { nx, ny, nz, nt };
  int i;
  site *s;

  QIO_Layout layout;
  QIO_Filesystem fs;
  QIO_Writer *outfile;
  int status;
  field_offset src = F_OFFSET(link[0]);
  gauge_file *gf;
  char *info;
  QIO_String *filexml;
  QIO_String *recxml;
  char default_file_xml[] = "<?xml version=\"1.0\" encoding=\"UTF-8\"?><title>MILC ILDG archival gauge configuration</title>";

  QIO_verbose(QIO_VERB_DEBUG);

  /* Set up the MILC layout for the replicated lattice */
  setup_layout2(reps);

  /* Build the layout structure */
  build_qio_layout2(&layout, reps);

  /* Define the I/O system  - Use the one in io_scidac.c */
  build_qio_filesystem(&fs);

  /* Make a dummy gauge file structure for MILC use */
  gf = setup_output_gauge_file();

  /* Set the filename in the gauge_file structure */
  gf->filename = filename;

  /* Open file for writing */
  filexml = QIO_string_create();
  QIO_string_set(filexml, default_file_xml);
  outfile = open_scidac_output(filename, volfmt, serpar, ildgstyle,
			       stringLFN, &layout, &fs, filexml);
  if(outfile == NULL)terminate(1);
  QIO_string_destroy(filexml);

  /* Create the QCDML string for this configuration */
  info = create_QCDML();
  recxml = QIO_string_create();
  QIO_string_set(recxml, info);

  /* Write the lattice field */
  status = write_F3_M_from_site2(outfile, recxml, src, LATDIM);
  if(status)terminate(1);

  /* Discard for now */
  QIO_string_destroy(recxml);
  
  /* Write information */
  if(volfmt == QIO_SINGLEFILE){
    node0_printf("Saved gauge configuration as a single binary file %s\n",
		 filename);
  }
  else if(volfmt == QIO_MULTIFILE){
    node0_printf("Saved gauge configuration as multifile to binary file %s\n",
	   filename);
  }
  else if(volfmt == QIO_PARTFILE){
    node0_printf("Saved gauge configuration in partition format to binary file %s\n",
	   filename);
  }
  else if(volfmt == QIO_PARTFILE_DIR){
    node0_printf("Saved gauge configuration in partition format to binary file %s\n",
	   filename);
  }

  node0_printf("Time stamp %s\n",gf->header->time_stamp);
  node0_printf("Checksums %x %x\n",
	       QIO_get_writer_last_checksuma(outfile),
	       QIO_get_writer_last_checksumb(outfile));

  /* Close the file */
  QIO_close_write(outfile);

  free_QCDML(info);
}

