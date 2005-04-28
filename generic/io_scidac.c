/********************** io_scidac.c **********************************/
/* MILC/QIO interface.  Produces a ILDG archivable file */
/* MIMD version 6 */
/* CD 11/2004
*/

#include "generic_includes.h"
#include <qio.h>
#include "../include/io_lat.h"
#include "../include/io_scidac.h"
#include <string.h>
#define LATDIM 4

/* Map QIO layout functions to MILC functions */

int qio_node_number(const int x[]){
  return node_number(x[0],x[1],x[2],x[3]);
}

int qio_node_index(const int x[]){
  return node_index(x[0],x[1],x[2],x[3]);
}

void qio_get_coords(int x[], int node, int index){
  if(node != this_node){
    printf("qio_get_coords: bad node number %d != this_node %d\n",
	   node,this_node);
    terminate(1);
  }
  x[0] = lattice[index].x;
  x[1] = lattice[index].y;
  x[2] = lattice[index].z;
  x[3] = lattice[index].t;
}

size_t qio_num_sites(int node){
  return num_sites(node);
}

void build_layout(QIO_Layout *layout){
  static int lattice_size[LATDIM];

  lattice_size[0] = nx;
  lattice_size[1] = ny;
  lattice_size[2] = nz;
  lattice_size[3] = nt;

  layout->node_number     = qio_node_number;
  layout->node_index      = qio_node_index;
  layout->get_coords      = qio_get_coords;
  layout->num_sites       = qio_num_sites;
  layout->latsize         = lattice_size;
  layout->latdim          = LATDIM;
  layout->volume          = volume;
  layout->sites_on_node   = sites_on_node;
  layout->this_node       = this_node;
  layout->number_of_nodes = number_of_nodes;
}

QIO_Writer *open_output(char *filename, int volfmt, QIO_Layout *layout,
			char *xml_write_file){
  QIO_String *xml_file_out;
  QIO_Writer *outfile;

  /* Create the file XML */
  xml_file_out = QIO_string_create();
  QIO_string_set(xml_file_out,xml_write_file);

  /* Open the file for writing */
  outfile = QIO_open_write(xml_file_out, filename, volfmt, layout, 0);
  if(outfile == NULL){
    printf("open_output(%d): QIO_open_write returned NULL\n",this_node);
    return NULL;
  }
  QIO_string_destroy(xml_file_out);
  return outfile;
}

QIO_Reader *open_input(char *filename, QIO_Layout *layout){
  QIO_String *xml_file_in;
  QIO_Reader *infile;
  char myname[] = "open_input";

  /* Create the file XML */
  xml_file_in = QIO_string_create();

  /* Open the file for reading */
  infile = QIO_open_read(xml_file_in, filename, layout, 0);
  if(infile == NULL){
    printf("%s(%d): QIO_open_read returns NULL.\n",myname,this_node);
    return NULL;
  }

  if(this_node==0){
    printf("Restoring binary SciDAC file %s\n",filename);
    printf("File info \n\"%s\"\n",QIO_string_ptr(xml_file_in));
  }

  QIO_string_destroy(xml_file_in);
  return infile;
}

void close_output(QIO_Writer *outfile)
{
  QIO_close_write(outfile);
}

void close_input(QIO_Reader *infile)
{
  QIO_close_read(infile);
}

/* Factory function for moving color matrices from site structure to
   single precision output */
void vget_F3_M_from_site(char *buf, size_t index, int count, void *arg)
{
  int dir;
  int i,j;
  /* Assume output lattice is single precision */
  fsu3_matrix *dest = (fsu3_matrix *)buf;
  /* arg contains pointer to field offset value */
  field_offset src = *((field_offset *)arg);
  site *s = &lattice[index];
  /* Source can be any precision */
  su3_matrix *src_mat = (su3_matrix *)F_PT(s,src);

  /* Copy, changing precision, if necessary */
  for(dir = 0; dir < count; dir++)
    for(i = 0; i < 3; i++)for(j = 0; j < 3; j++){
      dest[dir].e[i][j].real = src_mat[dir].e[i][j].real;
      dest[dir].e[i][j].imag = src_mat[dir].e[i][j].imag;
    }
}

/* Factory function for moving color matrices from temp to
   single precision output */
void vget_F3_M_from_field(char *buf, size_t index, int count, void *arg)
{
  int dir;
  int i,j;
  /* Assume output lattice is single precision */
  fsu3_matrix *dest = (fsu3_matrix *)buf;
  /* arg contains pointer to the temp su3_matrix array */
  su3_matrix *src = (su3_matrix *)arg;
  /* Source can be any precision */
  su3_matrix *src_mat = src + index * count;

  /* Copy, changing precision, if necessary */
  for(dir = 0; dir < count; dir++)
    for(i = 0; i < 3; i++)for(j = 0; j < 3; j++){
      dest[dir].e[i][j].real = src_mat[dir].e[i][j].real;
      dest[dir].e[i][j].imag = src_mat[dir].e[i][j].imag;
    }
}

int write_F3_M_from_site(QIO_Writer *outfile, char *xml_write_lattice,
	       field_offset src, int count){
  QIO_String *xml_record_out;
  int status;
  QIO_RecordInfo *rec_info;
  /* We assume output precision is single */
  char qdptype[] = "QDP_F3_ColorMatrix";
  char prec[] = "F";
  int datum_size = sizeof(fsu3_matrix);
  int word_size = sizeof(float);

  /* Create the record info for the field */
  rec_info = QIO_create_record_info(QIO_FIELD, qdptype, prec, 3,
				    0, datum_size, count);
  /* Create the record XML for the field */
  xml_record_out = QIO_string_create();
  QIO_string_set(xml_record_out,xml_write_lattice);

  /* Write the record for the field */
  status = QIO_write(outfile, rec_info, xml_record_out, vget_F3_M_from_site, 
		     count*datum_size, word_size, (void *)&src);
 if(status != QIO_SUCCESS)return 1;

  QIO_destroy_record_info(rec_info);
  QIO_string_destroy(xml_record_out);

  return 0;
}

int write_F3_M_from_field(QIO_Writer *outfile, char *xml_write_lattice,
	         su3_matrix *src, int count){
  QIO_String *xml_record_out;
  int status;
  QIO_RecordInfo *rec_info;
  /* We assume output precision is single */
  char qdptype[] = "QDP_F3_ColorMatrix";
  char prec[] = "F";
  int datum_size = sizeof(fsu3_matrix);
  int word_size = sizeof(float);

  /* Create the record info for the field */
  rec_info = QIO_create_record_info(QIO_FIELD, qdptype, prec, 3,
				    0, datum_size, count);
  /* Create the record XML for the field */
  xml_record_out = QIO_string_create();
  QIO_string_set(xml_record_out,xml_write_lattice);

  /* Write the record for the field */
  status = QIO_write(outfile, rec_info, xml_record_out, 
		     vget_F3_M_from_field, count*datum_size, word_size, 
		     (void *)src);
  if(status != QIO_SUCCESS)return 1;

  QIO_destroy_record_info(rec_info);
  QIO_string_destroy(xml_record_out);

  return 0;
}

/* Factory function for moving single precision gauge field from input
   to site structure */
void vput_F3_M_to_site(char *buf, size_t index, int count, void *arg)
{
  int dir;
  int i,j;
  /* Assume input lattice is single precision, 3 colors */
  fsu3_matrix *src = (fsu3_matrix *)buf;
  field_offset dest = *((field_offset *)arg);
  site *s = &lattice[index];
  /* Destination can be any precision */
  su3_matrix *dest_mat = (su3_matrix *)F_PT(s,dest);
  
  /* Copy, changing precision, if necessary */
  for (dir=0;dir<count;dir++)
    for(i = 0; i < 3; i++)for(j = 0; j < 3; j++){
      dest_mat[dir].e[i][j].real = src[dir].e[i][j].real;
      dest_mat[dir].e[i][j].imag = src[dir].e[i][j].imag;
    }
}

/* Factory function for moving single precision gauge field from input
   to field */
void vput_F3_M_to_field(char *buf, size_t index, int count, void *arg)
{
  int dir;
  int i,j;
  /* Assume input lattice is single precision, 3 colors */
  fsu3_matrix *src = (fsu3_matrix *)buf;
  su3_matrix *dest = (su3_matrix *)arg;
  /* Destination can be any precision */
  su3_matrix *dest_mat = dest + index * count;
  
  /* Copy, changing precision, if necessary */
  for (dir=0;dir<count;dir++)
    for(i = 0; i < 3; i++)for(j = 0; j < 3; j++){
      dest_mat[dir].e[i][j].real = src[dir].e[i][j].real;
      dest_mat[dir].e[i][j].imag = src[dir].e[i][j].imag;
    }
}

/* Read a set of color matrices */
int read_F3_M_to_site(QIO_Reader *infile, field_offset dest, int count)
{
  QIO_String *xml_record_in;
  QIO_RecordInfo rec_info;
  int status;
  /* We assume input precision is single */
  int datum_size = sizeof(fsu3_matrix);
  int word_size = sizeof(float);
  
  /* Read the field record */
  xml_record_in = QIO_string_create();
  status = QIO_read(infile, &rec_info, xml_record_in, 
		    vput_F3_M_to_site, datum_size*count, word_size, (void *)&dest);
  node0_printf("Record info \n\"%s\"\n",QIO_string_ptr(xml_record_in));
  if(status != QIO_SUCCESS)return 1;

  node0_printf("Checksums %x %x\n",
	       QIO_get_reader_last_checksuma(infile),
	       QIO_get_reader_last_checksumb(infile));

  QIO_string_destroy(xml_record_in);
  return 0;
}

/* Read a set of color matrices */
int read_F3_M_to_field(QIO_Reader *infile, su3_matrix *dest, int count)
{
  QIO_String *xml_record_in;
  QIO_RecordInfo rec_info;
  int status;
  /* We assume input precision is single */
  int datum_size = sizeof(fsu3_matrix);
  int word_size = sizeof(float);
  
  /* Read the field record */
  xml_record_in = QIO_string_create();
  status = QIO_read(infile, &rec_info, xml_record_in, 
		    vput_F3_M_to_field, datum_size*count, word_size, 
		    (void *)dest);
  node0_printf("Record info \n\"%s\"\n",QIO_string_ptr(xml_record_in));
  if(status != QIO_SUCCESS)return 1;

  node0_printf("Checksums %x %x\n",
	       QIO_get_reader_last_checksuma(infile),
	       QIO_get_reader_last_checksumb(infile));

  QIO_string_destroy(xml_record_in);
  return 0;
}

/* Save the single precision lattice in SciDAC format */
/* The QIO file is closed after writing the lattice */
gauge_file *save_scidac(char *filename, int volfmt){
  QIO_Layout layout;
  QIO_Writer *outfile;
  int status;
  field_offset src = F_OFFSET(link[0]);
  gauge_file *gf;
  char *qcdml;

  QIO_verbose(QIO_VERB_OFF);

  /* Build the layout structure */
  build_layout(&layout);

  /* Make a dummy gauge file structure for MILC use */
  gf = setup_output_gauge_file();

  /* Open file for writing */
  outfile = open_output(filename, volfmt, &layout,
			"MILC ILDG archival gauge configuration");
  if(outfile == NULL)terminate(1);

  /* Create the QCDML string for this configuration */
  qcdml = create_QCDML();

  /* Write the lattice field */
  status = write_F3_M_from_site(outfile, qcdml, src, LATDIM);
  if(status)terminate(1);
  
  /* Write information */
  if(volfmt == QIO_SINGLEFILE){
    node0_printf("Saved gauge configuration serially to binary file %s\n",
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

  node0_printf("Time stamp %s\n",gf->header->time_stamp);
  node0_printf("Checksums %x %x\n",
	       QIO_get_writer_last_checksuma(outfile),
	       QIO_get_writer_last_checksumb(outfile));

  /* Close the file */
  QIO_close_write(outfile);

  free_QCDML(qcdml);
  return gf;
}

/* The functions below constitute the API */

int read_lat_dim_scidac(char *filename, int *ndim, int dims[])
{
  QIO_Layout layout;
  int i;
  int *latsize;
  QIO_Reader *infile;

  QIO_verbose(QIO_VERB_REG);

  /* Build the layout structure */
  nx = 0; ny = 0; nz = 0; nt = 0;
  build_layout(&layout);
  /* Forces discovery */
  layout.latdim = 0;

  /* Get lattice dimensions from file */
  infile = open_input(filename, &layout);
  if(!infile)return 1;

  *ndim = QIO_get_reader_latdim(infile);
  latsize = QIO_get_reader_latsize(infile);

  for(i = 0; i < *ndim; i++)
    dims[i] = latsize[i];

  QIO_close_read(infile);

  return 0;
}

gauge_file *save_serial_scidac(char *filename){
  return save_scidac(filename,QIO_SINGLEFILE);
}

gauge_file *save_multifile_scidac(char *filename){
  return save_scidac(filename,QIO_MULTIFILE);
}

gauge_file *save_partition_scidac(char *filename){
  return save_scidac(filename,QIO_PARTFILE);
}

/* The QIO file is closed after reading the lattice */
gauge_file *restore_serial_scidac(char *filename){
  QIO_Layout layout;
  QIO_Reader *infile;
  int status;
  field_offset dest = F_OFFSET(link[0]);
  gauge_file *gf;

  QIO_verbose(QIO_VERB_OFF);

  /* Build the layout structure */
  build_layout(&layout);

  /* Make a dummy gauge file structure for MILC use */
  gf = setup_input_gauge_file(filename);

  /* Set the filename in the gauge_file structure */
  gf->filename = filename;

  /* Open file for reading */
  infile = open_input(filename, &layout);
  if(infile == NULL)terminate(1);

  /* Read the lattice field */
  status = read_F3_M_to_site(infile, dest, LATDIM);
  if(status)terminate(1);

  /* Close the file */
  QIO_close_read(infile);

  return gf;
}

/* Read color matrices in SciDAC format */
void restore_color_matrix_scidac_to_site(char *filename, 
				field_offset dest, int count){
  QIO_Layout layout;
  QIO_Reader *infile;
  int status;

  QIO_verbose(QIO_VERB_OFF);

  /* Build the layout structure */
  build_layout(&layout);

  /* Open file for reading */
  infile = open_input(filename, &layout);
  if(infile == NULL)terminate(1);

  /* Read the lattice field */
  status = read_F3_M_to_site(infile, dest, count);
  if(status)terminate(1);

  close_input(infile);
}

/* Read color matrices in SciDAC format */
void restore_color_matrix_scidac_to_field(char *filename, 
				su3_matrix *dest, int count){
  QIO_Layout layout;
  QIO_Reader *infile;
  int status;

  QIO_verbose(QIO_VERB_OFF);

  /* Build the layout structure */
  build_layout(&layout);

  /* Open file for reading */
  infile = open_input(filename, &layout);
  if(infile == NULL)terminate(1);

  /* Read the lattice field */
  status = read_F3_M_to_field(infile, dest, count);
  if(status)terminate(1);

  close_input(infile);
}

/* Write a set of color matrices in SciDAC format, taking data from the site
   structure */
void save_color_matrix_scidac_from_site(char *filename, char *filexml, 
      char *recxml, int volfmt,  field_offset src, int count)
{
  QIO_Layout layout;
  QIO_Writer *outfile;
  int status;

  QIO_verbose(QIO_VERB_OFF);

  /* Build the layout structure */
  build_layout(&layout);

  /* Open file for writing */
  outfile = open_output(filename, volfmt, &layout, filexml);
  if(outfile == NULL)terminate(1);

  /* Write the lattice field */
  status = write_F3_M_from_site(outfile, recxml, src, count);
  if(status)terminate(1);
  
  /* Write information */
  if(volfmt == QIO_SINGLEFILE){
    node0_printf("Saved KS matrix serially to binary file %s\n",
		 filename);
  }
  else if(volfmt == QIO_MULTIFILE){
    node0_printf("Saved KS matrix as multifile to binary file %s\n",
	   filename);
  }
  else if(volfmt == QIO_PARTFILE){
    node0_printf("Saved KS matrix in partition format to binary file %s\n",
	   filename);
  }

  node0_printf("Checksums %x %x\n",
	       QIO_get_writer_last_checksuma(outfile),
	       QIO_get_writer_last_checksumb(outfile));

  close_output(outfile);
}


/* Save a color matrix. */
void save_color_matrix_scidac_from_field(char *filename, 
        char *filexml, char *recxml, int volfmt, su3_matrix *src, int count)
{
  QIO_Layout layout;
  QIO_Writer *outfile;
  int status;

  QIO_verbose(QIO_VERB_OFF);

  /* Build the layout structure */
  build_layout(&layout);

  /* Open file for writing */
  outfile = open_output(filename, volfmt, &layout, filexml);
  if(outfile == NULL)terminate(1);

  /* Write the lattice field */
  status = write_F3_M_from_field(outfile, recxml, src, count);
  if(status)terminate(1);
  
  /* Write information */
  if(volfmt == QIO_SINGLEFILE){
    node0_printf("Saved KS matrix serially to binary file %s\n",
		 filename);
  }
  else if(volfmt == QIO_MULTIFILE){
    node0_printf("Saved KS matrix as multifile to binary file %s\n",
	   filename);
  }
  else if(volfmt == QIO_PARTFILE){
    node0_printf("Saved KS matrix in partition format to binary file %s\n",
	   filename);
  }

  node0_printf("Checksums %x %x\n",
	       QIO_get_writer_last_checksuma(outfile),
	       QIO_get_writer_last_checksumb(outfile));

  close_output(outfile);
}


