/*********************** io_scidac_types.c **************************/
/* Functions for reading and writing MILC fields through QIO        */
/* C. DeTar 6/7/07
*/

/* Conversion of MILC types between prevailing and fixed precision */

#include "generic_includes.h"
#include <qio.h>
#include "../include/io_lat.h"
#include "../include/io_scidac.h"

/* Real type */

static void
f2p_real(Real *dest, float *src){
  *dest = *src;
}

static void
p2f_real(float *dest, Real *src){
  *dest = *src;
}

#if 0
static void
d2p_real(Real *dest, double *src){
  *dest = *src;
}

static void
p2d_real(double *dest, Real *src){
  *dest = *src;
}
#endif

/* Complex type */

static void
f2p_complex(complex *dest, fcomplex *src){
  dest->real = src->real;
  dest->imag = src->imag;
}

static void
p2f_complex(fcomplex *dest, complex *src){
  dest->real = src->real;
  dest->imag = src->imag;
}

static void
d2p_complex(complex *dest, dcomplex *src){
  dest->real = src->real;
  dest->imag = src->imag;
}

static void
p2d_complex(dcomplex *dest, complex *src){
  dest->real = src->real;
  dest->imag = src->imag;
}

/* Color vector */

static void 
f2p_vec(su3_vector *dest, fsu3_vector *src){
  int i;
  
  for(i = 0; i < 3; i++){
    dest->c[i].real = src->c[i].real;
    dest->c[i].imag = src->c[i].imag;
  }
}

static void 
p2f_vec(fsu3_vector *dest, su3_vector *src){
  int i;
  
  for(i = 0; i < 3; i++){
    dest->c[i].real = src->c[i].real;
    dest->c[i].imag = src->c[i].imag;
  }
}

static void 
d2p_vec(su3_vector *dest, dsu3_vector *src){
  int i;
  
  for(i = 0; i < 3; i++){
    dest->c[i].real = src->c[i].real;
    dest->c[i].imag = src->c[i].imag;
  }
}

static void 
p2d_vec(dsu3_vector *dest, su3_vector *src){
  int i;
  
  for(i = 0; i < 3; i++){
    dest->c[i].real = src->c[i].real;
    dest->c[i].imag = src->c[i].imag;
  }
}

/* Wilson vector */

static void 
f2p_wvec(wilson_vector *dest, fwilson_vector *src){
  int i,j;
  
  for(j = 0; j < 4; j++)for(i = 0; i < 3; i++){
    dest->d[j].c[i].real = src->d[j].c[i].real;
    dest->d[j].c[i].imag = src->d[j].c[i].imag;
  }
}

static void 
p2f_wvec(fwilson_vector *dest, wilson_vector *src){
  int i,j;
  
  for(j = 0; j < 4; j++)for(i = 0; i < 3; i++){
    dest->d[j].c[i].real = src->d[j].c[i].real;
    dest->d[j].c[i].imag = src->d[j].c[i].imag;
  }
}

static void 
d2p_wvec(wilson_vector *dest, dwilson_vector *src){
  int i,j;
  
  for(j = 0; j < 4; j++)for(i = 0; i < 3; i++){
    dest->d[j].c[i].real = src->d[j].c[i].real;
    dest->d[j].c[i].imag = src->d[j].c[i].imag;
  }
}

static void 
p2d_wvec(dwilson_vector *dest, wilson_vector *src){
  int i,j;
  
  for(j = 0; j < 4; j++)for(i = 0; i < 3; i++){
    dest->d[j].c[i].real = src->d[j].c[i].real;
    dest->d[j].c[i].imag = src->d[j].c[i].imag;
  }
}

/* Color matrix */

static void 
f2p_mat(su3_matrix *dest, fsu3_matrix *src){
  int i,j;
  
  for(i = 0; i < 3; i++)for(j = 0; j < 3; j++){
    dest->e[i][j].real = src->e[i][j].real;
    dest->e[i][j].imag = src->e[i][j].imag;
  }
}

static void 
p2f_mat(fsu3_matrix *dest, su3_matrix *src){
  int i,j;
  
  for(i = 0; i < 3; i++)for(j = 0; j < 3; j++){
    dest->e[i][j].real = src->e[i][j].real;
    dest->e[i][j].imag = src->e[i][j].imag;
  }
}

static void 
d2p_mat(su3_matrix *dest, dsu3_matrix *src){
  int i,j;
  
  for(i = 0; i < 3; i++)for(j = 0; j < 3; j++){
    dest->e[i][j].real = src->e[i][j].real;
    dest->e[i][j].imag = src->e[i][j].imag;
  }
}

/* Compilation suppressed until we need it */
#if 0
static void 
p2d_mat(dsu3_matrix *dest, su3_matrix *src){
  int i,j;
  
  for(i = 0; i < 3; i++)for(j = 0; j < 3; j++){
    dest->e[i][j].real = src->e[i][j].real;
    dest->e[i][j].imag = src->e[i][j].imag;
  }
}
#endif

/* Factory function for moving data from MILC site structure to output
   buffer */

#define make_vget_from_site(P, C, T, FIXTYPE, VARTYPE, FUNC) \
void vget_##P##C##_##T##_from_site(char *buf, size_t index, int count, \
				   void *arg) \
{ \
  int i; \
  FIXTYPE *dest = (FIXTYPE *)buf; \
  /* arg contains pointer to field offset value */ \
  field_offset src_off = *((field_offset *)arg); \
  site *s = &lattice[index]; \
  VARTYPE *src = (VARTYPE *)F_PT(s,src_off); \
  for(i = 0; i < count; i++) \
    FUNC(dest+i, src+i); \
}

/* Factory function for moving data from a MILC field to the output
   buffer */

#define make_vget_from_field(P, C, T, FIXTYPE, VARTYPE, FUNC) \
void vget_##P##C##_##T##_from_field(char *buf, size_t index, int count, \
				    void *arg) \
{ \
  int i; \
  FIXTYPE *dest = (FIXTYPE *)buf; \
  /* arg contains pointer to the data array */ \
  VARTYPE *src_pt = (VARTYPE *)arg; \
  VARTYPE *src = src_pt + index * count; \
  for(i = 0; i < count; i++) \
    FUNC(dest+i, src+i); \
}

#define make_vget(P, C,  T, FIXTYPE, VARTYPE, FUNC) \
 make_vget_from_site(P, C,  T, FIXTYPE, VARTYPE, FUNC); \
 make_vget_from_field(P, C,  T, FIXTYPE, VARTYPE, FUNC);

/* Factory function for moving data from input buffer to MILC site 
   structure */

#define make_vput_to_site(P, C, T, FIXTYPE, VARTYPE, FUNC) \
void vput_##P##C##_##T##_to_site(char *buf, size_t index, int count, \
				 void *arg)			     \
{ \
  int i; \
  FIXTYPE *src = (FIXTYPE *)buf; \
  /* arg contains pointer to field offset value */ \
  field_offset dest_off = *((field_offset *)arg); \
  site *s = &lattice[index]; \
  VARTYPE *dest = (VARTYPE *)F_PT(s,dest_off); \
  for(i = 0; i < count; i++) \
    FUNC(dest+i, src+i); \
}

/* Factory function for moving data from input buffer to MILC field */

#define make_vput_to_field(P, C, T, FIXTYPE, VARTYPE, FUNC) \
void vput_##P##C##_##T##_to_field(char *buf, size_t index, int count, \
				 void *arg) \
{ \
  int i; \
  FIXTYPE *src = (FIXTYPE *)buf; \
  /* arg contains pointer to the data array */ \
  VARTYPE *dest_pt = (VARTYPE *)arg; \
  VARTYPE *dest = dest_pt + index * count; \
  for(i = 0; i < count; i++) \
    FUNC(dest+i, src+i); \
}

#define make_vput(P, C,  T, FIXTYPE, VARTYPE, FUNC) \
 make_vput_to_site(P, C,  T, FIXTYPE, VARTYPE, FUNC); \
 make_vput_to_field(P, C,  T, FIXTYPE, VARTYPE, FUNC);

/* Single precision */

make_vget(F,  , R, float,          Real,          p2f_real);
make_vget(F,  , C, fcomplex,       complex,       p2f_complex);
make_vget(F, 3, V, fsu3_vector,    su3_vector,    p2f_vec);
make_vget(F, 3, D, fwilson_vector, wilson_vector, p2f_wvec);
make_vget(F, 3, M, fsu3_matrix,    su3_matrix,    p2f_mat);

make_vput(F,  , R, float,          Real,          f2p_real);
make_vput(F,  , C, fcomplex,       complex,       f2p_complex);
make_vput(F, 3, V, fsu3_vector,    su3_vector,    f2p_vec);
make_vput(F, 3, D, fwilson_vector, wilson_vector, f2p_wvec);
make_vput(F, 3, M, fsu3_matrix,    su3_matrix,    f2p_mat);

/* Double precision */

make_vget(D,  , C, dcomplex,       complex,       p2d_complex);
make_vget(D, 3, V, dsu3_vector,    su3_vector,    p2d_vec);
make_vget(D, 3, D, dwilson_vector, wilson_vector, p2d_wvec);

make_vput(D,  , C, dcomplex,       complex,       d2p_complex);
make_vput(D, 3, V, dsu3_vector,    su3_vector,    d2p_vec);
make_vput(D, 3, D, dwilson_vector, wilson_vector, d2p_wvec);
make_vput(D, 3, M, dsu3_matrix,    su3_matrix,    d2p_mat);

/* Write MILC site structure data */

#define make_write_all_from_site(P, PSTRING, C, CVAL, SVAL, T, TYPESTRING, \
			     FIXTYPE, VARTYPE, MYREAL)		       \
int write_##P##C##_##T##_from_site(QIO_Writer *outfile, \
         QIO_String *xml_record_out, field_offset src, int count){ \
  int status; \
  QIO_RecordInfo *rec_info; \
  /* We assume output precision is single */ \
  char qdptype[] = TYPESTRING;	\
  char prec[] = PSTRING;  \
  int datum_size = sizeof(FIXTYPE); \
  int word_size = sizeof(MYREAL); \
 \
  /* Create the record info for the field */ \
  rec_info = QIO_create_record_info(QIO_FIELD, 0, 0, 0, qdptype, prec, CVAL, \
  				    SVAL, datum_size, count); \
 \
  /* Write the record for the field */ \
  status = QIO_write(outfile, rec_info, xml_record_out,  \
		     vget_##P##C##_##T##_from_site, \
		     count*datum_size, word_size, (void *)&src); \
 if(status != QIO_SUCCESS)return 1; \
 \
  QIO_destroy_record_info(rec_info); \
 \
  return 0; \
}

/* Write a time slice of MILC site structure data */

#define make_write_tslice_from_site(P, PSTRING, C, CVAL, SVAL, T, TYPESTRING, \
			     FIXTYPE, VARTYPE, MYREAL)		       \
int write_##P##C##_##T##_timeslice_from_site(QIO_Writer *outfile, \
         QIO_String *xml_record_out, field_offset src, int count, int t0){ \
  int status; \
  QIO_RecordInfo *rec_info; \
  /* We assume output precision is single */ \
  char qdptype[] = TYPESTRING;	\
  char prec[] = PSTRING;  \
  int datum_size = sizeof(FIXTYPE); \
  int word_size = sizeof(MYREAL); \
  int lower[4] = {0, 0, 0, t0}; \
  int upper[4] = {nx-1, ny-1, nz-1, t0}; \
 \
  /* Create the record info for the field */ \
  rec_info = QIO_create_record_info(QIO_HYPER, lower, upper, 4, qdptype, \
       prec, CVAL, SVAL, datum_size, count); \
 \
  /* Write the record for the field */ \
  status = QIO_write(outfile, rec_info, xml_record_out,  \
		     vget_##P##C##_##T##_from_site, \
		     count*datum_size, word_size, (void *)&src); \
 if(status != QIO_SUCCESS)return 1; \
 \
  QIO_destroy_record_info(rec_info); \
 \
  return 0; \
}

/* Write MILC field data */

#define make_write_all_from_field(P, PSTRING, C, CVAL, SVAL, T, TYPESTRING, \
 	 	                  FIXTYPE, VARTYPE, MYREAL) \
int write_##P##C##_##T##_from_field(QIO_Writer *outfile, \
          QIO_String *xml_record_out, VARTYPE *src, int count){ \
  int status; \
  QIO_RecordInfo *rec_info; \
  char qdptype[] = TYPESTRING; \
  char prec[] = PSTRING; \
  int datum_size = sizeof(FIXTYPE); \
  int word_size = sizeof(MYREAL); \
 \
  /* Create the record info for the field */ \
  rec_info = QIO_create_record_info(QIO_FIELD, 0, 0, 0, qdptype, prec, CVAL, \
				    SVAL, datum_size, count); \
 \
  /* Write the record for the field */ \
  status = QIO_write(outfile, rec_info, xml_record_out,  \
		     vget_##P##C##_##T##_from_field,  \
		     count*datum_size, word_size, (void *)src); \
  if(status != QIO_SUCCESS)return 1; \
 \
  QIO_destroy_record_info(rec_info); \
 \
  return 0; \
}

/* Write a time slice of MILC field data */

#define make_write_tslice_from_field(P, PSTRING, C, CVAL, SVAL, T, \
           TYPESTRING, FIXTYPE, VARTYPE, MYREAL) \
int write_##P##C##_##T##_timeslice_from_field(QIO_Writer *outfile, \
          QIO_String *xml_record_out, VARTYPE *src, int count, int t0){ \
  int status; \
  QIO_RecordInfo *rec_info; \
  char qdptype[] = TYPESTRING; \
  char prec[] = PSTRING; \
  int datum_size = sizeof(FIXTYPE); \
  int word_size = sizeof(MYREAL); \
  int lower[4] = {0, 0, 0, t0}; \
  int upper[4] = {nx-1, ny-1, nz-1, t0}; \
 \
  /* Create the record info for the field */ \
  rec_info = QIO_create_record_info(QIO_HYPER, lower, upper, 4, qdptype, \
               prec, CVAL, SVAL, datum_size, count); \
 \
  /* Write the record for the field */ \
  status = QIO_write(outfile, rec_info, xml_record_out,  \
		     vget_##P##C##_##T##_from_field,  \
		     count*datum_size, word_size, (void *)src); \
  if(status != QIO_SUCCESS)return 1; \
 \
  QIO_destroy_record_info(rec_info); \
 \
  return 0; \
}

#define make_write_all(P, PSTRING, C, CVAL, SVAL, T, TYPESTRING, \
		     FIXTYPE, VARTYPE, MYREAL) \
  make_write_all_from_site(P, PSTRING, C, CVAL, SVAL, T, TYPESTRING, \
  		       FIXTYPE, VARTYPE, MYREAL); \
  make_write_all_from_field(P, PSTRING, C, CVAL, SVAL, T, TYPESTRING, \
 		       FIXTYPE, VARTYPE, MYREAL);

#define make_write_tslice(P, PSTRING, C, CVAL, SVAL, T, TYPESTRING, \
		     FIXTYPE, VARTYPE, MYREAL) \
  make_write_tslice_from_site(P, PSTRING, C, CVAL, SVAL, T, TYPESTRING, \
  		       FIXTYPE, VARTYPE, MYREAL); \
  make_write_tslice_from_field(P, PSTRING, C, CVAL, SVAL, T, TYPESTRING, \
 		       FIXTYPE, VARTYPE, MYREAL);

/* Single precision */

make_write_all(F, "F",  , 0, 0, R, "QLA_F_Real", float, Real, float);
make_write_all(F, "F",  , 0, 0, C, "QLA_F_Complex", fcomplex, complex, float);
make_write_all(F, "F", 3, 3, 0, V, "USQCD_F3_ColorVector", fsu3_vector, su3_vector, float);
make_write_all(F, "F", 3, 3, 0, M, "USQCD_F3_ColorMatrix", fsu3_matrix, su3_matrix, float);
make_write_all(F, "F", 3, 3, 4, D, "USQCD_F3_DiracFermion", fwilson_vector, wilson_vector, float);

make_write_tslice(F, "F",  , 0, 0, R, "QLA_F_Real", float, Real, float);
make_write_tslice(F, "F",  , 0, 0, C, "QLA_F_Complex", fcomplex, complex, float);
make_write_tslice(F, "F", 3, 3, 0, V, "USQCD_F3_ColorVector", fsu3_vector, su3_vector, float);
make_write_tslice(F, "F", 3, 3, 4, D, "USQCD_F3_DiracFermion", fwilson_vector, wilson_vector, float);

/* Double precision */

make_write_all(D, "D",  , 0, 0, C, "QLA_D_Complex", dcomplex, complex, double);
make_write_all(D, "D", 3, 3, 0, V, "USQCD_D3_ColorVector", dsu3_vector, su3_vector, double);
make_write_all(D, "D", 3, 3, 4, D, "USQCD_D3_DiracFermion", dwilson_vector, wilson_vector, double);
make_write_tslice(D, "D",  , 0, 0, C, "QLA_D_Complex", dcomplex, complex, double);
make_write_tslice(D, "D", 3, 3, 0, V, "USQCD_D3_ColorVector", dsu3_vector, su3_vector, double);
make_write_tslice(D, "D", 3, 3, 4, D, "USQCD_D3_DiracFermion", dwilson_vector, wilson_vector, double);


/* Read MILC site structure data */

#define make_read_to_site(P, C, T, FIXTYPE, VARTYPE, MYREAL) \
int read_##P##C##_##T##_to_site(QIO_Reader *infile, QIO_String *xml_record_in, \
				field_offset dest, int count) \
{ \
  QIO_RecordInfo rec_info; \
  int status; \
  int datum_size = sizeof(FIXTYPE); \
  int word_size = sizeof(MYREAL); \
 \
  /* Read the field record */ \
  status = QIO_read(infile, &rec_info, xml_record_in,  \
		    vput_##P##C##_##T##_to_site, datum_size*count, \
		    word_size, (void *)&dest); \
  node0_printf("Record info \n\"%s\"\n",QIO_string_ptr(xml_record_in)); \
  if(status != QIO_SUCCESS)return 1; \
 \
  node0_printf("Checksums %x %x\n", \
	       QIO_get_reader_last_checksuma(infile), \
	       QIO_get_reader_last_checksumb(infile)); \
 \
  return 0; \
}

/* Read MILC field data */

#define make_read_to_field(P, C, T, FIXTYPE, VARTYPE, MYREAL) \
int read_##P##C##_##T##_to_field(QIO_Reader *infile, \
		 QIO_String *xml_record_in, VARTYPE *dest, int count) \
{ \
  QIO_RecordInfo rec_info; \
  int status; \
  int datum_size = sizeof(FIXTYPE); \
  int word_size = sizeof(MYREAL); \
 \
  /* Read the field record */ \
  status = QIO_read(infile, &rec_info, xml_record_in,  \
		    vput_##P##C##_##T##_to_field, datum_size*count, \
		    word_size, (void *)dest); \
  node0_printf("Record info \n\"%s\"\n",QIO_string_ptr(xml_record_in)); \
  if(status != QIO_SUCCESS)return 1; \
 \
  node0_printf("Checksums %x %x\n", \
	       QIO_get_reader_last_checksuma(infile), \
	       QIO_get_reader_last_checksumb(infile)); \
 \
  return 0; \
}

#define make_read(P, C, T, FIXTYPE, VARTYPE, MYREAL) \
  make_read_to_site(P, C, T, FIXTYPE, VARTYPE, MYREAL); \
  make_read_to_field(P, C, T, FIXTYPE, VARTYPE, MYREAL);

/* Single precision */
make_read(F,  , R, float, Real, float);
make_read(F,  , C, fcomplex, complex, float);
make_read(F, 3, V, fsu3_vector, su3_vector, float);
make_read(F, 3, M, fsu3_matrix, su3_matrix, float);
make_read(F, 3, D, fwilson_vector, wilson_vector, float);

/* Double precision */
make_read(D,  , C, dcomplex, complex, double);
make_read(D, 3, V, dsu3_vector, su3_vector, double);
make_read(D, 3, M, dsu3_matrix, su3_matrix, double);
make_read(D, 3, D, dwilson_vector, wilson_vector, double);

/* Factory function for moving random generator state from site structure to
   output */

void vget_S_from_site(char *buf, size_t index, int count, void *arg)
{
  char *dest = buf;
  /* arg contains pointer to field offset value */
  field_offset src = *((field_offset *)arg);
  site *s = &lattice[index];
  char *src_prn = (char *)F_PT(s,src);

  memcpy(dest, src_prn, sizeof(double_prn));
}

int write_S_from_site(QIO_Writer *outfile, QIO_String *xml_record_out,
		      field_offset src){
  int status;
  QIO_RecordInfo *rec_info;
  char qdptype[] = "MILC_RandomState";
  char prec[] = "";
  int datum_size = sizeof(double_prn);
  int word_size = sizeof(float);
  int count = 1;

  /* Create the record info for the field */
  rec_info = QIO_create_record_info(QIO_FIELD, 0, 0, 0, qdptype, prec, 0,
				    0, datum_size, count);
  /* Write the record for the field */
  status = QIO_write(outfile, rec_info, xml_record_out, vget_S_from_site, 
		     count*datum_size, word_size, (void *)&src);
  if(status != QIO_SUCCESS)return 1;

  QIO_destroy_record_info(rec_info);

  return 0;
}

/* Factory function for moving random number state from input
   to site structure */
void vput_S_to_site(char *buf, size_t index, int count, void *arg)
{
  /* Assume input lattice is single precision, 3 colors */
  char *src = buf;
  field_offset dest = *((field_offset *)arg);
  site *s = &lattice[index];
  char *dest_prn = (char *)F_PT(s,dest);
  
  memcpy(dest_prn, src, sizeof(double_prn));
}

/* Read random number state */
int read_S_to_site(QIO_Reader *infile,  QIO_String *xml_record_in,
		   field_offset dest)
{
  QIO_RecordInfo rec_info;
  int status;
  int count = 1;
  int datum_size = sizeof(double_prn);
  int word_size = sizeof(float);
  
  /* Read the field record */
  status = QIO_read(infile, &rec_info, xml_record_in, 
		    vput_S_to_site, datum_size*count, word_size, (void *)&dest);
  node0_printf("Record info \n\"%s\"\n",QIO_string_ptr(xml_record_in));
  if(status != QIO_SUCCESS)return 1;

  node0_printf("Checksums %x %x\n",
	       QIO_get_reader_last_checksuma(infile),
	       QIO_get_reader_last_checksumb(infile));

  return 0;
}

