/******** combine_files.c *****************/
/* MIMD version 7 */
/* Modifications ...

 * 3/25/12 CD Created  */

#include "file_combine_includes.h"
#include "qioxml_usqcd_ksprop_source.h"
#include "qioxml_usqcd_prop_source.h"
#include "../include/io_scidac_ks.h"
#include "../include/io_scidac_w.h"
#include "../include/io_scidac.h"

static int check_ks_color(QIO_String *recxml, int color){
  int status;
  int input_color;
  QIO_USQCDKSPropRecordInfo recinfo;
  char myname[] = "check_ks_color";

  status = QIO_decode_usqcd_ksproprecord_info(&recinfo, recxml);
  if(status != QIO_SUCCESS){
    node0_printf("%s: Can't decode the record info\n%s\n", myname, QIO_string_ptr(recxml));
    //    terminate(1);
    return 1;
  }
  input_color = QIO_get_usqcd_ksproprecord_color(&recinfo);
  if(color != input_color){
    node0_printf("%s(%d): Error: expected color %d got %d\n",
		 myname, this_node, color, input_color);
    return 1;
  }
  return 0;
}

static int check_wv_color_spin(QIO_String *recxml, int color, int spin){
  int status;
  int input_spin;
  int input_color;
  QIO_USQCDPropRecordInfo recinfo;
  char myname[] = "check_wv_color_spin";

  status = QIO_decode_usqcd_proprecord_info(&recinfo, recxml);
  if(status != QIO_SUCCESS){
    node0_printf("%s: Can't decode the record info\n%s\n", myname, QIO_string_ptr(recxml));
    //    terminate(1);
    return 1;
  }
  input_color = QIO_get_usqcd_proprecord_color(&recinfo);
  input_spin = QIO_get_usqcd_proprecord_spin(&recinfo);
  if(color != input_color || spin != input_spin){
    node0_printf("%s(%d): Error: expected color %d and spin %d got %d and %d\n",
		 myname, this_node, color, spin, input_color, input_spin);
    return 1;
  }
  return 0;
}

static void scalar_mult_add_latvec_field(su3_vector *a, su3_vector *b, Real s, su3_vector *c){
  int i;
  FORALLFIELDSITES(i){
    scalar_mult_add_su3_vector(a+i, b+i, s, c+i);
  }
}

static void scalar_mult_add_latwvec_field(wilson_vector *a, wilson_vector *b, Real s, wilson_vector *c){
  int i;
  FORALLFIELDSITES(i){
    scalar_mult_add_wvec(a+i, b+i, s, c+i);
  }
}

/* Read nfile files containing ncolor color vector fields, take the
   linear combination specified by coeff, and write the result as a
   file of the same format  */

static void combine_vector_field_files(int nfile, int ncolor, int t0, 
				       int startflag[],
				       char startfile[MAX_FILES][MAXFILENAME], 
				       Real coeff[], int saveflag, 
				       int savetype, char savefile[]){
  int file, color, volfmt, serpar, status;
  su3_vector *v_dst = create_v_field();
  su3_vector *v_src = create_v_field();
  QIO_String *xml_file, *xml_record;
  QIO_USQCDKSPropSourceFileInfo *sourcefile_info;
  QIO_Reader **infile = (QIO_Reader **)malloc(nfile*sizeof(QIO_Reader *));
  QIO_Writer *outfile;
  QIO_USQCDKSPropRecordInfo *recinfo;
  char fileinfo[] = "";
  char myname[] = "combine_vector_field_files";

  if(v_dst == NULL || v_src == NULL || infile == NULL){
    printf("(%s): No room\n", myname);
    terminate(1);
  }

  /* Open all the input files */
  for(file = 0; file < nfile; file++){
    if(startflag[file] == RELOAD_PARALLEL)
      serpar = QIO_PARALLEL;
    else
      serpar = QIO_SERIAL;

    xml_file = QIO_string_create();
    infile[file] = 
      r_open_scidac_file_xml(startfile[file], 
				       serpar, xml_file);
    QIO_string_destroy(xml_file);  /* Ignore for now */

    if(infile[file] == NULL){
      node0_printf("Failed to open %s for reading\n",startfile[file]);
      terminate(1);
    }
  }

  /* Open the output file */

  interpret_usqcd_ks_save_flag(&volfmt, &serpar, saveflag);
  xml_file = QIO_string_create();
  sourcefile_info = QIO_create_usqcd_kspropsourcefile_info(fileinfo);
  
  outfile = w_open_scidac_file(savefile, fileinfo, volfmt, serpar);
  QIO_string_destroy(xml_file);
  
  /* Loop over colors */

  for(color = 0; color < ncolor; color++){

    clear_v_field(v_dst);

    /* Read one color from each file and combine */
    for(file = 0; file < nfile; file++){
      xml_record = QIO_string_create();
      clear_v_field(v_src);
      status = read_ks_vector_scidac_xml(infile[file], v_src, 1, xml_record );
      if(status == QIO_SUCCESS){
	if(check_ks_color(xml_record, color))
	  terminate(1);
      }
      else if(status == QIO_EOF){
	node0_printf("Unexpected EOF on %s\n",startfile[file]);
	terminate(1);
      }
      QIO_string_destroy(xml_record);
      
      /* Accumulate */
      scalar_mult_add_latvec_field(v_dst, v_src, coeff[file], v_dst);
    }

    /* Write result for this color */
    /* Construct the record XML */
    xml_record = QIO_string_create();
    recinfo = QIO_create_usqcd_ksproprecord_c_info(color, "");
    QIO_encode_usqcd_ksproprecord_info(xml_record, recinfo);
    QIO_destroy_usqcd_ksproprecord_info(recinfo);
    if(t0 == ALL_T_SLICES){
      if(PRECISION==1)
	status = write_F3_V_from_field(outfile, xml_record, v_dst, 1 );
      else
	status = write_D3_V_from_field(outfile, xml_record, v_dst, 1 );
    } else {
      status = write_kspropsource_V_usqcd_xml(outfile, xml_record, v_dst, t0);
    }
    QIO_string_destroy(xml_record);
  }
    
  w_close_scidac_file(outfile);

  for(file = 0; file < nfile; file++){
    r_close_scidac_file(infile[file]);
  }
}

/* Read nfile files containing ncolor x nspin Dirac vector fields,
   take the linear combination specified by coeff, and write the
   result as a file of the same format.  The spin varies most rapidly. 
*/

static void combine_dirac_field_files(int nfile, int ncolor, int nspin, int t0,
				      int startflag[],
				      char startfile[MAX_FILES][MAXFILENAME], 
				      Real coeff[], int saveflag, 
				      int savetype, char savefile[]){
  int file, color, spin, volfmt, serpar, status;
  wilson_vector *wv_dst = create_wv_field();
  wilson_vector *wv_src = create_wv_field();
  QIO_String *xml_file, *xml_record;
  QIO_USQCDPropSourceFileInfo *sourcefile_info;
  QIO_Reader **infile = (QIO_Reader **)malloc(nfile*sizeof(QIO_Reader *));
  QIO_Writer *outfile;
  QIO_USQCDPropRecordInfo *recinfo;
  char fileinfo[] = "";
  char myname[] = "combine_dirac_field_files";

  if(wv_dst == NULL || wv_src == NULL || infile == NULL){
    printf("(%s): No room\n", myname);
    terminate(1);
  }

  /* Open all the input files */
  for(file = 0; file < nfile; file++){
    if(startflag[file] == RELOAD_PARALLEL)
      serpar = QIO_PARALLEL;
    else
      serpar = QIO_SERIAL;

    xml_file = QIO_string_create();
    infile[file] = 
      r_open_w_vector_scidac_file_xml(startfile[file], 
				      serpar, xml_file);
    QIO_string_destroy(xml_file);  /* Ignore for now */

    if(infile[file] == NULL){
      node0_printf("Failed to open %s for reading\n",startfile[file]);
      terminate(1);
    }
  }

  /* Open the output file */

  interpret_usqcd_w_save_flag(&volfmt, &serpar, saveflag);
  xml_file = QIO_string_create();
  sourcefile_info = QIO_create_usqcd_propsourcefile_info(fileinfo);
  
  outfile = w_open_w_vector_scidac_file(savefile, fileinfo,
					volfmt, serpar);
  QIO_string_destroy(xml_file);
  
  /* Loop over colors */

  for(color = 0; color < ncolor; color++)
    for(spin = 0; spin < nspin; spin++){
      
      clear_wv_field(wv_dst);
      
      /* Read one color and spin from each file and combine */
      for(file = 0; file < nfile; file++){
	xml_record = QIO_string_create();
	clear_wv_field(wv_src);
	status = read_w_vector_scidac_xml(infile[file], wv_src, 1, xml_record );
	if(status == QIO_SUCCESS){
	  if(check_wv_color_spin(xml_record, color, spin))
	    terminate(1);
	}
	else if(status == QIO_EOF){
	  node0_printf("Unexpected EOF on %s\n",startfile[file]);
	  terminate(1);
	}
	QIO_string_destroy(xml_record);
	
	/* Accumulate */
	scalar_mult_add_latwvec_field(wv_dst, wv_src, coeff[file], wv_dst);
      }
      
      /* Write result for this color */
      /* Construct the record XML */
      xml_record = QIO_string_create();
      recinfo = QIO_create_usqcd_proprecord_sc_info(spin, color, "");
      QIO_encode_usqcd_proprecord_info(xml_record, recinfo);
      QIO_destroy_usqcd_proprecord_info(recinfo);

      if(PRECISION==1) {
        if (t0 == ALL_T_SLICES)
	  status = write_F3_D_from_field(outfile, xml_record, wv_dst, 1 );
        else
	  status = write_F3_D_timeslice_from_field(outfile, xml_record, wv_dst, 1,t0);
      }

      else {
        if (t0 == ALL_T_SLICES)
	  status = write_D3_D_from_field(outfile, xml_record, wv_dst, 1 );
        else
	  status = write_D3_D_timeslice_from_field(outfile, xml_record, wv_dst, 1,t0);
      }

      QIO_string_destroy(xml_record);

    } /* color and spin */
  
  w_close_w_vector_scidac_file(outfile);
  
  for(file = 0; file < nfile; file++){
    r_close_w_vector_scidac_file(infile[file]);
  }
}

static void combine_vector_propagator_files(int nfile, int ncolor, int t0,
					    int startflag[],
					    char startfile[MAX_FILES][MAXFILENAME], 
					    Real coeff[], int saveflag, 
					    int savetype, char savefile[]){
}

static void combine_dirac_propagator_files(int nfile, int ncolor, int nspin, int t0,
					   int startflag[],
					    char startfile[MAX_FILES][MAXFILENAME], 
					    Real coeff[], int saveflag, 
					    int savetype, char savefile[]){
}


void combine_files(int nfile, int file_type, int ncolor, int nspin, int t0,
		   int startflag[], char startfile[MAX_FILES][MAXFILENAME], 
		   Real coeff[], int saveflag,
		   int savetype, char *savefile){
  
  if(file_type == VECTOR_FIELD_FILE)
    combine_vector_field_files(nfile, ncolor, t0, startflag, startfile, coeff,
			       saveflag, savetype, savefile);
  else if(file_type == DIRAC_FIELD_FILE)
    combine_dirac_field_files(nfile, ncolor, nspin, t0, startflag, startfile, coeff,
			      saveflag, savetype, savefile);
  else if(file_type == VECTOR_PROPAGATOR_FILE)
    combine_vector_propagator_files(nfile, ncolor, t0, startflag, startfile, coeff,
				    saveflag, savetype ,savefile);
  else if(file_type == DIRAC_PROPAGATOR_FILE)
    combine_dirac_propagator_files(nfile, ncolor, nspin, t0, startflag, startfile, coeff,
				   saveflag, savetype, savefile);
  else{
    node0_printf("combine_files: Unrecognized file type\n");
  }
}
