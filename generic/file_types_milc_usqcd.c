/******************** file_types_milc_usqcd.c *********************************/
/* MIMD version 7 */

/* For compilation with QIO */

#include "generic_includes.h"
#include <qio.h>

/********************************************************************/
/* Translate from MILC to USQCD file type */
int ks_prop_usqcd_to_milc(int usqcd_type){

  int milc_type = -1;

  switch(usqcd_type){
  case QIO_USQCDKSPROPFILETYPE_C1V3:
    milc_type = FILE_TYPE_KS_USQCD_C1V3;
    break;
  case QIO_USQCDKSPROPFILETYPE_VV_PAIRS:
    milc_type = FILE_TYPE_KS_USQCD_VV_PAIRS;
    break;
  case QIO_USQCDKSPROPFILETYPE_CV_PAIRS:
    milc_type = FILE_TYPE_KS_USQCD_CV_PAIRS;
    break;
  }
  return milc_type;
}

/********************************************************************/
/* Translate from MILC to USQCD file type */
int ks_prop_milc_to_usqcd(int milc_type){

  int usqcd_type = -1;

  switch(milc_type){
  case FILE_TYPE_KS_USQCD_C1V3:
    usqcd_type = QIO_USQCDKSPROPFILETYPE_C1V3;
    break;
  case FILE_TYPE_KS_USQCD_VV_PAIRS:
    usqcd_type = QIO_USQCDKSPROPFILETYPE_VV_PAIRS;
    break;
  case FILE_TYPE_KS_USQCD_CV_PAIRS:
    usqcd_type = QIO_USQCDKSPROPFILETYPE_CV_PAIRS;
    break;
  }
  return usqcd_type;
}

/********************************************************************/
/* Translate from MILC to USQCD file type */
int w_prop_milc_to_usqcd(int milc_type){

  int usqcd_type = -1;

  switch(milc_type){
  case FILE_TYPE_W_USQCD_C1D12:
    usqcd_type = QIO_USQCDPROPFILETYPE_C1D12;
    break;
  case FILE_TYPE_W_USQCD_DD_PAIRS:
    usqcd_type = QIO_USQCDPROPFILETYPE_DD_PAIRS;
    break;
  case FILE_TYPE_W_USQCD_CD_PAIRS:
    usqcd_type = QIO_USQCDPROPFILETYPE_CD_PAIRS;
    break;
  case FILE_TYPE_W_USQCD_LHPC:
    usqcd_type = QIO_USQCDPROPFILETYPE_LHPC;
    break;
  }
  return usqcd_type;
}

/********************************************************************/
/* Translate from MILC to USQCD file type */
int w_prop_usqcd_to_milc(int usqcd_type){

  int milc_type = -1;

  switch(usqcd_type){
  case QIO_USQCDPROPFILETYPE_C1D12:
    milc_type = FILE_TYPE_W_USQCD_C1D12;
    break;
  case QIO_USQCDPROPFILETYPE_DD_PAIRS:
    milc_type = FILE_TYPE_W_USQCD_DD_PAIRS;
    break;
  case QIO_USQCDPROPFILETYPE_CD_PAIRS:
    milc_type = FILE_TYPE_W_USQCD_CD_PAIRS;
    break;
  case QIO_USQCDPROPFILETYPE_LHPC:
    milc_type = FILE_TYPE_W_USQCD_LHPC;
    break;
  }
  return milc_type;
}

