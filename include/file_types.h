#ifndef FILE_TYPES_H
#define FILE_TYPES_H

/* LIME (SCIDAC/QIO)(all files) */
/* Files may be identified by the record XML containing site data
   sizes and counts */
#if HAVE_QIO
#include <qio.h>
#else
#define LIME_MAGIC_NO                0x456789ab /* decimal 1164413355
						   see lime_defs.h */
#endif

/* FNAL (all files) */
/* Files may be identified through the header which contains element
   sizes and counts */
#define IO_UNI_MAGIC                 0x71626434

/* Magic numbers for MILC gauge configuration files */

#define GAUGE_VERSION_NUMBER_V1      0xe7da  /* decimal 59354 Versions 1-4 */
#define GAUGE_VERSION_NUMBER         0x4e87  /* decimal 20103 Versions 5-7 */
#define GAUGE_VERSION_NUMBER_1996    0xd12a  /* decimal 53546 */
#define GAUGE_VERSION_NUMBER_ARCHIVE 0x42454749  /* 1111836489 decimal */

/* Magic numbers for MILC Wilson propagator file types (deprecated) */

#define W_PROP_VERSION_NUMBER        0x31ed /* 12781 decimal Versions 5-7 */
#define W_PROP_VERSION_NUMBER_1996   0xbca3 /* 48291 decimal */

/* Magic numbers for MILC KS propagator file types */

#define KSPROP_VERSION_NUMBER_V0     0x38339 /* 230201 decimal ca 2001 */
#define KSPROP_VERSION_NUMBER        0x5aa9  /* 23209 decimal ca June 2002 */

/* Coding for nonspecific types */

#define FILE_TYPE_UNKNOWN        -1
#define FILE_TYPE_FM              1
#define FILE_TYPE_LIME            2

/* Coding for gauge field types */

#define FILE_TYPE_GAUGE_V5       10
#define FILE_TYPE_GAUGE_FNAL     11
#define FILE_TYPE_GAUGE_ARCHIVE  12
#define FILE_TYPE_GAUGE_SCIDAC   13
#define N_GAUGE_TYPES             4

/* Coding for Wilson propagator types */

#define FILE_TYPE_W_FMPROP         20
#define FILE_TYPE_W_USQCD_C1D12    21
#define FILE_TYPE_W_USQCD_DD_PAIRS 22
#define FILE_TYPE_W_USQCD_CD_PAIRS 23   
#define FILE_TYPE_W_USQCD_LHPC     24   
#define N_WPROP_TYPES               5

/* For a Wilson propagator that was read and cached */
#define FILE_TYPE_W_STORE  29

/* Coding for staggered propagator types */

#define FILE_TYPE_KS_PROP           30
#define FILE_TYPE_KS_FMPROP         31
#define FILE_TYPE_KS_USQCD_C1V3     32
#define FILE_TYPE_KS_USQCD_VV_PAIRS 33
#define FILE_TYPE_KS_USQCD_CV_PAIRS 34
#define N_KSPROP_TYPES               5

/* For a KS propagator that was read and cached */
#define FILE_TYPE_KS_STORE  39

/* Table of broad file types distinguishable by their magic numbers */
/* To distinguish the FNAL and LIME subtypes requires further poking around */

typedef struct {
  int type;
  int32type magic_no;
} file_table;


#endif
