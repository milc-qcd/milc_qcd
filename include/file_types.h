#ifndef FILE_TYPES_H
#define FILE_TYPES_H

/* For tables of supported file types */
typedef struct {
  int type;
  int32type magic_no;
} file_type;

/* LIME (SCIDAC/QIO)(all files) */
/* Files may be identified by the record XML containing site data
   sizes and counts */
#ifndef HAVE_QIO
#define LIME_MAGIC_NO                0x456789ab /* decimal 1164413355
						   see lime_defs.h */
#endif

/* FNAL (all files) */
/* Files may be identified through the header which contains element
   sizes and counts */
#define IO_UNI_MAGIC                 0x71626434

/* Gauge configuration file types */

#define GAUGE_VERSION_NUMBER_V1      0xe7da  /* decimal 59354 Versions 1-4 */
#define GAUGE_VERSION_NUMBER         0x4e87  /* decimal 20103 Versions 5-7 */
#define GAUGE_VERSION_NUMBER_1996    0xd12a  /* decimal 53546 */
#define GAUGE_VERSION_NUMBER_ARCHIVE 0x42454749  /* 1111836489 decimal */

/* Wilson propagator file types */

#define W_PROP_VERSION_NUMBER        0x31ed /* 12781 decimal Versions 5-7 */
#define W_PROP_VERSION_NUMBER_1996   0xbca3 /* 48291 decimal */

/* KS propagator file types */

#define KSPROP_VERSION_NUMBER_V0     0x38339 /* 230201 decimal ca 2001 */
#define KSPROP_VERSION_NUMBER        0x5aa9  /* 23209 decimal ca June 2002 */

/* Tables */

/* Nonspecific types */

#define FILE_TYPE_FM              0
#define FILE_TYPE_LIME            1

#define FILE_TYPE_GAUGE_V1       10
#define FILE_TYPE_GAUGE_V5       11
#define FILE_TYPE_GAUGE_1996     12
#define FILE_TYPE_GAUGE_FNAL     13
#define FILE_TYPE_GAUGE_ARCHIVE  14
#define FILE_TYPE_GAUGE_SCIDAC   15
#define N_GAUGE_TYPES             6

#define FILE_TYPE_W_PROP           20
#define FILE_TYPE_W_PROP_1996      21
#define FILE_TYPE_W_FMPROP         22
#define FILE_TYPE_W_USQCD_C1D12    23
#define FILE_TYPE_W_USQCD_DD_PAIRS 24
#define FILE_TYPE_W_USQCD_CD_PAIRS 25   
#define FILE_TYPE_W_USQCD_LHPC     26   
#define N_WPROP_TYPES               7

#define FILE_TYPE_KSPROP     30
#define FILE_TYPE_KSFMPROP   31
#define FILE_TYPE_KSQIOPROP  32
#define N_KSPROP_TYPES        3

/* For a Wilson propagator that was read and cached */
#define FILE_TYPE_W_STORE  30

#endif
