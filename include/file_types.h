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
#define LIME_MAGIC_NO                0x456789ab /* decimal 1164413355
						   see lime_defs.h */

/* FNAL (all files) */
/* Files may be identified through the header which contains element
   sizes and counts */
#define IO_UNI_MAGIC                 0x71626434

/* Gauge configuration file types */

#define GAUGE_VERSION_NUMBER_V1      0xe7da  /* decimal 59354 Versions 1-4 */
#define GAUGE_VERSION_NUMBER         0x4e87  /* decimal 20103 Versions 5-7 */
#define GAUGE_VERSION_NUMBER_1996    0xd12a  /* decimal 53546 */
#define GAUGE_VERSION_NUMBER_FNAL    IO_UNI_MAGIC
#define GAUGE_VERSION_NUMBER_ARCHIVE 0x42454749  /* 1111836489 decimal */

/* Wilson propagator file types */

#define W_PROP_VERSION_NUMBER        0x31ed /* 12781 decimal Versions 5-7 */
#define W_PROP_VERSION_NUMBER_1996   0xbca3 /* 48291 decimal */
#define W_FMPROP_VERSION_NUMBER      IO_UNI_MAGIC

/* KS propagator file types */

#define KSPROP_VERSION_NUMBER_V0     0x38339 /* 230201 decimal ca 2001 */
#define KSPROP_VERSION_NUMBER        0x5aa9  /* 23209 decimal ca June 2002 */
#define KSFMPROP_VERSION_NUMBER      IO_UNI_MAGIC

/* Tables */

#define FILE_TYPE_GAUGE_V1       0
#define FILE_TYPE_GAUGE_V5       1
#define FILE_TYPE_GAUGE_1996     2
#define FILE_TYPE_GAUGE_FNAL     3
#define FILE_TYPE_GAUGE_ARCHIVE  4
#define FILE_TYPE_GAUGE_SCIDAC   5
#define N_GAUGE_TYPES            6

static file_type gauge_list[N_GAUGE_TYPES] =
  { {FILE_TYPE_GAUGE_V1,      GAUGE_VERSION_NUMBER_V1},
    {FILE_TYPE_GAUGE_V5,      GAUGE_VERSION_NUMBER},
    {FILE_TYPE_GAUGE_1996,    GAUGE_VERSION_NUMBER_1996},
    {FILE_TYPE_GAUGE_FNAL,    GAUGE_VERSION_NUMBER_FNAL},
    {FILE_TYPE_GAUGE_ARCHIVE, GAUGE_VERSION_NUMBER_ARCHIVE},
    {FILE_TYPE_GAUGE_SCIDAC,  LIME_MAGIC_NO}
  };

#define FILE_TYPE_W_PROP       0
#define FILE_TYPE_W_PROP_1996  1
#define FILE_TYPE_W_FMPROP     2
#define FILE_TYPE_W_QIOPROP    3
#define N_WPROP_TYPES          4

static file_type w_prop_list[N_WPROP_TYPES] =
  { {FILE_TYPE_W_PROP,       W_PROP_VERSION_NUMBER},
    {FILE_TYPE_W_PROP_1996,  W_PROP_VERSION_NUMBER_1996},
    {FILE_TYPE_W_FMPROP,     W_FMPROP_VERSION_NUMBER},
    {FILE_TYPE_W_QIOPROP,    LIME_MAGIC_NO}
  };

/* These are just an enumeration */
#define FILE_TYPE_KSPROP     0
#define FILE_TYPE_KSFMPROP   1
#define FILE_TYPE_KSQIOPROP  2
#define N_KSPROP_TYPES       3

static file_type ksprop_list[N_KSPROP_TYPES] =
  { {FILE_TYPE_KSPROP,       KSPROP_VERSION_NUMBER},
    {FILE_TYPE_KSFMPROP,     KSFMPROP_VERSION_NUMBER},
    {FILE_TYPE_KSQIOPROP,    LIME_MAGIC_NO}
  };

#endif
