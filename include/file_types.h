#ifndef FILE_TYPES_H
#define FILE_TYPES_H

/* For tables of supported file types */
typedef struct {
  int type;
  int32type magic_no;
} file_type;

/* LIME (SCIDAC/QIO) files (all fields) */
#define LIME_MAGIC_NO     1164413355 /* CHECK WITH lime_defs.h */

/* Gauge configuration file types */

#define GAUGE_VERSION_NUMBER_V1 59354   /* Versions 1-4 */
#define GAUGE_VERSION_NUMBER 20103  /* Versions 5 and 6 */
#define GAUGE_VERSION_NUMBER_1996 53546
#define GAUGE_VERSION_NUMBER_FNAL  0x71626434     /* field major order for FNAL
             project Nov. 2002 S.G. Note it is the same for different fields */
#define GAUGE_VERSION_NUMBER_ARCHIVE 1111836489

/* Wilson propagator file types */

#define W_PROP_VERSION_NUMBER 12781 /* Versions 5 and 6 */
#define W_PROP_VERSION_NUMBER_1996 48291

/* KS propagator file types */

#define KSPROP_VERSION_NUMBER_V0 230201 /* save_ksprop_ascii order, 2001 */
#define KSPROP_VERSION_NUMBER 23209     /* w_ascii_ks order, June 2002 */
#define KSFMPROP_VERSION_NUMBER  0x71626434     /* field major order for FNAL 
		project Nov. 2002 S.G.  */

/* These are just an enumeration */
#define FILE_TYPE_KSPROP 0
#define FILE_TYPE_KSFMPROP 1
#define FILE_TYPE_KSQIOPROP 2
#define N_KSPROP_TYPES 3

static file_type ksprop_list[N_KSPROP_TYPES] =
  { {FILE_TYPE_KSPROP,       KSPROP_VERSION_NUMBER},
    {FILE_TYPE_KSFMPROP,     KSFMPROP_VERSION_NUMBER},
    {FILE_TYPE_KSQIOPROP,    LIME_MAGIC_NO}
  };

#endif
