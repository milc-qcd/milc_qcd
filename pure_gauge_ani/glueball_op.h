/* Glueball operator shapes */

#define NGLUEBALL_BASIC_SHAPES 4
#define MAX_LOOP_LENGTH 8
#define MAX_ORIENTATION 48 
#define MAX_CLASS 48
#include "../include/dirs.h"

typedef struct
{
  char name[12];                  /* name of path */
  int length;                     /* length of path */
  int dir[MAX_LOOP_LENGTH];       /* list of directions in path */
} link_path;

/* Table of basic shapes.  If you change it, change
   NGLUEBALL_BASIC_SHAPES and possibly MAX_LOOP_LENGTH */

const link_path glueball_shape[NGLUEBALL_BASIC_SHAPES] =
{
  { "1x1",   4, XUP, YUP, XDOWN, YDOWN, NODIR, NODIR, NODIR, NODIR },
  { "1x2",   6, XUP, XUP, YUP,   XDOWN, XDOWN, YDOWN, NODIR, NODIR },
  { "1x1x1", 6, XUP, YUP, ZUP,   XDOWN, YDOWN, ZDOWN, NODIR, NODIR },
  { "E.B",   8, TUP, XUP, TDOWN, XDOWN, YUP, ZUP, YDOWN, ZDOWN }
};

/* Structure for holding a glueball operator to be measured */ 
typedef struct
{
  link_path path;               /* Definition of loop shape */
  int weight_r[N_CUBIC_IRREP];  /* Character sum for real part for each irrep */
  int weight_i[N_CUBIC_IRREP];  /* Character sum for imag part for each irrep */
  double_complex *trace;        /* Trace vs t */
} oriented_loop;

/* Structure for holding all orientations of a glueball operator */
typedef struct
{
  int norientations;                   /* Number of orientations */
  oriented_loop loop[MAX_ORIENTATION]; /* Oriented operators */
} loop_orientations;

/* Table of glueball operators to be measured */
/* This plays the role of "loop_table" in gauge_stuff.c" */
loop_orientations glueball_op[NGLUEBALL_BASIC_SHAPES];



