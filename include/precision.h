#ifndef _PRECISION_H
#define _PRECISION_H

/* generic floating point type.  Defaults to single precision. */
#ifndef MILC_PRECISION
#define MILC_PRECISION 1
#endif

#if (MILC_PRECISION == 2)
typedef double Real;
#else
typedef float Real;
#endif

#endif /* _PRECISION_H */
