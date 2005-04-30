#ifndef _PRECISION_H
#define _PRECISION_H

/* generic floating point type.  Defaults to single precision. */
#ifndef PRECISION
#define PRECISION 1
#endif

#if (PRECISION == 2)
typedef double Real;
#else
typedef float Real;
#endif

#endif /* _PRECISION_H */
