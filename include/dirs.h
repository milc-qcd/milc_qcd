#ifndef _DIRS_H
#define _DIRS_H

/* Directions, and a macro to give the opposite direction */
/*  These must go from 0 to 7 because they will be used to index an
    array. */
/* Also define NDIRS = number of directions */
#define XUP 0
#define YUP 1
#define ZUP 2
#define TUP 3
#define TDOWN 4
#define ZDOWN 5
#define YDOWN 6
#define XDOWN 7

#define NODIR -1  /* not a direction */

#define OPP_DIR(dir)	(7-(dir))	/* Opposite direction */
#define NDIRS 8				/* number of directions */

#endif /* _DIRS_H */
