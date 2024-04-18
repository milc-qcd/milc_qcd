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

/* defines for 3rd nearest neighbor (NAIK) stuff */
#define X3UP 8
#define Y3UP 9
#define Z3UP 10
#define T3UP 11
#define T3DOWN 12
#define Z3DOWN 13
#define Y3DOWN 14
#define X3DOWN 15

#define OPP_3_DIR(dir) (23-(dir))
#define DIR3(dir) ((dir)+8)
#define FORALL3UPDIR(dir) for(dir=X3UP; dir<=T3UP; dir++)

/* defines for 2nd nearest neighbor stuff */
/* Added by Hwancheol Jeong 4/2024 */
#define X2UP 16
#define Y2UP 17
#define Z2UP 18
#define T2UP 19
#define T2DOWN 20
#define Z2DOWN 21
#define Y2DOWN 22
#define X2DOWN 23

#define OPP_2_DIR(dir) (39-(dir))
#define DIR2(dir) ((dir)+16)
#define FORALL2UPDIR(dir) for(dir=X2UP; dir<=T2UP; dir++)

#endif /* _DIRS_H */
