/**************************** bsd_sum.c ****************************/
/* MIMD version 7 */
/*  This function calculates a simple checksum. It uses the same algorithm as
    the UNIX routine "sum".

    Usage:

    int   data[ndata]  ;     checksum = bsd_sum((char *) &data[0] , sizeof(int)*ndata)
    Real fdata[ndata]  ;    checksum = bsd_sum((char *) &fdata[0], sizeof(Real)*ndata) 

    
    This function produced the same checksum on:  IBM RISC/6000, cray T3D,
    and a paragon. 
    
    The checksum routine should be called BEFORE any byte reversing is done.

    This function was hacked out of the gnu distribution's sum.c code;
    see below.

 */


/* sum -- checksum and count the blocks in a file
   Copyright (C) 86, 89, 91, 95, 1996 Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.  */

/* Like BSD sum or SysV sum -r, except like SysV sum if -s option is given. */

/* Written by Kayvan Aghaiepour and David MacKenzie. */


/** required for the definition of int32type ***/

#include "generic_includes.h"


/* Right-rotate 32-bit integer variable C. */
#define ROTATE_RIGHT(c) if ((c) & 01) (c) = ((c) >>1) + 0x8000; else (c) >>= 1;


int32type bsd_sum (char *data,int32type total_bytes)
{
  register  int32type checksum = 0; /* The checksum mod 2^16. */
  int32type i ;

  for( i= 0 ; i < total_bytes ; ++i)
  {
    ROTATE_RIGHT (checksum);
    checksum += data[i] ;
    checksum &= 0xffff;	/* Keep it within bounds. */
  }

  return (int32type) checksum   ;
}

