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


#include <stdio.h>


/* Right-rotate 32-bit integer variable C. */
#define ROTATE_RIGHT(c) if ((c) & 01) (c) = ((c) >>1) + 0x8000; else (c) >>= 1;


/* Calculate and print the rotated checksum and the size in 1K blocks
   of file FILE, or of the standard input if FILE is "-".
   If PRINT_NAME is >1, print FILE next to the checksum and size.
   The checksum varies depending on sizeof(int).
   Return 0 if successful, -1 if an error occurs. */

int bsd_sum (char *data,int total_bytes)
{
  register FILE *fp;
  register unsigned long checksum = 0; /* The checksum mod 2^16. */
  int i ;

  for( i= 0 ; i < total_bytes ; ++i)
  {
    ROTATE_RIGHT (checksum);
    checksum += data[i] ;
    checksum &= 0xffff;	/* Keep it within bounds. */
  }


  return (int) checksum   ;
}

