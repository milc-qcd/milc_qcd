/*
 *  Convert from three spatial coordinates to a location
 *  in a three dimensional array.
 *
 */

#include"assert.h"

int where(int x, int y, int z ,int nx)
{
  int pt ;


  /*** range check the parameters *******/
  assert(x >= 0 && x < nx );
  assert(y >= 0 && y < nx );
  assert(z >= 0 && z < nx );

  return  x + nx*( y + nx*z)  ; 

}
