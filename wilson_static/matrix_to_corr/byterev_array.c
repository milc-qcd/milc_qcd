/*
   Byte reverse an array of integers.
   This is a minor modification of the 
   code byterev.c, written by  Doug Toussaint.

   Usage:
       int buff[dim-1] ; byte_rev_array(buf, dim)) ;

       Real buff[dim-1] ; byte_rev_array((int*)buf, dim)) ;

     I have not tried using "double"
*/


void byte_rev_array(int buf[], int words)
{
  register int i ;
  register int old,new;

  
  for(i=0;i<words;i++)
  {
    old = buf[i];
    new = old >> 24 & 0x000000ff;
    new |= old >> 8 & 0x0000ff00;
    new |= old << 8 & 0x00ff0000;
    new |= old << 24 & 0xff000000;
    buf[i] = new;
  }


}
