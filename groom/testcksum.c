#include <stdio.h>
#include <sys/types.h>
#include <errno.h>
#include <fcntl.h>
#include <string.h>
#include <unistd.h>

#ifdef SHORT32
typedef short type32;
typedef unsigned short u_type32;
#else
typedef int type32;
typedef unsigned int u_type32;
#endif

/******************************** byterevn.c ***************************/
/* MIMD version 6 */

/* For doing byte reversal on 32-bit words */

void byterevn(type32 w[], int n)
{
  register type32 old,newv;
  int j;
  
  for(j=0; j<n; j++)
    {
      old = w[j];
      newv = old >> 24 & 0x000000ff;
      newv |= old >> 8 & 0x0000ff00;
      newv |= old << 8 & 0x00ff0000;
      newv |= old << 24 & 0xff000000;
      w[j] = newv;
    }
} /* byterevn */

/************************************************************************/
int main(int argc, char* argv[]){
  int fd,fd2,oflg,oflg2;
  char filename[] = "testin";
  char filename2[] = "testout";
  u_type32 cksum;
  int rank;
  int i,nread;
  Real word;
  u_type32 *val;
  int byterevflag;

  oflg = O_RDONLY;
  oflg2 = O_WRONLY | O_CREAT;

  val = (u_type32 *)&word;

  if(argc < 3 ||
     (sscanf(argv[1],"%d",&nread) != 1) ||
     (sscanf(argv[2],"%d",&byterevflag) != 1))
    {
      printf("Usage %s <nwords> <byterevflg> (0 no byterev 1 byterev)\n",argv[0]);
      return 1;
    }

  if( (fd = open(filename, oflg, 00644)) == -1)
    {
      printf("error %d opening %s\n",errno,filename);
      return 1;
    }

  if( (fd2 = open(filename2, oflg2, 00644)) == -1)
    {
      printf("error %d opening %s\n",errno,filename2);
      return 1;
    }

  cksum = 0;
  rank = 0;
  
  for(i = 0; i < nread; i++)
    {
      if(read( fd, &word, 4 ) != 4)
	{
	  printf("error %d reading %s\n",errno,filename);
	  return 1;
	}

      if(i >= 64)
	{
	  if(write(fd2, &word, 4) != 4)
	    {
	      printf("error writing %s\n",filename2);
	      return 1;
	    }
	}

      if(byterevflag)
	{
	  byterevn((type32 *)&word,1);
	}

      printf("word = %f\n",word);

      cksum ^= (*val)<<rank | (*val)>>(32-rank);
      rank++;
    }

  printf("Checksum after %d values from %s is %x\n",
	 nread,filename,cksum);

  return 0;
}
