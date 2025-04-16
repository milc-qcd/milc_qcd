/*********************** io_grid_ks_eigen.c *************************/
/* MIMD version 7 */
/* For reading Grid eigenpack multifile eigenvectors using MPI I/O */

/* 11/24  C. DeTar code adapted from Grid parallelIO/BinaryIO.h */

/* The offset parameter specifies the location of the
   beginning of the data in the file.  This routine reads
   only the data and evaluates the checksum */

#include "generic_ks_includes.h"
#include <qio.h>
#include <qmp.h>
#define LOOPEND
#include "../include/loopend.h"
#include "../include/openmp_defs.h"
#include <zlib.h>
#include <assert.h>

#ifdef HAVE_MPI
#include <mpi.h>
#endif
/*------------------------------------------------------------------*/
/* Convert rank to coordinates */
static void lex_coords(int coords[], const int dim, const int size[], 
		       const size_t rank)
{
  int d;
  size_t r = rank;
  
  for(d = 0; d < dim; d++){
    coords[d] = r % size[d];
    r /= size[d];
  }
}

/*------------------------------------------------------------------*/
/* Convert coordinate to linear lexicographic PE rank (inverse of
   lex_coords) */

static size_t lex_rank(const int coords[], int dim, int size[])
{
  int d;
  size_t rank = coords[dim-1];

  for(d = dim-2; d >= 0; d--){
    rank = rank * size[d] + coords[d];
  }
  return rank;
}

/* Do byte reversal on n contiguous 32-bit words */
static void
byterevn32(uint32_t w[], size_t n)
{
  uint32_t old,newv;
  size_t j;

  assert(sizeof(uint32_t) == 4);
  
  for(j=0; j<n; j++)
    {
      old = w[j];
      newv = old >> 24 & 0x000000ff;
      newv |= old >> 8 & 0x0000ff00;
      newv |= old << 8 & 0x00ff0000;
      newv |= old << 24 & 0xff000000;
      w[j] = newv;
    }
}

/*------------------------------------------------------------------*/
/* Compute Grid-version of SciDac checksum on the odd-site data as
   read */
static void
grid_checksum(su3_vector *eigVec, uint32_t *suma, uint32_t *sumb){
  
  int latdim[4] = {nx, ny, nz, nt};
  int localdim[4];

  int squaresize[4];
  const int *mc = get_logical_coordinate();
  const int *nsquares = get_logical_dimensions();
  int d;
  FORALLUPDIR(d){
   localdim[d] = latdim[d]/nsquares[d];
  }

  localdim[0] /=2; /* Grid packs odd sites by dividing the x axis by 2 */

  int i;
  uint32_t csuma = 0, csumb = 0;
  FORODDFIELDSITES_OMP(i,reduction( ^ : csuma, csumb)){
    /* Here the index i is meant to specify the order of input, not
       the storage */
    int x[4];
    int coords[4];
    int64_t global_site;
    uint32_t site_crc = 0;
    QIO_Index k;

    uint32_t * site_buf = (uint32_t *)&eigVec[i];

    /* CRC of the su3_vector on this site */
    site_crc = crc32(0,(unsigned char *)site_buf, sizeof(dsu3_vector));

    /* Get Grid "coordinates" of this site */
    lex_coords( coords, 4, localdim, i - even_sites_on_node);
    for(int d = 0; d < 4; d++){
      coords[d] += localdim[d]*mc[d];
    }

    /* The global rank  of the site */
    global_site = lex_rank( coords, 4, latdim);

    // printf("%lu %5d coords %2d %2d %2d %2d", global_site,  i - even_sites_on_node, coords[0], coords[1], coords[2], coords[3]); dumpvec(&eigVec[i]);

    /* Then fix the byte order */
    byterevn64(site_buf, sizeof(dsu3_vector)/sizeof(double));

    /* Rotate and xor the hash */
    uint64_t gsite29 = global_site % 29;
    uint64_t gsite31 = global_site % 31;
    csuma ^= site_crc<<gsite29 | site_crc>>(32-gsite29);
    csumb ^= site_crc<<gsite31 | site_crc>>(32-gsite31);

  } END_LOOP_OMP;

  *suma = csuma;  *sumb = csumb;

  /* Combine results from all nodes */
  g_xor32(suma);
  g_xor32(sumb);
  
}

/*------------------------------------------------------------------*/
void
read_grid_eigenvector_data(su3_vector *eigVec, int typesize, char* file, off_t offset,
			   uint32_t *suma, uint32_t *sumb){
  char myname[] = "read_grid_eigenvector_data";

  if(typesize != 48){
    node0_printf("%s: Grid requires double.  Typesize %d != 48\n", myname, typesize);
    terminate(1);
  }

  const int latdim[4] = {nx, ny, nz, nt};
  int grid_reduced_dim[4] = {nx/2, ny, nz, nt};
  const int *nsquares = get_logical_dimensions();
  const int *mc = get_logical_coordinate();

  int squaresize[4];
  int coord_start[4];
  int d;
  FORALLUPDIR(d){
    squaresize[d] = grid_reduced_dim[d]/nsquares[d];
    coord_start[d] = mc[d]*squaresize[d];
  }

  double dtime = -dclock();

  if ( numnodes() > 1 ) {

#ifdef HAVE_MPI
    MPI_Datatype mpiObject;
    MPI_Datatype fileArray;
    MPI_Datatype localArray;
    MPI_Datatype mpiword;
    
    MPI_Offset disp = offset;
    MPI_File fh ;
    MPI_Status status;

    int numword = 6;  /* su3_vector */
  
    mpiword = MPI_DOUBLE; /* Grid requirement for eigenvectors */

    /* su3_vector in MPI phrasing */
    int ierr;
    ierr = MPI_Type_contiguous(numword, mpiword, &mpiObject);    assert(ierr==0);
    ierr = MPI_Type_commit(&mpiObject);
  
    /*  File global array data type */
    ierr=MPI_Type_create_subarray(4, grid_reduced_dim, squaresize, coord_start,
				  MPI_ORDER_FORTRAN, mpiObject, &fileArray);    assert(ierr==0);
    ierr=MPI_Type_commit(&fileArray);    assert(ierr==0);
    
    /*  File local array data type */
    int lStart[4] = {0, 0, 0, 0};
    ierr=MPI_Type_create_subarray(4, squaresize, squaresize, lStart,
				  MPI_ORDER_FORTRAN, mpiObject, &localArray);    assert(ierr==0);
    ierr=MPI_Type_commit(&localArray);    assert(ierr==0);

    /* Open the file and set the view */
    MPI_Comm *mpicm = (MPI_Comm *)mycomm();
    ierr=MPI_File_open(*mpicm, (char *) file, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);    assert(ierr==0);
    ierr=MPI_File_set_view(fh, disp, mpiObject, fileArray, "native", MPI_INFO_NULL);    assert(ierr==0);

    /* Do the parallel I/O */
    ierr=MPI_File_read_all(fh, (void *)&eigVec[even_sites_on_node], 1, localArray, &status);    assert(ierr==0);

    /* Done */
    MPI_File_close(&fh);
    MPI_Type_free(&fileArray);
    MPI_Type_free(&localArray);

#else

    /* MPI required here for multiple nodes */
    node0_printf("%s: ERROR: Must compile with MPI\n", myname);
    terminate(1);

#endif

  } else {

    /* For one node, use POSIX I/O */

    FILE *fin = fopen(file, "rb");
    fseek(fin, offset, SEEK_SET);
    size_t data_read = fread((char *)&eigVec[even_sites_on_node], typesize, odd_sites_on_node, fin);
    assert(ferror(fin) == 0);
    assert(data_read = odd_sites_on_node);
    fclose(fin);

  }
  
  grid_checksum(eigVec, suma, sumb);
  
  dtime += dclock();
  node0_printf("Time to read one eigenvector %f\n", dtime);
  
  g_sync();
  
} /*io_grid_ks_eigen.c */
