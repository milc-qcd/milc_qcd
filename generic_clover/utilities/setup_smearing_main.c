/*
 *  Set up the smearing functions that are
 *  required for the variational static calculation.
 *
 *  Currently this code produces a number of local
 *  smearing functions.
 *
 *
 *
 *
 */

#include<stdio.h>
#include<math.h>
#include<string.h>

/*** function declerations used in this code ****/
void titles (int nx, int ny, int nz, int, int );
void write_smear_wave(Real *data, int nx, int dim, char filename[]);

void fft_smearing_func(Real *array, int dim);

void create_local_oper(Real *src, int dim, char filename[]);
void create_wall_oper(Real *src, int dim, char filename[]);
void create_random_oper(Real *src, int *seed, int dim, char filename[]);
void create_expon_oper(Real *src, int nx, char filename[]);
void create_hydrogen(Real *src, int quantum_n, Real (*func)(Real),
int nx, char filename[]);



Real two_s_hydrogen_poly(Real x) ;
Real three_s_hydrogen_poly(Real x) ;
Real four_s_hydrogen_poly(Real x) ;
Real five_s_hydrogen_poly(Real x) ;


void create_2S_polyexpon_oper(Real *src, int nx, char filename[]);
void create_3S_polyexpon_oper(Real *src, int nx, char filename[]);

void create_gaussian_oper(Real *src, int nx, char filename[]) ; 
void create_milc_gaussian_oper(Real *src, int nx, char filename[]) ;

void create_squeeze_expon_oper(Real *src, int nx, char filename[]) ;
void create_2S_polyexpon_oper_x_squeeze(Real *src, int nx, char filename[]) ;
void create_3S_polyexpon_oper_x_squeeze(Real *src, int nx, char filename[]) ;

/** types of smearing functions that cen be created ****/
enum smearing_func { LOCAL , WALL , RANDOM , EXP , H_2S, H_3S, H_4S, H_5S ,
		       POLYEXP_2S , POLYEXP_3S  , GAUSS , GAUSS_MILC , EXP_SQUEEZE , 
		       POLYEXP_2S_SQUEEZE , POLYEXP_3S_SQUEEZE } ;

int main()
{
  int nx;
  int ny;
  int nz;
  int nxyz;

/**     const int nosmear = 3 ;   ***/
  const int nosmear = 1 ; 

  Real *wave_func ;
  int seed = -1020 ;
  int ismear ;
  char filename[80] ;
/**  int smear_type[] = { LOCAL, WALL , LOCAL, LOCAL };   ***/
/**     int smear_type[] = { RANDOM, EXP , RANDOM, WALL };   ***/ 
/**    int smear_type[] = { EXP , POLYEXP_2S  , POLYEXP_3S };  ****/
/**  int smear_type[] = { EXP , RANDOM };   ***/
/**  int smear_type[] = { POLYEXP_3S };    ****/
/** int smear_type[] = { GAUSS_MILC };  ***/
/**  int smear_type[] = { EXP_SQUEEZE , POLYEXP_2S_SQUEEZE , POLYEXP_3S_SQUEEZE } ;  ***/

/**  int smear_type[] = { GAUSS_MILC } ; **/
  int smear_type[] = { EXP } ;		       


  /********----------*********----------**************/

  /* Prompt for lattice size */
  ny = -1;
  while(ny <= 0){
    printf("Enter lattice dimension nx: ");
    if(scanf("%d",&nx) == 1){
      ny = nz = nx;
      nxyz = nx*ny*nz;
    }
  }

  /** titles ****/
  titles(nx, ny, nz, seed, nosmear);

  /** Reserve the required memory **********/
  if( (wave_func = (Real *)calloc( (size_t) nxyz, sizeof(Real) )  ) == NULL )
  {
    printf("ERROR: could not reserve  space for wave function");
    exit(1);
  }


  /** write to disk the local operator ****/
  create_local_oper(wave_func , nxyz, filename);  
  write_smear_wave(wave_func, nx,nxyz, filename);


  /*** loop over the smearing functions ****/
  for(ismear = 0 ; ismear < nosmear ; ++ismear )
  {

    /*** create the required smearing function *****/
    switch (smear_type[ismear] )
    {
    case LOCAL:
      create_local_oper(wave_func , nxyz, filename);
      break ;
    case WALL:
      create_wall_oper(wave_func , nxyz, filename);
      break ;
    case RANDOM:
      create_random_oper(wave_func , &seed, nxyz, filename);
      break ;
    case EXP:
      create_expon_oper(wave_func ,  nx, filename);
      break ;
    case H_2S:
      create_hydrogen(wave_func ,  2,&two_s_hydrogen_poly ,nx, filename);
      break ;
    case H_3S:
      create_hydrogen(wave_func ,  3,&three_s_hydrogen_poly ,nx, filename);
      break ;
    case H_4S:
      create_hydrogen(wave_func , 4,&four_s_hydrogen_poly ,nx, filename);
      break ;
    case H_5S:
      create_hydrogen(wave_func ,  5,&five_s_hydrogen_poly ,nx, filename);
      break ;
    case POLYEXP_2S:
      create_2S_polyexpon_oper(wave_func, nx, filename) ;
      break ;
    case POLYEXP_3S:
      create_3S_polyexpon_oper(wave_func, nx, filename) ;
      break ;
    case  GAUSS :
      create_gaussian_oper( wave_func  , nx, filename) ;
      break ;
    case  GAUSS_MILC :
      create_milc_gaussian_oper(wave_func  , nx, filename) ;
      break ;
    case EXP_SQUEEZE :
      create_squeeze_expon_oper(wave_func  , nx, filename) ;
      break ;
    case POLYEXP_2S_SQUEEZE : 
      create_2S_polyexpon_oper_x_squeeze(wave_func  , nx, filename) ;
      break ;
    case POLYEXP_3S_SQUEEZE : 
      create_3S_polyexpon_oper_x_squeeze(wave_func  , nx, filename) ;
      break ;
    default :
      printf("ERROR: Bad smearing function flag %d\n",ismear);
      exit(10);
  }

    write_smear_wave(wave_func, nx,nxyz, filename);
    printf("\n\n\n");

  } /*** end the loop over the smearing functions ***/


  /** free up the required memory ****/
  free(wave_func);


  return 0 ;
}














