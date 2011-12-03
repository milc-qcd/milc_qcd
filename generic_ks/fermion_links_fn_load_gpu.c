/**************** fermion_links_fn_load_milc_gpu.c **********************/
/* MILC Version 7 */

/* Originally Guochun Shi and Steve Gottlieb 2010 */

/* Entry points 

   load_fatlinks_gpu

*/

#include <test_util.h>
#include <blas_reference.h>
#include <quda.h>
#include <gauge_quda.h>
#include <llfat_quda.h>
#include <string.h>

#include "generic_ks_includes.h"
#include "../include/prefetch.h"
#include "../include/info.h"
#define FETCH_UP 1

#define LOOPEND
#include "../include/loopend.h"

extern int numnodes(void);
extern void initQuda(int dev);

//extern void comm_set_gridsize(int x, int y, int z, int t); // These declarations ought to be included
//extern void initDslashConstants(FullGauge, int, int); // Do I really need this?

static FullGauge cudaSiteLink;
static FullGauge cudaFatLink;
static FullStaple cudaStaple1;
static FullStaple cudaStaple;
static QudaGaugeParam gaugeParam;



static void comm_begin(const int gridsize[]){
  int i, total_num_processes = 1;

  for(i=0; i<4; ++i) total_num_processes *= gridsize[i];

  int size = numnodes();
  if( total_num_processes != size)
    errorQuda("Number of processes %d must match requested MPI volume %d", size, total_num_processes);

  // now initialize communication
  comm_set_gridsize(gridsize[0], gridsize[1], gridsize[2], gridsize[3]);
  comm_init();

  return;
}


// set the dimension of the sublattice on each cpu
static void set_local_dims(const int gridsize[], int dim[]){

  dim[0] = nx/gridsize[0];
  dim[1] = ny/gridsize[1];
  dim[2] = nz/gridsize[2];
  dim[3] = nt/gridsize[3];

  return;
}

// New function!
static void setLocalDims(int locdim[], const int latdim[], const int gridsize[]){

  int dir;
  for(dir=0; dir<4; ++dir){
    if(latdim[dir]%gridsize[dir] != 0){
	    errorQuda("latdim[%d] is not a multiple of gridsize[%d]\n", dir, dir);
	  }
    locdim[dir] = latdim[dir]/gridsize[dir];
  }
  return;
}


// work out the volumes of different slices
static void set_slice_volumes(const int dim[], int vol[]){
  int i, j;
  for(i=0; i<4; ++i){ 
    vol[i] = 1; 
    for(j=0; j<4; ++j){
      if(i != j){ vol[i] *= dim[j]; }
    }
  }
  return;
}



static int get_volume(const int dim[]){
  return dim[0]*dim[1]*dim[2]*dim[3];
}


// returns 1 if we should exchange even and odd sites
static int should_even_odd_exchange(const int dim[]){
  int i;
  int even_odd_exchange = 0;

  for(i=0; i<4; ++i){
    if(dim[i] % 2 == 1 && get_logical_coordinate()[i] % 2 == 1){
      even_odd_exchange = 1 - even_odd_exchange;
    }
  }
  return even_odd_exchange;
}



// Okay, I think this should do what I think
// At least this seems to be the format 
static void exchange_even_odd_links(const su3_matrix *src, int vol, su3_matrix *dst){

  const int half_vol = vol/2;
  const int odd_start_index = 4*half_vol;

  int i;
  for(i=0; i<odd_start_index; ++i){
    dst[i] = src[odd_start_index+i];
  }

  for(i=0; i<odd_start_index; ++i){
    dst[odd_start_index+i] = src[i];
  }
  return;
}






static void exchange_even_odd_links_in_place(su3_matrix *link, int vol){

  const int half_vol = vol/2;
  const int half_field_size = 4*half_vol*sizeof(su3_matrix);

  void *temp = (void*)malloc(half_field_size);

  memcpy(temp, link, half_vol*4*sizeof(su3_matrix));
  memcpy(link, ((char*)link)+half_field_size, half_field_size);
  memcpy(((char*)link)+half_field_size, temp, half_field_size);

  free(temp);
  return;
}




// put the field into a different format
static void convert_links_to_dir_major_format(su3_matrix *dst, const su3_matrix *src, int vol){
  int i, dir;
  for(dir=0; dir<4; ++dir){
    for(i=0; i<vol; ++i){
      dst[dir*vol + i] = src[i*4 + dir];
    }
  }
  return;
}


// inverse of the previous function
static void convert_links_to_site_major_format(su3_matrix *dst, const su3_matrix *src, int vol){
  int i, dir;
  for(i=0; i<vol; ++i){
    for(dir=0; dir<4; ++dir){
      dst[4*i + dir] = src[dir*vol + i];
    }
  }
  return;
}


// the preceeding functions can be subsumed into a single function call
//static void transpose(su3_matrix *dst, const su3_matrix *src, int nrows, int ncols)





// copies from a one dimensional array to two dimensional format
static void copy_links_from_1d_array_to_2d_array(su3_matrix **dst, const su3_matrix src[], int num_rows, int num_cols){
  int i, j;
  for(i=0; i<num_rows; ++i){
    for(j=0; j<num_cols; ++j){
      dst[i][j] = src[i*num_cols + j];
    }  
  }
  return;
}



static int get_max_2d_volume(int dim[]){
  int max_vol_2d = 0;
  int i, j, vol_2d;
  for(i=0; i<4; ++i){
    for(j=i+1; j<4; ++j){
		  vol_2d = dim[i]*dim[j];
      if(vol_2d > max_vol_2d){
        max_vol_2d = vol_2d;
      }
    }
  } // end loop over directions
  return max_vol_2d;
}



static size_t get_real_gauge_variable_size(const QudaGaugeParam *gauge_param){
  return (gaugeParam.cpu_prec == QUDA_DOUBLE_PRECISION) ? sizeof(double) : sizeof(float);
}



static void get_QUDA_links_from_field(su3_matrix **quda_link, int even_odd_exchange, const QudaGaugeParam *gauge_param, const su3_matrix *link){

  const int V = get_volume(gauge_param->X);
  const int field_size = 4*V*sizeof(su3_matrix);

  node0_printf("even_odd_exchange = %d\n", even_odd_exchange);
  node0_printf("V = %d\n", V);
  su3_matrix *site_major_link;


  if(even_odd_exchange){
    su3_matrix *templink = link;
    site_major_link = (su3_matrix*)malloc(field_size);
    exchange_even_odd_links(templink, V, site_major_link); 
  }else{
    site_major_link = link;
  }

  su3_matrix *dir_major_link  = (su3_matrix*)malloc(field_size); // need to free dir_major_link_field below
  convert_links_to_dir_major_format(dir_major_link, site_major_link, V);

  const int num_rows = 4;
  const int num_cols = V;
  copy_links_from_1d_array_to_2d_array(quda_link, dir_major_link, num_rows, num_cols);

  free(dir_major_link);
  if(even_odd_exchange){
    free(site_major_link);
  }
  return;
}




// I need some kind of function that initialises quda from inside MILC


void load_fatlinks_gpu(info_t *info, su3_matrix *fat, 
		       ks_component_paths *p, su3_matrix *link)
{
  const int *nsquares = get_logical_dimensions();
  int i; 
  int V, Vh;
  int Vs[4], Vsh[4];

  int fulldim[4] = {nx, ny, nz, nt};


  comm_begin(nsquares); // set the grid size, initialize communication
  
  set_local_dims(nsquares, gaugeParam.X); 
  V = get_volume(gaugeParam.X); 

  set_slice_volumes(gaugeParam.X, Vs);
  for(i=0; i<4; ++i){ Vsh[i] = Vs[i]/2; } // half slice volumes

  // This is essential!
#ifdef GPU_USE_12_RECON
  QudaReconstructType link_recon = QUDA_RECONSTRUCT_12;
#else
  QudaReconstructType link_recon = QUDA_RECONSTRUCT_NO;
#endif

#if (PRECISION==1)
  QudaPrecision  link_prec = QUDA_SINGLE_PRECISION;
  QudaPrecision  cpu_prec = QUDA_SINGLE_PRECISION;
#else
  QudaPrecision  link_prec = QUDA_DOUBLE_PRECISION;
  QudaPrecision  cpu_prec = QUDA_DOUBLE_PRECISION;
#endif

#ifdef GPU_FATLINK_SINGLE_PRECISION
  link_prec = QUDA_SINGLE_PRECISION;
#endif


  gaugeParam.cpu_prec = cpu_prec;
  gaugeParam.cuda_prec = link_prec;
  gaugeParam.reconstruct = link_recon;
  // Precision and reconstruction options set



  const size_t gSize = get_real_gauge_variable_size(&gaugeParam); 

  
  // allocate memory for qudalinks - This should be in it's own function
  su3_matrix *qudalink[4];
  int dir;
  for(dir=0; dir<4; ++dir){
    qudalink[dir] = (su3_matrix*)malloc(V*sizeof(su3_matrix));
    if(qudalink[dir] == NULL){
      fprintf(stderr, "Error: malloc failed for qudalink[%d]\n", dir);
      exit(1);
    }
  }


  int even_odd_exchange = should_even_odd_exchange(gaugeParam.X);  
  get_QUDA_links_from_field(qudalink, even_odd_exchange, &gaugeParam, link);

  int device=0;
  initQuda(device);

  const int Vh_2d_max = get_max_2d_volume(gaugeParam.X)/2;

  node0_printf("Vh_2d_max = %d\n", Vh_2d_max);


  // Need to change this, I think
  gaugeParam.site_ga_pad = gaugeParam.ga_pad = 3*(Vsh[0]+Vsh[1]+Vsh[2]+Vsh[3]) + 4*Vh_2d_max;
  gaugeParam.reconstruct = QUDA_RECONSTRUCT_NO;

  createLinkQuda(&cudaSiteLink, &gaugeParam);
  loadLinkToGPU(cudaSiteLink, (void**)qudalink, &gaugeParam);

  gaugeParam.staple_pad = 3*(Vsh[0] + Vsh[1] + Vsh[2]+ Vsh[3]);
  createStapleQuda(&cudaStaple, &gaugeParam);
  createStapleQuda(&cudaStaple1, &gaugeParam);

  gaugeParam.llfat_ga_pad = gaugeParam.ga_pad = Vsh[3]; // 3 => time
  gaugeParam.reconstruct = QUDA_RECONSTRUCT_NO;
  createLinkQuda(&cudaFatLink, &gaugeParam);


  double c[6];
  for(i=0; i<6; ++i){
   c[i] = p->act_path_coeff[i];
   node0_printf("c[%d] = %lf\n", i, c[i]);
  }

  
  // Perform the link fattening
  initCommonConstants(cudaSiteLink);
  llfat_init_cuda(&gaugeParam);

  llfat_cuda(cudaFatLink, cudaSiteLink, cudaStaple, cudaStaple1, &gaugeParam, c);
  //llfat_cuda(cudaFatLink, cudaSiteLink, cudaStaple, cudaStaple1, &gaugeParam, p->act_path_coeff);

  gaugeParam.ga_pad = gaugeParam.llfat_ga_pad;
  gaugeParam.reconstruct = QUDA_RECONSTRUCT_NO;
  storeLinkToCPU((void*)fat, &cudaFatLink, &gaugeParam);
 // loadLinkToGPU and storLinkToCPU have two different formats 
 // I need to address this first of all 
 // loadLinktoGPU -> load2DLinkToGPU, load1DLinkToGPU
 // storeLinkToGPU -> storeLinkToCPU 1D, storeLinkToCPUTD


  if(even_odd_exchange){
    exchange_even_odd_links_in_place(fat, V); 
  }

  // Free up memory:
  // This should be in its own function
  for(dir=0; dir<4; ++dir){
    free(qudalink[dir]);
    qudalink[dir] = NULL;
  }

  freeLinkQuda(&cudaSiteLink);
  freeLinkQuda(&cudaFatLink);
  freeStapleQuda(&cudaStaple);
  freeStapleQuda(&cudaStaple1);
  return;
}
