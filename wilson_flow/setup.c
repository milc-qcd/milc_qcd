/************************ setup.c ****************************/
/* Handles lattice and parameter setup                       */

#define IF_OK if(status==0)

/* definitions, files, and prototypes */
#include "wilson_flow_includes.h"
#include <string.h>
#include <unistd.h>
extern int gethostname (char *__name, size_t __len); // Should get this from unistd.h

/* Each node has a params structure for passing simulation parameters */
#include "params.h"
params par_buf;

int initial_set();

int
setup()
{
  int prompt;

  /* print banner, get initial parameters */
  prompt = initial_set();

  /* Initialize the layout functions, which decide where sites live */
  setup_layout();
  this_node = mynode();

  /* allocate space for lattice, set up coordinate fields */
  make_lattice();
  node0_printf("Made lattice\n"); fflush(stdout);

  /* set up neighbor pointers and comlink structures */
  make_nn_gathers();
  node0_printf("Made nn gathers\n"); fflush(stdout);

  /* initialize Runge-Kutta integrator for the flow */
  initialize_integrator();

  node0_printf("Finished setup\n"); fflush(stdout);
  return prompt;
}


/* SETUP ROUTINES */
int
initial_set()
{
  int prompt, status;

  if(mynode()==0){
    /* print banner */
    printf("Wilson/Symanzik Flow application\n");
    printf("MIMD version 7\n");
    printf("Machine = %s, with %d nodes\n", machine_type(), numnodes());
    gethostname(hostname, 128);
    printf("Host(0) = %s\n",hostname);
    printf("Username = %s\n", getenv("USER"));
    time_stamp("start");

    /* Print list of options selected */
    if(mynode()==0)printf("Options selected...\n");
    show_generic_opts();

    /* Read prompt type and lattice dimensions */
    status=get_prompt(stdin, &prompt);
    IF_OK status += get_i(stdin, prompt,"nx", &par_buf.nx );
    IF_OK status += get_i(stdin, prompt,"ny", &par_buf.ny );
    IF_OK status += get_i(stdin, prompt,"nz", &par_buf.nz );
    IF_OK status += get_i(stdin, prompt,"nt", &par_buf.nt );
#ifdef ANISOTROPY
    IF_OK status += get_f(stdin, prompt,"anisotropy", &par_buf.ani );
#endif

    if(status>0)
      par_buf.stopflag=1;
    else
      par_buf.stopflag=0;
  } /* end if(mynode()==0) */

  /* Node 0 broadcasts parameter buffer to all other nodes */
  broadcast_bytes((char *)&par_buf, sizeof(par_buf));

  if( par_buf.stopflag != 0 )
    normal_exit(0);

  /* Update global variables with parameters */
  nx=par_buf.nx;
  ny=par_buf.ny;
  nz=par_buf.nz;
  nt=par_buf.nt;
#ifdef ANISOTROPY
  ani=par_buf.ani;
#endif

  number_of_nodes = numnodes();
  volume=(size_t)nx*ny*nz*nt;

  return prompt;
}

/* Read in configuration specific parameters */
int
readin(int prompt)
{
  int status;

  /* On node zero, read parameters and send to all other nodes */
  if(this_node==0) {

    printf("\n\n");
    status=0;

    /* Identify the starting configuration */
    IF_OK status += ask_starting_lattice(stdin,  prompt, &(par_buf.startflag),
                                         par_buf.startfile);

    /* Get flow parameters */
    IF_OK {
      /* Get flow description */
      if (prompt==1)
        printf("enter 'wilson' or 'symanzik'\n");
      scanf("%s", par_buf.flow_description);
      printf("%s\n", par_buf.flow_description);

      /* Check flow description and set staple flag */
      if( strcmp("wilson", par_buf.flow_description) == 0 ) {
        par_buf.stapleflag = WILSON;
        printf("set staple to wilson\n");
      }
      else if( strcmp("symanzik", par_buf.flow_description) == 0 ) {
        par_buf.stapleflag = SYMANZIK;
        printf("set staple to symanzik\n");
      }
      else if( strcmp("zeuthen", par_buf.flow_description) == 0 ) {
        par_buf.stapleflag = ZEUTHEN;
        printf("set staple to zeuthen\n");
      }
      else {
        printf("Error: flow_description %s is invalid\n",
               par_buf.flow_description);
        status++;
      }
    } /*end: flow_description IF_OK */

    IF_OK status += get_i(stdin, prompt, "exp_order", &par_buf.exp_order);
    IF_OK status += get_f(stdin, prompt, "stepsize", &par_buf.stepsize);
#if GF_INTEGRATOR==INTEGRATOR_ADAPT_LUSCHER || \
    GF_INTEGRATOR==INTEGRATOR_ADAPT_CF3 || \
    GF_INTEGRATOR==INTEGRATOR_ADAPT_BS
    IF_OK status += get_f(stdin, prompt, "local_tol", &par_buf.local_tol);
#endif
    IF_OK status += get_f(stdin, prompt, "stoptime", &par_buf.stoptime);

    /* Determine what to do with the final configuration */
    IF_OK status += ask_ending_lattice(stdin,  prompt, &(par_buf.saveflag),
                                       par_buf.savefile );
    IF_OK status += ask_ildg_LFN(stdin, prompt, par_buf.saveflag,
                                 par_buf.stringLFN );

    if( status > 0)
      par_buf.stopflag=1;
    else
      par_buf.stopflag=0;

  } /* end if(this_node==0) */

  /* Node 0 broadcasts parameter buffer to all other nodes */
  broadcast_bytes((char *)&par_buf, sizeof(par_buf));

  if( par_buf.stopflag != 0 )
    return par_buf.stopflag;

  /* Update global variables with parameter buffer */
  startflag = par_buf.startflag;
  strcpy(startfile, par_buf.startfile);

  strcpy(flow_description, par_buf.flow_description);
  stapleflag = par_buf.stapleflag;
  exp_order = par_buf.exp_order;
  stepsize = par_buf.stepsize;
  stoptime = par_buf.stoptime;
#if GF_INTEGRATOR==INTEGRATOR_ADAPT_LUSCHER || \
    GF_INTEGRATOR==INTEGRATOR_ADAPT_CF3 || \
    GF_INTEGRATOR==INTEGRATOR_ADAPT_BS
  local_tol = par_buf.local_tol;
#endif

  saveflag = par_buf.saveflag;
  strcpy(savefile, par_buf.savefile);
  strcpy(stringLFN, par_buf.stringLFN);

  /* Load configuration (no phases in this program) */
  if( par_buf.startflag != CONTINUE )
    startlat_p = reload_lattice( par_buf.startflag, par_buf.startfile );

  return 0;
}

/* Initialize integrator depending on compile-time flags */
void
initialize_integrator()
{

#if GF_INTEGRATOR==INTEGRATOR_LUSCHER
  N_stages = 3;
  // Luscher coefficients in proper 2N-storage format
  A_2N[0] = 0;
  A_2N[1] = -17/32.;
  A_2N[2] = -32/27.;
  B_2N[0] = 1/4.;
  B_2N[1] = 8/9.;
  B_2N[2] = 3/4.;
  node0_printf("Integrator = INTEGRATOR_LUSCHER\n");
#elif GF_INTEGRATOR==INTEGRATOR_CK
  N_stages = 5;
  // Carpenter, Kennedy coefficients,
  // Technical Memorandum NASA-TM-109112
  A_2N[0] = 0;
  A_2N[1] = -567301805773/1357537059087.;
  A_2N[2] = -2404267990393/2016746695238.;
  A_2N[3] = -3550918686646/2091501179385.;
  A_2N[4] = -1275806237668/842570457699.;
  B_2N[0] = 1432997174477/9575080441755.;
  B_2N[1] = 5161836677717/13612068292357.;
  B_2N[2] = 1720146321549/2090206949498.;
  B_2N[3] = 3134564353537/4481467310338.;
  B_2N[4] = 2277821191437/14882151754819.;
  node0_printf("Integrator = INTEGRATOR_CK\n");
#elif GF_INTEGRATOR==INTEGRATOR_BBB
  N_stages = 6;
  // Berland, Bogey, Bailly coefficients,
  // Computers and Fluids 35 (2006) 1459-1463
  A_2N[0] = 0;
  A_2N[1] = -0.737101392796;
  A_2N[2] = -1.634740794341;
  A_2N[3] = -0.744739003780;
  A_2N[4] = -1.469897351522;
  A_2N[5] = -2.813971388035;
  B_2N[0] = 0.032918605146;
  B_2N[1] = 0.823256998200;
  B_2N[2] = 0.381530948900;
  B_2N[3] = 0.200092213184;
  B_2N[4] = 1.718581042715;
  B_2N[5] = 0.27;
  node0_printf("Integrator = INTEGRATOR_BBB\n");
#elif GF_INTEGRATOR==INTEGRATOR_CF3
  N_stages = 3;
  // Commutator-free 2N-storage with arbitrary coefficients

#ifdef READ_CF3_FROM_FILE
  FILE *fp;
  int st = 0;
  fp = fopen( "cf3_coeff.dat", "rt" );
  if( fp==NULL ) {
    node0_printf( "ERROR: cf3_coeff.dat file cannot be read\n" );
    terminate(1);
  }
#if (MILC_PRECISION==1)
  st += fscanf( fp, "%e", &(A_2N[0]) );
  st += fscanf( fp, "%e", &(A_2N[1]) );
  st += fscanf( fp, "%e", &(A_2N[2]) );
  st += fscanf( fp, "%e", &(B_2N[0]) );
  st += fscanf( fp, "%e", &(B_2N[1]) );
  st += fscanf( fp, "%e", &(B_2N[2]) );
#else
  st += fscanf( fp, "%le", &(A_2N[0]) );
  st += fscanf( fp, "%le", &(A_2N[1]) );
  st += fscanf( fp, "%le", &(A_2N[2]) );
  st += fscanf( fp, "%le", &(B_2N[0]) );
  st += fscanf( fp, "%le", &(B_2N[1]) );
  st += fscanf( fp, "%le", &(B_2N[2]) );
#endif
  fclose( fp );
  node0_printf( "%.16g\n", A_2N[0] );
  node0_printf( "%.16g\n", A_2N[1] );
  node0_printf( "%.16g\n", A_2N[2] );
  node0_printf( "%.16g\n", B_2N[0] );
  node0_printf( "%.16g\n", B_2N[1] );
  node0_printf( "%.16g\n", B_2N[2] );
#else
  // Williamson scheme 7
  A_2N[0] = 0;
  A_2N[1] = -5/9.;
  A_2N[2] = -153/128.;
  B_2N[0] = 1/3.;
  B_2N[1] = 15/16.;
  B_2N[2] = 8/15.;
  // optimized with Ralston procedure
/*A_2N[0] = 0;
  A_2N[1] = -0.637694471842202;
  A_2N[2] = -1.306647717737108;
  B_2N[0] = 0.457379997569388;
  B_2N[1] = 0.925296410920922;
  B_2N[2] = 0.393813594675071;*/
#endif
  node0_printf("Integrator = INTEGRATOR_CF3\n");
#elif GF_INTEGRATOR==INTEGRATOR_RKMK3
  N_stages = 3;
  p_order = 3;
  // Ralston coefficients
  // NOTE: indexing is shifted by -1 in C compared to usual Butcher tables
  a_RK[1][0] = 1/2.;
  a_RK[2][0] = 0; a_RK[2][1] = 3/4.;
  b_RK[0] = 2/9.; b_RK[1] = 1/3.; b_RK[2] = 4/9.;
  node0_printf("Integrator = INTEGRATOR_RKMK3\n");
#elif GF_INTEGRATOR==INTEGRATOR_RKMK4
  N_stages = 4;
  p_order = 4;
  // NOTE: indexing is shifted by -1 in C compared to usual Butcher tables
  a_RK[1][0] =  1/3.;
  a_RK[2][0] = -1/3.;
  a_RK[2][1] =  1;
  a_RK[3][0] =  1;
  a_RK[3][1] = -1;
  a_RK[3][2] =  1;
  b_RK[0] = 1/8.;
  b_RK[1] = 3/8.;
  b_RK[2] = 3/8.;
  b_RK[3] = 1/8.;
/* this set of coefficients gives higher truncation error
  a_RK[1][0] = 1/2.;
  a_RK[2][0] = 0; a_RK[2][1] = 1/2.;
  a_RK[3][0] = 0; a_RK[3][1] = 0; a_RK[3][2] = 1;
  b_RK[0] = 1/6.; b_RK[1] = 1/3.; b_RK[2] = 1/3.; b_RK[3] = 1/6.; */
  node0_printf("Integrator = INTEGRATOR_RKMK4\n");
#elif GF_INTEGRATOR==INTEGRATOR_RKMK5
  N_stages = 6;
  p_order = 5;
  // Butcher coefficients
  // NOTE: indexing is shifted by -1 in C compared to usual Butcher tables
  a_RK[1][0] = 1/4.;
  a_RK[2][0] = 1/8.; a_RK[2][1] = 1/8.;
  a_RK[3][0] = 0; a_RK[3][1] = -1/2.; a_RK[3][2] = 1;
  a_RK[4][0] = 3/16.; a_RK[4][1] = 0; a_RK[4][2] = 0; a_RK[4][3] = 9/16.;
  a_RK[5][0] = -3/7.; a_RK[5][1] = 2/7.; a_RK[5][2] = 12/7.;
  a_RK[5][3] = -12/7.; a_RK[5][4] = 8/7.;
  b_RK[0] = 7/90.; b_RK[1] = 0; b_RK[2] = 32/90.; b_RK[3] = 12/90.;
  b_RK[4] = 32/90.; b_RK[5] = 7/90.;
  // Dormand-Prince coefficients
  // NOTE: indexing is shifted by -1 in C compared to usual Butcher tables
/*  a_RK[1][0] = 1/5.;
  a_RK[2][0] = 3/40.; a_RK[2][1] = 9/40.;
  a_RK[3][0] = 44/45.; a_RK[3][1] = -56/15.; a_RK[3][2] = 32/9.;
  a_RK[4][0] = 19372/6561.; a_RK[4][1] = -25360/2187.;
  a_RK[4][2] = 64448/6561.; a_RK[4][3] = -212/729.;
  a_RK[5][0] = 9017/3168.; a_RK[5][1] = -355/33.; a_RK[5][2] = 46732/5247.;
  a_RK[5][3] = 49/176.; a_RK[5][4] = -5103/18656.;
  b_RK[0] = 35/384.; b_RK[1] = 0; b_RK[2] = 500/1113.; b_RK[3] = 125/192.;
  b_RK[4] = -2187/6784.; b_RK[5] = 11/84.;*/
  node0_printf("Integrator = INTEGRATOR_RKMK5\n");
#elif GF_INTEGRATOR==INTEGRATOR_RKMK8
  N_stages = 13;
  p_order = 8;
  // Dormand-Prince coefficients
  // NOTE: indexing is shifted by -1 in C compared to usual Butcher tables
  a_RK[1][0] = 1/18.;
  a_RK[2][0] = 1/48.; a_RK[2][1] = 1/16.;
  a_RK[3][0] = 1/32.; a_RK[3][1] = 0; a_RK[3][2] = 3/32.;
  a_RK[4][0] = 5/16.; a_RK[4][1] = 0; a_RK[4][2] = -75/64.; a_RK[4][3] = 75/64.;
  a_RK[5][0] = 3/80.; a_RK[5][1] = 0; a_RK[5][2] = 0;
  a_RK[5][3] = 3/16.; a_RK[5][4] = 3/20.;
  a_RK[6][0] = 29443841/614563906.; a_RK[6][1] = 0; a_RK[6][2] = 0;
  a_RK[6][3] = 77736538/692538347.; a_RK[6][4] = -28693883/1125000000.;
  a_RK[6][5] = 23124283/1800000000.;
  a_RK[7][0] = 16016141/946692911.; a_RK[7][1] = 0; a_RK[7][2] = 0;
  a_RK[7][3] = 61564180/158732637.; a_RK[7][4] = 22789713/633445777.;
  a_RK[7][5] = 545815736/2771057229.; a_RK[7][6] = -180193667/1043307555.;
  a_RK[8][0] = 39632708/573591083.; a_RK[8][1] = 0; a_RK[8][2] = 0;
  a_RK[8][3] = -433636366/683701615.; a_RK[8][4] = -421739975/2616292301.;
  a_RK[8][5] = 100302831/723423059.; a_RK[8][6] = 790204164/839813087.;
  a_RK[8][7] = 800635310/3783071287.;
  a_RK[9][0] = 246121993/1340847787.; a_RK[9][1] = 0; a_RK[9][2] = 0;
  a_RK[9][3] = -37695042795/15268766246.; a_RK[9][4] = -309121744/1061227803.;
  a_RK[9][5] = -12992083/490766935.; a_RK[9][6] = 6005943493/2108947869.;
  a_RK[9][7] = 393006217/1396673457.; a_RK[9][8] = 123872331/1001029789.;
  a_RK[10][0] = -1028468189/846180014.; a_RK[10][1] = 0; a_RK[10][2] = 0;
  a_RK[10][3] = 8478235783/508512852.; a_RK[10][4] = 1311729495/1432422823.;
  a_RK[10][5] = -10304129995/1701304382.; a_RK[10][6] = -48777925059/3047939560.;
  a_RK[10][7] = 15336726248/1032824649.; a_RK[10][8] = -45442868181/3398467696.;
  a_RK[10][9] = 3065993473/597172653.;
  a_RK[11][0] = 185892177/718116043.; a_RK[11][1] = 0; a_RK[11][2] = 0;
  a_RK[11][3] = -3185094517/667107341.; a_RK[11][4] = -477755414/1098053517.;
  a_RK[11][5] = -703635378/230739211.; a_RK[11][6] = 5731566787/1027545527.;
  a_RK[11][7] = 5232866602/850066563.; a_RK[11][8] = -4093664535/808688257.;
  a_RK[11][9] = 3962137247/1805957418.; a_RK[11][10] = 65686358/487910083.;
  a_RK[12][0] = 403863854/491063109.; a_RK[12][1] = 0; a_RK[12][2] = 0;
  a_RK[12][3] = -5068492393/434740067.; a_RK[12][4] = -411421997/543043805.;
  a_RK[12][5] = 652783627/914296604.; a_RK[12][6] = 11173962825/925320556.;
  a_RK[12][7] = -13158990841/6184727034.; a_RK[12][8] = 3936647629/1978049680.;
  a_RK[12][9] = -160528059/685178525.; a_RK[12][10] = 248638103/1413531060.;
  a_RK[12][11] = 0;
  b_RK[0] = 14005451/335480064.; b_RK[1] = 0; b_RK[2] = 0; b_RK[3] = 0;
  b_RK[4] = 0; b_RK[5] = -59238493/1068277825.; b_RK[6] = 181606767/758867731.;
  b_RK[7] = 561292985/797845732.; b_RK[8] = -1041891430/1371343529.;
  b_RK[9] = 760417239/1151165299.; b_RK[10] = 118820643/751138087.;
  b_RK[11] = -528747749/2220607170.; b_RK[12] = 1/4.;

  node0_printf("Integrator = INTEGRATOR_RKMK8\n");
#elif GF_INTEGRATOR==INTEGRATOR_ADAPT_LUSCHER
  N_stages = 3;
  // Luscher coefficients in proper 2N-storage format
  A_2N[0] = 0;
  A_2N[1] = -17/32.;
  A_2N[2] = -32/27.;
  B_2N[0] = 1/4.;
  B_2N[1] = 8/9.;
  B_2N[2] = 3/4.;
  Lambda[0] = -1;
  Lambda[1] = 2;
  Lambda[2] = 0;
  node0_printf("Integrator = INTEGRATOR_ADAPT_LUSCHER\n");
#elif GF_INTEGRATOR==INTEGRATOR_ADAPT_CF3
  N_stages = 3;
  // Commutator-free 2N-storage with arbitrary coefficients

#ifdef READ_ADPT_CF3_FROM_FILE
  FILE *fp;
  int st;
  fp = fopen( "cf3adpt_coeff.dat", "rt" );
  if( fp==NULL ) {
    node0_printf( "ERROR: cf3adpt_coeff.dat file cannot be read\n" );
    terminate(1);
  }
#if (MILC_PRECISION==1)
  st = fscanf( fp, "%e", &(A_2N[0]) );
  st += fscanf( fp, "%e", &(A_2N[1]) );
  st += fscanf( fp, "%e", &(A_2N[2]) );
  st += fscanf( fp, "%e", &(B_2N[0]) );
  st += fscanf( fp, "%e", &(B_2N[1]) );
  st += fscanf( fp, "%e", &(B_2N[2]) );
  st += fscanf( fp, "%e", &(Lambda[0]) );
  st += fscanf( fp, "%e", &(Lambda[1]) );
  st += fscanf( fp, "%e", &(Lambda[2]) );
#else
  st = fscanf( fp, "%le", &(A_2N[0]) );
  st += fscanf( fp, "%le", &(A_2N[1]) );
  st += fscanf( fp, "%le", &(A_2N[2]) );
  st += fscanf( fp, "%le", &(B_2N[0]) );
  st += fscanf( fp, "%le", &(B_2N[1]) );
  st += fscanf( fp, "%le", &(B_2N[2]) );
  st += fscanf( fp, "%le", &(Lambda[0]) );
  st += fscanf( fp, "%le", &(Lambda[1]) );
  st += fscanf( fp, "%le", &(Lambda[2]) );
#endif
  fclose( fp );
  node0_printf( "%.16g\n", A_2N[0] );
  node0_printf( "%.16g\n", A_2N[1] );
  node0_printf( "%.16g\n", A_2N[2] );
  node0_printf( "%.16g\n", B_2N[0] );
  node0_printf( "%.16g\n", B_2N[1] );
  node0_printf( "%.16g\n", B_2N[2] );
  node0_printf( "%.16g\n", Lambda[0] );
  node0_printf( "%.16g\n", Lambda[1] );
  node0_printf( "%.16g\n", Lambda[2] );
#else
  // Luscher coefficients -- default
  A_2N[0] = 0;
  A_2N[1] = -17/32.;
  A_2N[2] = -32/27.;
  B_2N[0] = 1/4.;
  B_2N[1] = 8/9.;
  B_2N[2] = 3/4.;
  Lambda[0] = -1;
  Lambda[1] = 2;
  Lambda[2] = 0;
#endif
  node0_printf("Integrator = INTEGRATOR_ADAPT_CF3\n");
#elif GF_INTEGRATOR==INTEGRATOR_ADAPT_BS
  // Bogacki-Shampine integrator based on Ralston coefficients
  // NOTE: there is a fourth stage that produces a lower order approximation
  N_stages = 3;
  // NOTE: indexing is shifted by -1 in C compared to usual Butcher tables
  a_RK[1][0] = 1/2.;
  a_RK[2][0] = 0; a_RK[2][1] = 3/4.;
  b_RK[0] = 2/9.; b_RK[1] = 1/3.; b_RK[2] = 4/9.;
  // additional fourth stage coefficients
  a_RK[3][0] = 7/24.; a_RK[3][1] = 1/4.;
  a_RK[3][2] = 1/3.; a_RK[3][3] = 1/8.;
  node0_printf("Integrator = INTEGRATOR_ADAPT_BS\n");
#endif

}
