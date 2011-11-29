/************************ control.c ******************************/
/* MIMD version 7 */
/* Main procedure for smoothing SU(3) instantons */
/* 2/19/98 if continuing from prev lattice, 
           allow total_sweeps to accumulate CD */

/* Uses APE or HYP blocking to achieve smoothing */

#define CONTROL

#include "smooth_inst_includes.h"
#include "lattice_qdp.h"
#ifdef HAVE_QIO
#include <qio.h>
#include "../include/io_scidac.h"
#endif

int main(int argc, char *argv[])
{
   int meascount;
   int todo;
   int prompt;
   double ssplaq;
   double stplaq;
   complex plp;
   double dtime;

   initialize_machine(&argc, &argv);

  /* Remap standard I/O */
  if(remap_stdio_from_args(argc, argv) == 1)terminate(1);

   g_sync();

   /* set up */
   prompt = setup();
   make_loop_table2();

   if(startflag != CONTINUE)total_sweeps = 0;

   /* loop over input sets */
   while ( readin(prompt) == 0 )
   {
      total_sweeps += sweeps;
      dtime = -dclock();

      /* perform smoothing sweeps, reunitarizing and measuring */

      meascount = 0;   /* number of measurements */
      for (todo=sweeps; todo > 0; --todo )
      {
         /* do one smoothing sweep */
         smooth();

         /* measure every "measinterval" trajectories */
         if ((todo%measinterval) == 0)
         {
            /* call plaquette measuring process */
            d_plaquette(&ssplaq, &stplaq);

            /* don't bother to */
            /* call the Polyakov loop measuring program */
            /* plp = ploop(); */
            plp = cmplx(99.9, 99.9);

            ++meascount;
            if (this_node==0)
            {
               /* Re(Polyakov) Im(Poyakov) cg_iters ss_plaq st_plaq */
               printf("GMES %e %e %e %e %e\n",
                      (double)plp.real, (double)plp.imag, 99.9,
                      (double)ssplaq, (double)stplaq);
            }

            fflush(stdout);
         }
      }       /* end loop over sweeps */


      /* perform gauge fixing after smoothing if requested */

      if ( fixflag == COULOMB_GAUGE_FIX )
      {
         if (this_node == 0)
         {
            printf("Fixing to Coulomb gauge\n");
         }
         gaugefix(TUP, (Real)1.5, 500, GAUGE_FIX_TOL);
      }
      else if (this_node == 0)
      {
         printf("GAUGE FIXING SKIPPED.\n");
      }

      /* measure the instanton charge on the fitted config */
      instanton_density();
	  
      if (this_node==0)
	{
	  printf("RUNNING COMPLETED\n");
	}
      
      dtime += dclock();
      if (this_node==0)
	{
	  printf("Time = %e seconds\n", dtime);
	}
      fflush(stdout);
      
      /* save FFdual */
      if(savetopoflag != FORGET) save_topo(topofile);

      /* save lattice if requested */
      if ( saveflag != FORGET )
      {
	save_lattice( saveflag, savefile, stringLFN );
      }
   }

   return 0;
}
