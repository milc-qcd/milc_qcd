/************* evolve_hqet_prop.c **************************/
/* MIMD version 6 */
#include "hqet_light_includes.h"

#ifdef DEBUGDEF
#include DEBUGDEF
#endif

/*   Evolve the HQET propagator forwards and backwards in time
 *   for the source time slice. This is just a driver routine.
 *
 *
 *  Subroutine arguments
 *    On input 
 *       tsrc :: the time slice that the source is on, this is where
 *               the HQET evolution starts from.
 *       v    :: pointer to the velocity
 *
 *    On output
 *       hqet_prop :: site structure pointer to su3_matrix, that contains the 
 *                    HQET propagator
 */



void evolve_hqet_forwards_and_backwards(field_offset hqet_prop, int tsrc, int v)
{
  int tcurrent = 2*nt ; /*** dummy variable ****/
  int tend ; 
  double t_start ; 

  t_start = dclock() ; 

  /*** generate the forward moving hqet propagator ***/
  tend = nt - 1  ; 
  generate_hqet_prop(hqet_prop, tsrc, tend, tcurrent, v, v) ;

  /*** generate the backward moving hqet propagator ***/
  tend = 0 ;
  generate_hqet_prop_back(hqet_prop, tsrc,tend ,tcurrent, v, v)  ; 


  IF_VERBOSE_ON(1)
    printf("evolve_hqet_forwards_and_backwards::Total time to generate HQET propagators = %g sec\n",dclock() - t_start) ;


}  /*** end of the function evolve_hqet_forwards_and_backwards   ***/





