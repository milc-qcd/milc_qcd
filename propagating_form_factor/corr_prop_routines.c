/************************ corr_prop_routines.c ************************/  
/* MIMD version 6 */
/*  Routines to set up and save the various types of   
 *  two and three point functions
 *
 *
 *
*/


#include "prop_form_includes.h"
#include "opertypes.h"
#include "corrlist.h"


#ifdef DEBUGDEF
#include "debug_form.h"
#endif


/*
 *  Set up the memory for the heavy-light correlator
 *
 *
 */

/* Modifications

   C. DeTar 4/21/98 (lean) Provision for storing values only for time slices on node

 */

void setup_w_meson_store()
{
  /* Provides for efficient storage of correlators -
     i.e. store only values for the time slices on this node */
  int x,y,z,t,i;
  int found;

  /* Allocate space for map */
  w_meson_my_t = (int *)malloc(nt*sizeof(int));
  if(w_meson_my_t == NULL)
    {
      printf("setup_w_meson_store: Can't malloc for w_meson_my_t\n");
      terminate(1);
    }

  /* Run through all sites and determine which time slices are on this node */
  w_meson_nstore = 0;
  for(t = 0; t < nt; t++)
    {
      found = 0;
      for(x=0;x<nx;x++)
	{
	  for(y=0;y<ny;y++)
	    {
	      for(z=0;z<nz;z++)
		if(node_number(x,y,z,t) == this_node)
		  {
		    w_meson_my_t[w_meson_nstore] = t;
		    w_meson_nstore++;
		    found = 1;
		    break;
		  }
	      if(found == 1)break;
	    }
	  if(found == 1)break;
	}
    }
  
  /* Allocate space for inverse map */
  w_meson_store_t = (int *)malloc(nt*sizeof(int));
  if(w_meson_store_t == NULL)
    {
      printf("setup_w_meson_store: Can't malloc for w_meson_store_t\n");
      terminate(1);
    }

  for(t = 0; t < nt; t++)
    w_meson_store_t[t] = 0;

  /* Make inverse map */
  for(i = 0; i < w_meson_nstore; i++)
    w_meson_store_t[w_meson_my_t[i]] = i;

  /* Allocate space for correlator */
  w_meson_corr = (complex *)malloc(nt*sizeof(complex));
  if(w_meson_corr == NULL)
    {
      printf("setup_w_meson_store: Can't malloc for w_meson_corr\n");
      terminate(1);
    }
}

void setup_HL3_corr(complex **corr, int *corr_dim, int *corr_stride )
{
  int i ; 
  int no_oper = MAX_THREEPT - 1 ; 
  int q_pt = no_q_values - 1 ;
  int p_pt = no_p_values - 1 ;
  int spect_pt = no_spectator - 1 , zonked_pt = no_zonked_light - 1 ;
  int seq_pt = no_sequential - 1 ; 
  /********************************************************************************/
                 
  *corr_stride = LIGHT_FORM_WHERE(w_meson_nstore,zonked_pt,seq_pt,spect_pt,q_pt,p_pt, no_oper )  ; 
  *corr_dim = three_pt_copies*(*corr_stride) ;



  /* allocate space for lattice, fill in parity, coordinates and index,*/
  node0_printf("Allocating %.1f MBytes per CPU for HL3 correlator\n",
	       (double)*corr_dim * sizeof(complex)/1e6);

  if( ( *corr = (complex *) calloc( (size_t) *corr_dim, sizeof(complex) )  ) == NULL )
  {
    printf("node[%d] ERROR: could not reserve buffer space for the H-L correlators\n",this_node);
    terminate(1);
  }


  for(i = 0 ; i  < *corr_dim ; ++i)
  {
    (*corr + i)->real = 0.0 ; 
    (*corr + i)->imag = 0.0 ; 
  }




}




/*
 *  Set up the memory for the heavy-heavy three point function.
 *
 *
 */

void setup_HH3_corr(complex **corr, int *corr_dim, int *corr_stride )
{
  int i ; 
  int no_oper = MAX_THREEPT - 1 ; 
  int q_pt = no_q_values - 1 ;
  int p_pt = no_p_values - 1 ;
  int spect_pt = no_spectator - 1 , zonked_pt = no_zonked_heavy - 1 ;
  int seq_pt = no_sequential - 1 ; 
  /************************************************************/

  *corr_stride = HEAVY_FORM_WHERE(w_meson_nstore,zonked_pt,seq_pt,spect_pt,q_pt,p_pt, no_oper )  ; 
  *corr_dim = three_pt_copies*(*corr_stride) ;
                 

  node0_printf("Allocating %.1f MBytes per CPU for HH3 correlator\n",
	       (double) *corr_dim * sizeof(complex)/1e6);

  if( ( *corr = (complex *) calloc( (size_t) *corr_dim, sizeof(complex) )  ) == NULL )
  {
    printf("node[%d] ERROR: could not reserve buffer space for the H-H correlators\n",this_node);
    terminate(1);
  }


  for(i = 0 ; i  < *corr_dim ; ++i)
  {
    (*corr + i)->real = 0.0 ; 
    (*corr + i)->imag = 0.0 ; 
  }




}



/*             LIGHT --- LIGHT 
 *  Set up the memory for the required two point
 *  functions.
 *
 */

void setup_LL2_corr(complex **two_corr, int *two_corr_dim )
{
  int i ; 
  int no_oper = MAX_TWOPT - 1 ; 
  int q_pt = no_k_values - 1 ;
  int k_spectator = no_spectator - 1  ;
  int k_zonked_light = no_zonked_light - 1 ;


  *two_corr_dim = LL_TWOPT_FORM_WHERE(w_meson_nstore,k_zonked_light,k_spectator,q_pt  ,no_oper ) ;

  node0_printf("Allocating %.1f MBytes per CPU for LL2 correlator\n",
	       (double) *two_corr_dim * sizeof(complex)/1e6);

  if( ( *two_corr = (complex *) calloc( (size_t) *two_corr_dim, sizeof(complex) )  ) == NULL )
  {
    printf("node[%d] ERROR: could not reserve buffer space for the two point correlators\n",this_node);
    terminate(1);
  }


  for(i = 0 ; i  < *two_corr_dim ; ++i)
  {
    (*two_corr + i)->real = 0.0 ; 
    (*two_corr + i)->imag = 0.0 ; 
  }

}




/*             HEAVY --- LIGHT 
 *  Set up the memory for the required two point
 *  functions.
 *
 */

void setup_HL2_corr(complex **two_corr, int *two_corr_dim )
{
 int i ; 
 int no_oper = MAX_TWOPT - 1 ; 
 int q_pt = no_k_values - 1 ;
 int k_spectator = no_spectator - 1  ;
 int k_zonked_heavy = no_zonked_heavy - 1 ;


  *two_corr_dim = HL_TWOPT_FORM_WHERE(w_meson_nstore,k_zonked_heavy,k_spectator,q_pt  ,no_oper ) ;

  node0_printf("Allocating %.1f MBytes per CPU for HL2 correlator\n",
	       (double) *two_corr_dim * sizeof(complex)/1e6);

  if( ( *two_corr = (complex *) calloc( (size_t) *two_corr_dim, sizeof(complex) )  ) == NULL )
  {
    printf("node[%d] ERROR: could not reserve buffer space for the two point correlators\n",this_node);
    terminate(1);
  }


  for(i = 0 ; i  < *two_corr_dim ; ++i)
  {
    (*two_corr + i)->real = 0.0 ; 
    (*two_corr + i)->imag = 0.0 ; 
  }

}






/*             HEAVY --- LIGHT 
 *  Set up the memory for the required two point
 *  functions. These two point function have the operators
 *  inserted at both the source and sink.
 */

void setup_HL2_corr_with_rotations(complex **two_corr, int *two_corr_dim, int *two_corr_stride )
{
 int i ; 
 int no_oper = MAX_TWOPT - 1 ; 
 int q_pt = no_k_values - 1 ;
 int k_spectator = no_spectator - 1  ;
 int k_zonked_heavy = no_zonked_heavy - 1 ;

 *two_corr_stride = HL_TWOPT_FORM_WHERE(w_meson_nstore,k_zonked_heavy,k_spectator,q_pt  ,no_oper ) ;
 *two_corr_dim = two_pt_copies * (*two_corr_stride) ;

  if( ( *two_corr = (complex *) calloc( (size_t) *two_corr_dim, sizeof(complex) )  ) == NULL )
  {
    printf("node[%d] ERROR: could not reserve buffer space for the two point correlators\n",this_node);
    terminate(1);
  }


  for(i = 0 ; i  < *two_corr_dim ; ++i)
  {
    (*two_corr + i)->real = 0.0 ; 
    (*two_corr + i)->imag = 0.0 ; 
  }

}


/*
 *  Assemble next set of t values
 *
 */

void prop_assemble_t(complex *corri,Real norm)
{
  int i,t;
  
  for(t = 0; t < nt; t++)
    w_meson_corr[t] = cmplx(0.,0.);

  for(i = 0; i < w_meson_nstore; i++)
    w_meson_corr[w_meson_my_t[i]] = corri[i];

  g_veccomplexsum(w_meson_corr,nt);

  for(t = 0; t < nt; t++)
    {
      w_meson_corr[t].real *= norm;
      w_meson_corr[t].imag *= norm;
    }
}

void prop_corr_write(FILE *fp,char *filename,
		     complex *corr,int corr_dim,Real norm)
{
  complex *corri;
  int i;

  for(i = 0, corri = corr; i < corr_dim/w_meson_nstore; 
      i++, corri += w_meson_nstore)
    {
      prop_assemble_t(corri,norm);
      IF_MASTER
	write_corr_item(fp,w_meson_corr,filename);
    }
}



/*
 *  A routine to finish up the calculation of the 
 *  heavy-light three point functions. The correlators
 *  are summed over the nodes, and the master node
 *  writes the result to a file.
 */

void finish_HL3_corr(complex *corr, int corr_dim, int copy_stride)
{
  FILE *fp ;
  Real norm;
  double time;

  time = -dclock();
  norm = 1.;

  /*** write the correlators to disk ****/
  if(saveflag_HL3 == SAVE_SERIAL)
    {
      IF_MASTER
	fp = write_prop_form_header(corr, filename_HL3, HEAVY_TO_LIGHT,
			no_zonked_light, MAX_THREEPT, three_pt_copies, 
			corr_dim*nt/w_meson_nstore) ;
      
      prop_corr_write(fp, filename_HL3, corr, corr_dim, norm);
      IF_MASTER
	write_prop_form_close(fp, filename_HL3);
    }

  else if(saveflag_HL3 == SAVE_ASCII)
    {
      dump_hl_prop_form_text(corr, filename_HL3, MAX_THREEPT , copy_stride, norm) ;
      time += dclock();
      node0_printf("Time to write HL3 %g sec\n",time);
    }
  else
    IF_MASTER
      printf("SKIPPING heavy_light form factor output as requested.\n");
}






/*
 *  A routine to finish up the calculation of the 
 *  heavy to heavy three point functions. The correlators
 *  are summed over the nodes, and the master node
 *  writes the result to a file.
 */

void finish_HH3_corr(complex *corr, int corr_dim, int copy_stride)
{
  FILE *fp ; 
  Real norm;
  double time;

  time = -dclock();
  norm = 1.;

  /*** write the correlators to disk ****/
  if(saveflag_HH3 == SAVE_SERIAL)
    {
      IF_MASTER
	fp = write_prop_form_header(corr, filename_HH3, HEAVY_TO_HEAVY,
			    no_zonked_heavy, MAX_THREEPT, three_pt_copies, 
			    corr_dim*nt/w_meson_nstore) ;
      
      prop_corr_write(fp, filename_HH3, corr, corr_dim, norm);
      IF_MASTER
	write_prop_form_close(fp, filename_HH3);
    }

  else if(saveflag_HH3 == SAVE_ASCII)
    {
      dump_hh_prop_form_text( corr  , filename_HH3, MAX_THREEPT , copy_stride, norm );
      time += dclock();
      node0_printf("Time to write HH3 %g sec\n",time);
    }
  else
    IF_MASTER
      printf("SKIPPING heavy_heavy form factor output as requested.\n");
}


/*
 *  A routine to finish up the calculation of the 
 *  two point functions
 *
 */

void finish_LL2_GG_corr(complex *corr, int corr_dim )
{
  FILE *fp ; 
  const int space_vol = nx*ny*nz ; 
  Real norm;
  double time;

  time = -dclock();

  norm = 1./space_vol;

  /*** write the correlators to disk ****/

  if(saveflag_LL2_GG == SAVE_SERIAL)
    {
      IF_MASTER
	fp = write_prop_twopt_header(corr, filename_LL2_GG  , LL_2PT,
			   no_spectator,no_zonked_light , MAX_TWOPT, 1, 
			   corr_dim*nt/w_meson_nstore )  ;
      
      prop_corr_write(fp, filename_LL2_GG, corr, corr_dim, 
		      norm);
      IF_MASTER
	write_prop_twopt_close(fp, filename_LL2_GG);
    }

  else if(saveflag_LL2_GG == SAVE_ASCII)
    {
      dump_ll_twopt_text(corr , filename_LL2_GG, MAX_TWOPT, norm ) ;
      time += dclock();
      node0_printf("Time to write LL2_GG %g sec\n",time);
    }
      
  else
    IF_MASTER
      printf("SKIPPING LL2_GG correlator output as requested.\n");
}


/*
 *  A routine to finish up the calculation of the 
 *  two point functions
 *
 */

void finish_HL2_GE_corr(complex *corr, int corr_dim )
{
  FILE *fp ; 
  const int space_vol = nx*ny*nz ; 
  Real norm;
  double time;

  time = -dclock();

  norm = 1./space_vol;

  /*** write the correlators to disk ****/
  if(saveflag_HL2_GE == SAVE_SERIAL)
    {
      IF_MASTER
	fp = write_prop_twopt_header(corr, filename_HL2_GE, HL_REL_2PT, 
			  no_spectator,no_zonked_heavy , MAX_TWOPT, 1, 
			  corr_dim*nt/w_meson_nstore   )  ;
      
      prop_corr_write(fp, filename_HL2_GE, corr, corr_dim, 
		      1./space_vol);
      IF_MASTER
	write_prop_twopt_close(fp, filename_HL2_GE);
    }
  else if(saveflag_HL2_GE == SAVE_ASCII)
    {
      dump_hl_twopt_text(corr,filename_HL2_GE, MAX_TWOPT, "GAUSS_SEQUENTIAL", norm  ); 
      time += dclock();
      node0_printf("Time to write HL2_GE %g sec\n",time);
    }
  else
    IF_MASTER
      printf("SKIPPING HL2_GE correlator output as requested.\n");
}



/*
 *  A routine to finish up the calculation of the 
 *  two point functions
 *
 */

void finish_HL2_GG_corr(complex *corr, int corr_dim )
{
  FILE *fp;
  const int space_vol = nx*ny*nz ; 
  Real norm;
  double time;

  time = -dclock();

  norm = 1./space_vol;

  /*** write the correlators to disk ****/
  if(saveflag_HL2_GG == SAVE_SERIAL)
    {
      IF_MASTER
	fp = write_prop_twopt_header(corr, filename_HL2_GG , HL_2PT_BAG,
			    no_spectator,no_zonked_heavy , MAX_TWOPT, 1, 
			    corr_dim*nt/w_meson_nstore  )  ;
      
      prop_corr_write(fp, filename_HL2_GG, corr, 
		      corr_dim, 1./space_vol);
      IF_MASTER
	write_prop_twopt_close(fp, filename_HL2_GG);
    }

  else if(saveflag_HL2_GG == SAVE_ASCII)
    {
      dump_hl_twopt_text(corr, filename_HL2_GG, MAX_TWOPT, "GAUSS_GAUSS", norm  ); 
      time += dclock();
      node0_printf("Time to write HL2_GG %g sec\n",time);
    }
  else
    IF_MASTER
      printf("SKIPPING HL2_GG correlator output as requested.\n");
}



/*
 *  A routine to finish up the calculation of the 
 *  two point functions. The smeared source, local sink
 *  correlators, with the O(a) current corrections.
 *
 *
 *
 */

void finish_HL2_GL_corr(complex *corr, int corr_dim , int copy_stride )
{
  FILE *fp;
  const int space_vol = nx*ny*nz ; 
  Real norm;
  double time;

  time = -dclock();

  norm = 1./space_vol;

  /*** write the correlators to disk ****/
  if(saveflag_HL2_GL == SAVE_SERIAL)
    {
      IF_MASTER
	fp = write_prop_twopt_header(corr,filename_HL2_GL,HL_LOCAL_SINK, 
		     no_spectator,no_zonked_heavy , MAX_TWOPT,two_pt_copies , 
		     corr_dim*nt/w_meson_nstore   )  ;
      
      prop_corr_write(fp, filename_HL2_GL, corr, 
		      corr_dim, 1./space_vol);
      IF_MASTER
	write_prop_twopt_close(fp, filename_HL2_GL);
    }

  else if(saveflag_HL2_GL == SAVE_ASCII)
    {
      dump_hl_twopt_text_oper_corrections(corr,filename_HL2_GL,MAX_TWOPT, copy_stride  ,"GAUSS_LOCAL", norm  ); 
      time += dclock();
      node0_printf("Time to write HL2_GL %g sec\n",time);
    }
  else
    IF_MASTER
      printf("SKIPPING HL2_GL correlator output as requested.\n");

}   /**** end of the finish_HL2_GL_corr function ***/




/*
 *
 *  Dump the heavy --> light two  point functions to 
 *  the screen.
 *
 *  This is a scalar routine.
 *
 */



void dump_hl_twopt_text(complex *corr, char *filename, int no_oper, char tag[], Real norm) 
{
  int where ;
  int t,zonk_pt,spect_pt,k_pt,oper_pt ;
  FILE *fp;



  if((fp = fopen(filename,"w")) == NULL)
    {
      printf("dump_hl_prop_form_text: Can't open %s for output\n",filename);
      terminate(1);
    }

  if(this_node==0)fprintf(fp,"HEAVY-LIGHT_TWO_POINT %s CORRELATORS START\n",tag);

  for(zonk_pt = 0 ; zonk_pt < no_zonked_heavy ; ++zonk_pt )
  {
    if(this_node==0)fprintf(fp,"zonked quark = %d kappa = %f\n",zonk_pt, kappa_zonked_heavy[zonk_pt]);
    for(spect_pt = 0 ; spect_pt < no_spectator ; ++spect_pt )
    {
      if(this_node==0)fprintf(fp,"spectator quark = %d kappa = %f\n",spect_pt,  kappa_spectator[spect_pt ] );
      for(k_pt = 0 ; k_pt < no_k_values ; ++k_pt ) 
      {
	if(this_node==0)fprintf(fp,"K--Momentum pointer = %d k = %d %d %d\n",k_pt,
	       k_momstore[k_pt][0] ,  k_momstore[k_pt][1] ,  k_momstore[k_pt][2] ); 
	  for(oper_pt = 0 ; oper_pt < no_oper ; ++oper_pt)
	  {
	    if(this_node==0)fprintf(fp,"Operator pointer = %d  name = %s\n",oper_pt, oper_name(oper_pt) );
	    where = HL_TWOPT_FORM_WHERE(0,zonk_pt ,spect_pt,k_pt, oper_pt) ;
	    prop_assemble_t(corr+where,norm);
	    if(this_node==0)for(t = 0 ; t < nt ; ++t)
	    {
	      fprintf(fp,"%d   %g %g\n",t , 
			   w_meson_corr[t].real, w_meson_corr[t].imag);
	    }

	  }
      }
    }
  }
  
  if(this_node==0)fprintf(fp,"HEAVY-LIGHT_TWO_POINT %s CORRELATORS END\n",tag);
  fclose(fp);


}  /*** end of the dump of the heavy-light 2pt functions ****/








/*
 *
 *  Dump the heavy --> light two  point functions to 
 *  the screen.
 *
 *  This is a scalar routine.
 *
 */



void dump_hl_twopt_text_oper_corrections(complex *corr, char *filename, int no_oper, int copy_stride, char tag[], Real norm) 
{
  int where ;
  int t,zonk_pt,spect_pt,k_pt,oper_pt ;
  int copy_pt ; 
  FILE *fp;



  if((fp = fopen(filename,"w")) == NULL)
    {
      printf("dump_hl_twopt_text_oper_corrections: Can't open %s for output\n",filename);
      terminate(1);
    }

  if(this_node==0)fprintf(fp,"HEAVY-LIGHT_TWO_POINT_OPER_CORRECTIONS %s CORRELATORS START\n",tag);

  for(zonk_pt = 0 ; zonk_pt < no_zonked_heavy ; ++zonk_pt )
  {
    if(this_node==0)fprintf(fp,"zonked quark = %d kappa = %f\n",zonk_pt, kappa_zonked_heavy[zonk_pt]);
    for(spect_pt = 0 ; spect_pt < no_spectator ; ++spect_pt )
    {
      if(this_node==0)fprintf(fp,"spectator quark = %d kappa = %f\n",spect_pt,  kappa_spectator[spect_pt ] );
      for(k_pt = 0 ; k_pt < no_k_values ; ++k_pt ) 
      {
	if(this_node==0)fprintf(fp,"K--Momentum pointer = %d k = %d %d %d\n",k_pt,
	       k_momstore[k_pt][0] ,  k_momstore[k_pt][1] ,  k_momstore[k_pt][2] ); 
	for(copy_pt = 0 ; copy_pt < two_pt_copies ; ++copy_pt) 
	  for(oper_pt = 0 ; oper_pt < no_oper ; ++oper_pt)
	  {
	    if(this_node==0)fprintf(fp,"Operator pointer = %d  name = %s\n",oper_pt, oper_name_and_corrections(oper_pt,copy_pt) );
	    where = copy_pt * copy_stride + HL_TWOPT_FORM_WHERE(0,zonk_pt ,spect_pt,k_pt, oper_pt) ;
	    prop_assemble_t(corr+where,norm);
	    if(this_node==0)for(t = 0 ; t < nt ; ++t)
	    {
	      fprintf(fp,"%d   %g %g\n",t , 
			   w_meson_corr[t].real, w_meson_corr[t].imag);
	    }

	  }
      }
    }
  }
  

  if(this_node==0)fprintf(fp,"HEAVY-LIGHT_TWO_POINT_OPER_CORRECTIONS  %s CORRELATORS END\n",tag);

  fclose(fp);

}  /*** end of the dump of the heavy-light 2pt functions ****/


/*
 *   Write the light-light two point function to the 
 *   screen.
*/



void dump_ll_twopt_text(complex *corr,  char *filename, int no_oper, Real norm)
{
  int where ;
  int t,zonk_pt,spect_pt,k_pt ,oper_pt ;
  FILE *fp;



  if((fp = fopen(filename,"w")) == NULL)
    {
      printf("dump_ll_twopt_text: Can't open %s for output\n",filename);
      terminate(1);
    }

  if(this_node==0)fprintf(fp,"LIGHT-LIGHT_TWO_POINT CORRELATORS START\n");

  for(zonk_pt = 0 ; zonk_pt < no_zonked_light ; ++zonk_pt )
  {
    if(this_node==0)fprintf(fp,"zonked quark = %d kappa = %f\n",zonk_pt, kappa_zonked_light[zonk_pt]);
    for(spect_pt = 0 ; spect_pt < no_spectator ; ++spect_pt )
    {
      if(this_node==0)fprintf(fp,"spectator quark = %d kappa = %f\n",spect_pt,  kappa_spectator[spect_pt ] );
      for(k_pt = 0 ; k_pt < no_k_values ; ++k_pt ) 
      {
	if(this_node==0)fprintf(fp,"K--Momentum pointer = %d k = %d %d %d\n",k_pt,
	       k_momstore[k_pt][0] ,  k_momstore[k_pt][1] ,  k_momstore[k_pt][2] ); 
	  for(oper_pt = 0 ; oper_pt < no_oper ; ++oper_pt)
	  {
	    if(this_node==0)fprintf(fp,"Operator pointer = %d  name = %s\n",oper_pt, oper_name(oper_pt) );
	    where = LL_TWOPT_FORM_WHERE(0,zonk_pt ,spect_pt,k_pt, oper_pt) ;
	    prop_assemble_t(corr+where,norm);
	    if(this_node==0)for(t = 0 ; t < nt ; ++t)
	    {
	      fprintf(fp,"%d   %g %g\n",t , 
			   w_meson_corr[t].real, w_meson_corr[t].imag);
	    }
	    
	  }
      }
    }
  }
  
  if(this_node==0)fprintf(fp,"LIGHT-LIGHT_TWO_POINT CORRELATORS END\n");
  
  fclose(fp);

}  /*** end of the dump of the light-light 2pt functions ****/


/* 
 *  Return the name of the two point function
 *  Subroutine arguments:
 *    On input
 *      n       :: number of dirac matrix in the operator
 *      copy_pt :: the type of operator insertion
 */

char *oper_name_and_corrections(int n , int copy_pt )
{
  static char name[MAXFILENAME] ; 
  
  assert( n >= 0 && n < NO_TWOPT_OPERS ); 
  assert( copy_pt >= 0 && copy_pt <  two_pt_copies  ) ; 


  switch( copy_pt ) 
  {
   case 0 :  
     sprintf(name,"%s",twopt_name[two_pt[n].oper])  ;
     break ;
   case 1 :  
     sprintf(name,"%s_3D_on_spect",  twopt_name[two_pt[n].oper])  ;
     break ;
   case 2 :  
     sprintf(name,"%s_3D_on_zonked",twopt_name[two_pt[n].oper]) ; 
     break ;
   case 3 :  
     sprintf(name,"%s_4D_on_spect",twopt_name[two_pt[n].oper]) ;
     break ;
   case 4 :  
     sprintf(name,"%s_4D_on_zonked",twopt_name[two_pt[n].oper] )  ;
     break ;
   default :
     printf("ERROR:oper_name_and_corrections: copy_pt = %d is out of range\n",copy_pt); 
     terminate(1) ;
     break ;
  }



  return &name[0] ;


}



/* 
 *  Return the name of the two point function
 *  when there are no insertions of O(a) corrections
 *  to the current.
 */

char *oper_name(int n )
{

  return oper_name_and_corrections(n,0) ;

}



/* 
 *  Return the name of the three point function
 *
 */



char *three_oper_name(int n , int copy_pt )
{
  static char name[MAXFILENAME] ; 


  assert( copy_pt >= 0 && copy_pt <  three_pt_copies  ) ; 
  assert( n >= 0 && n < MAX_THREEPT) ;

  switch( copy_pt ) 
  {
   case 0 :  
     sprintf(name,"%s",name_3pt[three_pt[n].oper])  ;
     break ;
   case 1 :  
     sprintf(name,"%s_3D_on_seq",name_3pt[three_pt[n].oper])  ;
     break ;
   case 2 :  
     sprintf(name,"%s_3D_on_zonked",name_3pt[three_pt[n].oper]) ; 
     break ;
   case 3 :  
     sprintf(name,"%s_4D_on_seq",name_3pt[three_pt[n].oper]) ;
     break ;
   case 4 :  
     sprintf(name,"%s_4D_on_zonked",name_3pt[three_pt[n].oper] )  ;
     break ;
   default :
     printf("ERROR:three_oper_name: copy_pt = %d is out of range\n",copy_pt); 
     terminate(1) ;
     break ;
  }


  return &name[0] ;


}






/*
 *
 *  Dump the heavy --> light three point function to the
 *  terminal screen.
 *
 *
 *
 */



void dump_hl_prop_form_text(complex *corr,  char *filename, int no_oper, int copy_stride, Real norm)
{
  int where ;
  int t,zonk_pt,spect_pt,q_pt,p_pt,oper_pt ;
  int seq_pt ; 
  int copy_pt ; 
  FILE *fp;



  if((fp = fopen(filename,"w")) == NULL)
    {
      printf("dump_hl_prop_form_text: Can't open %s for output\n",filename);
      terminate(1);
    }

  if(this_node==0)fprintf(fp,"HEAVY-LIGHT_THREE_POINT CORRELATORS START\n");

  for(zonk_pt = 0 ; zonk_pt < no_zonked_light ; ++zonk_pt )
  {
    if(this_node==0)fprintf(fp,"zonked quark = %d kappa = %f\n",zonk_pt, kappa_zonked_light[zonk_pt]);
    for(spect_pt = 0 ; spect_pt < no_spectator ; ++spect_pt )
    {
      if(this_node==0)fprintf(fp,"spectator quark = %d kappa = %f\n",spect_pt,  kappa_spectator[spect_pt ] );
	for(seq_pt = 0 ; seq_pt <no_sequential  ; ++seq_pt)
	{
	  if(this_node==0)fprintf(fp,"sequential quark = %d, kappa = %f\n",seq_pt,   kappa_sequential[seq_pt ] );
	  for(q_pt = 0 ; q_pt < no_q_values ; ++q_pt ) 
	  {
	    if(this_node==0)fprintf(fp,"Q--Momentum pointer = %d q = %d %d %d\n",q_pt,
		   q_momstore[q_pt][0] ,  q_momstore[q_pt][1] ,  q_momstore[q_pt][2] ); 
	    for(p_pt = 0 ; p_pt < no_p_values ; ++p_pt)
	    {
	    if(this_node==0)fprintf(fp,"P--Momentum pointer = %d p = %d %d %d\n",p_pt,
		   p_momstore[p_pt][0] ,  p_momstore[p_pt][1] ,  p_momstore[p_pt][2] ); 

	    for(copy_pt = 0 ; copy_pt < three_pt_copies  ; ++copy_pt)
	      for(oper_pt = 0 ; oper_pt < no_oper ; ++oper_pt)
	      {
		if(this_node==0)fprintf(fp,"Operator pointer = %d  name = %s\n",oper_pt, three_oper_name(oper_pt,copy_pt) );
		where = copy_pt*copy_stride + LIGHT_FORM_WHERE(0,zonk_pt,seq_pt,spect_pt,q_pt,p_pt,oper_pt) ;
		prop_assemble_t(corr+where,norm);
		if(this_node==0)for(t = 0 ; t < nt ; ++t)
		  {
		    fprintf(fp,"%d   %g %g\n",t , 
				 w_meson_corr[t].real, w_meson_corr[t].imag);
		  }
	      }
	    }
	  }
	}
    }
  }

  if(this_node==0)fprintf(fp,"HEAVY-LIGHT_THREE_POINT CORRELATORS END\n");

  fclose(fp);


}


/*
 *   Dump the heavy to heavy form three point functions
 *   to standard output
 *
 *   Subroutine arguments
 *    On input
 *     corr        :: complex array of the 3pt correlators
 *     no_oper     :: the number of operators
 *     copy_stride :: amout of data for one specific choice of ROTation
 */


void dump_hh_prop_form_text(complex *corr,  char *filename, int no_oper, int copy_stride, Real norm)
{
  int where ;
  int t,zonk_pt,spect_pt,q_pt,p_pt,oper_pt ;
  int seq_pt ; 
  Real total = 0 ; 
  int copy_pt ; 
  FILE *fp;


  if((fp = fopen(filename,"w")) == NULL)
    {
      printf("dump_hh_prop_form_text: Can't open %s for output\n",filename);
      terminate(1);
    }

  if(this_node==0)fprintf(fp,"HEAVY-HEAVY_THREE_POINT CORRELATORS START\n");

  for(zonk_pt = 0 ; zonk_pt < no_zonked_heavy ; ++zonk_pt )
  {
    if(this_node==0)fprintf(fp,"zonked quark = %d kappa = %f\n",zonk_pt, kappa_zonked_heavy[zonk_pt]);
    for(spect_pt = 0 ; spect_pt < no_spectator ; ++spect_pt )
    {
      if(this_node==0)fprintf(fp,"spectator quark = %d kappa = %f\n",spect_pt,  kappa_spectator[spect_pt ] );
	for(seq_pt = 0 ; seq_pt <no_sequential  ; ++seq_pt)
	{
	  if(this_node==0)fprintf(fp,"sequential quark = %d, kappa = %f\n",seq_pt,   kappa_sequential[seq_pt ] );
	  for(q_pt = 0 ; q_pt < no_q_values ; ++q_pt ) 
	  {
	    if(this_node==0)fprintf(fp,"Q--Momentum pointer = %d q = %d %d %d\n",q_pt,
		   q_momstore[q_pt][0] ,  q_momstore[q_pt][1] ,  q_momstore[q_pt][2] ); 

	    for(p_pt = 0 ; p_pt < no_p_values ; ++p_pt)
	    {
	    if(this_node==0)fprintf(fp,"P--Momentum pointer = %d p = %d %d %d\n",p_pt,
		   p_momstore[p_pt][0] ,  p_momstore[p_pt][1] ,  p_momstore[p_pt][2] ); 
	    for(copy_pt = 0 ; copy_pt < three_pt_copies  ; ++copy_pt)
	      for(oper_pt = 0 ; oper_pt < no_oper ; ++oper_pt)
	      {
		if(this_node==0)fprintf(fp,"Operator pointer = %d  name = %s\n",oper_pt, three_oper_name(oper_pt ,copy_pt) );
		where = copy_pt*copy_stride + HEAVY_FORM_WHERE(0,zonk_pt,seq_pt,spect_pt,q_pt,p_pt,oper_pt)  ;
		prop_assemble_t(corr+where,norm);
		if(this_node==0)for(t = 0 ; t < nt ; ++t)
		  {
		    fprintf(fp,"%d   %g %g\n",t , 
				 w_meson_corr[t].real, w_meson_corr[t].imag);
		  }

/**DEBUG***/	  if( oper_pt == 0 ) total += (corr+where)->real ;

	      }
	    
	    
	    }
	  }
	}
    }
  }


  if(this_node==0)fprintf(fp,"HEAVY-HEAVY_THREE_POINT CORRELATORS END\n");

  if( this_node == 0 ) printf("DEBUG total = %f\n",total) ; 

  fclose(fp);


} /*** end of dump_hh_prop_form  *****/
