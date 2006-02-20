/**************************** load_smearing.c ****************************/
/* MIMD version 7 */
/*
 *  This file contains a number of routines to load in 
 *  the smearing functions.
 *
 *
 */


#include "generic_form_includes.h"


/*
 *  Ths routine loads in a single smearing function,
 *  over the spatial lattice, from a disk file.
 *
 *  This is scalar code -- it should only be run from
 *  one node (the master).
 *
 *   The smearing function is stored in a binary file
 *
 *
 */

void load_scalar_smear(Real *data, int dim, char filename[])
{
  FILE *fp ;
  size_t nobj = (size_t) dim  ;
  int magic_number = 43241 ;
  enum rev_choices { byte_rev , do_nothing } ;
  int byte_rev_flag = do_nothing ; 
  char myname[] = "load_scalar_smear";

  /* Memory for some header information ***/
  size_t  three_object = (size_t) 3 ;
  int32type header_data[3] ;


  /** open the file ****/
  if( (fp = fopen(filename ,"rb")) == NULL )
  {
    printf("ERROR:%s: Could not open the file %s\n",myname,filename);
    terminate(1);
  }

  /*** read in the header information *********/
  if( fread(header_data,sizeof(int32type),three_object,fp) != three_object )
  {
    printf("There was an error during the reading of the smearing function HEADER \n");
    terminate(1);
  }
  
  /*** check the header information ****/

  if( header_data[0] != magic_number )
  {
    printf("Reading with byte reversal\n") ; 
    byte_rev_flag =  byte_rev ; 
  }

  if(  byte_rev_flag ==  byte_rev ) 
    byterevn( (int32type*) header_data , three_object ) ; 

  if( header_data[0] != magic_number )
  {
    printf("ERROR:%s: Mismatch of the magic numbers for the smearing function \n",myname);
    printf("file = %d code = %d\n",header_data[0],magic_number);
    terminate(1);
  }


  if( header_data[1] != dim )
  {
    printf("ERROR: %s: amount of data mismatch between file and code \n",myname);
    printf("file = %d code = %d \n",header_data[1], dim  );
    terminate(1);
  }


  if( header_data[2] != nx  )
  {
    printf("ERROR: %s: mismatch between code and file linear dimension\n",myname);
    printf("file = %d code = %d\n", header_data[2],nx);
    terminate(1);
  }


  /*** read the smearing function  ***/
  if( fread(data,sizeof(Real),nobj,fp) != nobj   )
  {
    printf("%s: There was an error during the reading of the smearing function \n",myname);
    terminate(1);
  }


  if(  byte_rev_flag ==  byte_rev ) 
    byterevn( (int32type*) data , nobj ) ; 



  /*** close the file ***/
  if( fclose(fp) != 0 )
  {
    printf("%s: There was an error during the closing of %s \n",myname,filename);
    terminate(1);
  }

  printf("%s: I have read a smearing function from the file %s\n",myname,filename);



}


/*
 *   Load a complex smearing function from disk
 *   into a position in the site structure
 *
 *   Subroutine arguments
 *      On entry
 *         filename :: the name of the file to read the smearing function from
 *      On return
 *         where_smear :: site structure pointer to complex smaering function
 *
 */



void load_smearing(field_offset where_smear, char filename[])
{
  Real *wave_func ;
  complex lbuf ;
  int nxyz = nx*ny*nz ;
  int x,y,z,t ;
  int newnode , currentnode ;
  int where ;
  int where_space ;
  int i ;
  double starttime,endtime ; 

  /****---------------------------------------*****/
  starttime = dclock() ; 

  /** Reserve space for the smearing function on the master node ***/
  if(mynode()==0)
  {
    if( (wave_func = (Real *) calloc( (size_t) nxyz, sizeof(Real) )  ) == NULL )
    {
      printf("ERROR: load_smearing: could not reserve buffer space for wave function on the master node \n");
      terminate(1);
    }

  }


  /** load in the smearing function into the master node buffer  **/
  if(mynode()==0)  
  {
    /** load the smearing function ****/
    load_scalar_smear(wave_func, nxyz, filename);
  }

  
  g_sync();
  currentnode = 0 ;

  /*** send the smearing function to the correct nodes **/
  for(t=0;t<nt;t++)
    for(z=0;z<nz;z++)
      for(y=0;y<ny;y++)
	for(x=0;x<nx;x++)
	{
	  where = x + nx*(y + ny*( z + nz*t)) ;
	  where_space = x + nx*(y + ny*z) ;


	  newnode=node_number(x,y,z,t);
	  if(newnode != currentnode)
	  {
	    g_sync();
	    currentnode=newnode;
	  }
 
	  /* Node 0 reads, and sends site to correct node */
	  if(this_node==0)
	  {
	    lbuf.real = *(wave_func + where_space ) ;
	    lbuf.imag = 0.0 ;

	    if(currentnode==0)
	    {			/* just copy links */
	      i = node_index(x,y,z,t);
	      *( (complex *)F_PT( &(lattice[i]), where_smear ) )=lbuf;
	    }
	    else 
	    {			/* send to correct node */
	      send_field((char *)&lbuf,sizeof(complex),currentnode);
	    }
	  }
 
	  /* The node which contains this site reads message */
	  else
	  { 
	    /* for all nodes other than node 0 */
	    if(this_node==currentnode)
	    {
	      get_field((char *)&lbuf,sizeof(complex),0);
	      i = node_index(x,y,z,t);
	      *( (complex *)F_PT( &(lattice[i]), where_smear ) )=lbuf;
	    }
	  }
	  
	} /** end the loop over the spatial lattice ****/




  /** free up the memory buffer on the master node ***/
  if(mynode()==0)  
    free(wave_func);

  endtime = dclock() - starttime ;

  node0_printf("Time to load smearing functions = %g sec \n",endtime );

}
