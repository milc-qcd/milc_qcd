/**************************** load_smear_routines.c ****************************/
/* MIMD version 6 */
/*
 *  This file contains a number of routines to load in 
 *  the smearing functions.
 *
 *
 */



#include "w_static_includes.h"

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
  const int magic_number = 43241 ;

  /* Memory for some header information ***/
  size_t  three_object = (size_t) 3 ;
  int header_data[3] ;


  /** open the file ****/
  if( (fp = fopen(filename ,"rb")) == NULL )
  {
    printf("Could not open the file %s\n",filename);
    terminate(1);
  }

  /*** read in the header information *********/
  if( fread(header_data,sizeof(int),three_object,fp) != three_object )
  {
    printf("There was an error during the reading of the smearing function HEADER \n");
    terminate(1);
  }
  
  /*** check the header information ****/

  if( header_data[0] != magic_number )
  {
    printf("ERROR::Mismatch of the magic numbers for the smearing function \n");
    printf("file = %d code = %d\n",header_data[0],magic_number);
    terminate(1);
  }


  if( header_data[1] != dim )
  {
    printf("ERROR: amount of data mismatch between file and code \n");
    printf("file = %d code = %d \n",header_data[1], dim  );
    terminate(1);
  }


  if( header_data[2] != nx  )
  {
    printf("ERROR mismatch between code and file linear dimension\n");
    printf("file = %d code = %d\n", header_data[2],nx);
    terminate(1);
  }


  /*** read the smearing function  ***/
  if( fread(data,sizeof(Real),nobj,fp) != nobj   )
  {
    printf("There was an error during the reading of the smearing function \n");
    terminate(1);
  }

  /*** close the file ***/
  if( fclose(fp) != 0 )
  {
    printf("There was an error during the closing of %s \n",filename);
    terminate(1);
  }

  printf("I have read a smearing function from the file %s\n",filename);



}


/*
 *  Read in all the smearing functions from disk into 
 *  memory.
 *
 *
 *
 *
 *
 */





void load_smearing()
{
  int aloop ;
  Real *wave_func ;
  complex lbuf ;
  int nxyz = nx*ny*nz ;
  int x,y,z,t ;
  int newnode , currentnode ;
  int where ;
  int where_space ;
  int i ;
  field_offset where_smear ;
  double starttime,endtime ; 

  /****---------------------------------------*****/
  starttime = dclock() ; 

  /** Reserve space for the smearing function on the master node ***/
  if(mynode()==0)
  {
    if( (wave_func = (Real *)calloc( (size_t) nxyz, sizeof(Real) )  ) == NULL )
    {
      printf("ERROR: could not reserve buffer space for wave function on the master node \n");
      terminate(1);
    }

  }


  /*** loop over the different smearing sources *****/
  for(aloop=0 ; aloop < nosmear ;++aloop)
  {

    /** load in the smearing function into the master node buffer  **/
    if(mynode()==0)  
    {

      /** load the smearing function ****/
      load_scalar_smear(wave_func, nxyz,smearfile_in[aloop]  );

    }


    /*** send the smearing function to the correct nodes **/
    currentnode = 0;
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
	      { /* just copy links */
		i = node_index(x,y,z,t);
		where_smear = F_OFFSET(smear_func[aloop] );
		*( (complex *)F_PT( &(lattice[i]), where_smear ) )=lbuf;
	      }
	      else 
	      {              /* send to correct node */
		send_field((char *)&lbuf,sizeof(complex),currentnode);
	      }
	    }
 
	    /* The node which contains this site reads message */
	    else
	    { 
	      /* for all nodes other than node 0 */
	      if(this_node==currentnode)
	      {
		get_field((char *)&lbuf,sizeof(complex));
		i = node_index(x,y,z,t);
		where_smear = F_OFFSET(smear_func[aloop] );
		*( (complex *)F_PT( &(lattice[i]), where_smear ) )=lbuf;
	      }
	    }
	  
	  } /** end the loop over the spatial lattice ****/

  } /*** end the loop over smearing functions ****/



  /** free up the memory buffer on the master node ***/
  if(mynode()==0)  
    free(wave_func);

  endtime = dclock() - starttime ;

  IF_MASTER
    printf("Time to load smearing functions = %g sec \n",endtime );
}


