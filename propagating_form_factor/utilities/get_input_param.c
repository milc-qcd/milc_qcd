/*
 *   Read the parameters for processing the matrix elements
 *
 *
 *
 */


#include "read_hl_form.h"
#include "prop_form_utilities.h"


/* Parameter format:
   
   three_point_select [fold | forward | backward ]
   SP <n> ZK <n> SQ <n> Q 1 P <n> OP 3 CP <n> WT <x>
   SP <n> ZK <n> SQ <n> Q 1 P <n> OP 3 CP <n> WT <x>
   ...

   two_point_recoil_select [fold | forward | backward ]
   SP <n> ZK <n> K <n> OP <n> CP <n> WT <x>
   SP <n> ZK <n> K <n> OP <n> CP <n> WT <x>
   ...

   two_point_sequential_select [fold | forward | backward ]
   SP <n> SQ <n> P <n> OP <n> CP <n> WT <x>
   SP <n> SQ <n> P <n> OP <n> CP <n> WT <x>
   ...

   filelist <filelistfilename>
   
   */

/* Sample format of <filelistfilename> = "HLfilelist"

   /work/u2219/prop_form/b560m01/0200/HL3_p0_0.1245.0200
      /work/u2219/prop_form/b560m01/0200/LL2_GG_0.1245.0200
      /work/u2219/prop_form/b560m01/0200/HL2_GE_0.1245.0200

   /work/u2219/prop_form/b560m01/0220/HL3_p0_0.1245.0220
      /work/u2219/prop_form/b560m01/0220/LL2_GG_0.1245.0220
      /work/u2219/prop_form/b560m01/0220/HL2_GE_0.1245.0220

      */

/*
 *  Read a single integer from the file, with the format
 *     tag 10
 */ 
int  load_single_i(char tag[], FILE *fp, char *buf, int *err )
{
  int ans ; 


  if( fscanf(fp,"%s",buf) != 1 )
  {
    /* Don't complain if EOF reached */
    if(!feof(fp))printf("ERROR in reading required tag %s\n",tag); 
    *err |= 1; 
    return 0;
  }

  if( strcmp(buf,tag ) != 0    )
  {
    *err |= 2;
    return 0;
  }

  if( fscanf(fp,"%d",&ans) != 1 )
  {
    printf("ERROR reading number after the tag \"%s\" \n",tag); 
    *err |= 1;
  }

  *err |= 0;
  return ans ; 

}

/*
 *  Read a single Real from the file, with the format
 *     tag 10
 */ 
Real  load_single_f(char tag[], FILE *fp, char *buf, int *err )
{

  Real ans ; 

  if( fscanf(fp,"%s",buf) != 1 )
  {
    printf("ERROR in reading required tag %s\n",tag); 
    *err |= 1;
    return 0;
  }

  if( strcmp(buf,tag ) != 0    )
  {
    *err |= 2;
    return 0;
  }


  if( fscanf(fp,"%lfHELP",&ans) != 1 )
  {
    printf("ERROR reading number after the tag \"%s\" \n",tag); 
    *err |= 1;
    return 0;
  }

  *err |= 0;
  return ans ; 

}

void read_input_param(
  char filename[80],           /* Name of parameter file */
  three_list *threept,         /* List of three-pt functions */
  two_list *twopt_recoil,      /* List of recoil meson propagators */
  two_list *twopt_sequential   /* List of incoming B meson propagators */
)
{
  FILE *fp,*fp_filelist;
  char buf[128] ; 
  char filelistfile[128];
  int quit,got,n;
  int err;
  
  /*** open the file containing the selection parameters  ******/
  if( (fp = fopen(filename ,"r")) == NULL )
    {
      printf("ERROR:read_input_param: Could not open the file %s\n",filename);
      exit(1);
    }
  
  if(fscanf(fp,"%s",buf) != 1)
    {
      printf("ERROR:read_input_param: Expecting data\n");
      exit(1);
    }

  if(strcmp(buf,"three_point_select") != 0)
    {
      printf("ERROR:read_input_param: got %s but wanted 'three_point_select'\n",buf);
      exit(1);
    }

  /* Get propagation direction */

  
  if(fscanf(fp,"%s",buf) != 1)
    {
      printf("ERROR:read_input_param: Expecting data\n");
      exit(1);
    }

  if(strcmp(buf,"forward") == 0)
    threept->forwback = FORWARD;
  else if (strcmp(buf,"backward") == 0)
    threept->forwback = BACKWARD;
  else if (strcmp(buf,"fold") == 0)
    threept->forwback = FOLD;
  else
    {
      printf("ERROR:read_input_param: got %s but wanted 'forward', 'backward' or 'fold'\n",buf);
      exit(1);
    }

  /* Read list of three_point_select parameters */
  err = 0;
  quit = 0;
  n = 0;
  
  while(!quit)
    {
      
      threept->select[n].spect  = load_single_i("SP",fp,buf,&err);
      if(err == 1)exit(1);
      if(err == 2)
	{
	  if(strcmp(buf,"two_point_recoil_select") == 0)break;
	  printf("ERROR:read_input_param: Unexpected tag %s\n",buf);
	  exit(1);
	}
      threept->select[n].zonked = load_single_i("ZK",fp,buf,&err);
      threept->select[n].seq    = load_single_i("SQ",fp,buf,&err);
      threept->select[n].q      = load_single_i("Q", fp,buf,&err);
      threept->select[n].p      = load_single_i("P", fp,buf,&err);
      threept->select[n].oper   = load_single_i("OP",fp,buf,&err);
      threept->select[n].copy   = load_single_i("CP",fp,buf,&err);
      threept->select[n].wt     = load_single_f("WT",fp,buf,&err);
      if(err == 1)exit(1);
      if(err == 2)
	{
	  printf("ERROR:read_input_param: Unexpected tag in three_point_select\n");
	  exit(1);
	}
      
      /**      printf("SP %d ZK %d SQ %d Q %d P %d OP %d CP %d\n",
	threept->select[n].spect,threept->select[n].zonked,
	threept->select[n].seq,threept->select[n].q,
	threept->select[n].p,threept->select[n].oper,
	threept->select[n].copy); **/
      n++;
      if(n >= MAX_SELECT_PER_FILE){
        printf("read_input_param: WARNING select array full at %d\n",n);
        break;
      }
    }
  
  threept->nselect = n;

  /* Get propagation direction for two_point_recoil */

  if(fscanf(fp,"%s",buf) != 1)
    {
      printf("ERROR:read_input_param: Expecting data\n");
      exit(1);
    }

  if(strcmp(buf,"forward") == 0)
    twopt_recoil->forwback = FORWARD;
  else if (strcmp(buf,"backward") == 0)
    twopt_recoil->forwback = BACKWARD;
  else if (strcmp(buf,"fold") == 0)
    twopt_recoil->forwback = FOLD;
  else
    {
      printf("ERROR:read_input_param: Expecting 'forward', 'backward' or 'fold'\n");
      exit(1);
    }

  /* Read list of two_point_recoil_select parameters */
  err = 0;
  quit = 0;
  n = 0;

  while(!quit)
    {
      
      twopt_recoil->select[n].spect   = load_single_i("SP",fp,buf,&err);
      if(err == 1)exit(1);
      if(err == 2)
	{
	  if(strcmp(buf,"two_point_sequential_select") == 0)break;
	  printf("ERROR:read_input_param: Unexpected tag %s\n",buf);
	  exit(1);
	}
      twopt_recoil->select[n].other   = load_single_i("ZK",fp,buf,&err);
      twopt_recoil->select[n].mom     = load_single_i("K", fp,buf,&err);
      twopt_recoil->select[n].oper    = load_single_i("OP",fp,buf,&err);
      twopt_recoil->select[n].copy    = load_single_i("CP",fp,buf,&err);
      twopt_recoil->select[n].wt      = load_single_f("WT",fp,buf,&err);
      if(err == 1)exit(1);
      if(err == 2)
	{
	  printf("ERROR:read_input_param: Unexpected tag in two_point_recoil_select\n");
	  exit(1);
	}

      n++;
      if(n >= MAX_SELECT_PER_FILE){
        printf("read_input_param: WARNING select array full at %d\n",n);
        break;
      }
    }

  twopt_recoil->nselect = n;
      
  /* Get propagation direction for two_point_sequential */

  if(fscanf(fp,"%s",buf) != 1)
    {
      printf("ERROR:read_input_param: Expecting data\n");
      exit(1);
    }

  if(strcmp(buf,"forward") == 0)
    twopt_sequential->forwback = FORWARD;
  else if (strcmp(buf,"backward") == 0)
    twopt_sequential->forwback = BACKWARD;
  else if (strcmp(buf,"fold") == 0)
    twopt_sequential->forwback = FOLD;
  else
    {
      printf("ERROR:read_input_param: Expecting 'forward', 'backward' or 'fold'\n");
      exit(1);
    }

  /* Read list of two_point_sequential_select parameters */
  err = 0;
  quit = 0;
  n = 0;

  while(!quit)
    {
      twopt_sequential->select[n].spect   = load_single_i("SP",fp,buf,&err);
      if(err == 1)exit(1);
      if(err == 2)
	{
	  if(strcmp(buf,"filelist") == 0)break;
	  printf("ERROR:read_input_param: wanted 'filelist' but got %s\n",buf);
	  exit(1);
	}
      twopt_sequential->select[n].other   = load_single_i("SQ",fp,buf,&err);
      twopt_sequential->select[n].mom     = load_single_i("P", fp,buf,&err);
      twopt_sequential->select[n].oper    = load_single_i("OP",fp,buf,&err);
      twopt_sequential->select[n].copy    = load_single_i("CP",fp,buf,&err);
      twopt_sequential->select[n].wt      = load_single_f("WT",fp,buf,&err);
      if(err == 1)exit(1);
      if(err == 2)
	{
	  printf("ERROR:read_input_param: Unexpected tag in two_point_sequential_select\n");
	  exit(1);
	}

      n++;
      if(n >= MAX_SELECT_PER_FILE){
        printf("read_input_param: WARNING select array full at %d\n",n);
        break;
      }
    }
  
  twopt_sequential->nselect = n;

  /* Get name of file containing file list */

  if(fscanf(fp,"%s",filelistfile) != 1)
    {
      printf("ERROR:read_input_param: Expecting name of filelist file\n");
      exit(1);
    }

  /* open the file containing the list of files */

  if( (fp_filelist = fopen(filelistfile ,"r")) == NULL )
    {
      printf("ERROR:read_input_param: Could not open the file %s\n",filelistfile);
      exit(1);
    }
  
  /* Read list of three and two point files */

  quit = 0;
  n = 0;

  while(!quit)
    {
      got = fscanf(fp_filelist,"%s%s%s",
		   threept->filename[n],
		   twopt_recoil->filename[n],
		   twopt_sequential->filename[n]);

      if(got<=0)break;
      if(got < 3)
	{
	  printf("ERROR:read_input_param: Expected 3 file names at n = %d\n",n);
	  exit(1);
	}

      n++;
      if(n >= MAX_NO_FILE){
        printf("read_input_param: WARNING file list array full at %d\n",n);
        break;
      }
    }

  threept->nfile = n;
  twopt_recoil->nfile = n;
  twopt_sequential->nfile = n;
      

  /*** close the files ****/
  if( fclose(fp_filelist) != 0 )
  {
    printf("ERROR:read_input_param:There was an error during the closing of %s \n",filelistfile);
    exit(1);
  }

  if( fclose(fp) != 0 )
  {
    printf("ERROR:read_input_param:There was an error during the closing of %s \n",filename);
    exit(1);
  }



}






