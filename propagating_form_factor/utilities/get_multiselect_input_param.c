/*
 *   Read the parameters 
 *
 *
 *
 */


#include "read_hl_form.h"
#include "prop_form_utilities.h"


/*
 *  Read a single integer from the file, with the format
 *     tag 10
 */ 
int  load_single_i(char tag[], FILE *fp, char *buf, int *err )
{
  int ans ; 


  if( fscanf(fp,"%s",buf) != 1 )
  {
    /* Don't complain if only the EOF */
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

void read_multiselect_twopt_param(
  int *nfile,             /* Number of files to write */
  twopt_oneselect *fp  /* two point selection parameter sets */
)
{
  char buf[128] ; 
  char filename[128];
  char oldfilename[128];
  int quit,n;
  int err;
  int first;
  int spect,other,mom,oper,copy;
  Real wt;

  /* Read list of two_point_recoil_select parameters */
  err = 0;
  quit = 0;
  *nfile = 0;
  n = 0;
  first = 1;

  while(!quit)
    {
      
      spect   = load_single_i("SP",stdin,buf,&err);
      if(err == 1)break;
      if(n>=MAX_SELECT_PER_FILE)
	{
	  printf("read_multiselect_twopt_param:  Selection array overflow\n");
	  exit(1);
	}
      if(*nfile >= MAX_NO_FILE)
	{
	  printf("read_multiselect_twopt_param:  Too many output files\n");
	  exit(1);
	}

      other   = load_single_i("ZK",stdin,buf,&err);
      mom     = load_single_i("K", stdin,buf,&err);
      oper    = load_single_i("OP",stdin,buf,&err);
      copy    = load_single_i("CP",stdin,buf,&err);
      wt      = load_single_f("WT",stdin,buf,&err);
      scanf("%s",filename);
      if(err == 1)exit(1);
      if(err == 2)
	{
	  printf("read_multiselect_twopt_param: Unexpected tag in two_point_select list\n");
	  exit(1);
	}

      /* If filename changed, start new selection list */
      if(!first && strcmp(oldfilename,filename) != 0)
	{
	  fp[*nfile].nselect = n;
	  n = 0;
	  (*nfile)++;
	}

      strcpy(oldfilename,filename);
      first = 0;

      fp[*nfile].select[n].spect  = spect;
      fp[*nfile].select[n].other  = other;
      fp[*nfile].select[n].mom    = mom;
      fp[*nfile].select[n].oper   = oper;
      fp[*nfile].select[n].copy   = copy;
      fp[*nfile].select[n].wt     = wt;
      strcpy(fp[*nfile].filename,filename);

      n++;
    }
  
  fp[*nfile].nselect = n;
  (*nfile)++;
}



void read_multiselect_form_param(
  int *nfile,             /* Number of files to write */
  threept_oneselect *fp   /* Selection parameter set for each file */
)
{
  char buf[128] ; 
  char filename[128];
  char oldfilename[128];
  int quit,n;
  int err;
  int first;
  int spect,zonked,seq,p,q,oper,copy;
  Real wt;
  

  /* Read list of three_point_select parameters */
  err = 0;
  quit = 0;
  *nfile = 0;
  n = 0;
  first = 1;
  
  while(!quit)
    {
      
      spect  = load_single_i("SP",stdin,buf,&err);
      if(err == 1)break;
      if(n>=MAX_SELECT_PER_FILE)
	{
	  printf("read_multiselect_form_param:  Selection array overflow\n");
	  exit(1);
	}
      if(*nfile >= MAX_NO_FILE)
	{
	  printf("read_multiselect_form_param:  Too many output files\n");
	  exit(1);
	}

      zonked = load_single_i("ZK",stdin,buf,&err);
      seq    = load_single_i("SQ",stdin,buf,&err);
      q      = load_single_i("Q", stdin,buf,&err);
      p      = load_single_i("P", stdin,buf,&err);
      oper   = load_single_i("OP",stdin,buf,&err);
      copy   = load_single_i("CP",stdin,buf,&err);
      wt     = load_single_f("WT",stdin,buf,&err);
      scanf("%s",filename);
      if(err == 1)exit(1);
      if(err == 2)
	{
	  printf("read_multiselect_form_param: Unexpected tag in three_point_select\n");
	  exit(1);
	}
      
      /* If filename changed, start new selection list */
      if(!first && strcmp(oldfilename,filename) != 0)
	{
	  fp[*nfile].nselect = n;
	  n = 0;
	  (*nfile)++;
	}

      strcpy(oldfilename,filename);
      first = 0;
  
/**      printf("SP %d ZK %d SQ %d Q %d P %d OP %d CP %d WT %f %s\n",
	     spect,zonked,seq,q,p,oper,copy,wt,filename); **/

      fp[*nfile].select[n].spect  = spect;
      fp[*nfile].select[n].zonked = zonked;
      fp[*nfile].select[n].seq    = seq;
      fp[*nfile].select[n].q      = q;
      fp[*nfile].select[n].p      = p;
      fp[*nfile].select[n].oper   = oper;
      fp[*nfile].select[n].copy   = copy;
      fp[*nfile].select[n].wt     = wt;
      strcpy(fp[*nfile].filename,filename);
       
      n++;
    }

  fp[*nfile].nselect = n;
  (*nfile)++;
}

