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

void read_select_twopt_param(
  int *nselect,            /* Number of selection sets */
  twopt_select *twoselect  /* two point selection parameter sets */
)
{
  char buf[128] ; 
  int quit,n;
  int err;


  /* Read list of two_point_recoil_select parameters */
  err = 0;
  quit = 0;
  n = 0;

  while(!quit)
    {
      
      twoselect[n].spect   = load_single_i("SP",stdin,buf,&err);
      if(err == 1)break;
      twoselect[n].other   = load_single_i("ZK",stdin,buf,&err);
      twoselect[n].mom     = load_single_i("K", stdin,buf,&err);
      twoselect[n].oper    = load_single_i("OP",stdin,buf,&err);
      twoselect[n].copy    = load_single_i("CP",stdin,buf,&err);
      twoselect[n].wt      = load_single_f("WT",stdin,buf,&err);
      if(err == 1)exit(1);
      if(err == 2)
	{
	  printf("ERROR:read_select_twopt_param: Unexpected tag in two_point_select list\n");
	  exit(1);
	}

      n++;
      if(n >= MAX_SELECT_PER_FILE){
        printf("read_select_twopt_param: WARNING select array full at %d\n",n);
        break;
      }
    }

  *nselect = n;
      
}



void read_select_form_param(
  int *nselect,                 /* Number of selection sets */
  threept_select *threeselect   /* Selection parameter sets */
)
{
  char buf[128] ; 
  int quit,n;
  int err;

  /* Read list of three_point_select parameters */
  err = 0;
  quit = 0;
  n = 0;
  
  while(!quit)
    {
      
      threeselect[n].spect  = load_single_i("SP",stdin,buf,&err);
      if(err == 1)break;
      threeselect[n].zonked = load_single_i("ZK",stdin,buf,&err);
      threeselect[n].seq    = load_single_i("SQ",stdin,buf,&err);
      threeselect[n].q      = load_single_i("Q", stdin,buf,&err);
      threeselect[n].p      = load_single_i("P", stdin,buf,&err);
      threeselect[n].oper   = load_single_i("OP",stdin,buf,&err);
      threeselect[n].copy   = load_single_i("CP",stdin,buf,&err);
      threeselect[n].wt     = load_single_f("WT",stdin,buf,&err);
      if(err == 1)exit(1);
      if(err == 2)
	{
	  printf("ERROR:read_select_form_param: Unexpected tag in three_point_select\n");
	  exit(1);
	}
      
      /**      printf("SP %d ZK %d SQ %d Q %d P %d OP %d CP %d WT %f\n",
	threeselect[n].spect,threeselect[n].zonked,
	threeselect[n].seq,threeselect[n].q,
	threeselect[n].p,threeselect[n].oper,
	threeselect[n].copy,threeselect[n].wt); **/

      n++;
      if(n >= MAX_SELECT_PER_FILE){
        printf("read_select_form_param: WARNING select array full at %d\n",n);
        break;
      }
    }
  
  *nselect = n;

}

