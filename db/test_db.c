#include "db.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>
#include "Wtimer.h"

unsigned long this_pid;

long unif(int a, int b)
{
  const double rmax = 2147483648. - 1;
  double f = a + (b-a)/rmax * lrand48();
  return(lround(f));
}

void set_vec_approx_normal(double *vec, int len)
{
  int j,k;
  for(j=0; j<len; ++j) vec[j] = drand48();
  for(k=0; k<11; ++k)
    for(j=0; j<12; ++j) vec[j] += drand48();
  for(j=0; j<12; ++j) vec[j] -= 6.0;
}

double model(int t, int nt, double E, double A)
{
  double f = exp(A-E*t) + exp(A-E*(nt-t));
  return(f);
}

void fake_data(int chan, int corr, int is_signal, char *json, int jlen, double *vec, int nt)
{
  double E = 0.5;
  double A = 1.0;
  double sigma = 0.02; // relative to unity
  double f = 0.1;
  double g = 1.25;
  switch(corr)
    {
    case 0:
      break;
    case 1:
    case 2:
      A *= f; // single smear
      sigma *= g;
      break;
    case 3:
      A *= f*f; // double smear
      sigma *= g*g;
      break;
    }
  if( !is_signal )
    {
      A *= 0.1;
      sigma *= 1.5;
    }
  if(!chan)
    {
      E += 0.1;
      A *= 3.0;
    }

  set_vec_approx_normal(vec,nt);

  int t;
  // noise
  for(t=0; t<nt; ++t) vec[t] = sigma * vec[t] * model(t,nt,E,A);
  // signal
  if(is_signal) for(t=0; t<nt; ++t) vec[t] += model(t,nt,E,A);
  // generate JSON array
  json[0] = '\0'; // reset string
  strcat(json,"[");
  for(t=0; t<nt; ++t)
    {
      int len = strlen(json);
      snprintf(json+len,jlen-1-len,"%.6e",vec[t]);
      if(t < nt-1) strcat(json,",");
    }
  strcat(json,"]");
}


int main(int argc, char *argv[])
{
  this_pid = getpid();

  int db_has_datetime = 1; // db name includes datetime
  int nt = 64; // lattice time extent

  sqlite3 *db;
  int ret;

  // timestamp: seconds since the epoch
  time_t utime = time(NULL);

  // build dbName
  const int max_dbName = 48;
  char dbName[max_dbName+1]; dbName[0]='\0';
  if(db_has_datetime)
    {
      // add current time to name
      struct tm *timeinfo = localtime(&utime);
      int off = strlen(dbName);
      strftime(dbName+off,max_dbName-off,"sample-%Y-%m-%d-%H-%M-%S.db",timeinfo);
    }
  else
    {
      strcat(dbName,"sample.db");
    }

  srand48(utime);

  // connect to db
  printf("%lu: connecting to %s\n",this_pid,dbName);
  ret = db_connect(&db,dbName);

  int wait_ms = 5000;

  // init tables if they do not exist
  Wtimer_t ct_timer; Wtimer_start(&ct_timer);
  ret = db_init_tables(db,2*wait_ms);
  printf("%lu ELAPSED-create-tables: %.3f\n",this_pid,Wtimer_stop(&ct_timer));
  if( ret != SQLITE_OK)
    {
      printf("%lu ERROR: cannot create tables. Quitting...\n",this_pid);
      exit(ret);
    }

  const int ncorr = 4;
  char *pi_name[] = { "pi_d_d_p000", "pi_d_1S_p000", "pi_1S_d_p000", "pi_1S_1S_p000", };
  char *ro_name[] = { "ro_d_d_p000", "ro_d_1S_p000", "ro_1S_d_p000", "ro_1S_1S_p000", };
  char *ckey[ncorr];
  const int do_pi = lrand48() % 2;
  if( do_pi ){
    int j;
    for(j=0; j<ncorr; ++j) ckey[j] = pi_name[j];
  } else {
    int j;
    for(j=0; j<ncorr; ++j) ckey[j] = ro_name[j];
  }

  // random jobid
  char jobid[16];
  snprintf(jobid,15,"fnal.%d",unif(100000,999999));

  // random series a-z
  char series[2]; series[1] = '\0';
  series[0] = 'a' + unif(0,25); 

  // random range of trajectories
  int traj_delta = 10;
  int traj_start = traj_delta * unif(2,20);
  int ntraj = 8;
  int traj_end = traj_start + ntraj * traj_delta;

  // vectors for data
  double *vec = static_cast(double*,malloc(nt*sizeof(double)));
  const int precision = 6;
  int jlen = (nt * (precision + 7) + 4)*sizeof(char);
  char *json_r = static_cast(char*,malloc(jlen));
  char *json_i = static_cast(char*,malloc(jlen));

  for(int c=0; c<ncorr; ++c)
    {

      // empty correltaor object
      db_correlator corr; corr.id = 0; corr.name = NULL; corr.metadata = NULL;

      // lookup ckey in db
      ret = db_query_correlator_by_name(db, ckey[c], &corr);
      if( corr.id < 1 )
	{
	  // ckey not found in db, so insert as a new correlator
	  corr.name = ckey[c];
	  corr.metadata = "{\"mass_quark\":\"0.0031\",\"mass_antiquark\":\"0.0031\"}";
	  ret = db_begin_transaction(db,wait_ms);
	  if(ret == SQLITE_OK) ret = db_insert_correlator(db,&corr);
	  if(ret == SQLITE_OK)
	    {
	      db_commit_transaction(db);
	    }
	  else
	    {
	      db_rollback_transaction(db);
	      printf("%lu ERROR: transaction failed inserting correlator c=%d. Quitting...\n",this_pid,c);
	      exit(ret);
	    }
	}
      //printf("%s id=%d\n",corr.name,corr.id);


      // insert many data items
      db_data data; data.id=0;data.correlator_id=0;data.series=NULL;data.trajectory=0;data.tsrc=0;data.jobid=NULL;data.timestamp=0;data.c_re=NULL;data.c_im=NULL;

      data.correlator_id = corr.id;
      data.jobid = jobid;
      data.timestamp = utime;
      data.series = series;

      // begin transaction
      Wtimer_t bt_timer; Wtimer_start(&bt_timer);
      ret = db_begin_transaction(db,wait_ms);
      printf("%lu ELAPSED-begin-transaction: %.3f\n",this_pid,Wtimer_stop(&bt_timer));
      if(ret == SQLITE_OK)
	{
	  Wtimer_t di_timer; Wtimer_start(&di_timer);
	  for(int traj=traj_start; traj<traj_end; traj+=traj_delta)
	    {
	      data.trajectory = traj;
	      for(int tsrc=0; tsrc<nt; tsrc+=16)
		{
		  data.tsrc = tsrc;
		  fake_data(do_pi,c,1,json_r,jlen,vec,nt);
		  fake_data(do_pi,c,0,json_i,jlen,vec,nt);
		  data.c_re = json_r;
		  data.c_im = json_i;
		  ret = db_update_data(db, &data);
		  if(ret != SQLITE_OK)
		    break;
		}
	    }
	  printf("%lu ELAPSED-insert-data: %.3f\n",this_pid,Wtimer_stop(&di_timer));
	}
      else
	{
	  printf("%lu WARNING: could not begin transaction inserting data c=%d. Skipping...\n",this_pid,c);
	  break;
	}
      if(ret == SQLITE_OK)
	{
	  // end transaction
	  db_commit_transaction(db);
	}
      else
	{
	  // rollback on error
	  db_rollback_transaction(db);
	  printf("%lu ERROR: transaction failed inserting data c=%d\n",this_pid,c);
	}
    }

  // disconnect from db
  ret = db_disconnect(db);

  // temp space
  free(json_i);
  free(json_r);
  free(vec);

  return(0);

}

