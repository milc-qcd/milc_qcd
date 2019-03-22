#include "db.h"
#include <stdio.h>
#include <time.h>

int main()
{
  const char *dbName = "sample.db";
  sqlite3 *db;
  int ret;

  // timestamp: seconds since the epoch
  time_t utime = time(NULL);

  // connect to db
  ret = db_connect(&db,dbName);

  // init tables if they do not exist
  ret = db_init_tables(db);

  // begin transaction
  int wait_ms = 50;
  db_begin_transaction(db,wait_ms);

  const int ncorr = 4;
  char* ckey[] = { "pi_d_d_p000", "pi_d_1S_p000", "pi_1S_d_p000", "pi_1S_1S_p000", };

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
	  ret = db_insert_correlator(db,&corr);
	}
      printf("%s id=%d\n",corr.name,corr.id);


      // insert many data items
      db_data data; data.id=0;data.correlator_id=0;data.series=NULL;data.trajectory=0;data.tsrc=0;data.jobid=NULL;data.timestamp=0;data.c_re=NULL;data.c_im=NULL;

      data.correlator_id = corr.id;
      data.jobid = "fnal.1234548";
      data.timestamp = utime;

      data.series = "a";
      for(int traj=100; traj<160; traj+=10)
	{
	  data.trajectory = traj;
	  for(int tsrc=0; tsrc<64; tsrc+=16)
	    {
	      data.tsrc = tsrc;
	      data.c_re = "[1.00e+00,5.00e-01,2.50e-01,1.25e-01,6.25e-02]";
	      data.c_im = "[-1.00e-15,1.25e-15,1.00e-15,-2.00e-15,1.50e-15]";
	      db_update_data(db, &data);
	    }
	}
    }
  // end transaction
  //db_rollback_transaction(db); exit(0);
  db_commit_transaction(db);

  // disconnect from db
  ret = db_disconnect(db);

  return(0);

}

