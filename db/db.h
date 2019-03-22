#ifndef __DB_H__
#define  __DB_H__

#include <stdlib.h>
#include <sqlite3.h>

#ifndef static_cast
#define static_cast(t,v) ((t)(v))
#endif

typedef struct db_correlator_t {
  int id;
  char *name;      // owned by application
  char *metadata;  // owned by application
} db_correlator;

typedef struct db_data_t {
  int id;
  int correlator_id;
  char *series;
  int trajectory;
  int tsrc; 
  char *jobid;
  unsigned long timestamp;
  char *c_re;
  char *c_im;
} db_data;

// prototypes

int db_connect(sqlite3 **db, const char* dbName);

int db_disconnect(sqlite3 *db);

int db_init_tables(sqlite3* db);

int db_begin_transaction(sqlite3* db, int wait_if_busy_ms);

int db_commit_transaction(sqlite3* db);

int db_rollback_transaction(sqlite3* db);

int db_correlator_free(db_correlator *corr);

int db_query_correlator_by_name(sqlite3 *db, const char *name, db_correlator *corr);

int db_insert_correlator(sqlite3 *db, db_correlator *corr);

int db_update_data(sqlite3 *db, db_data *data);

#endif
