#include "db.h"
#include <malloc.h>
#include <string.h>
#include <stdio.h>

int db_connect(sqlite3 **db, const char* dbName)
{
  int rc;
  rc = sqlite3_open(dbName,db);
  if( rc != SQLITE_OK )
    {
      fprintf(stderr,"db_connect:sqlite3_open ERROR %s\n",sqlite3_errmsg(*db));
      sqlite3_close(*db);
      return(rc);
    }
  return(rc);
}

int db_disconnect(sqlite3 *db)
{
  int rc;
  rc = sqlite3_close(db);
  if( rc != SQLITE_OK )
    {
      fprintf(stderr,"db_disconnect:sqlite3_close ERROR %s\n",sqlite3_errmsg(db));
    }
  return(rc);
}

int db_begin_transaction(sqlite3* db, int wait_if_busy_ms)
{
  int rc;
  char *err;
  // set busy timeout cumulatively waiting at least ms milliseconds
  rc = sqlite3_busy_timeout(db, wait_if_busy_ms);

  const char *sql = "BEGIN IMMEDIATE TRANSACTION";
  rc = sqlite3_exec(db, sql, NULL, NULL, &err);
  if( rc != SQLITE_OK )
    {
      fprintf(stderr, "db_begin_transaction: ERROR %s\n", err);
      sqlite3_free(err);
    }
return(rc);
}

int db_commit_transaction(sqlite3* db)
{
  int rc;
  char *err;
  const char *sql = "COMMIT TRANSACTION";
  rc = sqlite3_exec(db, sql, NULL, NULL, &err);
  if( rc != SQLITE_OK )
    {
      fprintf(stderr, "db_commit_transaction: ERROR %s\n", err);
      sqlite3_free(err);
    }
  sqlite3_busy_timeout(db, 0); // reset timeout
return(rc);
}

int db_rollback_transaction(sqlite3* db)
{
  int rc;
  char *err;
  const char *sql = "ROLLBACK TRANSACTION";
  rc = sqlite3_exec(db, sql, NULL, NULL, &err);
  if( rc != SQLITE_OK )
    {
      fprintf(stderr, "db_rollback_transaction: ERROR %s\n", err);
      sqlite3_free(err);
    }
  sqlite3_busy_timeout(db, 0); // reset timeout
return(rc);
}

int db_init_tables(sqlite3* db, int wait_if_busy_ms)
{
  int rc;
  char *err;
  
  rc = db_begin_transaction(db,wait_if_busy_ms);
  if( rc != SQLITE_OK )
    return(rc);

  const char *correlators = 
    "CREATE TABLE IF NOT EXISTS correlators ( -- correlator names and metadata\n\
       id INTEGER PRIMARY KEY NOT NULL,\n\
       name VARCHAR,\n\
       metadata VARCHAR -- JSON format\n\
       )";
  rc = sqlite3_exec(db, correlators, NULL, NULL, &err);
  if( rc != SQLITE_OK )
    {
      fprintf(stderr, "db_init_tables: %s\n", correlators);
      fprintf(stderr, "db_init_tables: ERROR %s\n", err);
      sqlite3_free(err);
      rc = db_rollback_transaction(db);
      return(rc);
    }

  const char *data =
    "CREATE TABLE IF NOT EXISTS data ( -- correlator data by configuration\n\
       id INTEGER PRIMARY KEY NOT NULL,\n\
       correlator_id INTEGER NOT NULL,\n\
       series VARCHAR,\n\
       trajectory INTEGER,\n\
       tsrc INTEGER,\n\
       jobid VARCHAR,\n\
       timestamp INTEGER, -- unix seconds since epoch\n\
       c_re VARCHAR, -- JSON encoded vector, real part\n\
       c_im VARCHAR, -- JSON imag part\n\
       FOREIGN KEY(correlator_id) REFERENCES correlators (id)\n\
       )";
  rc = sqlite3_exec(db, data, NULL, NULL, &err);
  if( rc != SQLITE_OK)
    {
      fprintf(stderr, "db_init_tables: %s\n", data);
      fprintf(stderr, "db_init_tables: ERROR %s\n", err);
      sqlite3_free(err);
      rc = db_rollback_transaction(db);
      return(rc);
    }

  // uniqueness constraint on a correlator name
  const char *cindx = "CREATE UNIQUE INDEX IF NOT EXISTS idx_name ON correlators (name)";
  rc = sqlite3_exec(db, cindx, NULL, NULL, &err);
  if( rc != SQLITE_OK)
    {
      fprintf(stderr, "db_init_tables: %s\n", cindx);
      fprintf(stderr, "db_init_tables: ERROR %s\n", err);
      sqlite3_free(err);
      rc = db_rollback_transaction(db);
      return(rc);
    }

  // uniqueness constraint on each datum
  const char *dindx = "CREATE UNIQUE INDEX IF NOT EXISTS idx_data ON data (correlator_id, series, trajectory, tsrc)";
  rc = sqlite3_exec(db, dindx, NULL, NULL, &err);
  if( rc != SQLITE_OK)
    {
      fprintf(stderr, "db_init_tables: %s\n", dindx);
      fprintf(stderr, "db_init_tables: ERROR %s\n", err);
      sqlite3_free(err);
      rc = db_rollback_transaction(db);
      return(rc);
    }

  rc = db_commit_transaction(db);
  return(rc);
}

int db_correlator_free(db_correlator *corr)
{
  if(corr->metadata) free(corr->metadata);
  if(corr->name) free(corr->name);
  corr->name = NULL;
  corr->metadata = NULL;
  return(0);
}

int db_query_correlator_by_name(sqlite3 *db, const char *name, db_correlator *corr)
{
  int rc;
  const char *query = "SELECT id,metadata FROM correlators WHERE name = ?1";
  sqlite3_stmt *stmt;
  rc = sqlite3_prepare_v2(db, query, strlen(query)+1, &stmt, NULL);
  if( rc != SQLITE_OK )
    {
      fprintf(stderr, "query: %s\n",query);
      fprintf(stderr, "db_query_correlator_by_name:sqlite3_prepare_v2 ERROR %s\n", sqlite3_errmsg(db));
      return(rc);
    }
  rc = sqlite3_bind_text(stmt, 1, name, -1, SQLITE_TRANSIENT);
  if( rc != SQLITE_OK )
    {
      fprintf(stderr, "db_query_correlator_by_name:sqlite3_bind_text ?1 ERROR %s\n", sqlite3_errmsg(db));
      return(rc);
    }
  rc = sqlite3_step(stmt);
  if( rc == SQLITE_ROW )
    {
      corr->id = sqlite3_column_int(stmt,0);
      corr->name = static_cast(char*,calloc(strlen(name)+1,sizeof(char)));
      strcpy(corr->name,name);
      const char *meta = static_cast(char*,sqlite3_column_text(stmt,1)); // WARNING: cast unsigned char* (UTF-8) to char* (ASCII)
      corr->metadata = static_cast(char*,calloc(strlen(meta)+1,sizeof(unsigned char)));
      strcpy(corr->metadata,meta);
    }
  else if(rc == SQLITE_DONE)
    {
      // not found
      corr->id = -1;
    }
  else
    {
      fprintf(stderr, "db_query_correlator_by_name:sqlite3_step ERROR %s\n", sqlite3_errmsg(db));
      return(rc);
    }
  rc = sqlite3_reset(stmt);
  if( rc != SQLITE_OK )
    {
      fprintf(stderr, "db_query_correlator_by_name:sqlite3_reset %s\n", sqlite3_errmsg(db));
      return(rc);
    }
  rc = sqlite3_finalize(stmt);
  if( rc != SQLITE_OK )
    {
      fprintf(stderr, "db_query_correlator_by_name:sqlite3_finalize %s\n", sqlite3_errmsg(db));
      return(rc);
    }
  return(rc);
}

int db_insert_correlator(sqlite3 *db, db_correlator *corr)
{
  int rc;
  const char *sql = "INSERT OR IGNORE INTO correlators(name,metadata) VALUES(?1,?2)";
  sqlite3_stmt *stmt;
  rc = sqlite3_prepare_v2(db, sql, strlen(sql)+1, &stmt, NULL);
  if( rc != SQLITE_OK)
    {
      fprintf(stderr, "sql: %s\n",sql);
      fprintf(stderr, "db_insert_correlator:sqlite3_prepare_v2 ERROR %s\n", sqlite3_errmsg(db));
      return(rc);
    }
  rc = sqlite3_bind_text(stmt, 1, corr->name, -1, SQLITE_TRANSIENT);
  if( rc != SQLITE_OK)
    {
      fprintf(stderr, "db_insert_correlator:sqlite3_bind_text ?1 ERROR %s\n", sqlite3_errmsg(db));
      return(rc);
    }
  rc = sqlite3_bind_text(stmt, 2, corr->metadata, -1, SQLITE_TRANSIENT);
  if( rc != SQLITE_OK)
    {
      fprintf(stderr, "db_insert_correlator:sqlite3_bind_text ?2 ERROR %s\n", sqlite3_errmsg(db));
      return(rc);
    }
  rc = sqlite3_step(stmt);
  if( rc != SQLITE_DONE )
    {
      fprintf(stderr, "db_insert_correlator:sqlite3_step ERROR %s\n", sqlite3_errmsg(db));
      return(rc);
    }
  rc = sqlite3_reset(stmt);
  if( rc != SQLITE_OK )
    {
      fprintf(stderr, "db_insert_correlator:sqlite3_reset %s\n", sqlite3_errmsg(db));
      return(rc);
    }
  rc = sqlite3_finalize(stmt);
  if( rc != SQLITE_OK )
    {
      fprintf(stderr, "db_insert_correlator:sqlite3_finalize %s\n", sqlite3_errmsg(db));
      return(rc);
    }


  // get the id from the db
  const char *query = "SELECT id FROM correlators WHERE name = ?1";
  rc = sqlite3_prepare_v2(db, query, strlen(query)+1, &stmt, NULL);
  if( rc != SQLITE_OK )
    {
      fprintf(stderr, "query: %s\n",query);
      fprintf(stderr, "db_insert_correlator:sqlite3_prepare_v2 ERROR %s\n", sqlite3_errmsg(db));
      return(rc);
    }
  rc = sqlite3_bind_text(stmt, 1, corr->name, -1, SQLITE_TRANSIENT);
  if( rc != SQLITE_OK )
    {
      fprintf(stderr, "db_insert_correlator:sqlite3_bind_text ?1 ERROR %s\n", sqlite3_errmsg(db));
      return(rc);
    }
  rc = sqlite3_step(stmt);
  if( rc != SQLITE_ROW )
    {
      fprintf(stderr, "db_insert_correlator:sqlite3_step ERROR %s\n", sqlite3_errmsg(db));
      return(rc);
    }
  corr->id = sqlite3_column_int(stmt,0);
  rc = sqlite3_reset(stmt);
  if( rc != SQLITE_OK )
    {
      fprintf(stderr, "db_insert_correlator:sqlite3_reset %s\n", sqlite3_errmsg(db));
      return(rc);
    }
  rc = sqlite3_finalize(stmt);
  if( rc != SQLITE_OK )
    {
      fprintf(stderr, "db_insert_correlator:sqlite3_finalize %s\n", sqlite3_errmsg(db));
      return(rc);
    }
  return(rc);
}

int db_update_data(sqlite3 *db, db_data *data)
{
  int rc;
  char *sql = "REPLACE INTO data (correlator_id,series,trajectory,tsrc,jobid,timestamp,c_re,c_im) VALUES (?1,?2,?3,?4,?5,?6,?7,?8)";
  sqlite3_stmt *stmt;
  rc = sqlite3_prepare_v2(db, sql, strlen(sql)+1, &stmt, NULL);
  if( rc != SQLITE_OK)
    {
      fprintf(stderr, "sql: %s\n",sql);
      fprintf(stderr, "db_update_data:sqlite3_prepare_v2 ERROR %s\n", sqlite3_errmsg(db));
      return(rc);
    }
  rc = sqlite3_bind_int(stmt, 1, data->correlator_id);
  if( rc != SQLITE_OK)
    {
      fprintf(stderr, "db_update_data:sqlite3_bind_int ?1 ERROR %s\n", sqlite3_errmsg(db));
      return(rc);
    }
  rc = sqlite3_bind_text(stmt, 2, data->series, -1, SQLITE_TRANSIENT);
  if( rc != SQLITE_OK)
    {
      fprintf(stderr, "db_update_data:sqlite3_bind_text ?2 ERROR %s\n", sqlite3_errmsg(db));
      return(rc);
    }
  rc = sqlite3_bind_int(stmt, 3, data->trajectory);
  if( rc != SQLITE_OK)
    {
      fprintf(stderr, "db_update_data:sqlite3_bind_int ?3 ERROR %s\n", sqlite3_errmsg(db));
      return(rc);
    }
  rc = sqlite3_bind_int(stmt, 4, data->tsrc);
  if( rc != SQLITE_OK)
    {
      fprintf(stderr, "db_update_data:sqlite3_bind_int ?4 ERROR %s\n", sqlite3_errmsg(db));
      return(rc);
    }
  rc = sqlite3_bind_text(stmt, 5, data->jobid, -1, SQLITE_TRANSIENT);
  if( rc != SQLITE_OK)
    {
      fprintf(stderr, "db_update_data:sqlite3_bind_text ?5 ERROR %s\n", sqlite3_errmsg(db));
      return(rc);
    }
  rc = sqlite3_bind_int64(stmt, 6, data->timestamp);
  if( rc != SQLITE_OK)
    {
      fprintf(stderr, "db_update_data:sqlite3_bind_int ?6 ERROR %s\n", sqlite3_errmsg(db));
      return(rc);
    }
  rc = sqlite3_bind_text(stmt, 7, data->c_re, -1, SQLITE_TRANSIENT);
  if( rc != SQLITE_OK)
    {
      fprintf(stderr, "db_update_data:sqlite3_bind_text ?7 ERROR %s\n", sqlite3_errmsg(db));
      return(rc);
    }
  rc = sqlite3_bind_text(stmt, 8, data->c_im, -1, SQLITE_TRANSIENT);
  if( rc != SQLITE_OK)
    {
      fprintf(stderr, "db_update_data:sqlite3_bind_text ?8 ERROR %s\n", sqlite3_errmsg(db));
      return(rc);
    }
  rc = sqlite3_step(stmt);
  if( rc != SQLITE_DONE )
    {
      fprintf(stderr, "db_update_data:sqlite3_step ERROR %s\n", sqlite3_errmsg(db));
      return(rc);
    }
  rc = sqlite3_reset(stmt);
  if( rc != SQLITE_OK )
    {
      fprintf(stderr, "db_update_data:sqlite3_reset %s\n", sqlite3_errmsg(db));
      return(rc);
    }
  rc = sqlite3_finalize(stmt);
  if( rc != SQLITE_OK )
    {
      fprintf(stderr, "db_update_data:sqlite3_finalize %s\n", sqlite3_errmsg(db));
      return(rc);
    }
  return(rc);
}






