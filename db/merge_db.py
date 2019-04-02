"""
"""
from __future__ import print_function
import sqlite3 as sql
import time
import argparse

def create_tables(db):
    """Create tables if they do not already exist."""
    db.execute('''
    CREATE TABLE IF NOT EXISTS correlators ( -- correlator names and metadata\n\
       id INTEGER PRIMARY KEY NOT NULL,\n\
       name VARCHAR,\n\
       metadata VARCHAR -- JSON format\n\
       )''')
    db.execute('''
    CREATE TABLE IF NOT EXISTS data ( -- correlator data by configuration\n\
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
       )''')
    db.execute('''CREATE UNIQUE INDEX IF NOT EXISTS idx_name ON correlators (name)''')
    db.execute('''CREATE UNIQUE INDEX IF NOT EXISTS idx_data ON data (correlator_id, series, trajectory, tsrc)''')
    db.commit()
    return

class MemDB:
    "An in-memory database to manage src <-> dst correlator relations"
    def close(self):
        self.db_mem.close()
        return
    def __init__(self,db_dst):
        "initialize in-memory table 'dst' to all correlators in destination db"
        self.db_dst = db_dst
        self.db_mem = sql.connect(':memory:')
        # in memory table 'dst' reflecting destination correlator table
        self.db_mem.execute('''
        CREATE TABLE dst (\n\
          id INTEGER PRIMARY KEY NOT NULL,\n\
          name VARCHAR\n\
          )''')
        self.db_mem.execute('''CREATE UNIQUE INDEX IF NOT EXISTS idx_dst ON dst (name)''')
        # select from db_dst 'correlators' table and insert into in-memory table 'dst'
        rows = self.db_dst.execute('''SELECT id, name FROM correlators''')
        self.db_mem.executemany('''INSERT INTO dst(id,name) VALUES(?1,?2)''',rows)
        self.db_mem.commit()
        return
    def update_correlators(self,db_src):
        "Add correlators from source to destination db. Map correlator_ids to maintain uniqueness."
        # for each new source db, drop then rebuild 'src' table and index
        self.db_mem.execute('''DROP INDEX IF EXISTS idx_src''')
        self.db_mem.execute('''DROP TABLE IF EXISTS src''')
        self.db_mem.execute('''
        CREATE TABLE src (\n\
          id INTEGER PRIMARY KEY NOT NULL,\n\
          name VARCHAR\n\
          )''')
        self.db_mem.execute('''CREATE UNIQUE INDEX idx_src ON src (name)''')
        # query source db and insert results into in-memory table 'src'
        rows = db_src.execute('''SELECT id,name FROM correlators''')
        self.db_mem.executemany('''INSERT INTO src(id,name) VALUES (?1,?2)''',rows)
        # drop then rebuild in-memory join table 'src_dst' relating id_src with id_dst
        self.db_mem.execute('''DROP TABLE IF EXISTS src_dst''')
        self.db_mem.execute('''
        CREATE TABLE src_dst AS
          SELECT s.id AS src_id, d.id AS dst_id FROM src s
          LEFT JOIN dst d ON s.name = d.name
        ''')
        # select source id, name from 'name_ids' where destination id is NULL
        rows = self.db_mem.execute('''SELECT s.name, r.src_id FROM src s INNER JOIN src_dst r ON s.id = r.src_id WHERE r.dst_id IS NULL''')
        mcur = self.db_mem.cursor()
        for row in rows:
            # add missing name to 'dst' (db will assign unique id)
            name = row[0]
            src_id = row[1]
            mcur.execute('''INSERT INTO dst(name) VALUES (?1)''',(name,))
            dst_id = mcur.execute('''SELECT id FROM dst WHERE name = ?''',(name,)).fetchone()[0]
            # fetch metadata
            metadata = db_src.execute('''SELECT metadata FROM correlators WHERE id = ?''',(src_id,)).fetchone()[0]
            # insert dst_db
            self.db_dst.execute('''INSERT INTO correlators(id,name,metadata) VALUES (?1,?2,?3)''',(dst_id,name,metadata))
            # update src_dst
            mcur.execute('''UPDATE src_dst SET dst_id = ?1 WHERE src_id = ?2''',(dst_id,src_id))
            pass
        self.db_mem.commit()
        self.db_dst.commit()
        return
    def update_data(self,db_src):
        "copy data table rows from source to destination, replacing any matching series,trajectory,tsrc entries"
        get = '''SELECT series,trajectory,tsrc,jobid,timestamp,c_re,c_im FROM data WHERE correlator_id = ? ORDER BY series,trajectory,tsrc'''
        rep = '''REPLACE INTO data (correlator_id,series,trajectory,tsrc,jobid,timestamp,c_re,c_im) VALUES (?1,?2,?3,?4,?5,?6,?7,?8)'''
        src_dst_list = self.db_mem.execute('''SELECT * FROM src_dst''')
        for src_id, dst_id in src_dst_list:
            rows = db_src.execute(get,(src_id,))
            for row in rows:
                out = [ dst_id, ]
                out.extend(row)
                self.db_dst.execute(rep,out)
                pass
            pass
        self.db_dst.commit()
        return
    pass

argparser = argparse.ArgumentParser(description='Merge sqlite3 database files.')
argparser.add_argument('-o','--output',action='store',default='sample-merged.db',help='Merged DB output file')
argparser.add_argument('dbList',action='store',nargs='+',help='List of input DBs')
args = argparser.parse_args()

dst_dbName = args.output

tstart = time.time()

# open destination merge db
dst_db = sql.connect(dst_dbName)

# create destination tables
create_tables(dst_db)

# in-memory db to coordinate merges
corrs = MemDB(dst_db)

for src_dbName in args.dbList:

    src_db = sql.connect(src_dbName)

    # merge correlator tables
    corrs.update_correlators(src_db)

    # merge correlator data
    corrs.update_data(src_db)

    src_db.close()
    pass
#

corrs.close()

dst_db.close()

elapsed = time.time()-tstart
print('merge time {:.3f} sec, {:.3f} sec/file'.format(elapsed,elapsed/len(args.dbList)))
