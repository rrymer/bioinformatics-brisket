# -*- coding: utf-8 -*-
"""
Created on Fri Apr  4 16:14:45 2014

@author: Richard Rymer

Uploads MultiHarmony results to a SQL database.  Generates an alignment index
for easier offline analysis, identifies and removes positions of >=90% gap,
generates another index.  Finally, a set of alignment "coordinates" is generated,
which maps the alignment position, to the sequence postion, to the sequence position
for a reference protein.
"""

import MySQLdb as mdb
import glob
import sys
import os
from Bio import AlignIO
import time
import logging

class Logger():
    logger=None
    def myLogger(self,path,name,level):
        print 'initializing logger {} located at {}'.format(name,path)
        if None == self.logger:
            self.logger=logging.getLogger(name)
            self.logger.setLevel(level)
            self.handler=logging.FileHandler('{}index_MH_data_{}.log'.format(path,name))
            formatter=logging.Formatter("%(message)s")
            self.handler.setFormatter(formatter)
            self.logger.addHandler(self.handler)
        return self.logger

def reinitialize_logging(path):
    reload(logging)
    main_loggers = Logger()
    main_logger = main_loggers.myLogger(basedir,'main',logging.DEBUG)
    return main_logger

def connect_to_mysqdb(host, db, un):
    connection = mdb.connect(host=host,db=db,user=un)
    return connection

def create_table(name):
    connection = connect_to_mysqdb('localhost', 'FamilyCpolymerases', 'sequser')
    cursor=connection.cursor()    
    sql=("CREATE TABLE `%s` (" +
      "`NCBI_GI` varchar(12) DEFAULT NULL," +
      "`Organism` varchar(200) DEFAULT NULL,"+
      "`PolGroup` varchar(45) DEFAULT NULL,"+
      "`Domains` varchar(45) DEFAULT NULL,"+
      "`Phylum` varchar(45) DEFAULT NULL,"+
      "`Class` varchar(45) DEFAULT NULL,"+
      "`Length` int(5) DEFAULT NULL,"+
      "`Sequence` longtext,"+
      "`Header` varchar(200) NOT NULL DEFAULT 'No entry',"+
      "`Accession_Number` varchar(8) DEFAULT 'X00000',"+
      "`OrganismID` varchar(12) DEFAULT NULL,"+
      "`Entry_Name` varchar(45) DEFAULT NULL,"+
      "`Lineage` varchar(90) DEFAULT NULL,"+
      "`Aligned_seq` varchar(10000) DEFAULT NULL,"+
      "`gfAligned_seq` varchar(10000) DEFAULT NULL,"+
      "`Aligned_seq_index` varchar(10000) DEFAULT '',"+
      "`Aligned_seq_coords` longtext,"+
      "`Indexing` int(11) DEFAULT '1',"+
      "`gfAligned_seq_index` varchar(10000) DEFAULT NULL,"+
      "`gfAligned_seq_coords` longtext,"+
      "`PolProfile` varchar(45) DEFAULT NULL,"+
      "PRIMARY KEY (`Header`)"+
      ") ENGINE=InnoDB DEFAULT CHARSET=latin1;")
    sql = sql % name
    cursor.execute(sql)
    connection.commit()
    connection.close()
   
def upload(dest_t,algnmt,pro):
    connection = connect_to_mysqdb('localhost', 'FamilyCpolymerases', 'sequser')
    cursor = connection.cursor()
    name = file_dict[algnmt][0]
    ref_t = file_dict[algnmt][1]
    algnmt = AlignIO.read(algnmt,'fasta')
    sql = "insert ignore into `{}` ({},{},{}) values('{}','{}','{}')"
    usql= "update `{}` a, `{}` b set a.{}=b.{} where a.{}=b.{}"
    for record in algnmt:
        if 'uref' in str(record.id):
            print 'inserting sequence {} into {}'.format(record.id,dest_t)
            main_logger.debug("insert statement:{}"
            .format(sql.format(dest_t,key_field, ID_field, algnmt_field,name,str(record.id)[-6:],str(record.seq))))
            try:
                cursor.execute(sql.format(dest_t,key_field, ID_field, algnmt_field,name,str(record.id)[-6:],str(record.seq)))
                connection.commit()
            except Exception, e:
                print e
    time.sleep(1)
    try:
        cursor.execute(usql.format(dest_t, ref_t, seq_field, seq_field, ID_field, ID_field))
        connection.commit()
    except Exception, e:
        print e
    time.sleep(1)
    connection.close()
    return index_alignment(dest_t,name)

def index_alignment(table,algnmt):
    connection = connect_to_mysqdb('localhost', 'FamilyCpolymerases', 'sequser')
    cursor = connection.cursor()
    sql="select length({}) from {} where {}='{}'"
    main_logger.debug('select alignment length sql query: {}'.format(sql.format(algnmt_field,table,key_field,algnmt)))    
    cursor.execute(sql.format(algnmt_field,table,key_field,algnmt))
    n=cursor.fetchall()[0]
    n=str(n[0])
    n=float(n)
    main_logger.info('alignment total length: {}'.format(str(n)))
    i=1
    print 'calculating alignment index for {}'.format(algnmt)
    sql1="update {} set {} = if(substring({},{},1) not like '-',concat({},' ',{}),concat({},' ','-')) where {}='{}'"
    sql2="update {} set {}=if(substring({},{},1) not like '-',{}+1,{}) where {}='{}'"
    main_logger.info=("update index sql statement: {}".format(sql1.format(table,index,algnmt_field,str(i),index,indexing_col,index,key_field,algnmt)))
    main_logger.info=("update indexing_col sql statement: {}".format(sql2.format(table,indexing_col,algnmt_field,str(i),indexing_col,indexing_col,key_field,algnmt)))
    while i<n:
        cursor.execute(sql1.format(table,index,algnmt_field,str(i),index,indexing_col,index,key_field,algnmt))
        cursor.execute(sql2.format(table,indexing_col,algnmt_field,str(i),indexing_col,indexing_col,key_field,algnmt))
        connection.commit()
        if i in range(1,int(n),500):
            print "{}".format(sql1.format(table,index,algnmt_field,str(i),index,indexing_col,index,key_field,algnmt))
            print "{}".format(sql2.format(table,indexing_col,algnmt_field,str(i),indexing_col,indexing_col,key_field,algnmt))
        i=i+1
    check_sql="select max(length({})) from {} where {}='{}'"
    cursor.execute(check_sql.format(index,table,key_field,algnmt))
    check_results = cursor.fetchall()[0]
    print "alignment indexing complete for file {}. {} positions indexed for an alignment {} characters long".format(algnmt,check_results,str(n))
    connection.close()
    return coords(table,algnmt)

def coords(table,algnmt):
    connection = connect_to_mysqdb('localhost', 'FamilyCpolymerases', 'sequser')
    cursor = connection.cursor()
    cursor.execute("update {} set {} = '' where {} = '{}'".format(table,coords_field,key_field,algnmt))
    connection.commit()
    cursor.execute("update {} set {} = 0 where {} = '{}'".format(table,indexing_col,key_field,algnmt))
    connection.commit()
    print "calculating alignment coordinates for {} in table {}".format(algnmt,table)
    sql="select avg(length({})) from {} where {}='{}'"
    print "calculate alignment length sql query:{}".format(sql.format(algnmt_field,table,key_field,algnmt))
    cursor.execute(sql.format(algnmt_field,table,key_field,algnmt))
    n=cursor.fetchall()[0]
    n=str(n[0])
    n=float(n)
    i=1
#    sql="update {} set {} = if(substring({},{},1) not like '-',concat({},' ',substring({},{},1),',',substring_index(substring_index({},' ',{}),' ',-1),',',{}),{}) where {}='{}'"
    sql="update {} set {} = concat({},' ',substring({},{},1),',',substring_index(substring_index({},' ',{}),' ',-1),',',{}) where {}='{}'"    
    main_logger.debug("update coords statement for file {}:{}"
    .format(algnmt,sql.format(table,coords_field,coords_field,algnmt_field,str(i),index,str(i),str(i+1),key_field,algnmt)))    
    print "calculating alignment coordinates for file {}".format(algnmt)     
    while i<n:
        cursor.execute(sql.format(table,coords_field,coords_field,algnmt_field,str(i),index,str(i+1),str(i),key_field,algnmt))
        i=i+1
        connection.commit()

if __name__ == '__main__':
##basic parameters
    refIDs={'DnaG':'P0ABS5','DnaC':'P0AEF0','DnaD':'P39787','DnaA':'P03004','DnaB':'P0ACB0'}
    basedir='/Users/pczrr/Documents/Work/Project-Folders/B-subtilis-DNA-Replication/Processing/'
    dir_dict = {'DnaA':'/Users/pczrr/Documents/Work/Project-Folders/B-subtilis-DNA-Replication/Processing/DnaA-alignment/Analysis',
        'DnaC':'/Users/pczrr/Documents/Work/Project-Folders/B-subtilis-DNA-Replication/Processing/DnaC-alignment/Analysis',
        'DnaG':'/Users/pczrr/Documents/Work/Project-Folders/B-subtilis-DNA-Replication/Processing/DnaG-alignment/Analysis',
        'DnaB':'/Users/pczrr/Documents/Work/Project-Folders/B-subtilis-DNA-Replication/Processing/DnaB-alignment/Analysis',
        }
    ID_field='Accession_Number'
    seq_field='Sequence'
    algnmt_field='Aligned_seq'
    key_field='Header'
    dest_t='mh_uref_alignments'
    index='Aligned_seq_index'
    indexing_col='Indexing'
    coords_field='Aligned_seq_coords'
    mhdata_t='MHData'
    main_logger = reinitialize_logging(basedir)
    file_dict={}
##main loop
    for pro,dir in dir_dict.iteritems():
        files=glob.glob(dir + '/*_aligned.fasta')
        for f in files:
            name = os.path.splitext(os.path.split(f)[1])[0]
            file_dict[f]=(name,pro+'s')
    main_logger.debug('___________________________________completed file list____________________________________{}'.format(file_dict))
    create_table(dest_t)
    for file,info in file_dict.iteritems():
        main_logger = reinitialize_logging(basedir)
        main_logger.info("processing: {} {}".format(file,str(info)))
        upload(dest_t,file,info)
    connection = connect_to_mysqdb('localhost', 'FamilyCpolymerases', 'sequser')
    cursor = connection.cursor()
    main_logger = reinitialize_logging(basedir)
    sql1="update {} a, {} b set a.{} = substring_index(substring_index(b.{},' ',a.Position+1),' ',-1) where a.DataName=b.{}"
    sql2="update {} a, {} b set a.{} = b.{} where a.DataName=b.{}"
    main_logger.debug("add coords sql statement:{}".format(sql1.format(mhdata_t,dest_t,'Data_coord', coords_field,key_field)))
    cursor.execute(sql1.format(mhdata_t,dest_t,'Data_coord',coords_field,key_field))
    connection.commit()
    cursor.execute(sql2.format(mhdata_t,dest_t,'uref_id',ID_field,key_field))
    connection.commit()
    connection.close()
    
