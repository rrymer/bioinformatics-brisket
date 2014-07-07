# -*- coding: utf-8 -*-
"""
Created on Fri Apr  4 16:14:45 2014

@author: Richard Rymer

Takes the raw output from PFAM for a given set of HMMs for a series of proteins,
and builds a sequence database with equivalent species represented for each
protein.
"""

import MySQLdb as mdb
import uniprot
import re
import time
import multiprocessing as mp
import glob
import itertools
import subprocess
import logging
import os

class Logger():
    logger= None
    def myLogger(self,path):
        if None == self.logger:
            self.logger=logging.getLogger(__name__)
            self.logger.setLevel(logging.DEBUG)
            handler=logging.FileHandler('{}/Build_DNARep_seqDB.log'.format(path))
            formatter=logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
            handler.setFormatter(formatter)
            self.logger.addHandler(handler)
        return self.logger

def find_and_test_data(path,logger):
    files = glob.iglob(path)
    for f in files:
        logger.debug(f)
    files = glob.iglob(path)
    for f in files:
        com="grep -c '>' {}".format(f)
        logger.debug('grep command for file {}:{}'.format(f,com))
        p = subprocess.Popen([com],shell=True,stdout=subprocess.PIPE)
        p.wait()
        count = p.stdout.readline()
        count=int(count)
        logger.info('processing %s sequences',str(count))
        if count < 1:
            print count
            logger.error('file in data dir %s is not a fasta file, skipped',f)
            return
        else:
            logger.info('processing sequence file %s',f)
            continue
    return combine_data(path,logger)

def convert_data(data,logger):        
    com = "convbioseq -i {} {} '{}'".format('fasta','tab',data)
    logger.debug('converting wrapped fasta file with {}'.format(com))
    p = subprocess.Popen([com],shell=True)   
    p.wait()
    return

def combine_data(path,logger):
    files = glob.iglob(path)
    for data in files:
        if debug:
            with open('{}/compiled_seqs.txt'.format(os.path.split(data)[0]),'a') as outfile:
                outfile.write('\n')
            data_size0 =os.path.getsize('{}/compiled_seqs.txt'.format(os.path.split(data)[0]))
            logger.debug('compiled seq starting size:{}'.format(data_size0))
        with open('{}/compiled_seqs.txt'.format(os.path.split(data)[0]),'a') as outfile:
            with open(data,'r') as infile:
                    outfile.write(infile.read())
        if debug:
            data_size1 =os.path.getsize('{}/compiled_seqs.txt'.format(os.path.split(data)[0]))
            logger.debug('compiled seq file: {}/compiled_seqs.txt'.format(os.path.split(data)[0]))
            logger.debug('data added to compiled file: {}'.format(str(data_size1-data_size0)))
            if data_size1 <= data_size0:
                logger.error('no data added from {} to compiled data file {}'.format(data,outfile))
    return

def split_every(n, iterable):
    i = iter(iterable)
    piece = list(itertools.islice(i, n))
    while piece:
        yield piece
        piece = list(itertools.islice(i, n))

def create_table(name,logger):
    connection = mdb.connect(host='localhost',db='FamilyCpolymerases',user='sequser')
    cursor=connection.cursor() 
    print 'assembling ' + name + ' sequence database'    
    sql=("CREATE TABLE `%s` (" +
      "`NCBI_GI` varchar(12) DEFAULT NULL," +
      "`Organism` varchar(200) DEFAULT NULL,"+
      "`PolGroup` varchar(45) DEFAULT NULL,"+
      "`Domains` varchar(45) DEFAULT NULL,"+
      "`Phylum` varchar(45) DEFAULT NULL,"+
      "`Class` varchar(45) DEFAULT NULL,"+
      "`Length` int(5) DEFAULT NULL,"+
      "`Sequence` varchar(3000) DEFAULT NULL,"+
      "`Header` varchar(200) NOT NULL DEFAULT 'No entry',"+
      "`Accession_Number` varchar(8) DEFAULT 'X00000',"+
      "`OrganismID` varchar(12) DEFAULT NULL,"+
      "`Entry_Name` varchar(45) DEFAULT NULL,"+
      "`Lineage` varchar(90) DEFAULT NULL,"+
      "`Aligned_seq` varchar(10000) DEFAULT NULL,"+
      "`gfAligned_seq` varchar(10000) DEFAULT NULL,"+
      "`Aligned_seq_index` varchar(10000) DEFAULT NULL,"+
      "`Indexing` int(11) DEFAULT '0',"+
      "`gfAligned_seq_index` varchar(10000) DEFAULT NULL,"+
      "`gfAligned_seq_coords` longtext,"+
      "`PolProfile` varchar(45) DEFAULT NULL,"+
      "PRIMARY KEY (`Header`)"+
      ") ENGINE=InnoDB DEFAULT CHARSET=latin1;")
    sql = sql % name
    logger.debug('create table statement for {}: {}'.format(name,sql))
    cursor.execute(sql)
    connection.commit()
    connection.close()

def check_tables(cursor):    
        cursor.execute("show tables")
        tables=[]
        for t in cursor:
            tables.append(t[0])
        return tables,cursor

def data_import(data_file,table,logger):
    logger.info('importing data from {} into table {}'.format(data_file,table))
    connection = mdb.connect(host='localhost',db='FamilyCpolymerases',user='sequser')
    cursor=connection.cursor()
    tables = check_tables(cursor)
    if table not in tables:
        try:
            connection = mdb.connect(host='localhost',db='FamilyCpolymerases',user='sequser')
            cursor=connection.cursor()
            logger.info('table does not exist, creating table %s',table)        
            create_table(table,logger)
            time.sleep(1)
            tables=check_tables(cursor)
            if table not in tables:
                logger.debug('why cannot it find the damn table!')
        except mdb.OperationalError:
            pass
#            logger.error('unable to create table, checking...')
#            tables = check_tables(cursor)
#            logger.error('tables in db: {}'.format(str(tables)))
#            logger.error('skipping to next table')
#            time.sleep(1)
#            return
    sql="load data local infile '{}' ignore into table `{}` fields terminated by '\t' (Header, Sequence)"
    logger.debug('load data statement for {}: {}'.format(table,sql.format(data_file,table)))
    cursor.execute(sql.format(data_file,table))
    connection.commit()
    connection.close()
    logger.info('data imported into table %s',table)
    return extract_accession_numbers(table,logger)
    
def extract_accession_numbers(table,logger):
    connection = mdb.connect(host='localhost',db='FamilyCpolymerases',user='sequser')
    cursor=connection.cursor() 
    cursor.execute("update " + table + " set Accession_Number=left(Header, 6)")
    logger.debug('extract ID sql statement for table {}:{}'.format(table,"update " + table + " set Accession_Number=left(Header, 6)"))
    connection.commit()
    logger.info('accession numbers extracted into Accession_Number for table %s',table)
    logger.debug('unique ID sql query for table {}:{}'.format(table,"select Header,Accession_Number from " + table + " group by Accession_Number"))
    cursor.execute("select Header,Accession_Number from " + table + " group by Accession_Number")
    entries=dict(cursor.fetchall())
    unique=[]
    for k,v in entries.iteritems():
        unique.append(k)
    unique = str(tuple(unique))
    logger.debug('unique accession numbers in table %s:%s',table,unique)
    cursor.execute("delete from " + table + " where Header not in " + unique)
    sql="select Accession_Number,count(*) from {} group by Accession_Number having count(*)>1"
    logger.debug('unique ID sql query for {}:{}'.format(table,sql.format(table)))
    cursor.execute(sql.format(table))        
    results=dict(cursor.fetchall())
    for k,v in results.iteritems():
        logger.info('duplicate IDs in {}: {}, {}'.format(table,str(k),str(v)))
    connection.commit()
    connection.close()
    logger.info('Accession numbers extracted and duplicates removed for members of %s', table)

def fetch_and_assign_taxids(seq,table):
    non_decimal = re.compile(r'[^\d.]+') 
    results=uniprot.fetch_uniprot_metadata(seq)
    time.sleep(1)
    connection = mdb.connect(host='localhost',db='FamilyCpolymerases',user='sequser')
    cursor=connection.cursor()
    for k,v in results.iteritems():  
        TaxID = non_decimal.sub('', results[k]['TaxID'])
        ID = results[k]['id']
        k=repr(k)
        sql="update {} set {}='{}' where Accession_Number={}"        
        cursor.execute(sql.format(table,'OrganismID',TaxID,k))
#        logger.debug('update TaxID sql statement for table {}:{}'.format(table,sql.format(table,'OrganismID',TaxID,k)))
        cursor.execute(sql.format(table,'Entry_Name',ID,k))
#        logger.debug('update ID sql statement for table {}:{}'.format(table,sql.format(table,'Entry_Name',ID,k)))
        connection.commit()
        time.sleep(1)
    connection.close()
    
def assemble_ids(table,logger):
    connection = mdb.connect(host='localhost',db='FamilyCpolymerases',user='sequser')
    cursor=connection.cursor()
    connection.commit()
    cursor.execute("select Accession_Number, '1' from " + table)
    results=dict(cursor.fetchall())
    connection.close()
    seqids =[]
    for k,v in results.iteritems():
        seqids.append(k)
    logger.info('Accession numbers assembled for TaxID retrieval for table %s',table)
    return seqids
    
def cull_update_db(table,logger):
    connection = mdb.connect(host='localhost',db='FamilyCpolymerases',user='sequser')
    cursor=connection.cursor()    
    logger.debug('delete sql statement for table {}:{}'.format(table,"delete from " + table + " where OrganismID not in (select OrganismID from RepPols)"))
    cursor.execute("delete from " + table + " where OrganismID not in (select OrganismID from RepPols)")
    connection.commit()
    IDs=['Organism','PolProfile','Phylum','Lineage','Class']
    sql='update {} a, RepPols b set a.{}=b.{} where a.OrganismID=b.OrganismID'
    for i in IDs:
        cursor.execute(sql.format(table,i,i))
        logger.debug('update tables statement for table {}:{}'.format(table,sql.format(table,i,i)))
    connection.commit()    
    connection.close()
    logger.info('Taxonomical data assigned for table %s',table)

def clean_up(path,logger):
    connection = mdb.connect(host='localhost',db='FamilyCpolymerases',user='sequser')
    cursor=connection.cursor()
    cursor.execute("SHOW TABLES")
    for table in cursor:
        if debug:
            sql="select Accession_Number,'' from {} where OrganismID is null"
            logger.debug('select null TaxIDs sql statement for table {}:{}'.format(table[0],sql.format(table[0])))
            cursor.execute(sql.format(table[0]))
            results=dict(cursor.fetchall())
            for k,v in results.iteritems():
                logger.debug('bad IDs in {}: {}'.format(table[0],str(k)))
        sql="delete from {} where OrganismID is null"        
        cursor.execute(sql.format(table[0]))
    missingIDs=[]
    cursor.execute("SHOW TABLES")
    for table in cursor:
        if 'RepPols' in table or 'Test' in table:
            continue
        sql="select OrganismID,Phylum from RepPols where OrganismID not in (select OrganismID from {} where OrganismID is not null) group by OrganismID"      
        logger.debug('select null TaxIDs sql statement for table {}:{}'.format(table[0],sql.format(table[0])))        
        cursor.execute(sql.format(table[0]))
        results=dict(cursor.fetchall())
        for k,v in results.iteritems():
            missingIDs.append(k + '\t' + v + '\t' + table[0])
    missingIDs=set(missingIDs)
    with open(path + 'missingIDs.txt','w') as outfile:
        for ID in missingIDs:
            outfile.write(ID+'\n')
    connection.commit()
    connection.close()
    
if __name__ == '__main__':
    debug=True
    loggers = Logger()
    basedir='/Users/pczrr/Documents/Work/Project-Folders/B-subtilis-DNA-Replication/Processing/'
    main_logger = loggers.myLogger(basedir)
    main_logger.info('############################ BEGIN DATABASE ASSEMBLY ################################')
    tables = ['DnaA','DnaB','DnaCI','DnaD']
    for t in tables:
        logger = loggers.myLogger('{}{}-alignment/'.format(basedir,t))
        logger.info('beginning construction of %s sequence table',t+'s')
        find_and_test_data('{}{}-alignment/Data/*'.format(basedir,t),logger)
#        files = glob.iglob('{}{}-alignment/Data/*.tab'.format(basedir,t))
        convert_data('{}{}-alignment/Data/{}'.format(basedir,t,'compiled_seqs.txt'),logger)
        data_import('{}{}-alignment/Data/{}'.format(basedir,t,'compiled_seqs.tab'),t+'s',logger)
        main_logger.info('Data imported for table %s',t+'s')
        seqs = assemble_ids(t + 's',logger)
        seqs = split_every(100,seqs)
        p0 = mp.Pool(8)
        for seq in seqs:
#            logger.debug('assigning organismIDs for %s IDs in table %s',len(seq),t+'s')
            p0.apply_async(fetch_and_assign_taxids,args=(seq,t + 's'))
        p0.close()
        p0.join()
        main_logger.info('TaxIDs and entry names assigned for table %s',t+'s')
        cull_update_db(t+'s',logger)
        main_logger.info('Completed assmebly of ' + t + ' sequence database.')
    clean_up(basedir,logger)
    main_logger.info('Completed assmebly of all sequence databases.')
    
