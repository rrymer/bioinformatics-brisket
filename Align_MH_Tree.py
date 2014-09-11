# -*- coding: utf-8 -*-
"""
Created on Thu Apr 24 16:21:47 2014

@author: Richard Rymer

Retrieves a balanced set of sequences from a sequence database for a given set 
of proteins, "balanced" refers to an approximately equal number of sequences for
each group in some group criterion, like phylum.

Appends that sequence set to a data set with all of the sequences of each protein
for each group level (e. g., proteobacteria DnaG sequences).

Makes an alignment with each sequence dataset.

Submits each dataset for analysis with MultiHarmony and tree construction with
RaXML.
"""

import MySQLdb as mdb
import glob
import subprocess as sp
import os
import multiprocessing as mp
from Bio.Align.Applications import MafftCommandline
#from Bio.Align.Applications import RaxmlCommandline #not implemented
import itertools
import logging
from cStringIO import StringIO
from Bio import AlignIO
import fileinput
import time
import numpy as np

###functions
class Logger():
    logger=None
    def myLogger(self,path,name,level):
        print 'initializing logger {} located at {}'.format(name,path)
        if None == self.logger:
            self.logger=logging.getLogger(name)
            self.logger.setLevel(level)
            self.handler=logging.FileHandler('{}align_and_tree_refseqonly_{}.log'.format(path,name))
            formatter=logging.Formatter("%(message)s")
            self.handler.setFormatter(formatter)
            self.logger.addHandler(self.handler)
        return self.logger

def reinitialize_logging(path):
    reload(logging)
    main_loggers = Logger()
    main_logger = main_loggers.myLogger(basedir,'main',logging.INFO)
    return main_logger

def connect_to_mysqdb(host, db, un):
    connection = mdb.connect(host=host,db=db,user=un)
    return connection

def worker(in_queue):
    item = in_queue.get()
    MakeTrees(item)

def process(fnclogger,table,ref_query,filenames):
    fnclogger.info("reference sequence sql query:{}".format(ref_query))
    control_groups = find_groups(fnclogger,table,control_field,query=ref_query)
    fnclogger.info("{} groups found in table {}:{}".format(control_field,table,str(control_groups)))
    f = ref_dexport(fnclogger,table,control_field,control_groups,ref_query)
    return f

def find_groups(fnclogger,table,group_field,query=''):
    connection = connect_to_mysqdb('localhost', 'FamilyCpolymerases', 'sequser')
    cursor=connection.cursor()
    query = 'where {}'.format(query)
    main_logger.info("searching table {} for a list of {}s.".format(table,group_field))
    sql="select {},count(*) from {} {} group by {}"
    fnclogger.debug=("determine counts sql query:{}".format(sql.format(group_field,table,query,group_field)))
    cursor.execute(sql.format(group_field,table,query,group_field))
    groups=dict(cursor.fetchall())
    main_logger.info("groups found in table {}:{}.".format(table,str(groups)))
    connection.close()
    return groups

def ref_dexport (fnclogger,table,group_field,groups,query):
    connection = connect_to_mysqdb('localhost', 'FamilyCpolymerases', 'sequser')
    cursor = connection.cursor()
    if query:
        query = 'and {}'.format(query)
    main_logger.info('exporting reference data')
    hsql="select concat(OrganismID, '|', PolProfile , '|', '{}', '|',Accession_Number), {} from {} where {} like '{}' {} {}" 
    lsql="order by rand() limit 10 into outfile '{}{}_Reference_data.fasta' fields terminated by '\n' lines starting by '>'"
    i=1
    for g,c in groups.iteritems():
       try:
           print i
           cursor.execute(hsql.format(g,seq_field,table, group_field,g,query,lsql
               .format(tmp_dir,table+str(i).zfill(3))))
           i=i+1
       except mdb.OperationalError, e:
           print e
           if os.path.isfile('{}{}_Reference_data.fasta'.format(tmp_dir,table+str(i).zfill(3))):
               print i
               i=i+1
               cursor.execute(hsql.format(g,seq_field,table, group_field,g,query,lsql
               .format(tmp_dir,table+str(i).zfill(3))))
               i=i+1
    fnclogger.info("assmebling reference file for {} {}".format(group_field,g))
    try:
       files = glob.glob(tmp_dir + table + '*' + '_Reference_data.fasta')
       with open(tmp_dir + table + '_Reference_data.fasta','a') as outFile:
           for f in files:
               with open(f,'r') as inFile:
                   outFile.write(inFile.read())
    except IOError:
       fnclogger.error("unable to assemble files for {}".format(table))
    connection.close()
    return tmp_dir + table + '_Reference_data.fasta'
        
def seq_dexport (fnclogger,table,group_field,groups,test_field,test_grp,query):
    connection = connect_to_mysqdb('localhost', 'FamilyCpolymerases', 'sequser')
    cursor = connection.cursor()
    main_logger.info("exporting primary data for table {}".format(table))
    query = 'and {}'.format(query)    
    hsql="select concat('Group:{}', '|', '{}','|',Accession_Number), {} from {} where {}='{}' {} {}"
    lsql="order by rand() limit 300 into outfile '/tmp/{}{}{}.fasta' fields terminated by '\n' lines starting by '>';"
    fnclogger.info("select group data statement: {}".format(hsql
    .format('::group::','tst_grp',seq_field,table,group_field,'::group::',query,lsql.format('::group::',table,seq_field))))
    search=[]
    for k,v in groups.iteritems():
        if int(v) >= 3:
            cursor.execute(hsql.format(k,test_grp,'count(*)',table,group_field,k,query,'group by {}'.format(control_field)))
            rslts=cursor.fetchall()
            try: 
                if int(rslts[0][1]) >= 3:
                    search.append(k)
                else:
                    fnclogger.info("subgroup {} {} {} is too small to analyze with only {} members".format(k,rslts,table,rslts[1]))
                    continue
            except IndexError, e:
                fnclogger.info("group {} {} is too small to analyze with only {} members".format(k,table,str(v)+str(e)))
                continue
#            search[k]=test_grp
        else:
            fnclogger.info("group {} {} is too small to analyze with only {} members".format(k,table,v))
            continue
    fnclogger.info("search list for {} {} = {}".format(table, test_grp, str(search)))
    for grp in search:
            fnclogger.info("exporting {} data for {} {} {} {} {}".format(seq_field,group_field,grp,test_field,test_grp,table))
            filenames['/tmp/{}{}{}.fasta'.format(grp+table,test_field+test_grp,seq_field)]=grp
            cursor.execute(hsql.format(grp,test_grp,seq_field,table,group_field,grp,query,lsql.format(grp+table,test_field+test_grp,seq_field)))
    connection.close()
    return hsql,filenames

def add_header( header, intFile ):
    print 'adding headers'
    for line in fileinput.input(intFile, inplace=True):
        if fileinput.isfirstline():
            print '\n',
        print line,
    for line in fileinput.input(intFile, inplace=True):
        if fileinput.isfirstline():
            print header,
        print line,

def create_table(name):
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
      "`Sequence` longtext,"+
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
    cursor.execute(sql)
    connection.commit()
    connection.close()

def MakeAlignments(seqs,name,path):   ##aligns exported data 
    if os.path.isfile(path + name + '_aligned.fasta') is False:
        in_file = seqs
        
#        try:
#            print file_table_dict
#        except KeyError,e:
#            print e
#            return
        for k,v in tables.iteritems():
            if k in seqs:
                file_table_dict[seqs]=v
            else:
                continue
        mafft_cline = MafftCommandline(input=in_file, auto=True, reorder=True)
        stdout, stderr = mafft_cline()
        algnmt = AlignIO.read(StringIO(stdout),'fasta')
        AlignIO.write(algnmt,path + name + '_aligned.fasta','fasta')
        connection = mdb.connect(host='localhost',db='FamilyCpolymerases',user='sequser')
        cursor=connection.cursor()
        create_table(name)
        sql = "insert ignore into `{}` ({},{},{}) values('{}','{}','{}')"
        usql= "update `{}` a, `{}` b set a.{}=b.{} where a.{}=b.{}"
        for record in algnmt:
            print 'inserting sequence {} into {}'.format(record.id,name)
            try:
                cursor.execute(sql.format(name,key_field, ID_field, seq_field,str(record.id),str(record.id)[-6:],str(record.seq)))
                connection.commit()
            except Exception, e:
                print e
        time.sleep(1)
        try:
            cursor.execute(usql.format(name, file_table_dict[seqs], taxid_field, taxid_field, ID_field, ID_field))
            connection.commit()
        except Exception, e:
            print e
        time.sleep(1)
        expsql="select {},{} from `{}` into outfile '{}' fields terminated by '{}' {}"
        print 'exporting fasta'
        cursor.execute(expsql.format(key_field, seq_field, name, path + name + '_aligned.fasta','\n',"lines starting by '>'"))
        print 'export phylip'
        cursor.execute(expsql.format(taxid_field, seq_field, name, path + name + '_aligned.phy','\t',''))
        cursor.execute("select max(length({})),count(*) from {}".format(seq_field,name,ID_field))
        params=cursor.fetchall()
        header = "{}\t{}".format(params[0][1],params[0][0])
        add_header(header,tmp_dir + name + '_aligned.phy')
#        cursor.execute("select count(*) from {} where {} like '%Reference%'".format(name,key_field))
#        ref_cnt=cursor.fetchall()
#        try:
#            cnt_dict[name + '_aligned']=(int(params[0][1])-int(ref_cnt[0][0]),int(ref_cnt[0][0]))
#        except Exception, e:
#            print e
        cursor.execute("drop table {}".format(name))
        connection.commit()
        time.sleep(10)
        connection.close()
    else:
        return

def fork(fasta,phylip,name):
    pass

def MakeTrees(nm,path,algnt):
#    raxml_cline=RaxmlCommandline(sequences=algnt,model=submod, name=nm)
#    stdout, stderr = raxml_cline()
    raxml_c = "raxml -T 2 -m PROTGAMMAAUTO -p 314159 -n {} -s {} -w {}".format(nm,algnt,path)
    p = sp.Popen([raxml_c],shell=True)
    p.wait()

def run_mh(algnt,name,out_arrays):
     print 'running MH with ' + name
     r=0
     t=0
     with open(algnt,'r') as a:    
         for line in a:
             if '>' in line:
                 if 'Reference' in line:
                     r = r+1
                 else:
                     t = t+1
             else:
                 continue
     local_array = np.array(['none',0,0,0,0,0,'a','b'])
#     groupsizes='{} {}'.format(cnt_dict[name][0],cnt_dict[name][1])
     groupsizes='{} {}'.format(str(t),str(r))
     print groupsizes
     p=sp.Popen(['perl', '/Users/pczrr/Google Drive/Coding/Scripts/mh_run.pl', algnt,'--groups', groupsizes],stdout=sp.PIPE)
     results = p.stdout
     for line in iter(results):
        if 'Output multi-Relief weights' in line:
            results.next()
            print line
            break
     while line:
        try:        
            tmp_line=results.next().split('\t')
            tmp_line.insert(0,name)
            del tmp_line[-1]
            tmp_array = np.array(tmp_line)
            local_array = np.vstack((local_array,tmp_array))
        except StopIteration:
            break
     with open('/tmp/out_arrays.txt','a') as outfile:
         outfile.write(local_array)
         
     main_logger.info('current arrays: {}'.format(str(out_arrays)))

def AssembleTreeFile(trees): ##assembles complete tree file
    with open(basedir+'/scratch/' + 'compiledtrees.txt', 'w') as outfile:     
        for t in trees:
            if os.path.getsize(t) <= 0:
                continue
            print 'appending ' + os.path.splitext(os.path.split(t)[1])[0] + ' to compiled tree file.'
            with open(t) as infile:        
                outfile.write(infile.read())

def process_tree():
    """
    Do something with trees
    """
    raise NotImplementedError
    #TODO

if __name__=='__main__':
##parameters and definitions
    reps=50
    exclude=False
    basedir='/Users/pczrr/Documents/Work/Project-Folders/B-subtilis-DNA-Replication/Processing/'
    DnaA_loggers=None
    DnaB_loggers=None
    DnaCI_loggers=None
    DnaG_loggers=None
    DnaA_logger=None
    DnaB_logger=None
    DnaCI_logger=None
    DnaG_logger=None
    hsql=None
    lsql=None
    test_files='/tmp/test*.txt'
    test_t='DnaGs'
    tmp_dir='/tmp/'
    cnt_dict = {}
    test= False
    table_loggers = {'DnaA':DnaA_loggers,'DnaB':DnaB_loggers,'DnaC':DnaCI_loggers,'DnaG':DnaG_loggers}
    logger_dict = {'DnaA':DnaA_logger,'DnaB':DnaB_logger,'DnaC':DnaCI_logger,'DnaG':DnaG_logger}
    control_field='Phylum'
    test_field='PolProfile'
    tables={'DnaA':'DnaAs','DnaB':'DnaBs','DnaC':'DnaCs','DnaG':'DnaGs'}
    done={}
    refIDs={'DnaG':'P0ABS5','DnaC':'P0AEF0','DnaD':'P39787','DnaA':'P03004','DnaB':'P0ACB0'}
    ID_field='Accession_Number'
    seq_field='Sequence'
    key_field='Header'
    taxid_field='OrganismID'
    tmp_table='Test'
    file_table_dict={}
    paths=[]
    out_arrays=[]
    algnts_file='/tmp/algnts.txt'
    seq_data_name='{}_Reference_data.fasta'
    data_array=np.array(['none',0,0,0,0,0,'a','b'])
    phyla=['Acidobacteria','Actinobacteria','Aquificae','Bacteroidetes','Caldiserica','Chlamydiae',
    'Chlorobi','Chloroflexi','Chrysiogenetes','Cyanobacteria','Deferribacteres','Deinococcus-Thermus','Dictyoglomi','Elusimicrobia',
    'Fibrobacteres','Firmicutes','Fusobacteria','Gemmatimonadetes','Ignavibacteria','Nitrospirae','Planctomycetes','Proteobacteria',
    'Spirochaetes','Synergistetes','Tenericutes','Thermodesulfobacteria','Thermotogae','Verrucomicrobia']
    dirs=['/Users/pczrr/Documents/Work/Project-Folders/B-subtilis-DNA-Replication/Processing/DnaA-alignment',
          '/Users/pczrr/Documents/Work/Project-Folders/B-subtilis-DNA-Replication/Processing/DnaB-alignment',
          '/Users/pczrr/Documents/Work/Project-Folders/B-subtilis-DNA-Replication/Processing/DnaC-alignment',
          '/Users/pczrr/Documents/Work/Project-Folders/B-subtilis-DNA-Replication/Processing/DnaG-alignment']
    cmp_tree_f='compiled_trees.txt'
##initialize logging
    main_logger = reinitialize_logging(basedir)    
##find and assess groups
    main_logger.info("BEGIN SEQUENCE SELECTION AND EXPORT")
    files=[]
    for k,v in tables.iteritems():
        if test == True:
            print 'test run'
            break
        else:
            filenames={}
            main_logger = reinitialize_logging(basedir) 
            table = v
            refID=refIDs[k]
            ref_query="PolGroup='ref{}'".format(str(k))
            workingdir="{}{}-alignment/".format(basedir,k)
            table_loggers[k] = Logger()        
            logger_dict[k] = table_loggers[k].myLogger(workingdir,k,logging.DEBUG)
            logger_dict[k].info('logger initialized for {}'.format(k))
            i = 1
            while i < reps:          
                f = process(logger_dict[k],table,ref_query,filenames)
                name = os.path.splitext(os.path.split(f)[1])[0]
                sp.Popen(['cp',f,'{}Scratch/{}{}.fasta'.format(workingdir,name,str(i).zfill(3))])
                print 'removing temp files'
                tmp_files = glob.glob(tmp_dir + '*.fasta')
                for f in tmp_files:
                    p = sp.Popen(['rm',f],stdout=sp.PIPE,stderr=sp.PIPE)
                    out,err = p.communicate()
                    print out,err
                i=i+1
            tmp_files=glob.glob('{}Scratch/{}.fasta'.format(workingdir,'*'))
            for f in tmp_files:            
                files.append(f)
    main_logger = reinitialize_logging(basedir) 
    if test == True:
        files = glob.glob(test_files)
        for f in files:
            file_table_dict[f]=test_t
        files = glob.glob(test_files)
    else:
#        files=[]
#        for k,v in tables.iteritems():
#            tmp_files = glob.glob('/tmp/' + seq_data_name.format(v))
#            files.append(tmp_files)
        pass
    if files:
        main_logger.info('data export completed')
        print files
    else:
        main_logger.error('data export was not successful')
    main_logger.info("BEGIN SEQUENCE ALIGNMENT")
    p1 = mp.Pool(8)
    paths=[]
    for f in files:
        main_logger.info('aligning file {}'.format(f))
        name = os.path.splitext(os.path.split(str(f))[1])[0]
        p = os.path.split(f)[0] +'/'
        paths.append(p)
        p1.apply_async(MakeAlignments, args=(f,name,p))
    p1.close()
    p1.join()  
#wait for all alignments to finish, then move on to making trees, after a little processing         
    if test == True:
        algnts = glob.glob(tmp_dir + 'test*_aligned.fasta')
    else:
        algnts = glob.glob(tmp_dir + '*_aligned.fasta')
    if algnts:
        main_logger.info("alignment complete")
    else:
        main_logger.error('sequence alignment was not successful')
    p4=mp.Pool(8)
    for a in algnts:
        name = os.path.splitext(os.path.split(a)[1])[0]
        run_mh(a,name,out_arrays)
#        p4.apply_async(run_mh, args=(a,name))
#        time.sleep(10)
    for a in out_arrays:
        data_array = np.vstack(data_array,a)
    print data_array
    np.savetxt(data_array,delimiter='\t')
    p4.close()
    p4.join()
    
    name = os.path.splitext(os.path.split(a)[1])[0]
    algnts=[]
    if paths:
        paths = set(paths) 
        for p in paths:
            tmp_files = glob.glob(p + '*_aligned.fasta')
            for f in tmp_files:
                algnts.append(f)
    else:
        with open(algnts_file,'r') as infile:
            for l in infile:
                algnts.append(l.rstrip())
    main_logger.info("BEGIN TREE CONSTRUCTION")      
    iters = itertools.chain(algnts, (None,)*8)
    p3 = mp.Pool(4)
    for a in enumerate(iters):
        if a[1]:
            name = os.path.splitext(os.path.split(a[1])[1])[0]
            path = os.path.split(a[1])[0]
            print 'making tree for {}'.format(name)
            p3.apply_async(MakeTrees,args=(name,path +'/',a[1]))
        else:
            break
    p3.close()  
    p3.join()
    #assemble trees into one file
    trees = glob.glob(tmp_dir + '*phyml_tree.txt')
    if trees:
        main_logger.info('tree construction complete')
    else:
        main_logger.error('tree construction was not successful')
    main_logger.info("BEGIN TREE COMPILATION")
    AssembleTreeFile(trees)
    process_tree()

                

    
