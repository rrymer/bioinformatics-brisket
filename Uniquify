# -*- coding: utf-8 -*-
"""
Created on Fri Apr  4 16:14:45 2014

@author: Richard Rymer

Trolls a database of sequences, and flags duplicate seqeunces, fragments, and sequence errors.
"""
##prequisites
import MySQLdb as mdb
import subprocess as sp
import os
import numpy as np
import itertools as it
from pylab import *
import glob
import logging
from Bio.Align.Applications import MafftCommandline
from Bio import pairwise2 as align
#import Alignments_to_DB as upload
from Bio.SubsMat import MatrixInfo as matlist
import operator

##functions
class Logger():
    logger=None
    def myLogger(self,path,name):
        if None == self.logger:
            self.logger=logging.getLogger(name)
            self.logger.setLevel(logging.DEBUG)
            self.handler=logging.FileHandler('{}uniquify_{}_seqDB.log'.format(path,name))
            formatter=logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
            self.handler.setFormatter(formatter)
            self.logger.addHandler(self.handler)
        return self.logger

def reinitialize_logging(path):
    reload(logging)
    main_loggers = Logger()
    main_logger = main_loggers.myLogger(basedir,'main')
    return main_logger

def connect_to_mysqdb(host, db, un):
    connection = mdb.connect(host=host,db=db,user=un)
    return connection

def find_groups(table,group_field,query=''):
    connection = connect_to_mysqdb('localhost', 'FamilyCpolymerases', 'sequser')
    cursor=connection.cursor()
    main_logger.info('searching table {} for a list of {}s.'.format(table,group_field))
    sql="select {},count(*) from {} {} group by {}"
    logger.debug=('determine counts sql query:{}'.format(sql.format(group_field,table,query,group_field)))
    cursor.execute(sql.format(group_field,table,query,group_field))
    groups=dict(cursor.fetchall())
    main_logger.info('groups found in table {}:{}.'.format(table,str(groups)))
    connection.close()
    return groups

def select_seq(table,ref_ID,groupby,*queries):
    connection = connect_to_mysqdb('localhost', 'FamilyCpolymerases', 'sequser')
    cursor = connection.cursor()
    logger.info("sql query for reference sequence:\n{}".format("select Accession_Number, Sequence from {} where Accession_Number = '{}'".format(table,ref_ID)))
    cursor.execute("select Accession_Number, Sequence from {} where Accession_Number = '{}'".format(table,ref_ID))
    ref_seqs=dict(cursor.fetchall())
    logger.info('executing argument queries: {}'.format(str(queries)))
    results_dict={}
    for q in queries:
        i=1
        logger.info("sql query for sample sequences:\n{}".format("select Accession_Number, Sequence from {} {} {}"
        .format(table,q,groupby)))
        cursor.execute("select Accession_Number, Sequence from {} {} {}".format(table,q,groupby))
        results=dict(cursor.fetchall())
        results_dict['q%s'%str(i)]=results
        i=i+1
    connection.close()
    return ref_seqs,results_dict

def pairwise(table,ref_seq,query_seq,query_ID):
    matrix = matlist.blosum62
    max_hits={}
    score_list=[]
    algnmt = align.align.globaldx(ref_seq,query_seq,matrix)
    for a in algnmt:
        score_list.append(a[2])
    top_hit = max(score_list)
    max_hits[query_ID]=top_hit
    return max(max_hits.iteritems(), key=operator.itemgetter(1))[1]     

def uniquify(table,ref_ID,cand_query,cand_id,uflag,dupflag,flag_col='PolGroup',cand_field='OrganismID',groupby=''):
    main_logger.info('uniquifying {} {}'.format(cand_field,cand_id))
    connection = connect_to_mysqdb('localhost', 'FamilyCpolymerases', 'sequser')
    cursor = connection.cursor()
    group_ref_seqs,cand_seqs=select_seq(table,ref_ID,cand_query,groupby)
    cand_seqs=cand_seqs['q1']
    cand_seq_scores={}
    for ID,seq in cand_seqs.iteritems():
        cand_seq_score = pairwise(table,group_ref_seqs[ref_ID],seq,ID)
        cand_seq_scores[ID]=(cand_seq_score)
    try:                         
        unqID=max(cand_seq_scores.iteritems(), key=operator.itemgetter(1))[0]
    except ValueError, e:
        main_logger.error('no good match sequences found for ID {}:{}'.format(k,cand_id,e))
    unq_sql = "update {} set {}='{}' where Accession_Number='{}'"
    dup_sql = "update {} set {}='{}' where Accession_Number not like '{}' and {} = '{}'"
    cursor.execute(unq_sql.format(table,flag_col,uflag,unqID))
    connection.commit()
    cursor.execute(dup_sql.format(table,flag_col,dupflag,unqID,cand_field,cand_id))
    connection.commit()
    connection.close()
    return

###main program
if __name__ == '__main__':
##parameters definitions##
    debug=True
    basedir='/Users/pczrr/Documents/Work/Project-Folders/B-subtilis-DNA-Replication/Processing/'
    tables={'DnaG':'P0ABS5','DnaCI':'P0AEF0','DnaD':'P39787'}
    done={'DnaA':'P03004','DnaB':'P0ACB0'}
    protein_names = {'DnaA':"where Entry_Name like '%DnaA%'",
    'DnaCI':"where (Entry_Name like '%DnaC%' or Entry_Name like '%DnaI%')",
    'DnaB':"where (Entry_Name like '%DnaB%' or Entry_Name like '%DnaC%')",
    'DnaD':"where (Entry_Name like '%DnaD%' or Entry_Name like '%DnaB%')",
    'DnaG':"where Entry_Name like '%DnaG%'"}
    ref_data_table='RepPols'    
    phyla=['Acidobacteria','Actinobacteria','Aquificae','Bacteroidetes','Caldiserica','Chlamydiae',
    'Chlorobi','Chloroflexi','Chrysiogenetes','Cyanobacteria','Deferribacteres','Deinococcus-Thermus','Dictyoglomi','Elusimicrobia',
    'Fibrobacteres','Firmicutes','Fusobacteria','Gemmatimonadetes','Ignavibacteria','Nitrospirae','Planctomycetes','Proteobacteria',
    'Spirochaetes','Synergistetes','Tenericutes','Thermodesulfobacteria','Thermotogae','Verrucomicrobia']
#    query = "where (Entry_Name like {})" ## use to restrict selection to "high-quality" sequences
    all_ref_seqs = {}
#   'DnaA': {'Spirochaetes': 'F4LK72', 'Verrucomicrobia': 'D5EKE9', 'Bacteroidetes': 'Q2S2T8', 'Fibrobacteres': 'C9RIX2',
#             'Elusimicrobia': 'B2KAM4', 'Chloroflexi': 'A9WAN1', 'Planctomycetes': 'D5SMD4', 'Actinobacteria': 'D7B601',
#             'Chlorobi': 'Q3AUF1', 'Thermodesulfobacteria': 'F8A7X1', 'Deferribacteres': 'D4H160', 'Fusobacteria': 'D1AIS7',
#             'Synergistetes': 'D1B7I9', 'Chrysiogenetes': 'E6W482', 'Cyanobacteria': 'Q2JQW9', 'Proteobacteria': 'P03004',
#             'Nitrospirae': 'D8P7E7', 'Chlamydiae': 'F4DII6', 'Dictyoglomi': 'B8DYG0', 'Tenericutes': 'A9NE65', 'Caldiserica': 'I0GH90',
#             'Ignavibacteria': 'I0AFF8', 'Thermotogae': 'H2J873', 'Firmicutes': 'C9LV38', 'Deinococcus-Thermus': 'G8NDC1',
#             'Aquificae': 'E8T2Q5', 'Gemmatimonadetes': 'C1A3G6', 'Acidobacteria': 'C1F9G6'},
#    'DnaB': {'Spirochaetes': 'F5YR56', 'Verrucomicrobia': 'D5EJK4', 'Bacteroidetes': 'Q2S6J9', 'Fibrobacteres': 'C9RMX7',
#              'Elusimicrobia': 'B2KBY0', 'Chloroflexi': 'B9KYY3', 'Planctomycetes': 'I0IDY3', 'Actinobacteria': 'H6N232',
#              'Chlorobi': 'A4SGH1', 'Thermodesulfobacteria': 'F8ABF1', 'Deferribacteres': 'E4TIR4', 'Fusobacteria': 'E3H9D4',
#              'Synergistetes': 'G7V663', 'Chrysiogenetes': 'E6W4V9', 'Tenericutes': 'A9NEP3', 'Gemmatimonadetes': 'C1A918',
#              'Nitrospirae': 'I0INN9', 'Chlamydiae': 'D6YSR9', 'Dictyoglomi': 'B8E0I7', 'Cyanobacteria': 'Q115S2',
#              'Caldiserica': 'I0GI72', 'Ignavibacteria': 'I0AHG6', 'Thermotogae': 'H2J6C2', 'Firmicutes': 'F6CKP2',
#              'Deinococcus-Thermus': 'Q1IWQ3', 'Aquificae': 'F0S203', 'Proteobacteria': 'P0ACB0', 'Acidobacteria': 'G2LJS9'},                    
#   'DnaD': {'Cyanobacteria': 'B7JZK9', 'Firmicutes': 'A3F2Y0', 'Chloroflexi': 'Q3Z916', 'Tenericutes': 'A9NF44',
#               'Actinobacteria': 'Q1B007'}, 
#   'DnaCI': {'Spirochaetes': 'B0SAC3', 'Chlamydiae': 'D6YT29', 
#           'Aquificae': 'B4U7X8', 'Chloroflexi': 'E8N2B4', 'Actinobacteria': 'C9YTH1', 'Chlorobi': 'Q3B1U7',
#           'Deferribacteres': 'D3PB17', 'Fusobacteria': 'D1AR15', 'Synergistetes': 'D5EFC2', 'Tenericutes': 'D3VRL5',
#           'Proteobacteria': 'P0AEF0', 'Nitrospirae': 'I0ILA8', 'Bacteroidetes': 'F9Z8M5', 'Cyanobacteria': 'Q7NE14',
#           'Thermotogae': 'A9BGC2', 'Firmicutes': 'B8D1A6', 'Deinococcus-Thermus': 'H8H2K4',
#           'Gemmatimonadetes': 'C1ADH4', 'Acidobacteria': 'G2LE44'},
#   'DnaG': {'Elusimicrobia': 'B2KBW9', 'Verrucomicrobia': 'D5EJD6', 'Bacteroidetes': 'F9Z3D2', 'Aquificae': 'O67465',
#            'Actinobacteria': 'Q9S1N4', 'Chlorobi': 'A1BE31', 'Deferribacteres': 'D4H0K5', 'Chrysiogenetes': 'E6W5Z0',
#            'Tenericutes': 'Q9PPZ6', 'Dictyoglomi': 'B8E0E5', 'Firmicutes': 'P47762', 'Deinococcus-Thermus': 'Q9RWR5',
#            'Proteobacteria': 'P0ABS5', 'Acidobacteria': 'G2LHR3', 'Spirochaetes': 'O83505', 'Nitrospirae': 'D8PBZ2',
#            'Chloroflexi': 'I0I4Q6', 'Planctomycetes': 'I0ICM2', 'Thermodesulfobacteria': 'F8A8B6', 'Fusobacteria': 'C7NC71',
#            'Synergistetes': 'D5EDU0', 'Gemmatimonadetes': 'C1A4A8', 'Chlamydiae': 'Q9Z6W4', 'Cyanobacteria': 'A5GKN4',
#            'Caldiserica': 'I0GJM1', 'Ignavibacteria': 'I0AGC8', 'Thermotogae': 'C5CDI3', 'Fibrobacteres': 'C9RLX5'}
#               }
    groupby="group by `OrganismID`"
    uniquify_lvl='OrganismID'
    group_lvl='Phylum'
##initialize logging
    main_logger = reinitialize_logging(basedir)
##identify reference sequences
    if all_ref_seqs:
        main_logger.info("Pre-defined reference sequence pairs:{}".format(str(all_ref_seqs)))
    for k,v in tables.iteritems():
        main_logger = reinitialize_logging(basedir)
        main_logger.info("BEGIN SEQUENCE SELECTION")
        workingdir="/Users/pczrr/Documents/Work/Project-Folders/B-subtilis-DNA-Replication/Processing/{}-alignment/".format(k)
        table_loggers = Logger()        
        logger = table_loggers.myLogger(workingdir,k)
        main_logger.info("{} specific log in {}".format(k,workingdir))
        main_logger.info("examining {} sequences for good reference sequences".format(k))
        try:
            table_ref_seqs = dict(all_ref_seqs[k])
            logger.info('found previous reference sequence dictionary: {}'.format(str(table_ref_seqs)))
        except KeyError, e:
            logger.error('did not find previous dictionary: {}'.format(e))
            table_ref_seqs = {}
        name_query=protein_names[k]
        table = k + 's'
        groups = find_groups(table,group_lvl,query=name_query)
        for p in phyla:
                 try:
                     table_ref_seqs[p]
                     main_logger.error('{} is already in reference seq dictionary for table {} with value: {}'
                     .format(p,table,str(table_ref_seqs[p])))
                     search = False
                 except KeyError:
                    main_logger.error('adding {} to search list for table {}'.format(p,table))
                    search = True
                 try:
                     groups[p]
                     cand_query=name_query + " and {}='{}'".format(group_lvl,p)
                     main_logger.info('group {} found in table {}, restricting to annotated sequences for reference search'.format(p,table))
                 except KeyError:
                     cand_query="where {}='{}'".format(group_lvl,p)
                     main_logger.info('group {} not found in table {}, not restricting reference search.'.format(p,table))
                     #logger.debug("candidate query:{}".format(cand_query))
#                 else:
#                     main_logger.error('something else went wrong in table {} for {} {}'.format(table,group_column,p))
#                     continue
                 if search:
                     logger.info("candidate query:{}".format(cand_query))
                     group_ref_seqs,cand_seqs=select_seq(table,v,groupby,cand_query)
                     cand_seqs=cand_seqs['q1']
                     cand_seq_scores={}
                     for ID,seq in cand_seqs.iteritems():
                         cand_seq_score = pairwise(table,group_ref_seqs[v],seq,ID)
#                         logger.debug('reference sequence found for group {} in table {}:{}'.format(g,table[0],ref_seqs))
                         cand_seq_scores[ID]=(cand_seq_score)
                     try:                         
                         table_ref_seqs[p]=max(cand_seq_scores.iteritems(), key=operator.itemgetter(1))[0]
                     except ValueError, e:
                         main_logger.error('no {} sequences found for {}:{}'.format(k,p,e))
                         continue
                     logger.info('current {} reference sequences: {}'.format(k,str(table_ref_seqs)))
        main_logger.info('final reference sequences for {}:{}'.format(k,str(table_ref_seqs)))
        all_ref_seqs[k]=table_ref_seqs
    main_logger.info('Complete reference sequence dictionary:\n{}'.format(str(all_ref_seqs)))
    with open('/Users/pczrr/Documents/Work/Project-Folders/B-subtilis-DNA-Replication/Processing/reference_seqs.txt','w') as outfile:
        for k,v in all_ref_seqs.iteritems():
            for i,j in v.iteritems():
                outfile.write('{} \t {} \t {} \n'.format(k,i,j))
    main_logger.info('Completed export of reference sequence dictionary.')
##flag unique sequences
    main_logger.info('begin uniquification')
    for k,v in tables.iteritems():
        main_logger = reinitialize_logging(basedir)
        table = k + 's'
        workingdir="/Users/pczrr/Documents/Work/Project-Folders/B-subtilis-DNA-Replication/Processing/{}-alignment/".format(k)
        table_loggers = Logger()        
        logger = table_loggers.myLogger(workingdir,k)
        groups = find_groups(table,group_lvl)
        for p in phyla:
            try:
                print groups[p]
                search =True
            except KeyError:
                main_logger.error('no sequences for {} {} in table {}'.format(group_lvl,p,table))
                search = False
            if search == True:
                group_query = "where {}='{}'".format(group_lvl,p)
                main_logger.info('searching for {}s from phylum {} in table {}'.format(uniquify_lvl,p,table))
                org_groups = find_groups(table,uniquify_lvl,query=group_query)
                for org,cnt in org_groups.iteritems():
                    logger.info('uniquifying {} IDs with {}s {} in table {}'.format(cnt,uniquify_lvl,table,org))
                    ref_ID = all_ref_seqs[k][p]
                    cand_query="where {}='{}'".format(uniquify_lvl,org)
                    uniquify(table,ref_ID,cand_query,org,'ref'+k,'dup'+k,flag_col='PolGroup',cand_field=uniquify_lvl,groupby='')
                    logger.info('{} {} uniquified'.format(uniquify_lvl,org))
            main_logger.info('completed uniquification of phylum {} in table {}'.format(p,table))
        main_logger.info('completed uniquification of of table {}'.format(table))

