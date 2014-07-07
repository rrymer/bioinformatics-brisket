# -*- coding: utf-8 -*-
"""
Created on Thu Jun 26 13:25:24 2014

@author: Richard Rymer

Connects to a database of sequences, and retrieves sequences for a given set of
organisms of a set of proteins.

Carries out all unique pairwise alignments for each protein, normalizes the a
lignment scores, and stores the results in a dictionary.

Generates a delimited text file suitable for importing into cytoscape or similar
as a network file.
"""
from Bio import pairwise2 as align
from Bio.SubsMat import MatrixInfo as matlist
import multiprocessing as mp
import MySQLdb as mdb
import MySQLdb.cursors as cursors

def connect_to_mysqldb(host, db, un, cursor_to_use = cursors.SSDictCursor):
    """
    takes the connection parameters for connecting to a MySQL database of sequences, and 
    returns connection and cursor objects through MySQLdb
    """
    connection = mdb.connect(host=host,db=db,user=un,cursorclass=cursor_to_use)
    cursor = connection.cursor()
    return connection,cursor

def close_connection(connection):
    """
    closes the argument connection
    """
    connection.close()

def get_groups(table,group_field,query=''):
    """
    returns the hits for group_field in a given table, table.
    """
    sql="select {},count(*) from {} {} group by {}"
    groups = run_query(sql.format(group_field,table,query,group_field))
    return groups

def run_query(query):
    """
    runs the argument sql query, and returns the results in a tuple of dictionaries format
    """
    connection, cursor = connect_to_mysqldb(host,db,un)
    cursor.execute(query)
    results = cursor.fetchall()
    close_connection(connection)
    return results

def get_taxa(n=100, table= 'RepPols', query = ''):
    """
    returns a list of NCBI taxids of length n, default 100.
    Other options include selection of a table with table = '<table name>' and
    an SQL 'where' statement with query = '<where ...>'
    """
    results = run_query( 'select {},{} from {} {} group by {} order by rand() limit {}'.format(taxid_field,taxon_name_field,table,query,taxid_field,n) )
    taxa = {}
    for r in results:
        taxa[r[taxid_field]]=r[taxon_name_field]
    return taxa

def get_table(protein_name):
    """
    returns the table associated with a protein name in an SQL database
    """
    return protein_name_table_map[protein_name]

def get_unique_reference(protein_name):
    """
    returns the uniquifier reference string for a protein named protein_name
    """
    return protein_uniquifier_references[protein_name]
    

def get_seqs(taxids,protein_name):
    """
    takes a list of taxonomic ids and the name of a protein, and returns the amino acid sequences along with identifier information
    """
    table = get_table(protein_name)
    return run_query("select {},{},{} from {} where {} = '{}' and {} in {}"
            .format(taxid_field,seq_field,id_field,table,uniquifier_field,get_unique_reference(protein_name),taxid_field,taxids))
        
def do_pairwise(ref_seq,query_seq):
    """
    takes two single letter amino acid sequences, aligns them, and returns an alignment score 
    based on the BLOSUM100 score matrix
    """
    matrix = matlist.blosum100
    score_list=[]
    algnmt = align.align.globalds(ref_seq,query_seq,matrix,-2,-1)
    for a in algnmt:
        score_list.append(a[2])
#    print max(score_list)
    return max(score_list)

def do_mp_pairwise(seq1,seq2,mp_scores_dict,id_field):
    """
    takes two single letter amino acid sequences, aligns them, and returns an alignment score 
    based on the BLOSUM100 score matrix, multithread safe, requires a shared dict object
    """
    matrix = matlist.blosum100
    score_list=[]
    algnmt = align.align.globalds(seq1.get_sequence(),seq2.get_sequence(),matrix,-2,-1)
    for a in algnmt:
        score_list.append(a[2])
#    print max(score_list)
    mp_scores_dict[(seq1.get_identifier(id_field),seq2.get_identifier(id_field))] = max(score_list)

class Sequence(object):
    def __init__(self,sequence_record):
        """
        initializes a sequence record object with a sequence record from a database in dictionary format
        """
        self.__dict__.update(sequence_record)
    
    def get_sequence(self):
        """
        returns the aa sequence for the record
        """
        return self.__getattribute__(seq_field)
    
    def get_identifier(self,identifier):
        """
        returns the value of the attribute named in the identifer
        """
        return self.__getattribute__(identifier)

    def __str__(self):
        print 'Sequence record object of length', len(self.get_identifier(seq_field))
        for object in dir(self):
            print object

class Protein(object):
    def __init__(self,name,**identifiers):
        """
        Initializes a Protein object with name name, and an arbitrary number of identifiers
        identifiers is a dictionary of field:value mappings, for example, a protein may be
        associated with a particular table in a database, in which case the entry would be
        table:<table name>
        """
        self.name = name
        self.__dict__.update(identifiers)
        
    def get_name(self):
        """
        returns the name of the protein associated with the object
        """
        return self.name
            
    def get_all_taxa_info(self):
        """
        Return a dict of taxa associated with this protein
        """
        raise NotImplementedError
        return self.taxa
    
    def pull_sequences(self,taxids):
        """
        retrieve sequences for the protein from a data base or file
        """
        return get_seqs(taxids,self.name)
    
    def get_sequences(self,Targets):
        """
        assemble sequences into sequence objects, and return a list of the sequence
        object records
        """
        self.sequences = []
        for sequence in self.pull_sequences(self.get_taxids_by_group(group_field)):
            self.sequences.append(Sequence(sequence))
        return self.sequences
    
    def get_unique_pairs(self):
        """
        determine the unique pairs of sequences in the data set, and return a dictionary
        mapping of the sequence pair as a tuple in the dict keys
        """
        sequences = self.get_sequences()
        copy = sequences[:]
        pairs = {}
        for sequence1 in sequences:
            for sequence2 in copy:
                if sequence1 is sequence2:
                    continue
                try:
                    pairs[(sequence1,sequence2)]
                except KeyError:
                    pairs[(sequence1,sequence2)] = True
                    pairs[(sequence2,sequence1)] = True
        self.pairs = pairs
        return self.pairs                    

    def mp_score_pairs(self,id_field):
        """
        conduct pairwise sequence alignments for all unique pairs, and return a score dict
        mapping sequence pair to alignment score, multithreaded version
        """
        p = mp.Pool(8)
        man = mp.Manager()
        scores_dict = man.dict()
        for pair in self.get_unique_pairs():
            p.apply_async(do_mp_pairwise, args=(pair[0],pair[1],scores_dict,id_field))
        p.close()
        p.join()
        return self.calculate_scores(scores_dict)
                
    def score_pairs(self,id_field):
        """
        conduct pairwise sequence alignments for all unique pairs, and return a score dict
        mapping sequence pair to alignment score, serial version
        """
        scores_dict = {}
        for pair in self.get_unique_pairs():
            scores_dict[(pair[0].get_identifier(id_field),pair[1].get_identifier(id_field))] = do_pairwise(pair[0].get_identifier(seq_field),pair[1].get_identifier(seq_field))
        return self.calculate_scores(scores_dict)
        
    def calculate_scores(self,scores_dict):
        """
        normalize scores to the min and max of the data, assign scores as a dictionary object
        attribute
        """        
        scores = {}
        max_score = max(scores_dict.values())
        min_score = min(scores_dict.values())
        for k in scores_dict.keys():
            v = scores_dict[k]
            score = ((v - min_score)/(max_score - min_score))
            if score >= cutoff:
                scores[k] = score
        del scores_dict
        self.scores = scores

    def get_scores(self):
        """
        return a dictionary of sequences pairs and scores
        """        
        return self.scores

    def make_edge_tuples(self):
        """
        assemble a list of tuples containing the ids of the sequence pairs, the associated
        protein name, and the score value after normalization.  Store this list as an 
        attribute of the object
        """
        tuples = []
        for k,v in self.get_scores().iteritems():
            tuples.append((k[0], k[1], self.name, v))
        self.edge_tuples = tuples
    
    def get_edge_tuples(self):
        """
        return a list of tuples, where each tuple contains the data for a single edge.
        Format: (source node, target node, protein name, normalized alignment score)
        """
        return self.edge_tuples        

class Targets(object):
    def __init__(self, group_field, taxid_field, n, default_table):
        self.group_field = group_field
        self.taxid_field = taxid_field
        self.n = n
        self.default_table = default_table
        pass
    
    def pull_base_taxids_by_group(self):
        """
        returns a list of n taxids for which there is least 1 sequence of the protein
        available in the defined group_field, where group_field is some field in a database.
        n is an int.
        """
        results = get_groups(self.default_table, self.group_field)
        taxa_dict = {}
        for result in results:
            taxa_dict = dict(taxa_dict.items() + get_taxa(self.n, self.default_table, query = "where {} = '{}' ".format(self.group_field,result[self.group_field])).items())
        self.taxids = tuple(taxa_dict.keys())
        self.taxa = taxa_dict
    
    def get_taxids(self,name):
        taxa_dict = {}
        for result in results:
            taxa_dict = dict(taxa_dict.items() + get_taxa(self.n, self.default_table, query = "where {} in '{}' ".format(self.taxid_field,self.taxids)).items())
        self.taxids = tuple(taxa_dict.keys())
        self.taxa = taxa_dict
        
if __name__ == '__main__':
##Defs and params
    host = 'localhost'
    db = 'FamilyCpolymerases'
    un = 'sequser'
    id_field = 'Accession_Number'
    seq_field = 'Sequence'
    taxid_field = 'OrganismID'
    taxon_name_field = 'Organism'
    uniquifier_field = 'PolGroup'
    group_field = 'Phylum'
    n = 1
    cutoff = 0.0
    finished = {}
    protein_name_table_map={'PC':'RepPols','DnaA':'DnaAs','DnaG':'DnaGs','DnaC':'DnaCs','DnaD':'DnaDs','DnaB':'DnaBs',
                       'E1':'RepPols','E2':'RepPols','E3':'RepPols'}
    protein_uniquifier_references = {'DnaG':'refDnaG','DnaC':'refDnaC','DnaD':'refDnaD','DnaA':'refDnaA','DnaB':'refDnaB',
                       'PC':'PC','E1':'E1','E2':'E2','E3':'E3'}
    debug = False
    outfile = '/tmp/test_data.txt'
##initialize output file:
    with open(outfile,'w') as fout:
        fout.write('{}\t{}\t{}\t{}\n'.format('source','target','protein','alignment_score'))
## MAIN ##
    if debug:
        import pdb
        pdb.set_trace()
        protein_name_table_map = {'PC':'RepPols'}
#    taxa_id_names_dict = get_taxa()
#    taxids = tuple(taxa_id_names_dict.keys())
    proteins = []
    #iterate through a dictionary of protein names and associated sequence tables, and assemble is list of Protein objects
    for name,seq_table in protein_name_table_map.iteritems(): 
        proteins.append(Protein(name,table = seq_table))
    #iterate over Protein objects, calculate, and then write out edges to a tab-separated text file
    for protein in proteins:
        print 'calculating scores for protein',protein.get_name()
        protein.mp_score_pairs(taxid_field)
        protein.make_edge_tuples()
        with open(outfile,'a') as fout:
            for edge_tuple in protein.get_edge_tuples():
                fout.write('{}\t{}\t{}\t{}\n'.format(protein.get_all_taxa_info()[edge_tuple[0]],protein.get_all_taxa_info()[edge_tuple[1]],edge_tuple[2],edge_tuple[3]))
        
