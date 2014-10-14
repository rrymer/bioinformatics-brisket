# -*- coding: utf-8 -*-
"""
Created on Sat Oct 11 20:40:31 2014

@author: Richard Rymer

taxon name retriever unittests
"""

import unittest
import name_retriever as nr
import get_genome
from Bio import SeqIO

class RetrieverTests(unittest.TestCase):
    def setUp(self):
        self.dbs = ['pubmed', 'protein', 'nuccore', 'nucleotide', 'nucgss', 
        'nucest', 'structure', 'genome', 'assembly', 'genomeprj', 'bioproject',
        'biosample', 'blastdbinfo', 'books', 'cdd', 'clinvar', 'clone', 'gap',
        'gapplus', 'dbvar', 'epigenomics', 'gene', 'gds', 'geoprofiles',
        'homologene', 'medgen', 'journals', 'mesh', 'ncbisearch', 'nlmcatalog',
        'omim', 'orgtrack', 'pmc', 'popset', 'probe', 'proteinclusters',
        'pcassay', 'biosystems', 'pccompound', 'pcsubstance', 'pubmedhealth',
        'seqannot', 'snp', 'sra', 'taxonomy', 'toolkit', 'toolkitall',
        'toolkitbook', 'unigene', 'gencoll', 'gtr']
        nr.Entrez.email = 'richard.rymer@rrcc.edu'
        self.handle = nr.get_taxid_entry_from_entrez('316385')
        self.test_name = 'Escherichia coli str. K-12 substr. DH10B'
    
    def test_00_name_retriever_get_taxon_name(self):
        """
        Verifies that the Entrez taxonomy database can be communicated with,
        and a search of taxid 316385 return E. coli
        """
        self.assertEqual(nr.get_taxon_name(nr.parse_xml(self.handle)),self.test_name)
    
    def test_01_get_genome_gets_correct_gi(self):
        self.assertEqual(get_genome.find_genome_gi(self.test_name,restr='refseq'),'170079663')
    
    def test_02_get_genome_returns_correct_object(self):
        self.assertIsInstance(get_genome.get_genome_seq(self.test_name),SeqIO.SeqRecord)
    
    def test_03_get_genome_returns_correct_genome(self):
        self.assertEqual(get_genome.get_genome_seq(self.test_name).name,'gi|170079663|ref|NC_010473.1|')
    
    def test_04_get_genome_returns_full_genome(self):
        self.assertEqual(len(get_genome.get_genome_seq(self.test_name).seq._data),4686137)
        
        