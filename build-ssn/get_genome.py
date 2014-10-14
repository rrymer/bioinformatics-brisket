# -*- coding: utf-8 -*-
"""
Created on Sat Oct 11 21:52:32 2014

@author: Richard Rymer

Retrieves the full genome from Entrez nucleotide
"""

from urllib2 import HTTPError
from Bio import SeqIO, Entrez

class InputError(Exception):
    pass

class NoMatchingSequence(Warning):
    pass

def find_genome_gi(ORG_NAME,restr=''):
    """
    Takes an organism name, and returns a gi ID for one genome associated with
    the organism
    """
    records = Entrez.read(Entrez.esearch(db="nucleotide",term=ORG_NAME + "[ORGN] complete genome " + restr))
    try:
        return records['IdList'][0]
    except IndexError:
        return

def get_genome_seq(ORG_NAME):
    """
    Takes an organism name, and returns a SeqIO record object, or raises a warning,
    and returns None if no sequence was found
    """
    id = find_genome_gi(ORG_NAME,restr='refseq')
    if not id:
        id = find_genome_gi(ORG_NAME)
    if not id:
        raise NoMatchingSequence, 'No sequence was found for organism: ' + ORG_NAME
        return
    try:
        seq_record = SeqIO.read(Entrez.efetch(db="nucleotide", id=id,rettype="fasta", retmode="text"),\
                                'fasta')                        
        if len(seq_record.seq._data) < 1:
            raise NoMatchingSequence, 'No sequence was found for organism: ' + ORG_NAME
            return
    except HTTPError:
        raise NoMatchingSequence, 'No sequence was found for organism: ' + ORG_NAME
        return
    
    return seq_record