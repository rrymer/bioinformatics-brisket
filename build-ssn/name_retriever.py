# -*- coding: utf-8 -*-
"""
Created on Sat Oct 11 20:13:34 2014

@author: Richard Rymer

Functions to retrieve the formal name for an organism given a taxid.
"""
from Bio import Entrez
import xml.etree.ElementTree as ET

def get_taxid_entry_from_entrez(ID):
    """
    Retrieves the full record for the taxon with taxid ID, returns the 
    record handle
    """
    return Entrez.efetch(db='taxonomy', id=ID)

def parse_xml(XML):
    """
    Takes a stream of XML data as a file-like object, and return a parsed object
    """
    return ET.parse(XML)

def get_taxon_name(parsedXMLobject):
    """
    Takes a parsed XML object, and returns the taxon name
    """
    return parsedXMLobject.getroot()[0][1].text
