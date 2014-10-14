# -*- coding: utf-8 -*-
"""
Created on Sat Oct 11 20:13:34 2014

@author: Richard Rymer

Returns the formal name for an organism given a taxid.
"""
from Bio import Entrez
import xml.etree.ElementTree as ET

def _retrieve_taxid_entry_from_entrez(ID):
    """
    Private function
    Retrieves the full record for the taxon with taxid ID, returns the 
    record handle
    """
    return Entrez.efetch(db='taxonomy', id=ID)

def _parse_xml(XML):
    """
    Private function
    Takes a stream of XML data as a file-like object, and return a parsed object
    """
    return ET.parse(XML)

def _retrieve_taxon_name(parsedXMLobject):
    """
    Private function
    Takes a parsed XML object, and extracts and returns the taxon name from the XML object
    """
    return parsedXMLobject.getroot()[0][1].text

def get_taxon_name(ID):
    """
    Takes a Taxonomical ID, and returns the taxon name
    """
    return _retrieve_taxon_name(_parse_xml(_retrieve_taxid_entry_from_entrez(ID)))
