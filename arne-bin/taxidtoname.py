#!/bin/env python3
from Bio import Entrez
import sys

def get_common_name(taxid):
    Entrez.email = "your.email@example.com"  # Always tell NCBI who you are
    handle = Entrez.efetch(db="Taxonomy", id=taxid, retmode="xml")
    records = Entrez.read(handle)
    return records[0]["ScientificName"]

taxid = sys.argv[1]
print(get_common_name(taxid))