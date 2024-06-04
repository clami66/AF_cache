#!/bin/env python3
from Bio import SeqIO
import gzip
import argparse

def extract_entries(filename):
    entries = []
    found=set()
    #print ("Test",filename)
    with gzip.open(filename, "rt") as handle:
        for record in SeqIO.parse(handle, "swiss"):
            #print ("ID:",record.annotations["accessions"])
            #print ("Name:",record.annotations["gene_name"])
            #print ("Organism:",record.annotations["organism"])
            #print ("TaxID",record.annotations["ncbi_taxid"])
            #print ("HostID",record.annotations["host_ncbi_taxid"])
            #print  ((len(record.annotations["ncbi_taxid"])) , (len(record.annotations["host_ncbi_taxid"])))
            try:   
                if ((len(record.annotations["ncbi_taxid"])>0) and (len(record.annotations["host_ncbi_taxid"])>0)):
                    #entries.append(record)
                    
                    for i in record.annotations["ncbi_taxid"]:
                        for j in record.annotations["host_ncbi_taxid"]:
                            tmp=str(i)+","+str(j)
                            if tmp not in found:
                                #entries.append(record)
                                found.add(tmp)
                                print (f"{j},{i}")
                    
            except:
                continue


    return entries


parser = argparse.ArgumentParser()
parser.add_argument("filename", help="Path to the input file")
args = parser.parse_args()

filename = args.filename
    #filename = "/proj/berzelius-2021-29/users/x_arnel/herpes/uniprot_trembl.dat.gz"  # replace with your filename
print ('"Host NCBI taxid","Virus NCBI taxid"')
entries = extract_entries(filename)
#for entry in entries:
#    print(entry.annotations["host_ncbi_taxid"])