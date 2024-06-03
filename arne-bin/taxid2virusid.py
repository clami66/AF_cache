#!/bin/env python3

import pickle
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import argparse
from sys import argv

parser = argparse.ArgumentParser(description="Replace Host ID with Virus IDS")
parser.add_argument("--infile", help = "Path to the directory containing the taxdata ",default="/proj/beyondfold/apps/colabfold_databases/seqid2taxid.pkl")
parser.add_argument("--csv",nargs="+",help="CSV files with changes",required=True)
parser.add_argument("--out",nargs="+",help="Output files",required=True)
parser.add_argument("--duplicate",action="store_true",help="Allow duplicate entries",default=False)
args, unknownargs = parser.parse_known_args()


if len(args.csv) != len(args.out):
    raise ValueError("Number of CSV files and output files must be the same")
    sys.exit(1)

with open(args.infile, "rb") as file:
    taxdata = pickle.load(file)
df_taxdata = pd.DataFrame.from_dict(taxdata, orient='index', columns=['TaxID'])


for i in range(len(args.csv)):
    df=df_taxdata.copy()
    #del taxdata
    links = pd.read_csv(args.csv[i])
    if args.duplicate:
        newlinks=pd.DataFrame({"Virus NCBI taxid":[],"Host NCBI taxid":[]}).astype("int64")
        for j in links["Virus NCBI taxid"].unique() :
            l=links.loc[links["Virus NCBI taxid"]==j]["Host NCBI taxid"].to_list()
            string=",".join(str(x) for x in l)
            newlinks=newlinks.append({"Virus NCBI taxid":j,"Host NCBI taxid":string},ignore_index=True)            
            #newlinks=newlinks.append({"Virus NCBI taxid":j,"Host NCBI taxid":l},ignore_index=True)
        links=newlinks
        print (links)    
    for index, row in links.iterrows():
        print (index,row)
        df.loc[df.TaxID==row["Virus NCBI taxid"],"TaxID"]=row["Host NCBI taxid"]
    taxdata_dict = df.to_dict("dict")["TaxID"]
    testfile=args.out[i]+".csv"
    df.to_csv(testfile)
    with open(args.out[i], "wb") as file:
        pickle.dump(taxdata_dict, file)
    del taxdata_dict
#
#df=df_taxdata.copy()    
#family = pd.read_csv("GenusLinks.csv")
#for index, row in family.iterrows():
#    df.loc[df.TaxID==row["Virus NCBI taxid"],"TaxID"]=row["Host NCBI taxid"]
#taxdata_dict = df.to_dict("dict")["TaxID"]
#with open("genus/seqid2taxid.pkl", "wb") as file:
#    pickle.dump(taxdata_dict, file)
#del taxdata_dict
#
#
#df=df_taxdata.copy()
#family= pd.read_csv("OrderLinks.csv")
#for index, row in family.iterrows():
#    df.loc[df.TaxID==row["Virus NCBI taxid"],"TaxID"]=row["Host NCBI taxid"]
#taxdata_dict = df.to_dict("dict")["TaxID"]
#with open("order/seqid2taxid.pkl", "wb") as file:
#    pickle.dump(taxdata_dict, file)
#del taxdata_dict
#
#df=df_taxdata.copy()  
#family = pd.read_csv("SubFamilyLinks.csv")
#for index, row in family.iterrows():
#    df.loc[df.TaxID==row["Virus NCBI taxid"],"TaxID"]=row["Host NCBI taxid"]
#taxdata_dict = df.to_dict("dict")["TaxID"]
#with open("subfamily/seqid2taxid.pkl", "wb") as file:
#    pickle.dump(taxdata_dict, file)
#del taxdata_dict
#           
#
