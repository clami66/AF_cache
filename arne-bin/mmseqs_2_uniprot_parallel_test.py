import re
import os
import glob
import pickle
import string
import base64
import hashlib
import argparse
from sys import argv
from pathlib import Path
from functools import partial
from multiprocessing import Pool, TimeoutError
from Bio import AlignIO
import pandas as pd

seqid_to_id = dict()

def get_accession2id(db_path, taxid=True):
    if taxid:
        pkl_path = f"{db_path}/seqid2taxid.pkl"
    else:
        pkl_path = f"{db_path}/seqid2repid.pkl"

    if not os.path.exists(pkl_path):
        print("Generating accession to taxid dictionary")
        make_accession2id(db_path)
    
    print("Loading accession to taxid dictionary")
    with open(pkl_path, "rb") as input_file:
        return pickle.load(input_file)

def make_accession2id(db_path):
    print("Generating accession2taxid pickle file")
    seqid_to_taxid_dic = {}
    seqid_to_repid_dic = {}
    header_files = glob.glob(f"{db_path}/*_h.tsv")
    ptax = re.compile("TaxID=([0-9]+)")
    prep = re.compile("RepID=[A-Z0-9]+[_]([A-Z0-9]+)")
    
    for header_file in header_files:
        print(f"Looking for taxids in {header_file}")
        with open(header_file) as headers:
            for header in headers:
                match_taxid = ptax.search(header)
                match_repid = prep.search(header)
                accession = header.split()[1]
                if match_taxid:
                    taxid = match_taxid.group(1)
                    seqid_to_taxid_dic[accession] = int(taxid)
                if match_repid:
                    repid = match_repid.group(1)
                    seqid_to_repid_dic[accession] = repid
    
    with open(seqid_to_taxid, 'wb') as f:
      pickle.dump(seqid_to_taxid_dic, f, protocol=4)
    with open(seqid_to_repid, 'wb') as f:
      pickle.dump(seqid_to_repid_dic, f, protocol=4)

def convert_alignment(in_alignment, out_dir, taxid,duplicate=False,shortname=False):
    seqids = set()
    tolower = str.maketrans('', '', string.ascii_lowercase)
    print(f"Opening {in_alignment}")
    with open(in_alignment) as aln:
        a3m_data = iter(aln.readlines())
    print(f"{in_alignment} loaded")
    # always write first (target) sequence to file
    target_header = next(a3m_data)
    target_seq = next(a3m_data)
    
    # prepare the directory structure as AF wants it
    target_id = target_header.strip().strip(">")
    print(f"Target ID: {target_id}")
    Path(args.out_dir, target_id, "msas", "A").mkdir(parents=True, exist_ok=True)
    pseudo_uniprot = open(f"{args.out_dir}/{target_id}/msas/A/uniprot_hits.a3m", "w")
    pseudo_uniprot.write(target_header)
    pseudo_uniprot.write(target_seq)
    print(f"Starting conversion: {target_id}")
    #duplicateentries=False
    memid=None
    if (duplicate):     
        for line in a3m_data:
            if line.startswith(">"):
                if memid:
                    for i in range(len(memid)):
                        pseudo_uniprot.write(uniprotid[i])
                        pseudo_uniprot.write(sequence)   
                memid = []
                seqid = line.split()[0].strip(">") # gets accession ID
                tempid=[]
                #print ("Testing1 seqid",seqid,len(seqid_to_id))
                longseqid="UniRef100_"+seqid
                shortseqid=re.sub("UniRef100_", "", seqid)
                if  (longseqid not in seqids) and  (shortseqid not in seqids) and (seqid not in seqids): # avoids duplicates
                    #print ("Adding seqid",seqid)
                    seqids.add(seqid)
                    seqids.add(longseqid)
                    seqids.add(shortseqid)
                else:
                    continue
                #print ("Testing2 seqid",seqid)                  
                if seqid in seqid_to_id:
                    try:
                        memid = seqid_to_id[seqid].split(",")
                    except:
                        memid.append(str(seqid_to_id[seqid]))
                    #memid = seqid_to_id[seqid]
                    tempid=memid.copy()
                    #print ("updating seqid1",seqid,memid)
                elif longseqid in seqid_to_id:
                    try:
                        memid = seqid_to_id[longseqid].split(",")
                    except:
                        memid.append(str(seqid_to_id[longseqid]))
                    tempid=memid.copy()
                    #print ("updating seqid2",longseqid,memid)
                elif shortseqid in seqid_to_id: 
                    try:
                        memid = seqid_to_id[shortseqid].split(",")
                    except:
                        memid.append(str(seqid_to_id[shortseqid]))
                    tempid=memid.copy()
                    #print ("updating seqid3",shortseqid,memid)
                #else:
                #    print ("no match:",seqid,memid,len(seqid_to_id))
                if taxid:
                    if (memid):
                        j=0
                        for i in memid:
                            memid[j] = base64.urlsafe_b64encode(hashlib.md5(str(i).encode('utf-8')).digest()).decode("utf-8").replace("_", "").replace("-", "")[:5].upper()
                            j+=1                    
                    else:
                        memid=["NOSPECIEFOUND"]      
                        tempid=["NOSPECIEFOUND"]      
                if memid:
                    uniprotid=[]
                    for i in range(len(memid)):
                        m=memid[i]
                        t=tempid[i]   
                        seqid_alpha=re.sub(r"[^a-zA-Z0-9]", '', seqid)+str(i)
                        if shortname:
                            seqid_alpha=base64.urlsafe_b64encode(hashlib.md5(str(seqid_alpha).encode('utf-8')).digest()).decode("utf-8").replace("_", "").replace("-", "")[:7].lower()
                        #pseudo_uniprot.write(f">tr|{seqid_alpha}|{seqid_alpha}_{memid}/1-{len(target_seq)} Seq:{seqid} Temp:{tempid} mem:{memid}\n")
                        #pseudo_uniprot.write(f">tr|{seqid_alpha}|{seqid_alpha}_{memid}/1-{len(target_seq)} Memid:{tempid}\n")
                        uniprotid.append(f">tr|{seqid_alpha}|{seqid_alpha}_{m}/1-{len(target_seq)} TaxID:{t} Seqid:{seqid}\n")
                sequence=''
            elif memid:
                #pseudo_uniprot.write(line.translate(tolower))
                sequence+=line.translate(tolower)
        if memid:
            for i in range(len(memid)):
                pseudo_uniprot.write(uniprotid[i])
                pseudo_uniprot.write(sequence)    
    else:
        for line in a3m_data:
            if line.startswith(">"):
                memid = None
                seqid = line.split()[0].strip(">") # gets accession ID
                tempid=""
                #print ("Testing1 seqid",seqid,len(seqid_to_id))
                longseqid="UniRef100_"+seqid
                shortseqid=re.sub("UniRef100_", "", seqid)
                if  (longseqid not in seqids) and  (shortseqid not in seqids) and (seqid not in seqids): # avoids duplicates
                    #print ("Adding seqid",seqid)
                    seqids.add(seqid)
                    seqids.add(longseqid)
                    seqids.add(shortseqid)
                else:
                    continue
                #print ("Testing2 seqid",seqid)                  
                if seqid in seqid_to_id:
                    memid = seqid_to_id[seqid]
                    tempid=memid
                    #print ("updating seqid1",seqid,memid)
                elif longseqid in seqid_to_id:
                    memid = seqid_to_id[longseqid]
                    tempid=memid
                    #print ("updating seqid2",seqid,memid)
                elif shortseqid in seqid_to_id: 
                    memid = seqid_to_id[shortseqid]
                    tempid=memid
                    #print ("updating seqid3",seqid,memid)
                #else:
                #    print ("no match:",seqid,memid,len(seqid_to_id))
                if taxid:
                    if (memid):
                        memid = base64.urlsafe_b64encode(hashlib.md5(str(memid).encode('utf-8')).digest()).decode("utf-8").replace("_", "").replace("-", "")[:5].upper()
                    else:
                        memid="NOSPECIEFOUND"      
                if memid:
                    seqid_alpha = re.sub(r"[^a-zA-Z0-9]", '', seqid)
                    if shortname:
                        seqid_alpha=base64.urlsafe_b64encode(hashlib.md5(str(seqid_alpha).encode('utf-8')).digest()).decode("utf-8").replace("_", "").replace("-", "")[:7].lower()
                    #pseudo_uniprot.write(f">tr|{seqid_alpha}|{seqid_alpha}_{memid}/1-{len(target_seq)} Seq:{seqid} Temp:{tempid} mem:{memid}\n")
                    pseudo_uniprot.write(f">tr|{seqid_alpha}|{seqid_alpha}_{memid}/1-{len(target_seq)} Taxid:{tempid} Seqid:{seqid}\n")
      
            elif memid:
                pseudo_uniprot.write(line.translate(tolower))



                
    pseudo_uniprot.close()

    input_handle  = open(f"{out_dir}/{target_id}/msas/A/uniprot_hits.a3m", "r")
    output_handle = open(f"{out_dir}/{target_id}/msas/A/uniprot_hits.sto", "w")

    alignments = AlignIO.parse(input_handle, "fasta")
    AlignIO.write(alignments, output_handle, "stockholm")

    output_handle.close()
    input_handle.close()
    
    # bfd symlink
    src = os.path.abspath(in_alignment)
    dst = f"{os.path.abspath(out_dir)}/{target_id}/msas/A/bfd_uniref_hits.a3m"

    if os.path.exists(dst):
        os.remove(dst)
    os.symlink(src, dst)


    
def main(args):
    #print ("ARGS",args)
    taxid = False if args.repid else True
    #seqid_to_id = get_accession2id(args.db_path, taxid=taxid)
    print ("SEQID:",taxid,len(seqid_to_id))
    a3ms = glob.glob(f"{args.in_path}/*.a3m")
    #print(f"Converting alignments: {a3ms}")
    with Pool(processes=int(args.n_cpu)) as pool:
        pool.map(partial(convert_alignment, out_dir=args.out_dir, taxid=taxid, 
                         duplicate=args.duplicate,shortname=args.shortname), a3ms)
    #for aln in a3ms:
    #    print(f"Converting alignments: {aln}")
    #    convert_alignment(aln,args.out_dir,seqid_to_id,taxid)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Converts an a3m alignment (e.g. MMseqs2 alignment) to a Stockholm alignment \
                                 that can be used for MSA pairing in AlphaFold-multimer runs")
    parser.add_argument("in_path", help = "Path to the a3m alignment file directory")
    parser.add_argument("out_dir", help = "Path to output directory")
    parser.add_argument("db_path", help = "Path to ColabFold databases")
    parser.add_argument("--n_cpu", default=8, type=int)
    parser.add_argument("--taxid", action="store_true", default=True, help = "If tax IDs should be used in the pairing procedure")
    parser.add_argument("--repid", action="store_true", default=False, help = "If mnemonic species IDs should be used in the pairing procedure")
    parser.add_argument("--duplicate", action="store_true", default=False, help = "If we should add duplicate species to the uniprot file")
    parser.add_argument("--shortname", action="store_true", default=False, help = "To replace long names in the databasw with 7char names - necessary for pairing")
    args = parser.parse_args()
    taxid = False if args.repid else True    
    seqid_to_id = get_accession2id(args.db_path, taxid=taxid)
    #print ("tt",len(seqid_to_id))
    main(args)
