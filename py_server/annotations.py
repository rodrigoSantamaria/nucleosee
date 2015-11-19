# -*- coding: utf-8 -*-
"""
Different methods to read annotations
@author: rodri
"""

#%% -------------------------- ADDITIONAL ANNOTATIONS ----------------------

def gff():
    f=open("/Users/rodri/Documents/investigacion/IBFG/visualization/spombe/gff/schizosaccharomyces_pombe.I.gff3")
    import csv
    reader=csv.DictReader(f, delimiter="\t")
    next(reader)
    reader.fieldnames=["chromosome", "source", "type", "start", "end", "xx", "sense", "xx", "id"]
    data=[]
    for row in reader:
        data.append(row)
    return data

            

#%%
def go():
    f=open("/Users/rodri/Documents/investigacion/IBFG/visualization/spombe/goa/gene_association.pombase")
    lines=f.readlines()
    data=[]
    for l in lines:
        if(l.startswith("!")==False):
            vals=l.split("\t")
            data.append({"gene_id":vals[1],"gene_name":vals[2], "go_id":vals[4], "go_type":vals[8], "gene_desc":vals[9]})
    return data


#%%
def goa():
    f=open("/Users/rodri/Documents/investigacion/IBFG/visualization/spombe/goa/gene_association.pombase")
    lines=f.readlines()
    data=[]
    for l in lines:
        if(l.startswith("!")==False):
            vals=l.split("\t")
            data.append({"gene_id":vals[1],"gene_name":vals[2], "go_id":vals[4], "go_type":vals[8], "gene_desc":vals[9]})
    return data
    

#%%
def fasta(ch):
    f=open("/Users/rodri/Documents/investigacion/IBFG/FASTA/chromosome"+(str)(ch)+".fasta")
    reader=f.readlines()
    reader=reader[1:]
    seq=""
    for i in xrange(len(reader)):
        seq=seq+reader[i].replace("\n", "")
    return seq
    
