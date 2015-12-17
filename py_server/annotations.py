# -*- coding: utf-8 -*-
"""
Different methods to read annotations
@author: rodri
"""

#%% -------------------------- ADDITIONAL ANNOTATIONS ----------------------

#%%
def gff(filename="genomes/annotations/spombe/gff/schizosaccharomyces_pombe.I.gff3"):
    f=open(filename)
    import csv
    cad=f.readline()
    skip=0
    while(cad.startswith("#")==True):
        cad=f.readline()
        skip+=1
    tam=len(f.readlines())+1
    f.seek(0)
    reader=csv.DictReader(f, delimiter="\t")
    reader.fieldnames=["chromosome", "source", "type", "start", "end", "xx", "sense", "xx", "id"]
    
    import numpy as np
    data=np.empty(tam,dtype=[("chromosome", "a2"),("type", "a40"), ("start", "i8"), ("end", "i8"), ("sense", "a1"), ("id", "a200")])
    for i in range(skip):
        next(reader)
    for i in xrange(tam):
        row=reader.next()
        tow=data[i]
        #tow["start"]=(int)(row["start"])
        #tow["end"]=(int)(row["end"])
        tow["start"]=(row["start"])
        tow["end"]=(row["end"])
        tow["chromosome"]=row["chromosome"]
        tow["type"]=row["type"]
        tow["sense"]=row["sense"]
        tow["id"]=row["id"]
    return data
          

#%%  Loads the gene ontology into a dic where keys are go_ids and values are just go names by now
def go(filename="genomes/annotations/go/gene_ontology_ext.obo"):
    f=open(filename)
    lines=f.readlines()
    data={}
    i=0
    while i<len(lines):
        if(lines[i].startswith("[Term]")):
            go_id=lines[i+1].replace("id:", "").strip()
            go_name=lines[i+2].replace("name:", "").strip()
            data[go_id]=go_name
        i+=1
    return data


#%% Loads the gene ontology annottion of a given species into an array which only stores gene_id, gene_name, go_id, go_type and gene_desc
def goa(filename="genomes/annotations/spombe/goa/gene_association.pombase"):
    f=open(filename)
    lines=f.readlines()
    data=[]
    for l in lines:
        if(l.startswith("!")==False):
            vals=l.split("\t")
            data.append({"gene_id":vals[1],"gene_name":vals[2], "go_id":vals[4], "go_type":vals[8], "gene_desc":vals[9]})
    return data
    

#%% Loads the fasta sequence of a given file (by now only working for S pombe files)
def fasta(ch):
    f=open("genomes/annotations/spombe/fasta/chromosome"+(str)(ch)+".fasta")
    reader=f.readlines()
    reader=reader[1:]
    seq=""
    for i in xrange(len(reader)):
        seq=seq+reader[i].replace("\n", "")
    return seq
    
#%%
def annotate(mm, dataGFF, types=["any"], ws=1000):
    import numpy as np
    print len(dataGFF)
    data2=dataGFF
    if(types[0]!="any"):
        wanted_set = set(types)  # Much faster look up than with lists, for larger lists
        @np.vectorize
        def selected(elmt): return elmt in wanted_set  # Or: selected = numpy.vectorize(wanted_set.__contains__)
        data2=dataGFF[selected(dataGFF["type"])]
    print len(data2)
    em={} #enriched (i.e. detailed) matches, including for each the thigs found at GFF
    interval=ws*0.5
    import re
    for x in mm:
        sel=data2[(data2["start"]>x-interval) & (data2["end"]<x+interval)]
        if(len(sel)>0):     
            em[x]=[]
            for s in sel:
                gid=s["id"] 
                gid=re.sub(".*SP", "SP", gid)
                gid=re.sub(":.*", "", gid)
                gid=re.sub(";.*$", "", gid)
                gid=re.sub(".1$", "", gid)
                em[x].append({"type":s["type"], "id":gid, "start":s["start"], "end":s["end"], "sense":s["sense"]})
    return em
