# -*- coding: utf-8 -*-
"""
Different methods to read annotations
@author: rodri
"""

#%% -------------------------- ADDITIONAL ANNOTATIONS ----------------------

#%%
def gff(filename="genomes/annotations/spombe/gff/schizosaccharomyces_pombe.I.gff3"):
    import time
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
    for i in range(tam):
        row=reader.next()
        data[i]["start"]=(int)(row["start"])
        data[i]["end"]=(int)(row["end"])
        data[i]["chromosome"]=row["chromosome"]
        data[i]["type"]=row["type"]
        data[i]["sense"]=row["sense"]
        data[i]["id"]=row["id"]

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
    

#%% Loads the fasta sequence of a given file (by now only working forS pombe files)
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
    
#%%
def annotateGO(em, dataGOA, discard=['GO:0003674','GO:0005575','GO:0008150']):
    ego={}
    for x in dataGOA:
        if((x["go_id"] in discard) ==False):
            for y in em.iterkeys():
                if(len(em[y])>0 and em[y][0]["id"]==x["gene_id"]):
                    #print "{} -> {}({})".format(x["gene_id"], x["go_id"], x["go_type"])
                    if(x["go_id"] in ego.keys()):
                        ego[x["go_id"]].add(x["gene_id"])
                    else :
                        ego[x["go_id"]]=set()
                        ego[x["go_id"]].add(x["gene_id"])
    return ego

# Check GO terms with names in the matches
def annotateGOnames(em, dataGOA, dataGO):
    ego=annotateGO(em, dataGOA)
    egon={}
    for k in ego.keys():
        try:
            egon[k]={"name": dataGO[k], "genes": (list)(ego[k])} #'genes' must be a list to be JSON serializable
        except KeyError:
            print 'Key {} not found'.format(k)
    return egon

#%% Fisher's enrichment
#According to https://pypi.python.org/pypi/fisher/0.1.4
# gis is a set of genes of interest by id
# th is the threshold
def enrichmentFisher(gis, dataGOA, th=0.01, correction="none"):
    #0) Prepare sets    
    # Retrieve a dict where k=go id and value=set of genes
    goids=[x["go_id"] for x in dataGOA]
    goterms={}
    for x in goids:
        goterms[x]=set()
    for x in dataGOA:
        goterms[x["go_id"]].add(x["gene_id"])
    # Compute universe genes
    unigenes=set()
    for k in goterms.keys():
            unigenes |= set(goterms[k])
    unigenes=unigenes-gis
    uni=len(unigenes)
    sel=len(gis)
    #
    #1) Fisher's test
    from fisher import pvalue
    pvals={}
    for k in goterms.keys():
    
        selgo=len(gis.intersection(goterms[k])) #number of gis in the term
        selnogo=sel-selgo
    
        unigo=len(goterms[k])-selgo #number of non-gis in the term
        uninogo=uni-unigo
    
        p = pvalue(unigo, uninogo, selgo, selnogo)
        pvals[k]={"pval":p.two_tail, "ngis":selgo, "ngo":len(goterms[k])}
        
    #2) Multiple hypotheses correction
    if correction=="bonferroni":
        th=th/len(goterms.keys())
    if correction=="fwer" or correction=="fdr":#sort first
        pvalso=[]
        for key, value in sorted(ego.iteritems(), key=lambda (k,v): (v["pval"],k)):
            pvalso.append(value["pval"])
        if correction=="fwer":
            for k in range(1,len(pvalso)):
                if pvalso[k] > th/(len(pvalso)-k+1):
                    th=pvalso[k-1]
                    break
        if correction=="fdr":
            for k in range(0,len(pvalso)):
                if pvalso[k] <= th*(k+1)/len(pvalso):
                    th=pvalso[k-1]
                    break
        
    print "th is {}".format(th)
    # and filter out terms
    pvalsf={}
    for p in pvals.keys():
        pv=pvals[p]
        if(pv["pval"]<th):
            pvalsf[p]=pv

    return pvalsf
