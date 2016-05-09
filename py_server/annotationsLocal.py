# -*- coding: utf-8 -*-
"""
Searching for coincidences with Quique's perfect nucleosomes
@author: rodri
"""
import sys
sys.path.append("/Users/rodri/WebstormProjects/seqview/py_server")
import analysisLocal as al
import annotations as ann
import suffixSearch as ss
import motifSearch as ms

import numpy as np
import time
import os
os.chdir("/Users/rodri/WebstormProjects/seqview/py_server")
#------------------------------------
#%%40M data test (3 whole pombe chromosomes)
#------------------------------------
""" read wig 2s
    compute sizes 2s
    formatting 20s
    TOTAL      47s
    So, normalize + discretize is about 20s, not bad. The formatting part should be checked
    Search in one chromosome is about 1s, not bad
"""
#f="/Users/rodri/WebstormProjects/seqview/py_server/genomes/rodri/23479_h90_wlt_mean.wig"
f="/Users/rodri/WebstormProjects/seqview/py_server/genomes/rodri/Mei3h_center_wl-peque2.wig"
ws=30
res=al.preprocess(filename=f, windowSize=ws)

res["mean"]
len(res["seq"])
max(res["seq"])
ds=res["dseq"]
text=''.join(ds)+"$"
t=res["bwt"]

#%%
dataGFF=ann.gff("genomes/annotations/spombe/gff/schizosaccharomyces_pombe.III.gff3")
#%%
em=ann.annotate([61000], dataGFF, ["any"], ws=1500)

em
#%%
t0=time.clock()
match=ss.bwMatchingV8(text, "aba", t["bwt"], t["firstOccurrence"],t["suffixArray"],t["checkpoints"], d=0)
len(match)
print "search {}s".format(time.clock()-t0)#1.4s


#%%     ABCBA(2)
t0=time.clock()
match=ss.bwMatchingV8(text, "abcba", t["bwt"], t["firstOccurrence"],t["suffixArray"],t["checkpoints"], d=0)
len(match)
print "search {}s".format(time.clock()-t0)#1.4s

#%%

#%% Large perfect nucleosomes of moderate height
match=ss.bwMatchingV8(text, "abcba"*2, t["bwt"], t["firstOccurrence"],t["suffixArray"],t["checkpoints"], d=0)
len(match)
#match
#%%
match=ss.bwMatchingV8(text, "abcba"*3, t["bwt"], t["firstOccurrence"],t["suffixArray"],t["checkpoints"], d=0)
len(match)




#%% 3 consecutive perfect nucleosomes
match=ss.bwMatchingV8(text, "abcba"*3, t["bwt"], t["firstOccurrence"],t["suffixArray"],t["checkpoints"], d=2)
len(match)
#%% 6 consecutive almost perfect nucleosomes
match=ss.bwMatchingV8(text, "abcba"*6, t["bwt"], t["firstOccurrence"],t["suffixArray"],t["checkpoints"], d=3)
len(match)
#%% 8 consecutive almost perfect nucleosomes
match=ss.bwMatchingV8(text, "abcba"*8, t["bwt"], t["firstOccurrence"],t["suffixArray"],t["checkpoints"], d=5)
len(match)
#%% Nucleosome depleted large zones (.9K nucleotides)
match=ss.bwMatchingV8(text, "a"*30, t["bwt"], t["firstOccurrence"],t["suffixArray"],t["checkpoints"], d=0)
match
#%% Areas of indefinition of up to 600 nucleotides
match=ss.bwMatchingV8(text, "b"*20, t["bwt"], t["firstOccurrence"],t["suffixArray"],t["checkpoints"], d=0)
len(match)

#%% Areas with an NDR plus 3 nucleosomes
pattern=''.join(["aaaaa","abcba"*3])
match=ss.bwMatchingV8(text, pattern, t["bwt"], t["firstOccurrence"],t["suffixArray"],t["checkpoints"], d=1)
match







#%%
#remove close overlaps
t0=time.clock()
mm=[x*ws+74 for x in match]
mm2=[]
for i in range(len(mm)-1):
    if(mm[i+1]-mm[i]>=(3)*ws):
        mm2.append(mm[i])
mm2.append(match[-1])
mm=mm2

len(mm)
print "overlap removal {}s".format(time.clock()-t0)#1.4s

#%% -------------------------- ADDITIONAL ANNOTATIONS ----------------------

def gff0(filename="/Users/rodri/Documents/investigacion/IBFG/visualization/spombe/gff/schizosaccharomyces_pombe.I.gff3"):
    #filename="/Users/rodri/Documents/investigacion/IBFG/visualization/spombe/gff/schizosaccharomyces_pombe.I.gff3"
    f=open(filename)
    import csv
    reader=csv.DictReader(f, delimiter="\t")
    next(reader)
    reader.fieldnames=["chromosome", "source", "type", "start", "end", "xx", "sense", "xx", "id"]
    data=[]
    for row in reader:
        row["start"]=(int)(row["start"])
        row["end"]=(int)(row["end"])
        data.append(row)
    return data
    

#%%
def gff(filename="/Users/rodri/Documents/investigacion/IBFG/visualization/spombe/gff/schizosaccharomyces_pombe.I.gff3"):
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
        tow=data[i]
        tow["start"]=(int)(row["start"])
        tow["end"]=(int)(row["end"])
        tow["chromosome"]=row["chromosome"]
        tow["type"]=row["type"]
        tow["sense"]=row["sense"]
        tow["id"]=row["id"]
    return data
          

#%%  Loads the gene ontology into a dic where keys are go_ids and values are just go names by now
def go(filename="/Users/rodri/Desktop/vorontoParsing/VorontoGOremap/gene_ontology_ext.obo"):
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
def goa(filename="/Users/rodri/Documents/investigacion/IBFG/visualization/spombe/goa/gene_association.pombase"):
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
    f=open("/Users/rodri/Documents/investigacion/IBFG/FASTA/chromosome"+(str)(ch)+".fasta")
    reader=f.readlines()
    reader=reader[1:]
    seq=""
    for i in xrange(len(reader)):
        seq=seq+reader[i].replace("\n", "")
    return seq
#%%
import time    
t0=time.clock()    
dataGFF=gff()
dataGO=go()
dataGOA=goa()
dataFASTA=fasta(1)
print "annotation read {}s".format(time.clock()-t0)#1.4s
        
        
#%% Annotates each position in mm with the corresponding annotations found in dataGFF
# Only annotations of the types specified will be returned, or all if types=["any"]
# Annotations returned ar only the id and the type of annotation
# NOTE: by now it is responsability of the caller to control that the positions in mm and the annotations in dataGFF correspond to the same chromosome
#NOTE: with >1k matches and any type takes more than 1s
def annotate(mm, dataGFF, types=["any"], ws=1000):
    import numpy as np
    data2=dataGFF
    if(types[0]!="any"):
        wanted_set = set(types)  # Much faster look up than with lists, for larger lists
        @np.vectorize
        def selected(elmt): return elmt in wanted_set  # Or: selected = numpy.vectorize(wanted_set.__contains__)
        data2=dataGFF[selected(dataGFF["type"])]

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
t0=time.clock()
em=annotate(mm, dataGFF, ["gene", "CDS", "exon"], ws=500)
print "annotat takes {}s".format(time.clock()-t0)#2.5s for 1500 matches (any type)

#%% Check GO terms related to the genes in the em structure respect to dataGOA
    #Terms in discard are not considered (by default CC, BP y MF roots)
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
            egon[k]={"name": dataGO[k], "genes":ego[k]}
        except KeyError:
            print 'Key {} not found'.format(k)
    return egon

agon=annotateGOnames(em, dataGOA, dataGO)
#%%

agoH=filter(lambda x:len(ego[x])>5, ego.keys())
filter(lambda x:len(egon[x])>5, egon.keys())

#%%
#According to https://pypi.python.org/pypi/fisher/0.1.4
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

# Genes of interest given some positions and after annotation (em)
gis=set()
for x in em.keys():
    for y in em[x]:
        gis.add(y["id"])
    
ego=enrichmentFisher(list(gis), dataGOA, 0.01, "fdr")
for k in ego.keys():
    ego[k]["go_name"]=dataGO[k]
print len(ego),'enriched terms'

#%%
for key in sorted(ego.iterkeys()):
    print "%s: %s" % (key, ego[key])
#%% FWER
pvalso=[]
for key, value in sorted(ego.iteritems(), key=lambda (k,v): (v["pval"],k)):
    pvalso.append(value["pval"])
pvalso
for k in range(1,len(pvalso)):
    if pvalso[k] > 0.01/(len(pvalso)-k+1):
        print pvalso[k-1]
        break
pvalso
#%% --------------------- MOTIF -------------------------
#Given a list of positions mm and a sequence length slength, the method
#uses a Gibbs sampler to check which are the best krange-motifs for the correponding
#sequences in dataFASTA
#NOTE: the method is slow (1min for 30matches) and requires being slower for more precision (more iterations)
#NOTE: the method might bias towards lower kmers
def motif(mm, slength, dataFASTA, krange=range(8,9)):
    import time
    t0=time.clock()
    mots=[]
    for m in mm:
        mots.append(dataFASTA[m:m+slength].upper())
    print "Time in getting seqs {}".format(time.clock()-t0)
    t0=time.clock()
    bs=-1
    for k in krange:
        print '----------{}------------'.format(k)
        res=ms.gibbsSampler(mots,k,100)
        sc=ms.score(res)
        if(bs==-1 or sc<bs):
            bs=sc
            bm=res            
        for i in range(200):
            res=ms.gibbsSampler(mots,k,100)
            sc=ms.score(res)
            if(sc<bs):
                bs=sc
                bm=res
                print "{}) mejora: {}\t{}".format(i,ms.consensus(bm), bs)
    print "El algoritmo tarda {}".format(time.clock()-t0)
    return {'motifs':bm, 'score':bs, 'consensus':ms.consensus(bm)}
    
motif(mm, len(pattern)*ws, dataFASTA, range(8,9)) #about 1 min per k range
