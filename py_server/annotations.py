# -*- coding: utf-8 -*-
"""
Different methods to read annotations
@author: rodri
"""

#%% -------------------------- ADDITIONAL ANNOTATIONS ----------------------

#%%
def gff(filename="genomes/annotations/spombe/gff/schizosaccharomyces_pombe.III.gff3"):
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
    import re


    data=np.empty(tam,dtype=[("chromosome", "a2"),("type", "a40"), ("start", "i8"), ("end", "i8"), ("sense", "a1"), ("id", "a50"), ("name", "a50")])
    for i in range(skip):
        next(reader)
    for i in range(tam):
        row=reader.next()
        data[i]["start"]=(int)(row["start"])
        data[i]["end"]=(int)(row["end"])
        data[i]["chromosome"]=row["chromosome"]
        data[i]["type"]=row["type"]
        data[i]["sense"]=row["sense"]
        gid=re.sub("gene:", "", re.sub(";Name.*","",re.sub(".*ID=", "", row["id"])))
        data[i]["id"]=gid
        name=re.sub(";.*","",re.sub(".*ID=.*;Name=", "", row["id"]))
        data[i]["name"]=name
    return data

#%%  Loads the gene ontology into a dic where keys are go_ids and values are just go names by now
def go(filename="genomes/annotations/go/go-basic.obo"):
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
#it takes only 0.125s for the pombe genomw (the three chromosomes)
def fasta(ch):
    f=open("genomes/annotations/spombe/fasta/"+(str)(ch)+".fasta")
    reader=f.readlines()
    reader=reader[1:]
    return "".join(reader).replace("\n","")
    
#import time
#t0=time.clock()
#for x in ["chromosome1","chromosome2","chromosome3"]:
#    tal=fasta(x)
#print(time.clock()-t0)


#%%
def annotate(mm, dataGFF, types=["any"], ws=1000, align="center"):
    import numpy as np
    data2=dataGFF
    if(types[0]!="any"):
        wanted_set = set(types)  # Much faster look up than with lists, for larger lists
        @np.vectorize
        def selected(elmt): return elmt in wanted_set  # Or: selected = numpy.vectorize(wanted_set.__contains__)
        data2=dataGFF[selected(dataGFF["type"])]
    em={} #enriched (i.e. detailed) matches, including for each the thigs found at GFF
    interval=ws*0.5
    
    for x in mm:
        if(align=="left"):
            s1=x
            e1=x+ws
        else:
            e1=x+interval
            s1=x-interval
        #sel=data2[(data2["start"]<s1) & (data2["end"]>e1)] #search window fully inside the annotated interval
        #sel=data2[(data2["start"]<s1) & (data2["end"]>e1)] #annotated interval fully inside the search window
        sel=data2[((data2["end"]>s1) & (data2["end"]<e1)) | ((data2["start"]>s1) & (data2["start"]<e1)) | ((data2["start"]<s1) & (data2["end"]>e1))] #intersecting (more time expensive and not sure it makes a difference)
        if(len(sel)>0):     
            em[x]=[]
            for s in sel:
                m[x].append({"type":s["type"], "id":s["id"], "name":s["name"], "start":s["start"], "end":s["end"], "sense":s["sense"]})
    return em    
#%%
#import time
#t0=time.clock()
#gis=np.random.randint(0,2.5e6,3284)
#tal=annotate(gis, dataGFF,types=["gene"], ws=4000)
#print("it tiook",(time.clock()-t0))
#%%
#for x in dataGFF:
#    if(x["id"].find("SPCC757.15")>=0):
#        print "{}\t{}, {}".format(x["type"],x["start"],x["end"])
#for x in tal[61000]:
#    print "{}\t{}, {}".format(x["id"],x["type"],x["start"])
#    
#%%
#em is just a list of gene ids
#dataGOA is a table with GOA data as retrieved by goa()
#discard is a list of GO terms you don't want to see (basically the level 0's)    
def annotateGO(em, dataGOA, discard=['GO:0003674','GO:0005575','GO:0008150']):
    ego={}
    for x in dataGOA:
        if((x["go_id"] in discard) ==False):
            for g in em:
                if(g==x["gene_id"]):
                    #print "{} -> {}({})".format(x["gene_id"], x["go_id"], x["go_type"])
                    if(x["go_id"] in ego.keys()):
                        ego[x["go_id"]].add(x["gene_id"])
                    else:
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
            print('Key',k,' not found')
    return egon

#%% Fisher's enrichment
#According to https://pypi.python.org/pypi/fisher/0.1.4
# gis is a set of genes of interest by id
# th is the threshold
# minGO minimum number of genes in a GO term to be considered
# maxGO maximum numer of genes in a GO term to be considered
def enrichmentFisher(gis, dataGOA, th=0.01, correction="none", minGO=2, maxGO=500):
    #0) Prepare sets    
    # Retrieve a dict where k=go id and value=set of genes
    print("enrichment fisher")
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
    print("universe created")
    print(sel,"gis on a universe of",uni)
    #
    #1) Fisher's test
    from fisher import pvalue #at least for mac os must be downloaded from here: https://pypi.python.org/pypi/fisher/0.1.4
    pvals={}
    for k in goterms.keys():
    
        gisInTerm=gis.intersection(goterms[k])
        selgo=len(gisInTerm) #gi and go
        selnogo=sel-selgo #gi not go
    
        unigo=len(goterms[k])-selgo #number of non-gis in the term
        uninogo=uni-len(gis)-len(goterms[k])+selgo
        if(unigo>=minGO and unigo<=maxGO and selgo>0):
            p = pvalue(selgo, selnogo, unigo, uninogo)
            #pvals[k]={"pval":p.right_tail, "ngis":selgo, "ngo":len(goterms[k])}
            pvals[k]={"pval":p.right_tail, "ngis":selgo, "ngo":len(goterms[k]), "gis":list(gisInTerm)}
    print("fisher test finished with ", len(pvals)," terms enriched")
        
    #2) Multiple hypotheses correction
    if correction=="bonferroni":
        th=th/len(goterms.keys())
    if correction=="fwer" or correction=="fdr":#sort first
        pvalso=[]
        def keyL(k):
            return k[1]["pval"],k[0]
        if correction=="fwer":
            for key, value in sorted(pvals.items(), key=keyL):
                pvalso.append(value["pval"])
            for k in range(1,len(pvalso)):
                if pvalso[k] > th/(len(pvalso)-k+1):
                    th=pvalso[k-1]
                    break
        if correction=="fdr":
            for key, value in sorted(pvals.items(), key=keyL, reverse=True):
                    pvalso.append(value["pval"])
            for k in range(0,len(pvalso)):
                #print(pvalso[k],"   ",th*(k+1)/len(pvalso))
                if pvalso[k] <= th*(k+1)/len(pvalso):
                    th=pvalso[k-1]
                    break
        print("multiple hypotheses correction finished with ",len(pvalso),"{} terms")
    print("th is ",th)
    # and filter out terms
    pvalsf={}
    for p in pvals.keys():
        pv=pvals[p]
        if(pv["pval"]<th):
            pvalsf[p]=pv
    print("number of enriched terms is",len(pvalsf))

    return pvalsf
#%%
#tal=enrichmentFisher(set(gis),dataGOA, 0.01, "fdr")
#fisher.pvalue(2,294,255,4584)
#tal=enrichmentFisher(set(gis),dataGOA,0.01,"fdr")
#%%
#def keyL(k):
#    return k[1]["pval"],k[0]
#sorted(pvals.items(), key=keyL)    
##%%
#pvalso=[]
#for key, value in sorted(pvals.items(), key=keyL):
#        pvalso.append(value["pval"])
#        
##%%
#path="/Users/rodri/WebstormProjects/seqview/py_server/genomes/annotations/spombe/goa/gene_association.pombase"
#enrichmentFisher(set(["SPBC3B8.06"]), goa(path), th=0.01, correction="fdr")
#
##%%
#dataGOA=dg
#gis=set(["SPBC3B8.06"])
#
##%%
#print "enrichment fisher"
#goids=[x["go_id"] for x in dataGOA]
#goterms={}
#for x in goids:
#    goterms[x]=set()
#for x in dataGOA:
#    goterms[x["go_id"]].add(x["gene_id"])
## Compute universe genes
#unigenes=set()
#for k in goterms.keys():
#        unigenes |= set(goterms[k])
#unigenes=unigenes-gis
#uni=len(unigenes)
#sel=len(gis)
#print "universe created"
#    #
