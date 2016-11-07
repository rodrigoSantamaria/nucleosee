# -*- coding: utf-8 -*-
"""
Different methods to read annotations
@author: rodri
"""

#%% -------------------------- ADDITIONAL ANNOTATIONS ----------------------

#%%
"""
Returns a np table ("single") or a dictionary with np tables (structure="multilple") 
with genome functional annotations.
NOTE: for larger GFF files (e.g. dmelanogaster, 600MB, this is not performing at all)
"""
def gff(filename="genomes/annotations/spombe/gff/schizosaccharomyces_pombe.III.gff3", structure="single"):
    import csv
    import time
    
    print("reading GFF csv...")
    f=open(filename)
    regions=["gene","exon","ncRNA_gene", "tRNA_gene", "snRNA_gene", "snoRNA_gene", "rRNA_gene", "three_prime_UTR", "five_prime_UTR"]
    fieldnames=["seqid", "source", "type", "start", "end", "score", "sense", "phase", "attributes"]
    reader=csv.DictReader(f, fieldnames=fieldnames, delimiter="\t")
    
    t0=time.clock()
    entries=[]
    
    print("filtering out comments...")
    for row in reader:
        if(row['seqid'].startswith("#") or not (row["type"] in regions)):   #case with comments between entries. NOTE: tam will be miscalculated in these cases
                continue
        entries.append(row)
        
    print("time in filtering out entries", (time.clock()-t0))
    t0=time.clock()
 
    import numpy as np
    import re
    print("populating annotations...")
    sc=("cerevisiae" in filename)
    
    data=np.empty(len(entries),dtype=[("chromosome", "a20"),("type", "a30"), ("start", "i8"), ("end", "i8"), ("sense", "a1"), ("id", "a40"), ("name", "a40")])
                
    for i in xrange(len(entries)):
        row=entries.pop()

        data[i]["start"]=(int)(row["start"])
        data[i]["end"]=(int)(row["end"])
        
        seqid=row["seqid"]
        try:
            rs=re.search("[XIV]+$",seqid).start() #roman numbers, happen in yeasts :s
            if(rs>=0):
                an=roman2arabic(seqid[rs:])
                if(an<10):
                    an="0"+(str)(an)
                else:
                    an=(str)(an)
                seqid=seqid[:rs]+an
        except:
            pass
        data[i]["chromosome"]=seqid
        
        data[i]["type"]=row["type"]
        data[i]["sense"]=row["sense"]
        if(sc==True):#SGD are 'flexible' about standars... grrr
            gid=re.sub("^.*ID=", "",re.sub(";.*$","",re.sub("^.*SGD:", "", row["attributes"])))
            name=re.sub("^.*ID=", "", re.sub(";.*$","",re.sub("^.*gene=", "", row["attributes"])))
        else:
            gid=re.sub("gene:", "", re.sub(";Name.*$","",re.sub("^.*ID=", "", row["attributes"])))
            name=re.sub(";.*$","",re.sub("^.*ID=.*;Name=", "", row["attributes"]))
            
        
        data[i]["id"]=gid
        data[i]["name"]=name            
    if(structure=="single"):
        return data
    else:
        d={}
        wanted_set = set(data["chromosome"])  # Much faster look up than with lists, for larger lists
        for k in wanted_set:
            @np.vectorize
            def selected(elmt): return elmt in k  # Or: selected = numpy.vectorize(wanted_set.__contains__)
            d[k]=data[selected(data["chromosome"])]

        return d
#import time
#t0=time.clock()
#dataGFF=gff("genomes/annotations/dmelanogaster/dmel-all-no-analysis-r6.12.gff", "multiple")    

##dataGFF=gff("genomes/annotations/scerevisiae/gff/saccharomyces_cerevisiae.gff", "multiple")
###dataGFF["chr01"][100]
#print("it took",(time.clock()-t0))

#%%
"""
Wraps annotations.gff to return a dictionary, with keys as provided, for the
corresponding organism.

This is highly heterogeneous: some organisms have all their chromosomes in a single
gff, some have one gff per chromosome, chromosome names does not mach with wig files, etc..

By now we are dealing it with this multiplexer function and by now ONLY for S pombe and S cerevisiae

Another option is to force gff files to a given format.

org - organism for the GFF (currently suppurted sce, spo, dme)
"""

def gffData(org="Schizosaccharomyces pombe", tracks=[]):    
    import re
    import os
    path=(org[0]+org[org.find(" ")+1:]).lower()    
    print path
    d={}
    if(org=="Schizosaccharomyces pombe"):
        for k in tracks:
            roman="I"
            if(k.find("3")>=0 or k.find("III")>=0):
                roman+="II"
            elif(k.find("2")>=0 or k.find("II")>=0):
                roman+="I"
            ret="genomes/annotations/spombe/gff/"+"schizosaccharomyces_pombe."+roman+".gff3"
            d[k]=gff(ret, "single")
    if(org=="Saccharomyces cerevisiae"):
        data=gff("genomes/annotations/"+path+"/gff/saccharomyces_cerevisiae.gff", "multilpe")
        print("gff read")
        for k in data.keys():
            k0=k
            k=k.replace("chromosome", "").replace("chr", "") 
            k=k.upper()
            if("I" in k or "X" in k or "V" in k):
                k=roman2arabic(k)
            k2="chr"
            if(k<10):
                k2+="0"
            k2+=(str)(k)
            d[k2]=data[k0]
    if(org=="Drosophila melanogaster"):
        contents=os.listdir("genomes/annotations/"+path)
        for c in contents:
            if(re.search(".gff$", c)):
                filename=c
        data=gff("genomes/annotations/"+path+"/"+filename, "multiple")
        d=data
        print("gff read")
        
    return d
    
#import time
#t0=time.clock()
#dataGFF2=gffData("Drosophila melanogaster", ["X","Y","2R"])
#print("GFF takes", (time.clock()-t0))#takes 1.5 min (4M entries filterest to 150K)

#%%
"""
Roman to arabic number conversion from 1 to 39 (no L, D, C, M contempled)
"""
def roman2arabic(roman):
    roman=roman.upper()
    arab=0
    for i in range(len(roman)):
        if(roman[i]=="X"):
            arab+=10
        if(roman[i]=="V"):
            arab+=5
        if(roman[i]=="I"):
            if(i<(len(roman)-1) and roman[i+1]!="I"):
                arab-=1
            else:
                arab+=1
    return arab
#roman2arabic("xxxviii")
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


#%% Loads the gene ontology annotation of a given species into an array which only stores gene_id, gene_name, go_id, go_type and gene_desc
def goa(org="Schizosaccharomyces pombe"):
    path=org[0]+org[org.find(" ")+1:]
    path=path.lower()
    if(path=="scerevisiae"):
        f=open("genomes/annotations/"+path+"/goa/gene_association.sgd")
    if(path=="spombe"):
        f=open("genomes/annotations/"+path+"/goa/gene_association.pombase")
    if(path=="dmelanogaster"):
        f=open("genomes/annotations/"+path+"/gene_association.fb")
    if(path=="mmusculus"):
        f=open("genomes/annotations/"+path+"/goa/gene_association.mgi")
        
    lines=f.readlines()
    data=[]
    for l in lines:
        if(l.startswith("!")==False):
            vals=l.split("\t")
            data.append({"gene_id":vals[1],"gene_name":vals[2], "go_id":vals[4], "evidence":vals[6], "go_type":vals[8], "gene_desc":vals[9]})
    return data
    

#%% Loads the fasta sequence of a given file (by now only working forS pombe files)
#it takes only 0.125s for the pombe genomw (the three chromosomes)
def fasta(ch, org="Schizosaccharomyces pombe"):
    path=org[0]+org[org.find(" ")+1:]
    path=path.lower()
    try:
        f=open("genomes/annotations/"+path+"/fasta/"+(str)(ch)+".fasta")
    except:  
        f=open("genomes/annotations/"+path+"/fasta/"+(str)(ch)+".fsa")
    reader=f.readlines()
    reader=reader[1:]
    return "".join(reader).replace("\n","")
  
#dataFASTA=fasta("chr01", "Saccharomyces cerevisiae")  
#import time
#t0=time.clock()
#for x in ["chromosome1","chromosome2","chromosome3"]:
#    tal=fasta(x)
#print(time.clock()-t0)


#%% Returns the annotations on areas around positions in mm
#dataGFF gff data where to search annotations (as returned by gff())
#types   specific lifst of annotations to retrieve (such as 'gene' or 'CDS') (default ["any"])
#ws      interval size from mm where to search for annotations (or a list of numbers for the same length than mm)
#align   wheter mm is taken as the start ("left" or the center of the range
#intersect  whether the range must be fully inside of the annotations ("hard") ur just have a non-null intersection ("soft")
#TODO: something broken here?
def annotate(mm, dataGFF, types=["any"], winSize=1000, align="center", intersect="soft"):
    import numpy as np
    data2=dataGFF
    if(types[0]!="any"):
        wanted_set = set(types)  # Much faster look up than with lists, for larger lists
        @np.vectorize
        def selected(elmt): return elmt in wanted_set  # Or: selected = numpy.vectorize(wanted_set.__contains__)
        data2=dataGFF[selected(dataGFF["type"])]
    em={} #enriched (i.e. detailed) matches, including for each the thigs found at GFF
    
    if(type(winSize)==int):
        ws=winSize
        interval=ws*0.5
    
    for i in range(len(mm)):
        x=mm[i]
        if(type(winSize)==list):
            ws=winSize[i]
            interval=ws*0.5
        
        if(align=="left"):
            s1=x
            e1=x+ws
        else:
            e1=x+interval
            s1=x-interval
        if(intersect=="hard"):
            sel=data2[(data2["start"]<=s1) & (data2["end"]>=e1)] #annotated interval fully inside the search window
        else:    
            sel=data2[((data2["end"]>=s1) & (data2["end"]<=e1)) | ((data2["start"]>=s1) & (data2["start"]<=e1)) | ((data2["start"]<=s1) & (data2["end"]>=e1))] #intersecting (more time expensive and not sure it makes a difference)

        if(len(sel)>0):     
            em[x]=[]
            for s in sel:
                em[x].append({"t":s["type"], "id":s["id"], "n":s["name"], "s":s["start"], "e":s["end"], "ss":s["sense"]})
                #em[x].append({"type":s["type"], "id":s["id"], "name":s["name"], "start":s["start"], "end":s["end"], "sense":s["sense"]})
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
def searchGene(text, dataGFF, types=["gene"]):
    import numpy as np
    import string
    data=dataGFF
    result=[]
    if(types[0]!="any"):
        wanted_set = set(types)  # Much faster look up than with lists, for larger lists
        @np.vectorize
        def selected(elmt): return elmt in wanted_set  # Or: selected = numpy.vectorize(wanted_set.__contains__)
        data=dataGFF[selected(dataGFF["type"])]
    for x in data:
      if(string.find(x["id"],text)>=0 or string.find(x["name"],text)>=0):
          result.append({"start":x["start"], "end":x["end"]})
    return result
#searchGene("raf1", dataGFF)    
#%%
#Searches a given text in go descriptions, then searches for genes annotated
#with the corresponding go terms and returns their locations as start-end
def searchGO(text, dataGOA, dataGO, dataGFF):
    import numpy as np
    import string
    result=[]
    for k in dataGO.keys():
      if(string.find(dataGO[k],text)>=0):
          result.append(k)
    result=set(result)
    res2=[]
    for x in dataGOA:
        if x["go_id"] in result:
            res2.append(x["gene_name"])
    res2=set(res2)
    
    wanted_set = set(["gene"])  # Much faster look up than with lists, for larger lists
    @np.vectorize
    def selected(elmt): return elmt in wanted_set  # Or: selected = numpy.vectorize(wanted_set.__contains__)
    data=dataGFF[selected(dataGFF["type"])]

    res3=[]
    for x in data:
      if x["name"] in res2:
         res3.append({"start":x["start"], "end":x["end"]})
    return res3

#result=searchGO("ribosomal", dataGOA, dataGO, dataGFF)    
#result
#
#import numpy as np
#wanted_set = set(result)  # Much faster look up than with lists, for larger lists
#@np.vectorize
#def selected(elmt): return elmt in wanted_set  # Or: selected = numpy.vectorize(wanted_set.__contains__)
#data=dataGOA[selected(dataGOA["go_id"])]

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
#%%
#goids=[x["go_id"] for x in dataGOA]
#goterms={}
#for x in goids:
#    goterms[x]=set()
#for x in dataGOA:
#    goterms[x["go_id"]].add(x["gene_id"])
    
#%% Fisher's enrichment
#According to https://pypi.python.org/pypi/fisher/0.1.4
# gis is a set of genes of interest by id
# th is the threshold
# minGO minimum number of genes in a GO term to be considered
# maxGO maximum numer of genes in a GO term to be considered
#discard TODO: indicate terms that we want to discard (typically inferred electronically annotations)
def enrichmentFisher(gis, dataGOA, th=0.01, correction="none", minGO=5, maxGO=500, discard=["IEA"]):
    #0) Prepare sets    
    # Retrieve a dict where k=go id and value=set of genes
    print("enrichment fisher")
    goids=[x["go_id"] for x in dataGOA]
    goterms={}
    for x in goids:
        goterms[x]=set()
    for x in dataGOA:
        if((x["evidence"] in discard) == False):
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

        if(len(goterms[k])>=minGO and len(goterms[k])<=maxGO and selgo>0):
            p = pvalue(selgo, selnogo, unigo, uninogo)
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
#tal=enrichmentFisher(set(gis),dataGOA,0.01,"fdr")