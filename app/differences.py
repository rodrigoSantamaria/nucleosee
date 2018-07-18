# -*- coding: utf-8 -*-
"""
Agnostic difference searches
@author: rodri
"""
#%%
def loadData(dataName="Test", track="None", clear="false"):
    import time
    import os
    import annotations as ann
    import pickle #i tried cPickle but is way slower!
    
    t00=time.clock()
    
    
    #0) Get info from the track file
    #basePath=os.path.join(app.config['UPLOAD_FOLDER'],user)
    basePath="./genomes"
    f=open(os.path.join(basePath,"tracks.txt"))
    for l in f.readlines():
        if(l.split("\t")[0].strip()==dataName):
            break

    chars=l.split("\t") 
    picklePaths=chars[1].strip().split(",")
    if(len(picklePaths)>1):
        picklePath=chars[0].strip()+".pic"
    else:
        picklePath=picklePaths[0]
   
    organism=chars[4].strip()
    
    #----Load the corresponding pickle
    data={}
    
    t0=time.clock()
    print("---------------- FILE: "+picklePath+" ------------------")
    f=open(os.path.join(basePath,picklePath))
    pdata=pickle.load(f)
    print("KEYS", pdata.keys())
    
    filenames=filter(lambda x:x.endswith(".wig") or x.endswith(".bw"),pdata.keys())
    if(len(filenames)>0):
        filename=filenames[0]
    else: #batch data
        filename="processed"
 
    if track=="None":
        track=sorted(pdata[filename]["maximum"].keys())[0]
    if(len(pdata[filename]["gff"].keys())!=len(pdata[filename]["seq"].keys())):
            print("One or more chromosome tracks do not match with GFF names:")      
    
   
    print("load data takes",(time.clock()-t0))
    
    #--------------- COMPILE BATCH
    import analysis
    data=analysis.batchCompile(pdata, filenames, filename,organism)
    #session[user][dataName]=data#trying to load several datasets at once
    print('file preprocess takes ',(time.clock()-t00),"s")
    #print("USER KEYS:", session[user].keys())
    
    dbp=data["batch"]["processed"]
    
    return dbp["dseq"][track]
#    return jsonify(seq=dbp["res"][track], 
#                fullLength=dbp["fullLength"],
#                maximum=(float)(dbp["maximum"][track]),
#                minimum=(float)(dbp["minimum"][track]), 
#                mean=(float)(dbp["mean"][track]), 
#                stdev=(float)(dbp["stdev"][track]), 
#                dseq=dbp["dseq"][track], 
#                chromosomes=sorted(dbp["maximum"].keys()),
#                bins=dbp["bins"][track],
#                windowSize=data["windowSize"])
#%%
def agnosticSearch(l=300, v=200, ws=30, seq1=[], seq2=[]):
    import time
    import numpy as np

    ll=l/ws
    vv=v/ws
    
    n=min(len(seq1), len(seq2))#TODO: workaround to avoid errors, should raise an exception if sizes differ
    t0=time.clock()
    difs=np.zeros(n)
    contv=0
    contl=0

    for i in xrange(n):
        if(seq1[i]!=seq2[i]):
            difs=1
    res=[]
    for i in xrange(n-ll):
        if(sum(difs[i,i+ll])>vv):
            res.append(i*ws)
    print(time.clock()-t0)
    return res

#%%
wt=loadData("h-972", track="chromosome3")
hta=loadData("DHTA1", track="chromosome3")

#%%
res=agnosticSearch(2000,900,30,wt, hta)
#%%
#Bruteforce differences -- good enough!
#def agnosticSearch(window=300, difwindow=200):
#    ws=30
#    l=300
#    v=250
#    ll=l/ws
#    vv=v/ws
#    
#    import time
#    t0=time.clock()
#    difs=[]
#    contv=0
#    contl=0
#    for i in range(len(hta)):
#        contl+=1
#        if(hta[i]!=wt[i]):
#            contv+=1
#            if(contv>vv):
#                #print(contv,contl)
#                difs.append([(i-ll*ws), "".join(wt[i-ll:i]), "".join(hta[i-ll:i])])
#        if(contl>ll):
#            contv=0
#            contl=0
#    print(time.clock()-t0)
#    len(difs)
        
