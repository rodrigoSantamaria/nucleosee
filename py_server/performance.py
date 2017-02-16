# -*- coding: utf-8 -*-
"""
Time performance tests

@author: rodri
"""
import os
import time
import string
import pickle #i tried cPickle but is way slower!

import sys

#Our methods
import suffixSearch as ss
import motifSearch as ms
import annotations as ann
import analysis as ana

import helpers
import re



#%%S pombe - preprocess
numBins=3
windowSize=30

#These are visualization or client/server parameters that we can leave to default in local console
maxSize=100000
stdev=3
track="None"
dataName="None"
interpolation="mean"
basePath="."

#These are parameters yet to implement, by now defaults
organism="Saccharomyces cerevisiae"

ana.batchPreprocessLocal(filenames, outfile, dataName, windowSize, numBins, maxSize, stdev, track, organism, interpolation, basePath)
print("Data successfuly preprocessed")
    
#%% S pombe searches LOAD
#ifile="/Users/rodri/WebstormProjects/seqview/py_server/genomes/MN-Damien-WT-1_S1_wlt_meanc3.0ws30nb3imeanorgSchizosaccharomyces_pombe.pic"
print("Opening",ifile)
f=open(os.path.join(ifile))
pdata=pickle.load(f)
fnames=[]
print(pdata.keys())
for k in pdata.keys():
    if(k.find(".")>-1):
        fnames.append(k)
if(len(fnames)==0):
    print("ERROR: No wig associated to pic data")
    exit
print fnames
data=ana.batchCompile(pdata, fnames, fnames[0],"")
print"DATA LOADED"

#%% ------------------- D EFFECT ----------------
ifile="/Users/rodri/WebstormProjects/seqview/py_server/genomes/MN-Damien-WT-1_S1_wlt_meanc3.0ws30nb3imeanorgSchizosaccharomyces_pombe.pic"
pattern=""
psize=5
geo="none"
d=0
print("Opening",ifile,"for pattern", pattern,"(",d,") restricted to",geo)
times=[]
for d in range(4):
    print("D is", d)
    tt=[]
    for i in range(20):
        pattern=""
        for j in range(psize):
            pattern=pattern+letters[random.randint(0,2)]
        print pattern

        t0=time.clock()
        ret=ana.searchLocal(data=data, pattern=pattern, d=d, geo=geo, intersect="soft", softMutations="false")
        tt.append(time.clock()-t0)
    times.append(np.mean(tt))

print"SEARCH FINISHED"
#%%
#t0=time.clock()
#ret=ana.searchLocal(data=data, pattern=pattern, d=d, geo=geo, intersect="soft", softMutations="false")
#tt.append(t0)

#%% xkcd plotting fun!
import matplotlib.pyplot as plt
import pylab
#plt.xkcd() #coolio

fig=plt.figure()
#fig.set_figwidth(12)
ax=fig.add_subplot(1,1,1)
ax.plot(range(4), times, "b-")
#pylab.ylim([-4,4])
ax.spines["top"].set_color("none")
ax.spines["right"].set_color("none")
plt.tick_params(top=False, right=False)
plt.xticks(np.arange(0,4,1))
plt.ylabel("time (s)")
plt.xlabel("allowed variations (v)")

fig.show()

#%% ------------------- PATTERN SIZE EFFECT ----------------

letters=["a","b","c"]
pattern=""
geo="none"
d=0
times=[]
for k in range(2,21):
    tt=[]
    for i in range(20):
        pattern=""
        for j in range(k):
            pattern=pattern+letters[random.randint(0,2)]
        print pattern
        t0=time.clock()
        ret=ana.searchLocal(data=data, pattern=pattern, d=0, geo=geo, intersect="soft", softMutations="false")
        tt.append(time.clock()-t0)
    times.append(np.mean(tt))

print"SEARCH FINISHED"

#%% xkcd plotting fun!
import matplotlib.pyplot as plt
import pylab
#plt.xkcd() #coolio

fig=plt.figure()
#fig.set_figwidth(12)
ax=fig.add_subplot(1,1,1)
ax.plot(range(2,len(times)+2), times, "b-")
#pylab.ylim([-4,4])
ax.spines["top"].set_color("none")
ax.spines["right"].set_color("none")
plt.tick_params(top=False, right=False)
#plt.xticks(np.arange(0,4,1))
plt.ylabel("time (s)")
plt.xlabel("pattern size")

fig.show()


#%% ------------------- GENOME SIZE EFFECT ----------------

genomes=[
         "/Users/rodri/WebstormProjects/seqview/py_server/genomes/GSM585199_FC14703_YE_MNase_Lane_1_eland_resultc3.0ws30nb3imeanorgSaccharomyces_cerevisiae.pic",
         "/Users/rodri/WebstormProjects/seqview/py_server/genomes/MN-Damien-WT-1_S1_wlt_meanc3.0ws30nb3imeanorgSchizosaccharomyces_pombe.pic",
         "/Users/rodri/WebstormProjects/seqview/py_server/genomes/whitek1_depth_wl_trimmed_PEc3.0ws30nb3imeanorgCandida_albicans_WO.pic",
         "/Users/rodri/WebstormProjects/seqview/py_server/genomes/WormGSM2417781_N2_Pr_neg_normc3.0ws30nb3imeanorgCaenorhabditis_elegans.pic",
         "/Users/rodri/WebstormProjects/seqview/py_server/genomes/GSE36212_Dmel_Combined_+c3.0ws30nb3imeanorgDrosophila_melanogaster.pic",
         "/Users/rodri/WebstormProjects/seqview/py_server/genomes/GSE81436_Cr_DJc3.0ws30nb3imeanorgOryza_sativa_jap.pic"]
         
names=["cerevisiae","pombe","albicans","worm", "fly", "rice"]
pattern="abcba"
letters=["a","b","c"]
geo="none"
d=0
times=[]
sizes=[]
numch=[]
for ifile in genomes:
    
    print("Opening",ifile)
    f=open(os.path.join(ifile))
    pdata=pickle.load(f)
    fnames=[]
    print(pdata.keys())
    for k in pdata.keys():
        if(k.find(".")>-1):
            fnames.append(k)
    if(len(fnames)==0):
        print("ERROR: No wig associated to pic data")
        exit
    print fnames
    data=ana.batchCompile(pdata, fnames, fnames[0],"")
    print"DATA LOADED"
    
    s=0
    c=0
    dbps=data["batch"]["processed"]["seq"]
    for k in dbps.keys():
        s+=len(dbps[k])
        c+=1
        
    print("GENOME SIZE:", s)
    sizes.append(s)
    numch.append(c)

    tt=[]
    for i in range(20):
        pattern=''.join([letters[random.randint(0,2)] for x in range(5)])
        t0=time.clock()
        ret=ana.searchLocal(data=data, pattern=pattern, d=0, geo="none", intersect="soft", softMutations="false")
        tt.append(time.clock()-t0)
    times.append(np.mean(tt))

print"SEARCH FINISHED"
#Candida and pombe are genomes which very precise values and annotations?

#%% xkcd plotting fun!
import matplotlib.pyplot as plt
import pylab
#plt.xkcd() #coolio

fig=plt.figure()
start=0
#fig.set_figwidth(12)
ax=fig.add_subplot(1,1,1)
ax.plot(sizes[start:], times[start:], "b-")
ax.plot(sizes[start:], times[start:], "bo")
for i in range(start, len(names)):
    ax.text(sizes[i]+10000000,times[i]-0.01,names[i])
#pylab.ylim([-4,4])
ax.spines["top"].set_color("none")
ax.spines["right"].set_color("none")
plt.tick_params(top=False, right=False)
#plt.xticks(np.arange(0,4,1))
plt.ylabel("time (s)")
plt.xlabel("genome size (hundreds of Mbps)")

fig.show()