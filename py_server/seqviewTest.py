# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 18:36:20 2017

@author: rodri
"""
import os
os.chdir("/Users/rodri/WebstormProjects/seqview/py_server")
import suffixSearch as ss
import analysis as an
import pickle


#%%

#pdata=pickle.load(open("/Users/rodri/WebstormProjects/seqview/py_server/genomes/h-972_Rep2_depth_wl_trimmed_PE-chonlyc3.0ws30nb3imeanorgSchizosaccharomyces_pombe.pic"))
pdata=pickle.load(open("/Users/rodri/WebstormProjects/seqview/py_server/genomes/972h.pic"))
#%%
filenames=filter(lambda x:x.endswith(".wig") or x.endswith(".bw"),pdata.keys())
#data=an.batchCompile(pdata, filenames, 'h-972_Rep1_depth_wl_trimmed_PE-chonly.wig',"Schizosaccharomyces_pombe")
data=an.batchCompile(pdata, filenames, 'processed',"Schizosaccharomyces_pombe")
dbp=data["batch"]["processed"]

k="chromosome1"
pattern="abcba"
#%%
t=dbp["bwt"][k]
res=ss.bwMatchingV7(pattern, t["bwt"], t["firstOccurrence"], t["suffixArray"], t["checkpoints"], k=100000)

#%%
t=dbp["bwt"][k]
text="".join(data["batch"]["processed"]["dseq"][k])
res=ss.bwMatchingV8(text, pattern, t["bwt"], t["firstOccurrence"], t["suffixArray"], t["checkpoints"], 1000, 0)
