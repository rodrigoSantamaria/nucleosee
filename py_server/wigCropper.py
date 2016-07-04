# -*- coding: utf-8 -*-
"""
This script crops a wig file so max values cannot exceed a given value based
on standard deviations
@author: rodri
"""
import helpers
import numpy as np
#%%
path="/Users/rodri/WebstormProjects/seqview/py_server/genomes/jpiriz/mono-H3K9_me2_norm_center_wlt.wig"
seq=helpers.readWig(path)

#%%
f=open("/Users/rodri/WebstormProjects/seqview/py_server/genomes/jpiriz/mono-H3K9_me2_norm_center_wlt_cropped2.wig", "w")
for k in seq.keys():
    sk=seq[k]
    tal=np.clip(sk, 0, 2)
    f.write("track type=wiggle_0 name="+k+" description=\"H3K9_me2_MAP cropped values >2\"\n")
    f.write("fixedStep chrom="+k+" start=1 step=1\n")
    for n in tal:
        f.write((str)(n)+"\n")
f.close()
