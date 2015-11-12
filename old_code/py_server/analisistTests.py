# -*- coding: utf-8 -*-
"""
Created on Mon Jan 26 11:02:37 2015

@author: rodri
"""

"""
Dada una secuencia de números reales (seq) calcula su media m y desviación estándar sd.
Para cada intervalo (tamaño windowSize) calcula la media del intervalo im y le
asigna un valor alfanumérico (de entre numBins valores) dependiendo del valor
de im respecto a m/sd
Por ejemplo, si seq=[0,0,34,22] y windowSize=2 y numBins=5 --> [a,d]
"""
def discretize(seq=[0,0,0], windowSize=2,numBins=5):
    import numpy as np
    alphabetTotal=['a','b','c','d','e', 'f', 'g','h','i','j','k','l','m','n','o','p','q','r','s','t']
    alphabet=alphabetTotal[0:numBins]   
    #windowSize=int(request.args.get("windowSize"))
    #numBins=int(request.args.get("numBins"))
    #seq=""+request.args.get("seq")
    #seq= [float(x) for x in seq.split(",")]
    dseq=[]
    sm=np.mean(seq)
    ssd=np.std(seq)
    for i in range(0, len(seq),windowSize):
        im=np.mean(seq[i:i+windowSize])
        dseq.append(alphabet[max(0,min(len(alphabet)/2+int(np.round((im-sm)/ssd)),numBins-1))])
    print len(dseq)
    #return jsonify(result=dseq)
    return dseq
      
#%%
seq=[1,7,7,7,0,2,3,4,5,6,0,0,0]
discretize(seq, 2)
#%%
import time
import string
t0=time.clock()
#f=open("/Users/rodri/Documents/investigacion/IBFG/nucleosomas/dwtMini2.wig")
f=open("/Users/rodri/Documents/investigacion/IBFG/nucleosomas/Mei3h_raw_ch1.wig")
seq=f.readlines()
del seq[0:2]
seq=[float(string.replace(x, "\n", "")) for x in seq]
print "time for reading: {}".format(time.clock()-t0)#9s for 1ch of pombe (17M)
#%%
t0=time.clock()
seqn=normalize(seq)
print "time for normalization: {}".format(time.clock()-t0)#15s for 1ch of pombe

import string
seq=[float(string.replace(x, "\n", "")) for x in seq]

ds=discretize(seq,100)

#%%
#%% READ WIG
def readWig(path="/Users/rodri/Documents/investigacion/IBFG/nucleosomas/Mei3h_center.wig"):
    import numpy
    import time
    t0=time.clock()
    f=open(path)
    seq=f.readlines()
    print '{} s in reading'.format(time.clock()-t0) #5 secs
    t0=time.clock()
    chsize=[]
    cont=0
    for i in range(len(seq)):
        s=seq[i]
        if s[0]=='t' and i>0:#new chromosome
         chsize.append(cont-1)
         cont=0
        else:
            cont=cont+1
    chsize.append(cont-1)
    print '{} s in computing sizes'.format(time.clock()-t0) #39 seqs, go to numpy.array
    t0=time.clock()
    ch=[]
    cont=0
    for i in chsize:
        print i
        cont=cont+2
        chi=numpy.empty(i,dtype=float)
        for j in range(0,i-1):
            chi[j]=round(float(seq[cont+j]),2)
        cont=cont+i
        ch.append(chi)
    print '{} s in formatting'.format(time.clock()-t0) #39 seqs, go to numpy.array
    return ch

#%%
def preprocess(path="/Users/rodri/WebstormProjects/untitled/py/genomes/dwtMini2.wig", windowSize=100, numBins=5, maxSize=100000):
    import numpy as np

    #0) read
    print 'reading...'
    #seq=readWig(str(request.args.get("path")))
    seq=readWig(path)
    seq=seq[0] #(TODO: by now, only first chromosome)
    #1) normalize 
    print 'normalizing...'
    m=np.mean(seq)
    sd=np.std(seq)
    nseq=[(x-m)/sd for x in seq]
    #2) discretize
    print 'discretizing...'
    #windowSize=int(request.args.get("windowSize"))
    #numBins=int(request.args.get("numBins"))
    #maxSize=int(request.args.get("maxSize"))
    dseq=discretize(seq, windowSize, numBins)
    print 'done!'
    #TODO: broken pipe returning the object if whole 
            #Although it is not always happening, it's clearly not sensible
            #to send the full sequence. We sample instead 100K values 
    #return jsonify(result=[nseq[x] for x in range(0,len(seq),max(1,len(seq)/maxSize))], fullLength=len(nseq), maximum=np.max(nseq), minimum=np.min(nseq), mean=m, sdev=sd, dseq=dseq)
    return {"result":[seq[x] for x in range(0,len(seq),max(1,len(seq)/maxSize))], "fullLength":len(nseq), "maximum":np.max(seq), "minimum":np.min(seq), "mean":m, "sdev":sd, "dseq":dseq}
    
#%%
res=preprocess()
res["result"][26271:26421]
res["mean"]