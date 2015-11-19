# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 11:39:03 2015

@author: rodri
"""
#%%
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
    alphabet=alphabetTotal[:numBins]   
    dseq=[]
    maximo=max(seq)
    minimo=min(seq)
    for i in xrange(0, len(seq)-windowSize+1,windowSize):
        im=np.mean(seq[i:i+windowSize])
        dseq.append(alphabet[(int)(np.round((numBins-1)*(im-minimo)/(maximo-minimo)))])
    print len(dseq)
    return dseq
      

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
    for i in xrange(len(seq)):
        s=seq[i]
        if s[0]=='t' and i>0:#new chromosome
         chsize.append(cont-1)
         print 'Chromosome with size {}'.format(cont-1)
         cont=0
        else:
            cont=cont+1
    chsize.append(cont-1) 
    print 'Chromosome with size {}'.format(cont-1)
    print '{} s in computing sizes'.format(time.clock()-t0) #39 seqs, go to numpy.array
    t0=time.clock()
    ch=[]
    cont=0
    for i in chsize:
        print i
        cont=cont+2
        chi=numpy.empty(i,dtype=float)
        for j in xrange(0,i-1):
            chi[j]=round(float(seq[cont+j]),2)
        cont=cont+i
        ch.append(chi)
    print '{} s in formatting'.format(time.clock()-t0) #39 seqs, go to numpy.array
    return ch

#%%
def preprocess(path="/Users/rodri/WebstormProjects/untitled/py/genomes/dwtMini2.wig", windowSize=100, numBins=5, maxSize=100000, track=0, onlyPositive=True):
    import numpy as np

    #0) read
    print 'reading...'
    #seq=readWig(str(request.args.get("path")))
    seq=readWig(path)
    seq=seq[track]
    print 'size of track {} is {}'.format(track, len(seq))
    
    #1) normalize 
    print 'normalizing...'
    m=np.mean(seq)
    sd=np.std(seq)
    nseq=[(x-m)/sd for x in seq]
    print 'mean: {}\tsd: {}'.format(m, sd)

    #2) discretize
    print 'discretizing...'
    #windowSize=int(request.args.get("windowSize"))
    #numBins=int(request.args.get("numBins"))
    #maxSize=int(request.args.get("maxSize"))
    if onlyPositive:
        seq=[min(max(0,x),m+2*sd) for x in seq]
        
    dseq=discretize(seq, windowSize, numBins)
    print 'done!'
    #TODO: broken pipe returning the object if whole 
            #Although it is not always happening, it's clearly not sensible
            #to send the full sequence. We sample instead 100K values 
    #return jsonify(result=[nseq[x] for x in range(0,len(seq),max(1,len(seq)/maxSize))], fullLength=len(nseq), maximum=np.max(nseq), minimum=np.min(nseq), mean=m, sdev=sd, dseq=dseq)
    #return {"result":[seq[x] for x in range(0,len(seq),max(1,len(seq)/maxSize))], "fullLength":len(nseq), "maximum":np.max(seq), "minimum":np.min(seq), "mean":m, "sdev":sd, "dseq":dseq}
    step=max(1,len(seq)/maxSize)
    return {"result":[np.mean(seq[x:x+step]) for x in range(0,len(seq),step)], "step":step, "fullLength":len(nseq), "maximum":np.max(seq), "minimum":np.min(seq), "mean":m, "sdev":sd, "dseq":dseq}
#%%
def search(pattern="", d=0):
    global data
    
    print pattern
    pattern=convertString(pattern)
    print pattern
    t=data["bwt"]
    print t["firstOccurrence"]
    if(False in [x in t["firstOccurrence"] for x in set(pattern)]):
        return "There are characters in pattern that do not correspond to the sequence characters: {}".format(t["firstOccurrence"].keys())
    else:
        t0=time.clock()
        match=ss.bwMatchingV8("".join(data["dseq"]), pattern, t["bwt"], t["firstOccurrence"],t["suffixArray"],t["checkpoints"],1000, d)
        print "Search takes {}".format((time.clock()-t0))
        return (str)(match)#return "hola {} {}".format(pattern, d)
