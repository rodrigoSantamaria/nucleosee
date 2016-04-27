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
    for i in range(0, len(seq)-windowSize+1,windowSize):
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
    print '{} s in reading'.format(time.clock()-t0) #2 secs
    t0=time.clock()
    chsize=[]
    cont=0
    for s in seq[1:]:
        if s.startswith("track"):#new chromosome
         chsize.append(cont-1)
         print 'Chromosome with size {}'.format(cont-1)
         cont=0
        else:
            cont=cont+1
    chsize.append(cont-1) 
    print 'Chromosome with size {}'.format(cont-1)
    print '{} s in computing sizes'.format(time.clock()-t0) #5 seqs
    t0=time.clock()
    ch=[]
    cont=0
    for i in chsize:
        print i
        cont=cont+2
        chi=numpy.empty(i,dtype=float)
        for j in range(0,i-1):
            chi[j]=float(seq[cont+j])
        #chi.round(2)
        cont=cont+i
        ch.append(chi)
    print '{} s in formatting'.format((time.clock()-t0)) #6 seqs
    return ch
#%%

path="/Users/rodri/Documents/investigacion/IBFG/nucleosomas/quique/23479_h90_wlt_mean.wig"
#tal=readWig(path)
tal=preprocess(path, track="chromosome3")
#
##%%
#t0=time.clock()
#f=open(path)
#seq=f.readlines()
#print '{} s in reading'.format(time.clock()-t0) 
##%%
#import pandas
#import time
#t0=time.clock()
#seq=pandas.read_table(path, skiprows=2, names=["level"])
#print '{} s in reading with pandas'.format(time.clock()-t0) 
##searching for different tracks
#t0=time.clock()
#seq[seq.level.apply(lambda x: "track" in (str)(x))] 
#print '{} s in searching tracks with pandas'.format(time.clock()-t0) 
#%%
def preprocess(filename="dwtMini2.wig", windowSize=100, numBins=5, maxSize=100000, track=0):
#    global data
#    global session
    import numpy as np
    import suffixSearch as ss
    import annotations as ann
    import time
    t00=time.clock()
    
    #0) read
    print 'reading...'
    t0=time.clock()
#    path=os.path.join(app.config['UPLOAD_FOLDER'],user,str(request.args.get("filename")))
    path=filename
    #track=int(request.args.get("track"))
    seq=readWig(path)
    seq=seq[track] #(TODO: by now, only first chromosome)
    print 'done! in {}s'.format((time.clock()-t0))
    #1) normalize 
    print 'computing statistics...'
    t0=time.clock()
    m=np.mean(seq)
    sd=np.std(seq)
    maximum=np.max(seq)
    minimum=np.min(seq)
    #adding upper limit    
    upperlim=m+2*sd#avoid outliers? testing
    seq=[max(0,min(x,upperlim)) for x in seq]
    print 'max seq mod {}'.format(np.max(seq))
    #nseq=[(x-m)/sd for x in seq]
    print '{}\t{}'.format(m,sd)
    print 'done! in {}s'.format((time.clock()-t0))
    #2) discretize
    print 'discretizing...'
    t0=time.clock()
    #windowSize=int(request.args.get("windowSize"))
    #numBins=int(request.args.get("numBins"))
    #maxSize=int(request.args.get("maxSize"))
    dseq=discretize(seq, windowSize, numBins)
    print 'done! in {}s'.format((time.clock()-t0))
    print 'bwt...'
    t0=time.clock()
    t=ss.bwt(''.join(dseq)+"$")
    print 'done! in {}s'.format((time.clock()-t0))
    print 'squeezing...'
    t0=time.clock()
    step=max(1,len(seq)/maxSize)
    print 'step is {}'.format(step)
    res=[round(seq[x],2) for x in range(0,len(seq),step)]
    #res=[seq[i:(i+step)]/step for i in range(0,len(seq)-step,step)]
    print 'done! in {}s'.format((time.clock()-t0))
    print 'loading annotations...'
    t0=time.clock()
    dataGFF=ann.gff() #TODO: careful - select the right chr gff
    t0=time.clock()
    dataGO=ann.go()
    dataGOA=ann.goa()
    #dataFASTA=fasta(1)
    print 'done! in {}s'.format((time.clock()-t0))
    
    print '--- Preprocessing took {}s ---'.format((time.clock()-t00))

#    data={"seq":res, "fullLength":len(seq), "maximum":maximum, "minimum":minimum,
#          "mean":m, "stdev":sd, "dseq":dseq, "bwt":t,
#          "gff":dataGFF, "go":dataGO, "goa":dataGOA}
#    session[user]=data
    return {"seq": res, "fullLength":len(seq), "maximum":maximum, 
    "minimum":minimum, "mean":m, "stdev":sd, "dseq":dseq,
    "bwt":t, "gff":dataGFF, "go":dataGO, "goa":dataGOA}


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






