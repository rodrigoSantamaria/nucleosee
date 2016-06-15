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
def discretize(seq, windowSize, minimo, maximo, numBins=5, percentile=True):
    import numpy as np
    #import sys
    #sys.path.append("/Users/rodri/WebstormProjects/seqview/py_server")
    import helpers
    
    alphabetTotal=['a','b','c','d','e', 'f', 'g','h','i','j','k','l','m','n','o','p','q','r','s','t']
    alphabet=alphabetTotal[:numBins]   
    dseq=[]
    
    sseq=helpers.rolling_window(seq,windowSize)
    #sseq=np.split(np.array(seq[:windowSize*(len(seq)/windowSize)]), len(seq)/windowSize)
    mseq=np.mean(sseq, axis=1, keepdims=True)
    
    if(percentile==True):
        mseq=np.array(mseq);
        pers=[0]
        for i in range(1,numBins+1):
            p=np.percentile(mseq,(100.0/numBins)*i)
            print("percentile",(100.0/numBins)*i,"=",p)
            pers.append(p)
        digseq=np.digitize(mseq,pers)
        for s in digseq:
           # print(s)
            dseq.append(alphabet[min(s-1,len(alphabet)-1)])           
    else:
        factor=(numBins-1.0)/float(maximo-minimo)
        for im in mseq:
            dseq.append(alphabet[(int)(factor*(im-minimo))])
    return dseq 

#%% READ WIG
def readWig(path="/Users/rodri/Documents/investigacion/IBFG/nucleosomas/Mei3h_center.wig"):
    import numpy
    import time
    import re
    
    t0=time.clock()
    f=open(path)
    seq=f.readlines()
    print ((time.clock()-t0),' s in reading') 
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
    print((time.clock()-t0),' s in computing sizes')
    t0=time.clock()
    ch={}
    cont=0
    name=re.sub("\n", "", re.sub(" .*$", "", re.sub("^.*chrom=", "", seq[cont+1])))
    ch[name]=[]
    print("name is ", name)
    for i in chsize:
        print(i)
        cont=cont+2
        #chi=numpy.array(seq[cont:cont+i-1],float)
        chi=numpy.array(seq[cont:cont+i-1], dtype=numpy.float16) #gives issues with jsonify but is much more memory efficient
        ch[name].append(chi)
        cont=cont+i
        if(cont<len(seq)):
            name=re.sub("\n", "",re.sub(" .*$", "", re.sub("^.*chrom=", "", seq[cont])))
            print("name is ", name)
            ch[name]=[]
    for k in ch.keys():
        ch[k]=ch[k][0]
    print ((time.clock()-t0),' s in formatting')
    
    return ch
 #%%
def readWig0(path="/Users/rodri/Documents/investigacion/IBFG/nucleosomas/Mei3h_center.wig"):
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

#path="/Users/rodri/Documents/investigacion/IBFG/nucleosomas/quique/23479_h90_wlt_mean.wig"
#tal=readWig(path)
#tal=preprocess(path, track="chromosome3")
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
def preprocess(filename="dwtMini2.wig", windowSize=100, numBins=5, maxSize=100000, stdev=3, track="None", recharge="False"):
    import numpy as np
    import time

    import sys
    sys.path.append("/Users/rodri/WebstormProjects/seqview/py_server")
    import suffixSearch as ss
    import helpers
    
    t00=time.clock()
    #0) read
    print('reading...')
    t0=time.clock()
#    basePath=os.path.join(app.config['UPLOAD_FOLDER'],user)
#    filename=str(request.args.get("filename"))
#    path=os.path.join(basePath,filename)
#    track=request.args.get("track")
#    forceReload=request.args.get("recharge")
#    windowSize=int(request.args.get("windowSize"))
#    numBins=int(request.args.get("numBins"))
#    maxSize=int(request.args.get("maxSize"))
#    stdev=float(request.args.get("stdev"))
#    picklePath=os.path.join(basePath,re.sub(r"\..*$", ".pic", filename))
#    savePickle=False
    recharge=bool(recharge)
    path=filename
    
    m={}
    sd={}
    maximum={}
    minimum={}
    dseq={}#discretized (alphanumeric) sequence
    res={}#reduced sequence
    t={} #burrows-wheeler transform
    dataGFF={}
    dataFASTA={}
    seqd={}
    
#    if os.path.isfile(picklePath) and forceReload=="False":
#        print("Pickle exists!!!, recharge=", forceReload)
#        f=open(picklePath)
#        datap=pickle.load(f)
#        genome=datap["seq"]
#        t=datap["bwt"]
#        dseq=datap["dseq"]
#        maximum=datap["maximum"]
#        minimum=datap["minimum"]
#        m=datap["mean"]
#        sd=datap["stdev"]
#    else:
    genome=readWig(path)
    savePickle=True

    if track=="None":
        track=sorted(genome.keys())[0]
    
    print("load wig takes",(time.clock()-t0))
   
    #for each separate chromosome: TODO: include ds and bwt into pickle!
    for k in genome.keys():
        tk=time.clock()
        seq=genome[k]
        print("ch", k, "with length", len(seq))
        if(savePickle):
            #1) normalize 
            t0=time.clock()
            
            m[k]=np.mean(seq, dtype=float)
            sd[k]=np.std(seq, dtype=float)
            upperlim=m[k]+stdev*sd[k]#avoid outliers? testing
            seq=np.clip(seq,0,upperlim)
        
            m[k]=np.mean(seq, dtype=float)
            sd[k]=np.std(seq, dtype=float)
            maximum[k]=np.max(seq)
            minimum[k]=np.min(seq)
            print('\\tstats in ',(time.clock()-t0), "s")
        
        
            #2) discretize
            t0=time.clock()
            dseq[k]=discretize(seq, windowSize, minimum[k], maximum[k], numBins)
            print('\\tdiscretize in',(time.clock()-t0),' s')
        
            t0=time.clock()
            t[k]=ss.bwt(''.join(dseq[k])+"$")
            print('\tbwt in ',(time.clock()-t0),'s')
    
        t0=time.clock()
        res[k]=list(np.mean(helpers.rolling_window(seq, max(1,len(seq)/maxSize)),-1, dtype=float)) #maybe round?  
        print('\tsampling in',(time.clock()-t0),'s')
    
        t0=time.clock()
        #dataGFF[k]=ann.gff(helpers.gffPath(ch=k))
        #dataFASTA[k]=ann.fasta(k)
        print('\ttime in annotations (GFF and FASTA):',(time.clock()-t0),'s')
        seqd[k]=seq
        print("processing",k,"takes",(time.clock()-tk))
        
    #3) annotations
    t0=time.clock()
    #dataGO=ann.go()
    #dataGOA=ann.goa()
    print("done! ... GO annotations takes",(time.clock()-t0))
    
    #data={"seq":seqd, "fullLength":len(seq), "maximum":maximum, "minimum":minimum,
    #      "mean":m, "stdev":sd, "dseq":dseq, "bwt":t, "gff":dataGFF, "res":res,
    #      "go":dataGO, "goa":dataGOA, "fasta":dataFASTA}
    #session[user]=data
    savePickle=False
    if(savePickle):
        tpickle=time.clock()
        print("pickle path:", picklePath)
        datap={"seq":seqd, "dseq":dseq, "bwt":t, "maximum":maximum, "minimum":minimum,
          "mean":m, "stdev":sd}
    
        f=open(picklePath, 'w')
        pickle.dump(datap,f)
        print("serialize data takes",(time.clock()-tpickle))
        #NOTE: should remove the .wig here to avoid double memory space?

    print('whole preprocess takes ',(time.clock()-t00),"s")
    print('returning track',track)
    print('max values are',maximum)
    
#    return {"seq":res[track], "fullLength":len(genome[track]), "maximum":(float)(maximum[track]), "minimum":(float)(minimum[track]), "mean":(float)(m[track]), "stdev":(float)(sd[track]), "dseq":dseq[track], "chromosomes":sorted(genome.keys())}
    return {"seq":res, "bwt":t, "dseq":dseq, "chromosomes":sorted(genome.keys())}
#%%
def preprocess0(filename="dwtMini2.wig", windowSize=100, numBins=5, maxSize=100000, track=0):
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
    dseq=discretize(seq, windowSize, minimum, maximum, numBins, True)
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








