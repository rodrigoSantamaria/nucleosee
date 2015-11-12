#----------------- BURROWS-WHEELER TRANSFORM
#Returns a dict with the transform (key 'bwt') and also the first occurrences
#and the suffix array
def bwt(text, numCols=1000):
    import numpy
    import time
    t=time.clock()
    numCols=min(numCols,len(text))
    #Build cyclic rotations
    cr=numpy.empty([len(text)],dtype=[("pos", "i4"),("text", "a"+(str)(numCols)), ("initial", "a1"), ("final", "a1")]) #we now must change to numpy.array!
    for i in xrange(len(text)-numCols):
         p0=i-1
         if(i-1<0):
             p0=len(text)-i-1
         cr[i]=(i, text[i:i+numCols],text[i],text[p0])
    for i in xrange(numCols):
         cr[i+len(text)-numCols]=(i+len(text)-numCols, text[len(text)-numCols+i : len(text)]+ text[0:i], text[len(text)-numCols+i], text[len(text)-numCols+i-1])
    print 'time in building M: {}'.format(time.clock()-t)
    t=time.clock()
    cr=numpy.sort(cr, order="text") #best method is default, quicksort. It is still the bottleneck in building
    print 'time in sorting M: {}'.format(time.clock()-t)
    t=time.clock()
    print 'time in taking first and last: {}'.format(time.clock()-t)
    t=time.clock()
    fo={}
    for x in set(cr["initial"]):
        fo[x]=numpy.searchsorted(cr["initial"], x)
    print 'time in getting firstOccurrences: {}'.format(time.clock()-t)
    cp=checkpoints(cr["final"],numCols)
    print 'time in getting checkpoints: {}'.format(time.clock()-t)
    return {"bwt": cr["final"], "firstOccurrence":fo, "suffixArray":cr["pos"], "checkpoints":cp}
#text="panamabananas$"
#pattern="ana"
#t=bwt(text, len(pattern))
#t

#%% ------------ CHECKPOINTS
#Returns a dict with the checkpoints for each symbol in text. 
#Ancillary function for suffix array searches (where text must be the bwt)
def checkpoints(text, k=1000):
    counts={}
    checks={}
    symbols=set(text)
    for s in symbols:
        checks[s]=[]
        counts[s]=0
        checks[s].append(0)
    for i in xrange(len(text)):
        counts[text[i]]+=1    
        if((i+1) % k ==0):
            for s in symbols:
                checks[s].append(counts[s])
    return checks
#text="smnpbnnaaaaa$a"
#cp=checkpoints(text,5)
#cp

#%% ------------- COUNT ---------------------
#Returns the number of occurrences of symbol in text[:pos], using for it 
#the checkpoints of text.
#Ancillary function for suffix array searches.
def count(symbol, pos, text, checkpoints, k=1000):
    if(pos==len(text)-1):
        pos-=1
    for i in xrange(0,len(text),k):
        if(pos-i<0):
            break
    pos0=max(0,i-1)/k
    count=checkpoints[symbol][pos0]
    for i in xrange(pos0*k, pos):
        if(text[i]==symbol):
            count+=1
    return count
#text="smnpbnnaaaaa$a"
#cp=checkpoints(text,5)
#count("a", 13, text, cp, 5)




#%% ------------ SUFFIX ARRAY SEARCH (version 7.0) ---------------
#Searches for EXACT coincidences of in a text characterized by a given 
#suffix array (sa), fist occurrences (fo), bwt (cf) and checkpoints
#All these structures can be obtained by the ancillary functions above.

#We separate data structure construction (sa, fo,cf, checkpoints) from search
#for performance
def bwMatchingV7(pattern, cf, fo, sa, checkpoints, k=1000):
    top=0
    bottom=len(cf)-1
    while top <= bottom:
        if(len(pattern)>0):
            symbol=pattern[len(pattern)-1]
            pattern=pattern[:len(pattern)-1]
            fos=fo[symbol]
            top=fos + count(symbol, top, cf, checkpoints, k)
            bottom=fos + count(symbol, bottom+1, cf, checkpoints, k) - 1
        else:
            r=sa[top:bottom+1]
            return r
    return []




#%% --------- SUFFIX ARRAY SEARCH (version 8.0) -------------------
#Searches for a pattern in a text with up to d mutations (d<3 on large seqs or performance issues)
#It requires the bwt (cf), the first occurrences (fo) and the suffix array (sa)
#from text (provided by ancillary funtion bwt)
#It also requires checkpoints from bwt (provided by ancillary function checkpoints)

#Performance:
#For 1Mbps genomes (V. cholerae) it takes 1) less than 10s to build the data structures
# 2) Less than 0.1s to search for exact matches
# 3) less than 1s to search for up to 1 mutation matches
# 4) about 5s to search for up to 2 mutation mathces
# 5) Less than 1G RAM consumption
#Testing now with human (ch1, 230M) it takes 4min to build the structures and <1G RAM

def bwMatchingV8(text, pattern, cf, fo, sa, checkpoints, k=1000, d=0):
    import numpy
    import time
    step=(int)(round((float)(len(pattern))/(d+1)))
    print "{} {} {}".format(pattern, d, step)
    matches=[]
    for i in range(0, len(pattern), step):
        #1) seed definition
        seed=pattern[i:i+step]
        print "seed: {}".format(seed)
        #2) seed detection
        t0=time.clock()
        result=bwMatchingV7(seed, cf, fo, sa, checkpoints, k)
        print "seed exact matchig takes {}s".format(time.clock()-t0)
        t0=time.clock()
        #print "\tfound in {}".format(result)
        #3) seed extension
        #a) matrix reconstruction
        cr=numpy.empty(len(result),dtype=[("pos", "i4"),("text", "a"+(str)(len(pattern))), ("mismatches", "i4")])
        cont=0
        for r in result:
             if(r-i>=0):
                 cr[cont]=(r, text[(r-i):(r+len(pattern)-i)],0)
                 if((r+len(pattern)-i)>=len(text)):
                     cr[cont]["text"]+=text[:(len(text)-r-len(pattern)-i+1)]
             else:
                 cr[cont]=(r, text[len(text)+r-i:]+text[:r+len(seed)], 0)
             cont+=1
        print "matrix reconstruction takes {}s and has {} candidates".format(time.clock()-t0, len(cr))
        t0=time.clock()
        #print "CYCLIC ROTATIONS: {}".format(cr)
        #b) approximate searh
        for m in xrange(len(pattern)):
            #print "candidates left: {}".format(len(cr))
            symbol=pattern[m]
            for j in xrange(len(cr)):
                t=cr[j]
                if(t["text"][m]!=symbol):
                    t["mismatches"]+=1
            cr=cr[numpy.where(cr["mismatches"]<=d)]
        parray=cr["pos"]-i
        parray=parray[numpy.where(parray>=0)]
        matches.append(parray)
        print "approx search takes {}s".format(time.clock()-t0)
    ret=set()
    for i in xrange(len(matches)):
        ret=ret.union(matches[i])
    ret=list(ret)
    ret.sort()
    return ret


