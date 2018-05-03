# -*- coding: utf-8 -*-
"""
Search methods based on the Burrows-Wheeler Transform.

@author: Rodrigo Santamaría (rodri@usal.es). Universidad de Salamanca
            http://vis.usal.es/rodrigo

License: -GPL3.0 with authorship attribution (extension 7.b) -

    
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  
    
    If not, see <https://www.gnu.org/licenses/gpl.txt>; applying 7.b extension:
    Requiring preservation of specified reasonable legal notices or
    author attributions in that material or in the Appropriate Legal
    Notices displayed by works containing it;   
"""

#----------------- BURROWS-WHEELER TRANSFORM
#Returns a dict with the transform (key 'bwt') and also the first occurrences
#and the suffix array
def bwt(text, numCols=1000):
    import numpy
    import time
    t=time.clock()
    numCols=min(numCols,len(text))
    #Build cyclic rotations
    cr=numpy.empty([len(text)],dtype=[("pos", "i4"),("text", "U"+(str)(numCols)), ("initial", "U1"), ("final", "U1")]) #we now must change to numpy.array!
    #cr=numpy.empty([len(text)],dtype=[("pos", "i4"),("text", a), ("initial", a1), ("final", a1)]) #we now must change to numpy.array!
    for i in range(len(text)-numCols):
         p0=i-1
         if(i-1<0):
             p0=len(text)-i-1
         cr[i]=(i, text[i:i+numCols],text[i],text[p0])
    for i in range(numCols):
         cr[i+len(text)-numCols]=(i+len(text)-numCols, text[len(text)-numCols+i : len(text)]+ text[0:i], text[len(text)-numCols+i], text[len(text)-numCols+i-1])
    print('time in building M:',(time.clock()-t))
    #t=time.clock()
    #cr2=numpy.sort(cr, order="text") #best method is default, quicksort. It is still the bottleneck in building
    #print('time in sorting M:',(time.clock()-t))
    
    t=time.clock()
    #cr=numpy.sort(cr, order="text") #best method is default, quicksort. It is still the bottleneck in building
    cr.sort(order="text") #best method is default, quicksort. It is still the bottleneck in building (a bit faster than the prevous one)
    #sorted(cr["text"]) #maybe a bit faster but not with numpy
    print('time in sorting M:',(time.clock()-t))
    
    t=time.clock()
    print('time in taking first and last:',(time.clock()-t))
    t=time.clock()
    fo={}
    for x in set(cr["initial"]):
        fo[x]=numpy.searchsorted(cr["initial"], x)
    print('time in getting firstOccurrences:',(time.clock()-t))
    cp=checkpoints(cr["final"],numCols)
    #cp=checkpoints((str)(cr["final"]),numCols)
    print('time in getting checkpoints:',(time.clock()-t))
    return {"bwt": cr["final"], "firstOccurrence":fo, "suffixArray":cr["pos"], "checkpoints":cp}
#%%
#text="panamabananas$"
#pattern="ana"
#t=bwt(text, len(pattern))
##t=bwt("".join(dtp["dseq"]), 100)
#res=bwMatchingV7("pa", t["bwt"], t["firstOccurrence"], t["suffixArray"], t["checkpoints"], k=1000)

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
    for i in range(len(text)):
        counts[text[i]]+=1    
        if((i+1) % k ==0):
            for s in symbols:
                checks[s].append(counts[s])
    return checks
#%% TODO: this method is a bit faster but somehow this is giving a list index out of range in count
#def checkpoints1(text, k=1000):
#    checks={}
#    counts={}
#    symbols=set(text)
#    for s in symbols:
#        checks[s]=[]
#        counts[s]=0 
#        checks[s].append(0)
#    for i in range(0,len(text),k):
#        if(i+k<=len(text)):
#            tt=text[i:i+k]
#            for s in symbols:
#                counts[s]+=tt.count(s)
#                checks[s].append(counts[s])
#    return checks
#import time
#text="smnpbnnaaaaa$a"
#
#t0=time.clock()
#cp0=checkpoints0(text,5)
#print((time.clock()-t0))
#t0=time.clock()
#cp=checkpoints(text,5)
#print((time.clock()-t0))

#%% ------------- COUNT ---------------------
#Returns the number of occurrences of symbol in text[:pos], using for it 
#the checkpoints of text.
#Ancillary function for suffix array searches.
def count(symbol, pos, text, checkpoints, k=1000):
    if(pos==len(text)-1):
        pos-=1
    for i in range(0,len(text),k):
        if(pos-i<0):
            break
    pos0=int(max(0,i-1)/k)
    print(symbol, pos0, i, k)
    count=checkpoints[symbol][pos0]
    for i in range(pos0*k, pos):
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
            try:
                fos=fo[symbol]
            except:#in the rare case that the symbol is not in the whole sequence
                print("rare case; ", symbol)
                return []
            top=fos + count(symbol, top, cf, checkpoints, k)
            bottom=fos + count(symbol, bottom+1, cf, checkpoints, k) - 1
        else:
            r=sa[top:bottom+1]
            return r
    return []
#%%
#
#import pickle
#f=open("/home/rodri/workspace/nucleosee/py_server/genomes/h972.pic", "rb")
#data=pickle.load(f, encoding="bytes")
#t=data["processed"]["bwt"]["chromosome1"]
#res=bwMatchingV7(b"c", t["bwt"], t["firstOccurrence"], t["suffixArray"], t["checkpoints"], k=1000)

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
    t00=time.clock()
    step=(int)(round((float)(len(pattern))/(d+1)))
    matches=[]
    for i in range(0, len(pattern), step):
        #1) seed definition
        seed=pattern[i:i+step]
        #2) seed detection
        t0=time.clock()
        result=bwMatchingV7(seed, cf, fo, sa, checkpoints, k)
        print("seed exact matchig takes",(time.clock()-t0))
        t0=time.clock()
        #print "\tfound in {}".format(result)
        #3) seed extension
        #a) matrix reconstruction
        #cr=numpy.empty(len(result),dtype=[("pos", "i4"),("text", "a"+(str)(len(pattern))), ("mismatches", "i4")])
        cr=numpy.empty(len(result),dtype=[("pos", "i4"),("text", "U"+(str)(len(pattern))), ("mismatches", "i4")]) #Python 3
        cont=0
        for r in result:
             if(r-i>=0):
                 cr[cont]=(r, text[(r-i):(r+len(pattern)-i)],0)
                 if((r+len(pattern)-i)>=len(text)):
                     cr[cont]["text"]+=text[:(len(pattern)-(len(text)-r)-i+1)]
             else:
                 cr[cont]=(r, text[len(text)+r-i:]+text[:len(pattern)+r-i], 0)
             cont+=1
        print("matrix reconstruction takes", (time.clock()-t0), "s and has",len(cr)," candidates")
        t0=time.clock()
        #print "CYCLIC ROTATIONS: {}".format(cr)
        #b) approximate searh
        for m in range(len(pattern)):
            #print "candidates left: {}".format(len(cr))
            symbol=pattern[m]
            for j in range(len(cr)):
                t=cr[j]
                if(t["text"][m]!=symbol): #con a*5+abcba(2) falla aqui, se sale de t["text"]=ccaaaa con m=6
                    t["mismatches"]+=1
            cr=cr[numpy.where(cr["mismatches"]<=d)]
        parray=cr["pos"]-i
        parray=parray[numpy.where(parray>=0)]
        matches.append(parray)
        print("approx search takes",(time.clock()-t0),"s")
    ret=set()
    for i in range(len(matches)):
        ret=ret.union(matches[i])
    ret=list(ret)
    ret.sort()
    print("BWT search takes ",(time.clock()-t00),"s")
    return ret

#%%
#res=bwMatchingV8("".join(dtp["dseq"]), "bbbbb", t["bwt"], t["firstOccurrence"], t["suffixArray"], t["checkpoints"])

