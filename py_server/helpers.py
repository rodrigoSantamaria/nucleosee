# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 16:31:59 2015

Ancillary methods for python server

@author: rodri
"""

#%%
"""
This script crops a wig file so max values cannot exceed a given value based
on standard deviations
@author: rodri
"""
    
    # --------------------- INTERNAL METHODS -----------------
#%% -----------  READ WIG --------------
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
    
#tal=readWig("/Users/rodri/WebstormProjects/seqview/py_server/genomes/jpiriz/dwtMini2.wig")
#tal=readWig("/Users/rodri/WebstormProjects/seqview/py_server/genomes/jpiriz/23479_h90_wlt_mean.wig")
#seq=tal["chromosome1"][0]
#np.mean(seq)
#%% ------------------ DISCRETIZATION -------------------
"""Given a numerical sequence seq, this method binarizes based on the average 
and standard deviations on windows of size windowSize. Binarzation is done
in categories a to z (z the larger), as many as detailed by numBins

Percentile - if true, percentiles are used based in numBins, instead of just
   a division of the range between min and maximum (default true)"""

#NOTE: maybe a good idea to optimize this method is to use numerical bins instead of letters
def discretize(seq, windowSize, minimo, maximo, numBins=5, percentile=True):
    import numpy as np
    alphabetTotal=['a','b','c','d','e', 'f', 'g','h','i','j','k','l','m','n','o','p','q','r','s','t']
    alphabet=alphabetTotal[:numBins]   
    dseq=[]
    factor=(numBins-1.0)/float(maximo-minimo)
    
    sseq=rolling_window(seq,windowSize)
    #sseq=np.split(np.array(seq[:windowSize*(len(seq)/windowSize)]), len(seq)/windowSize)
    mseq=np.mean(sseq, axis=1, keepdims=True)
    
    if(percentile==True):
        mseq=np.array(mseq);
        pers=[0]
        for i in range(1,numBins+1):
            p=np.percentile(mseq,(100.0/numBins)*i)
            print("percentile",(100.0/numBins)*i,"=",p)
            pers.append(p)
        bins=pers
        digseq=np.digitize(mseq,pers)
        for s in digseq:
           # print(s)
            dseq.append(alphabet[min(s-1,len(alphabet)-1)])           
    else:
        for im in mseq:
            dseq.append(alphabet[(int)(factor*(im-minimo))])
        bins=[]
        for i in range(numBins):
            bins.append(factor*i)
    return {'dseq':dseq, 'bins':bins} 
    
#0.20s faster for 5M    
#def discretize(seq, windowSize, minimo, maximo, numBins=5):
#    import numpy as np
#    #alphabetTotal=['a','b','c','d','e', 'f', 'g','h','i','j','k','l','m','n','o','p','q','r','s','t']
#    #alphabet=alphabetTotal[:numBins]   
#    #alphabet=np.arange(numBins)
#    factor=(numBins-1.0)/float(maximo-minimo)
#    
#    sseq=helpers.rolling_window(seq,windowSize)
#    mseq=np.mean(sseq, axis=1, keepdims=True)
#    dseq=np.empty(len(mseq),dtype=int)
#    #for im in mseq:
#    for i in xrange(len(mseq)):
#        #dseq.append(alphabet[(int)(factor*(im-minimo))])
#        #dseq.append((int)(factor*(im-minimo)))
#        dseq[i]=factor*(mseq[i]-minimo)
#    return dseq 
#%%
#import time
#import numpy as np
#t0=time.clock()   
#tal=discretize(seq,30, np.min(seq), np.max(seq))
#print '{}'.format((time.clock()-t0))
##%%
#t0=time.clock()   
#tal0=discretize0(seq,30, np.min(seq), np.max(seq))
#print '{}'.format((time.clock()-t0))
#%%
#wig="/Users/rodri/Documents/investigacion/IBFG/quique/mono-H3K9_me2_norm_center_wlt.wig"
#ch=a.readWig(wig)
#f=open(wig.replace(".wig", "")+"-cropped.wig", "w")
#for k in ch.keys():
#    print("computing mean and values")
#    chi=ch[k]
#    th=np.mean(chi)+numSD*np.std(chi)
#    for j in xrange(len(chi)):
#        chi[j]=str(np.min([th,(float)(chi[j])]))
#    print("writing header")
#    f.write("track type=wiggle_0 name=cropped_wig description=\"cropped to "+str(numSD)+" standard deviations\"")
#    f.write("fixedStep chrom="+str(k)+" start=1 step=1")
#    print("writing lines")
#    f.writelines([chi])
#   # return
#f.close()
#%%    


#Superslow for interaction: 10 seqs of 600 nucleotides takes 14s
#def align(seqs, method="clustalw"):
#    import os
#    from Bio.Alphabet import generic_dna
#    from Bio.Seq import Seq
#    from Bio.SeqRecord import SeqRecord
#    
#    sr=[]
#    for k in seqs.keys():
#        sr.append(SeqRecord(Seq(seqs[k], generic_dna), id=(str)(k)))
#    
#    from Bio import SeqIO
#    output=open("unaligned.fasta", "w")
#    SeqIO.write(sr,output, "fasta")
#    output.close()
#
#    import subprocess
#    p=subprocess.Popen(["/usr/bin/env", "t_coffee","-quiet -method "+method+" -infile unaligned.fasta -outfile aligned.aln -output=fasta"], bufsize=-1, cwd=u'/Users/rodri/WebstormProjects/seqview/py_server')
#    #p=subprocess.Popen(["/usr/bin/env", "t_coffee","-quiet -output "+method+" -infile mitDNAprimates.fasta -outfile aligned.aln"], bufsize=-1, cwd=u'/Users/rodri/WebstormProjects/seqview/py_server')
#    #p.wait()
#    p.communicate()
#   
#    lines=open("aligned.aln").readlines()
#    aln={}
#    k=""
#    for l in lines:
#        l=l.replace("\n","")
#        if(l[0]=='>'):
#            k=l.replace(">","").replace("  <unknown description>","")
#            aln[k]=""
#        else:
#            aln[k]+=l
#    #os.remove("aligned.aln")
#    os.remove("unaligned.fasta")
#    os.remove("unaligned.dnd")
#    return aln

#%%
#import time
#t0=time.clock()
#t=align({},"clustalw")
#print("Clustal takes ",(time.clock()-t0))
##t0=time.clock()
##t=align({},"mafft")
##print("Clustal takes ",(time.clock()-t0))
#t0=time.clock()
#t=align({},"tcoffee")
#print("Clustal takes ",(time.clock()-t0))
##t0=time.clock()
##t=align({},"kalign")
##print("Clustal takes ",(time.clock()-t0))
#    
#%%
#Evolved from http://www.rigtorp.se/2011/01/01/rolling-statistics-numpy.html
# (however that solution get overlapping windows)
#a must be a 1D array, does not contemplate overlap so far
def rolling_window(a, window):
    import numpy as np
    shape = (a.shape[0]/window, window)
    strides = (a.strides[0]*window, a.strides[0])
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)

#%%
def parse(text):
    chunks = ['']

    for character in text:
        if character.isalpha():
            if chunks[-1].isalpha():   # If the last chunk is already a number
                chunks[-1] += character  # Add onto that number
            else:
                chunks.append(character) # Start a new number chunk
        if character.isdigit():
            if chunks[-1].isdigit():   # If the last chunk is already a number
                chunks[-1] += character  # Add onto that number
            else:
                chunks.append(character) # Start a new number chunk
        elif character in '+*':
            chunks.append(character)  # This doesn't account for `1 ++ 2`.
    return chunks[1:]

"""
Parses a text which may contain + and * operations (no parenthesis allowed)
"""   
def convertString(text):
    print("Initial text is ",text)
    s=parse(text)
    if(len(s)==1):
        return s[0]
    #solve *
    s2=[]
    i=0
    while(i<len(s)):
        if(i<len(s)-1 and s[i+1]=='*'):
            try:
                int(s[i+2])
                s2.append(s[i]*int(s[i+2]))
                i+=2
            except ValueError:
                s2.append(s[i+2]*int(s[i]))
        elif(s[i]!="*"):
            s2.append(s[i])
        i+=1
    #solve +
    if(len(s2)==1):
       return s2[0]
    s3=""
    i=0
    while i<len(s2)-1:
        if(s2[i+1]=='+'):
            if(s3==""):
                s3+=(s2[i]+s2[i+2])
            else:
                s3+=s2[i+2]
        #else:
            #s3+=(s2[i])
         #   i+=1
        i+=1
    print("Search string is ",s3)
    return s3
    
#convertString("abcba+a*5+abcba")
#%%    
"""
Return the right gff for a given species/chromosome.
This is highly heterogeneous: some organisms have all their chromosomes in a single
gff, some have one gff per chromosome, etc.
By now we are dealing it with this multiplexer function and by now ONLY for S pombe

Another option is to force gff files to a given format.
"""
def gffPath(organism="Schizosaccharomyces pombe", ch="chromosome1"):
    ret="schizosaccharomyces_pombe.I.gff3" #only pombe by now
    if(organism=="Schizosaccharomyces pombe"):
        roman="I"
        if(ch.find("3")>=0 or ch.find("III")>=0):
            roman+="II"
        elif(ch.find("2")>=0 or ch.find("II")>=0):
            roman+="I"
        ret="schizosaccharomyces_pombe."+roman+".gff3"
    return "genomes/annotations/spombe/gff/"+ret
#%%
    
#tal=readWig()
#seq=tal["chromosome2"][0]
##%%
#import numpy as np
#import time
#t0=time.clock()
#np.mean(seq[:30])
#print 'mean takes {}s'.format((time.clock()-t0))

