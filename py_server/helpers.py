# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 16:31:59 2015

Ancillary methods for python server

@author: rodri
"""

#Superslow for interaction: 10 seqs of 600 nucleotides takes 14s
def align(seqs, method="clustalw"):
    import os
    from Bio.Alphabet import generic_dna
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    
    sr=[]
    for k in seqs.keys():
        sr.append(SeqRecord(Seq(seqs[k], generic_dna), id=(str)(k)))
    
    from Bio import SeqIO
    output=open("unaligned.fasta", "w")
    SeqIO.write(sr,output, "fasta")
    output.close()

    import subprocess
    p=subprocess.Popen(["/usr/bin/env", "t_coffee","-quiet -method "+method+" -infile unaligned.fasta -outfile aligned.aln -output=fasta"], bufsize=-1, cwd=u'/Users/rodri/WebstormProjects/seqview/py_server')
    #p=subprocess.Popen(["/usr/bin/env", "t_coffee","-quiet -output "+method+" -infile mitDNAprimates.fasta -outfile aligned.aln"], bufsize=-1, cwd=u'/Users/rodri/WebstormProjects/seqview/py_server')
    #p.wait()
    p.communicate()
   
    lines=open("aligned.aln").readlines()
    aln={}
    k=""
    for l in lines:
        l=l.replace("\n","")
        if(l[0]=='>'):
            k=l.replace(">","").replace("  <unknown description>","")
            aln[k]=""
        else:
            aln[k]+=l
    #os.remove("aligned.aln")
    os.remove("unaligned.fasta")
    os.remove("unaligned.dnd")
    return aln

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


