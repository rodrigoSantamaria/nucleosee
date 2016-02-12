# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 16:31:59 2015

Ancillary methods for python server

@author: rodri
"""

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
    print "s2 is {}".format(s2)
    #solve +
    if(len(s2)==1):
       return s2[0]
    s3=[]
    i=0
    while i<len(s2)-1:
        if(s2[i+1]=='+'):
            s3.append(s2[i]+s2[i+2])
        elif(s2[i]!="+"):
            s3.append(s2[i])
            i+=1
        i+=1
    return s3[0]
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

