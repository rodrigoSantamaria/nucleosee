# -*- coding: utf-8 -*-
"""
SEARCH METHOD FOR GENOME BROWSER
@author: rodri
"""

#%% Perfil de los motivos (frecuencias o probabiliades)
def profile(motifs):
    freqs={'A':[], 'C':[], 'G':[],'T':[],'N':[]}
    for j in xrange(len(motifs[0])): #para cada columna
        count={'A':0, 'C':0, 'T':0, 'G':0, 'N':0}
        for m in motifs:    #para cada motivo
            count[m[j]]+=1
        for k in count.keys():
            freqs[k].append(round((float)(count[k])/len(motifs),1))
    return freqs



#%% ------------------- CONSENSUS ---------------------
#Retorna la secuencia de consenso de una lista de secuencias
#This version makes key control
def consensus0(seqs):
    import operator
    match=[]
    for letter in xrange(len(seqs[0])):
        freq=dict(A=0,G=0,T=0,C=0)
        for s in seqs:
            if(s[letter] in freq.keys()):
                freq[s[letter]]+=1
            else:
                freq[s[letter]]=1
        match.append(max(freq.iteritems(), key=operator.itemgetter(1))[0])
    return reduce(lambda a,b: '{0}{1}'.format(a,b), match)  #HIC SUNT DRACONES: más info sobre las funciones lambda: http://www.secnetix.de/olli/Python/lambda_functions.hawk

#This version requires that non-nucleotides are turned to N (2x speed)
def consensus(seqs):
    import operator
    match=[]
    for j in xrange(len(seqs[0])):
        freq=dict(A=0,G=0,T=0,C=0,N=0)
        for s in seqs:
            freq[s[j]]+=1
        match.append(max(freq.iteritems(), key=operator.itemgetter(1))[0])
    return reduce(lambda a,b: '{0}{1}'.format(a,b), match)  #HIC SUNT DRACONES: más info sobre las funciones lambda: http://www.secnetix.de/olli/Python/lambda_functions.hawk

#%% Puntuación del grupo de motivos, a parte las probabilidades
def score0(motifs):
    s=0
    con=consensus(motifs)
    for m in motifs:
        for j in xrange(len(m)):
            if(con[j]!=m[j]):
                s+=1
    return s
    
#a bit faster version integrating consensus     -> not really
def score(motifs):
    import numpy as np
    import operator
    s=0
    consensus=[]
    motifs2=np.matrix([list(x) for x in motifs])
    for j in xrange(motifs2.shape[1]):
        freq=dict(A=0,G=0,T=0,C=0,N=0)
        col=motifs2[:,j].tostring()
        for k in freq.keys():
            freq[k]=col.count(k)
        conj=max(freq.iteritems(), key=operator.itemgetter(1))[0]#the most frequent letter
        consensus.append(conj)
        s+=len(motifs)-freq[conj]
    return {'score':s, 'consensus':consensus}

# a much faster version with previously computed occurrence profiles
def scoreInc(motifs, profiles):
    import operator
    s=0
    consensus=[]
    for j in xrange(len(motifs[0])):
        freq=dict(A=0,G=0,T=0,C=0,N=0)
        for k in freq.keys():
            freq[k]=profiles[k][j]
        conj=max(freq.iteritems(), key=operator.itemgetter(1))[0]#the most frequent letter
        consensus.append(conj)
        s+=len(motifs)-freq[conj]
    return {'score':s, 'consensus':consensus}   
    
#%% PROBABILITY
import numpy as np
def pr(m,p):
    prod=p[m[0]][0]
    for i in xrange(1, len(m)):
        prod*=p[m[i]][i]
    return prod

def pr1(m,p): #with np.prod (slower)
    prod=np.empty(len(m))
    for i in xrange(len(m)):
        prod[i]=p[m[i]][i]
    return np.prod(prod)

def pr2(m,p): #with np.prod and np.choose (requires transformations, slower) p must be an array -not dict- with the profiles for A,G,C,T resp.)
    #kmer2=m.replace("T", "D").replace("G", "B")
    #ch=[ord(x)-65 for x in kmer2]
    prob=np.prod(np.choose(m, p))
    return prob

#%% MOST PROBABLE PROFILE
def profileMostProbableKmer0(text, k, perfil):
    bestPr=-1
    bestSeq=""
    for i in xrange(len(text)-k+1): #para cada posible k-mer en text
        kmer=text[i:(i+k)]
        newPr=pr(kmer, perfil)
        if(newPr>bestPr):
            bestPr=newPr
            bestSeq=kmer
    return bestSeq
 
#pr internalized   (not a real gain in time though) 
#incremental detection of kmers
def profileMostProbableKmer(text, k, perfil):
    #1) initial kmer
    kmer=text[:k]
    bestSeq=kmer
    bestPr=pr(kmer,perfil)
        
    for i in xrange(k+1, len(text)):
        kmer=kmer[1:]+text[i]
        newPr=pr(kmer, perfil)
        if(newPr>bestPr):
            bestPr=newPr
            bestSeq=kmer
    return bestSeq


# testing a bit faster version by np.choose/np.prod -> not really, it's slower!
def profileMostProbableKmer2(text, k, perfil):
    bestPr=-1
    bestSeq=""
    perfil2=[perfil["A"],  perfil["G"], perfil["C"], perfil["T"]]
    text2=text.replace("T", "D").replace("G", "B")
    text2=[ord(x)-65 for x in text2]
    
    import time
    t0=time.time()
    for i in xrange(len(text2)-k+1): #para cada posible k-mer en text
        kmer=text2[i:(i+k)]
        newPr=pr2(kmer, perfil2)
        if(newPr>bestPr):
            bestPr=newPr
            bestSeq=kmer
    print("internal loop takes",(time.time()-t0))
    ret="".join([chr(x+65) for x in bestSeq]).replace("D", "T").replace("B", "G")
    return ret
#%%
#import time
#motifs=["AAGCTTCACCGGCGCAGTCATTCTCATAATCGCCCACGGACTCACATCCTCATTACTATTCTGCCTAGCAAACTCAAACTACGAACGCACTCACAGTCGCATCATAATCCTCTCTCAAGGACTTCAAACTCTACTCCCACTAATAGCTTTTTGATGACTTCTAGCAAGCCTCGCTAACCTCGCCTTACCCCCCACTATTAACCTACTGGGAGAACTCTCTGTGCTAGTAACCACGTTCTCCTGATCAAATATCACTCTCCTACTTACAGGACTCAACATACTAGTCACAGCCCTATACTCCCTCTACATATTTACCACAACACAATGGGGCTCACTCACCCACCACATTAACAACATAAAACCCTCATTCACACGAGAAAACACCCTCATGTTCATACACCTATCCCCCATTCTCCTCCTATCCCTCAACCCCGACATCATTACCGGGTTTTCCTCTTGTAAATATAGTTTAACCAAAACATCAGATTGTGAATCTGACAACAGAGGCTTACGACCCCTTATTTACCGAGAAAGCTCACAAGAACTGCTAACTCATGCCCCCATGTCTAACAACATGGCTTTCTCAACTTTTAAAGGATAACAGCTATCCATTGGTCTTAGGCCCCAAAAATTTTGGTGCAACTCCAAATAAAAGTAATAACCATGCACACTACTATAACCACCCTAACCCTGACTTCCCTAATTCCCCCCATCCTTACCACCCTCGTTAACCCTAACAAAAAAAACTCATACCCCCATTATGTAAAATCCATTGTCGCATCCACCTTTATTATCAGTCTCTTCCCCACAACAATATTCATGTGCCTAGACCAAGAAGTTATTATCTCGAACTGACACTGAGCCACAACCCAAACAACCCAGCTCTCCCTAAGCTT",
#"AAGCTTCACCGGCGCAGTCATTCTTATAATCGCCCACGGACTTACATCCTCATTACTATTCTGCCTAGCAAACTCAAATTACGAACGCACCCACAGTCGCATCATAATTCTCTCCCAAGGACTTCAAACTCTACTCCCACTAATAGCCTTTTGATGACTCCTAGCAAGCCTCGCCAACCTCGCCCTACCCCCCACCATTAATCTCCTAGGAGAACTCTCCGTGCTAGTAACCTCATTCTCCTGATCAAATACTACCCTCCTACTCACAGGATTCAACATACTAATTACAGCCCTGTACTCCCTCTACATGTTTACCACAACACAATGAGGCTCACTCACCCACCACATTAATAACATAAAACCCTCATTCACACGAGAAAACACTCTCATATTTATACACCTATCCCCCATCCTCCTCCTATCCCTCAATCCTGATATTATCACTGGATTCACCTCCTGTAAATATAGTTTAACCAAAACATCAGATTGTGAATCTGACAACAGAGGCTCACGACCCCTTATTTACCGAGAAAGCTTATAAGAACTGCTAATTCATATCCCCATGCCTAACAACATGGCTTTCTCAACTTTTAAAGGATAACAGCCATCCGTTGGTCTTAGGCCCCAAAAATTTTGGTGCAACTCCAAATAAAAGTAATAACCATGTATACTACCATAACCACCTTAACCCTAACTCCCTTAATTCTCCCCATCCTCACCACCCTCATTAACCCTAACAAAAAAAACTCATATCCCCATTATGTGAAATCCATTATCGCGTCCACCTTTATCATTAGCCTTTTCCCCACAACAATATTCATATGCCTAGACCAAGAAGCTATTATCTCAAACTGGCACTGAGCAACAACCCAAACAACCCAGCTCTCCCTAAGCTT",
#"AAGCTTCACCGGCGCAATTATCCTCATAATCGCCCACGGACTTACATCCTCATTATTATCCTGCCTAGCAAACTCAAATTATGAACGCACCCACAGTCGCATCATAATTCTCTCCCAAGGACTTCAAACTCTACTCCCACTAATAGCCTTTTGATGACTCCTGGCAAGCCTCGCTAACCTCGCCCTACCCCCTACCATTAATCTCCTAGGGGAACTCTCCGTGCTAGTAACCTCATTCTCCTGATCAAATACCACTCTCCTACTCACAGGATTCAACATACTAATCACAGCCCTGTACTCCCTCTACATGTTTACCACAACACAATGAGGCTCACTCACCCACCACATTAATAGCATAAAGCCCTCATTCACACGAGAAAACACTCTCATATTTTTACACCTATCCCCCATCCTCCTTCTATCCCTCAATCCTGATATCATCACTGGATTCACCTCCTGTAAATATAGTTTAACCAAAACATCAGATTGTGAATCTGACAACAGAGGCTCACGACCCCTTATTTACCGAGAAAGCTTATAAGAACTGCTAACTCGTATTCCCATGCCTAACAACATGGCTTTCTCAACTTTTAAAGGATAACAGTTATCCATTGGTCTTAGGCCCCAAAAATTTTGGTGCAACTCCAAATAAAAGTAATAACCATGTATGCTACCATAACCACCTTAGCCCTAACTTCCTTAATTCCCCCCATCCTCGGCGCCCTCATTAACCCTAACAAAAAAAACTCATACCCCCATTACGTGAAATCCATTATCGCATCCACCTTTATCATTAGCCTTTTCCCCACAACAATATTCATATGCCTAGACCAAGAAACTATTATCTCGAACTGACACTGAGCAACAACCCAAACAACCCAACTCTCCCTGAGCTT",
#"AAGCTTCACCGGCGCAGTTGTTCTTATAATTGCCCACGGACTTACATCATCATTATTATTCTGCCTAGCAAACTCAAACTACGAACGAACCCACAGCCGCATCATAATTCTCTCTCAAGGACTCCAAACCCTACTCCCACTAATAGCCCTTTGATGACTTCTGGCAAGCCTCGCCAACCTCGCCTTACCCCCCACCATTAACCTACTAGGAGAGCTCTCCGTACTAGTAACCACATTCTCCTGATCAAACACCACCCTTTTACTTACAGGATCTAACATACTAATTACAGCCCTGTACTCCCTTTATATATTTACCACAACACAATGAGGCCCACTCACACACCACATCACCAACATAAAACCCTCATTTACACGAGAAAACATCCTCATATTCATGCACCTATCCCCCATCCTCCTCCTATCCCTCAACCCCGATATTATCACCGGGTTTACCTCCTGTAAATATAGTTTAACCAAAACATCAGATTGTGAATCTGATAACAGAGGCTCACAACCCCTTATTTACCGAGAAAGCTCGTAAGAGCTGCTAACTCATACCCCCGTGCTTAACAACATGGCTTTCTCAACTTTTAAAGGATAACAGCTATCCATTGGTCTTAGGACCCAAAAATTTTGGTGCAACTCCAAATAAAAGTAATAACTATGTACGCTACCATAACCACCTTAGCCCTAACTTCCTTAATTCCCCCTATCCTTACCACCTTCATCAATCCTAACAAAAAAAGCTCATACCCCCATTACGTAAAATCTATCGTCGCATCCACCTTTATCATCAGCCTCTTCCCCACAACAATATTTCTATGCCTAGACCAAGAAGCTATTATCTCAAGCTGACACTGAGCAACAACCCAAACAATTCAACTCTCCCTAAGCTT",
#"AAGCTTCACCGGCGCAACCACCCTCATGATTGCCCACGGACTCACATCCTCCCTACTATTCTGCCTAGCAAACTCAAACTACGAACGAACCCACAGCCGCATCATAATCCTCTCTCAAGGCCTTCAAACTCTACTCCCCCTGATAGCCCTCTGATGACTTCTAGCAAGCCTCACTAACCTTGCCCTACCACCCACCATCAACCTACTAGGAGAGCTCTCCGTACTAATAGCCATATTCTCTTGATCTAACATCACCATCCTACTAACAGGACTCAACATATTAATCACAGCCCTATACTCCCTCTACATATTCATCACAACACAACGAGGCACACCCTCACACCACATCAACAACATAAAACCCTCTTCCACACGTGAAAACACCCTCATGCTCATACACCTATTCCCTATCCTCCTCCTATCCCTCAACCCCAGCATCATTGCTGGGCTCACCTACTGTAAATATAGTTTAACCAAAACATTAGATTGTGAATCTAATAATAGGGCCCCACAACCCCTTATTTACCGAGAAAGCTCACAAGAACTGCTAACTCCCAGCCCCATGTATAACAACATGGCTTTCTCGACTTTTAAAGGATAACAGCTATCCCTTGGTCTTAGGACCCAAAAATTTTGGTGCAACTCCAAATAAAAGTAACAGCCATGTTCACCACCATAACCGCCCTCACCTTAACTTCCCTAATCCCCCCCATTACCGCTACCCTCATCAACCCCAACAAAAAGAACTCATACCCCCACTATGTAAAAACGGCTATCGCATCCGCCTTTACTATCAGCCTTATCCCAACAACAATATTTATCTGCCTAGGGCAAGAAACCATCATCACAAACTGATGTTGAACAACCACCCAAACACTGCAACTCTCACTAAGCTT",
#"AAGCTTTACAGGTGCAACCGTCCTCATAATCGCCCACGGACTAACCTCTTCCCTGCTATTCTGCCTTGCAAACTCAAACTACGAACGAACTCACAGCCGCATCATAATCCTATCTCGAGGGCTCCAAGCCTTACTCCCACTGATAGCCTTCTGATGACTCGCAGCAAGCCTCGCTAACCTCGCCCTACCCCCCACTATTAACCTCCTAGGTGAACTCTTCGTACTAATGGCCTCCTTCTCCTGGGCAAACACTACTATTACACTCACCGGGCTCAACGTACTAATCACGGCCCTATACTCTCTTTACATATTTATCATAACACAACGAGGCACACTTACACACCACATTAAAAACATAAAACCCTCACTCACACGAGAAAACATATTAATACTTATGCACCTCTTCCCCCTCCTCCTCCTAACCCTCAACCCTAACATCATTACTGGCTTTACTCCCTGTAAACATAGTTTAATCAAAACATTAGATTGTGAATCTAACAATAGAGGCTCGAAACCTCTTGCTTACCGAGAAAGCCCACAAGAACTGCTAACTCACTATCCCATGTATAACAACATGGCTTTCTCAACTTTTAAAGGATAACAGCTATCCATTGGTCTTAGGACCCAAAAATTTTGGTGCAACTCCAAATAAAAGTAATAGCAATGTACACCACCATAGCCATTCTAACGCTAACCTCCCTAATTCCCCCCATTACAGCCACCCTTATTAACCCCAATAAAAAGAACTTATACCCGCACTACGTAAAAATGACCATTGCCTCTACCTTTATAATCAGCCTATTTCCCACAATAATATTCATGTGCACAGACCAAGAAACCATTATTTCAAACTGACACTGAACTGCAACCCAAACGCTAGAACTCTCCCTAAGCTT"]
#motifs=["AAGCTTCAC","AAGGTTCAC","AAGATACAC","ATGCTTCAC","AAGCTTCGG","AAGCTTCTT"]
#t0=time.time()
#p=profileLaplace(motifs)
#print("Laplace took", (time.time()-t0))
#seq="AAGCTTCACCGGCGCAGTCATTCTCATAATCGCCCACGGACTCACATCCTCATTACTATTCTGCCTAGCAAACTCAAACTACGAACGCACTCACAGTCGCATCATAATCCTCTCTCAAGGACTTCAAACTCTACTCCCACTAATAGCTTTTTGATGACTTCTAGCAAGCCTCGCTAACCTCGCCTTACCCCCCACTATTAACCTACTGGGAGAACTCTCTGTGCTAGTAACCACGTTCTCCTGATCAAATATCACTCTCCTACTTACAGGACTCAACATACTAGTCACAGCCCTATACTCCCTCTACATATTTACCACAACACAATGGGGCTCACTCACCCACCACATTAACAACATAAAACCCTCATTCACACGAGAAAACACCCTCATGTTCATACACCTATCCCCCATTCTCCTCCTATCCCTCAACCCCGACATCATTACCGGGTTTTCCTCTTGTAAATATAGTTTAACCAAAACATCAGATTGTGAATCTGACAACAGAGGCTTACGACCCCTTATTTACCGAGAAAGCTCACAAGAACTGCTAACTCATGCCCCCATGTAAAACAACATGGCTTTCTCAACTTTTAAAGGATAACAGCTATCCATTGGTCTTAGGCCCCAAAAATTTTGGTGCAACTCCAAATAAAAGTAATAACCATGCACACTACTATAACCACCCTAACCCTGACTTCCCTAATTCCCCCCATCCTTACCACCCTCGTTAACCCTAACAAAAAAAACTCATACCCCCATTATGTAAAATCTATTGTCGCATCCACCTTTATTATCAGTTTCTTCCCCACAACAATATTCATGTGCCTAGACCAAGAAGTTATTATCTCGAACTGACACTGAGCCACAACCCAAACAACCCAGCTCTCCCTAAGCTT"
#t0=time.time()
#mpk=profileMostProbableKmer(seq, 9,p)
#print(mpk) # AAGCTTCAC
#print("MPK took", (time.time()-t0))
#t0=time.time()
#mpk=profileMostProbableKmer2(seq, 9,p)
#print(mpk) # AAGCTTCAC
#print("MPK2 took", (time.time()-t0))


#%% Profile Laplace
def profileLaplace0(motifs):
    freqs={'A':[], 'C':[], 'G':[],'T':[], 'N':[]}
    den=(float)(len(motifs)*2)
    for j in xrange(len(motifs[0])): #para cada columna
        count={'A':1, 'C':1, 'T':1, 'G':1, 'N':1}
        for m in motifs:    #para cada motivo
            count[m[j]]+=1
        for k in count.keys():
            freqs[k].append(count[k]/den)
    return freqs
    
# Profile Laplace (returning also raw occurences)
def profileLaplace(motifs):
    freqs={'A':[], 'C':[], 'G':[],'T':[], 'N':[]}
    ocs={'A':[], 'C':[], 'G':[],'T':[], 'N':[]}
    den=(float)(len(motifs)*2)
    for j in xrange(len(motifs[0])): #para cada columna
        count={'A':1, 'C':1, 'T':1, 'G':1, 'N':1}
        for m in motifs:    #para cada motivo
            count[m[j]]+=1
        for k in count.keys():
            freqs[k].append(count[k]/den)
            ocs[k].append(count[k]-1)
    return {"frequencies":freqs, "occurrences":ocs}


#Incremental version to save time: saves the for m in motifs loop
#we just take the old profile an substituthe the ith row by the profile of the
#new motif m
def profileLaplaceInc(old_ocs, old_m, m, num_motifs):
    freqs={'A':[], 'C':[], 'G':[],'T':[], 'N':[]}
    den=(float)(num_motifs*2)
    
    ocs=old_ocs
    for j in xrange(len(m)):
        ocs[m[j]][j]+=1
        ocs[old_m[j]][j]-=1
        for k in ocs.keys():
            freqs[k].append((ocs[k][j]+1)/den)
    return {"frequencies":freqs, "occurrences":ocs}

#import time
##motifs=["AAGCTTCACCGGCGCAGTCATTCTCATAATCGCCCACGGACTCACATCCTCATTACTATTCTGCCTAGCAAACTCAAACTACGAACGCACTCACAGTCGCATCATAATCCTCTCTCAAGGACTTCAAACTCTACTCCCACTAATAGCTTTTTGATGACTTCTAGCAAGCCTCGCTAACCTCGCCTTACCCCCCACTATTAACCTACTGGGAGAACTCTCTGTGCTAGTAACCACGTTCTCCTGATCAAATATCACTCTCCTACTTACAGGACTCAACATACTAGTCACAGCCCTATACTCCCTCTACATATTTACCACAACACAATGGGGCTCACTCACCCACCACATTAACAACATAAAACCCTCATTCACACGAGAAAACACCCTCATGTTCATACACCTATCCCCCATTCTCCTCCTATCCCTCAACCCCGACATCATTACCGGGTTTTCCTCTTGTAAATATAGTTTAACCAAAACATCAGATTGTGAATCTGACAACAGAGGCTTACGACCCCTTATTTACCGAGAAAGCTCACAAGAACTGCTAACTCATGCCCCCATGTCTAACAACATGGCTTTCTCAACTTTTAAAGGATAACAGCTATCCATTGGTCTTAGGCCCCAAAAATTTTGGTGCAACTCCAAATAAAAGTAATAACCATGCACACTACTATAACCACCCTAACCCTGACTTCCCTAATTCCCCCCATCCTTACCACCCTCGTTAACCCTAACAAAAAAAACTCATACCCCCATTATGTAAAATCCATTGTCGCATCCACCTTTATTATCAGTCTCTTCCCCACAACAATATTCATGTGCCTAGACCAAGAAGTTATTATCTCGAACTGACACTGAGCCACAACCCAAACAACCCAGCTCTCCCTAAGCTT",
##"AAGCTTCACCGGCGCAGTCATTCTTATAATCGCCCACGGACTTACATCCTCATTACTATTCTGCCTAGCAAACTCAAATTACGAACGCACCCACAGTCGCATCATAATTCTCTCCCAAGGACTTCAAACTCTACTCCCACTAATAGCCTTTTGATGACTCCTAGCAAGCCTCGCCAACCTCGCCCTACCCCCCACCATTAATCTCCTAGGAGAACTCTCCGTGCTAGTAACCTCATTCTCCTGATCAAATACTACCCTCCTACTCACAGGATTCAACATACTAATTACAGCCCTGTACTCCCTCTACATGTTTACCACAACACAATGAGGCTCACTCACCCACCACATTAATAACATAAAACCCTCATTCACACGAGAAAACACTCTCATATTTATACACCTATCCCCCATCCTCCTCCTATCCCTCAATCCTGATATTATCACTGGATTCACCTCCTGTAAATATAGTTTAACCAAAACATCAGATTGTGAATCTGACAACAGAGGCTCACGACCCCTTATTTACCGAGAAAGCTTATAAGAACTGCTAATTCATATCCCCATGCCTAACAACATGGCTTTCTCAACTTTTAAAGGATAACAGCCATCCGTTGGTCTTAGGCCCCAAAAATTTTGGTGCAACTCCAAATAAAAGTAATAACCATGTATACTACCATAACCACCTTAACCCTAACTCCCTTAATTCTCCCCATCCTCACCACCCTCATTAACCCTAACAAAAAAAACTCATATCCCCATTATGTGAAATCCATTATCGCGTCCACCTTTATCATTAGCCTTTTCCCCACAACAATATTCATATGCCTAGACCAAGAAGCTATTATCTCAAACTGGCACTGAGCAACAACCCAAACAACCCAGCTCTCCCTAAGCTT",
##"AAGCTTCACCGGCGCAATTATCCTCATAATCGCCCACGGACTTACATCCTCATTATTATCCTGCCTAGCAAACTCAAATTATGAACGCACCCACAGTCGCATCATAATTCTCTCCCAAGGACTTCAAACTCTACTCCCACTAATAGCCTTTTGATGACTCCTGGCAAGCCTCGCTAACCTCGCCCTACCCCCTACCATTAATCTCCTAGGGGAACTCTCCGTGCTAGTAACCTCATTCTCCTGATCAAATACCACTCTCCTACTCACAGGATTCAACATACTAATCACAGCCCTGTACTCCCTCTACATGTTTACCACAACACAATGAGGCTCACTCACCCACCACATTAATAGCATAAAGCCCTCATTCACACGAGAAAACACTCTCATATTTTTACACCTATCCCCCATCCTCCTTCTATCCCTCAATCCTGATATCATCACTGGATTCACCTCCTGTAAATATAGTTTAACCAAAACATCAGATTGTGAATCTGACAACAGAGGCTCACGACCCCTTATTTACCGAGAAAGCTTATAAGAACTGCTAACTCGTATTCCCATGCCTAACAACATGGCTTTCTCAACTTTTAAAGGATAACAGTTATCCATTGGTCTTAGGCCCCAAAAATTTTGGTGCAACTCCAAATAAAAGTAATAACCATGTATGCTACCATAACCACCTTAGCCCTAACTTCCTTAATTCCCCCCATCCTCGGCGCCCTCATTAACCCTAACAAAAAAAACTCATACCCCCATTACGTGAAATCCATTATCGCATCCACCTTTATCATTAGCCTTTTCCCCACAACAATATTCATATGCCTAGACCAAGAAACTATTATCTCGAACTGACACTGAGCAACAACCCAAACAACCCAACTCTCCCTGAGCTT",
##"AAGCTTCACCGGCGCAGTTGTTCTTATAATTGCCCACGGACTTACATCATCATTATTATTCTGCCTAGCAAACTCAAACTACGAACGAACCCACAGCCGCATCATAATTCTCTCTCAAGGACTCCAAACCCTACTCCCACTAATAGCCCTTTGATGACTTCTGGCAAGCCTCGCCAACCTCGCCTTACCCCCCACCATTAACCTACTAGGAGAGCTCTCCGTACTAGTAACCACATTCTCCTGATCAAACACCACCCTTTTACTTACAGGATCTAACATACTAATTACAGCCCTGTACTCCCTTTATATATTTACCACAACACAATGAGGCCCACTCACACACCACATCACCAACATAAAACCCTCATTTACACGAGAAAACATCCTCATATTCATGCACCTATCCCCCATCCTCCTCCTATCCCTCAACCCCGATATTATCACCGGGTTTACCTCCTGTAAATATAGTTTAACCAAAACATCAGATTGTGAATCTGATAACAGAGGCTCACAACCCCTTATTTACCGAGAAAGCTCGTAAGAGCTGCTAACTCATACCCCCGTGCTTAACAACATGGCTTTCTCAACTTTTAAAGGATAACAGCTATCCATTGGTCTTAGGACCCAAAAATTTTGGTGCAACTCCAAATAAAAGTAATAACTATGTACGCTACCATAACCACCTTAGCCCTAACTTCCTTAATTCCCCCTATCCTTACCACCTTCATCAATCCTAACAAAAAAAGCTCATACCCCCATTACGTAAAATCTATCGTCGCATCCACCTTTATCATCAGCCTCTTCCCCACAACAATATTTCTATGCCTAGACCAAGAAGCTATTATCTCAAGCTGACACTGAGCAACAACCCAAACAATTCAACTCTCCCTAAGCTT",
##"AAGCTTCACCGGCGCAACCACCCTCATGATTGCCCACGGACTCACATCCTCCCTACTATTCTGCCTAGCAAACTCAAACTACGAACGAACCCACAGCCGCATCATAATCCTCTCTCAAGGCCTTCAAACTCTACTCCCCCTGATAGCCCTCTGATGACTTCTAGCAAGCCTCACTAACCTTGCCCTACCACCCACCATCAACCTACTAGGAGAGCTCTCCGTACTAATAGCCATATTCTCTTGATCTAACATCACCATCCTACTAACAGGACTCAACATATTAATCACAGCCCTATACTCCCTCTACATATTCATCACAACACAACGAGGCACACCCTCACACCACATCAACAACATAAAACCCTCTTCCACACGTGAAAACACCCTCATGCTCATACACCTATTCCCTATCCTCCTCCTATCCCTCAACCCCAGCATCATTGCTGGGCTCACCTACTGTAAATATAGTTTAACCAAAACATTAGATTGTGAATCTAATAATAGGGCCCCACAACCCCTTATTTACCGAGAAAGCTCACAAGAACTGCTAACTCCCAGCCCCATGTATAACAACATGGCTTTCTCGACTTTTAAAGGATAACAGCTATCCCTTGGTCTTAGGACCCAAAAATTTTGGTGCAACTCCAAATAAAAGTAACAGCCATGTTCACCACCATAACCGCCCTCACCTTAACTTCCCTAATCCCCCCCATTACCGCTACCCTCATCAACCCCAACAAAAAGAACTCATACCCCCACTATGTAAAAACGGCTATCGCATCCGCCTTTACTATCAGCCTTATCCCAACAACAATATTTATCTGCCTAGGGCAAGAAACCATCATCACAAACTGATGTTGAACAACCACCCAAACACTGCAACTCTCACTAAGCTT",
##"AAGCTTTACAGGTGCAACCGTCCTCATAATCGCCCACGGACTAACCTCTTCCCTGCTATTCTGCCTTGCAAACTCAAACTACGAACGAACTCACAGCCGCATCATAATCCTATCTCGAGGGCTCCAAGCCTTACTCCCACTGATAGCCTTCTGATGACTCGCAGCAAGCCTCGCTAACCTCGCCCTACCCCCCACTATTAACCTCCTAGGTGAACTCTTCGTACTAATGGCCTCCTTCTCCTGGGCAAACACTACTATTACACTCACCGGGCTCAACGTACTAATCACGGCCCTATACTCTCTTTACATATTTATCATAACACAACGAGGCACACTTACACACCACATTAAAAACATAAAACCCTCACTCACACGAGAAAACATATTAATACTTATGCACCTCTTCCCCCTCCTCCTCCTAACCCTCAACCCTAACATCATTACTGGCTTTACTCCCTGTAAACATAGTTTAATCAAAACATTAGATTGTGAATCTAACAATAGAGGCTCGAAACCTCTTGCTTACCGAGAAAGCCCACAAGAACTGCTAACTCACTATCCCATGTATAACAACATGGCTTTCTCAACTTTTAAAGGATAACAGCTATCCATTGGTCTTAGGACCCAAAAATTTTGGTGCAACTCCAAATAAAAGTAATAGCAATGTACACCACCATAGCCATTCTAACGCTAACCTCCCTAATTCCCCCCATTACAGCCACCCTTATTAACCCCAATAAAAAGAACTTATACCCGCACTACGTAAAAATGACCATTGCCTCTACCTTTATAATCAGCCTATTTCCCACAATAATATTCATGTGCACAGACCAAGAAACCATTATTTCAAACTGACACTGAACTGCAACCCAAACGCTAGAACTCTCCCTAAGCTT"]
#motifs=["AAGCTTCAC","AAGGTTCAC","AAGATACAC","ATGCTTCAC","AAGCTTCGG","AAGCTTCTT"]
#t0=time.time()
#p=profileLaplace(motifs)
#print("Laplace took", (time.time()-t0))
#print(p["occurrences"])
#
#t0=time.time()
#pi=profileLaplaceInc(p["occurrences"], motifs[2], "TTTTTTTTT")
#print("LaplaceInc took", (time.time()-t0))
#print(pi["occurrences"])

#%% GIBBS SAMPLER
#N is the number of iterations we want for the optimization
#TODO: change dna to numpy matrices for better performance?
def gibbsSampler0(dna, k, N):
    import time
    t0=time.time()
    t=len(dna)
    #0) Random selection of initial motifs
    import random
    random.seed()
    motivos=[]
    for i in xrange(t):
        ri=random.randint(0,len(dna[i])-k) 
        motivos.append(dna[i][ri:ri+k])
    bestMotifs= [x[:k] for x in dna]
    sb=score(bestMotifs)
    t0=time.time()
    #1) iteration
    for j in xrange(N):
        i=random.randint(0,t-1)
        del motivos[i]
        
        prof=profileLaplace0(motivos)
        motivos.insert(i, profileMostProbableKmer(dna[i], k, prof)) #perfil más probable en la iésima cadena de Dna
        temp=score(motivos)
        s=temp['score']
        
        if(s<sb):
            bestMotifs=motivos
            sb=s            
    print("Iterations took", (time.time()-t0))
    print("best motif is ",consensus(bestMotifs), "with score",sb)
    return {"motifs":bestMotifs, "score":sb}

#%%    
def gibbsSampler(dna, k, N):
    import time
    t0=time.time()
    t=len(dna)
    
    #0) Random selection of initial motifs
    t00=time.time()
    import random
    random.seed()
    motivos=[]
    for i in xrange(t):
        ri=random.randint(0,len(dna[i])-k) 
        motivos.append(dna[i][ri:ri+k])
    bestMotifs= [x[:k] for x in dna]
    temp=score(bestMotifs)
    sb=temp["score"]
    
    #1.0) First iteration (for incremental implementation of Laplace)
    i=random.randint(0,t-1)
    oldmot=motivos[i]
    del motivos[i]
    antprof=profileLaplace(motivos)
    newmot=profileMostProbableKmer(dna[i], k, antprof["frequencies"])
    motivos.insert(i, newmot) #perfil más probable en la iésima cadena de Dna
    temp=scoreInc(motivos, antprof["frequencies"]) #incremental score (faster)
    s=temp['score']
    if(s<sb):
        bestMotifs=motivos
        sb=s            
    
    #1) iterations
    for j in xrange(1,N):
        i=random.randint(0,t-1)

        oldmot=motivos[i]#NEW
        del motivos[i]
        
        #ti=time.time()
        #prof=profileLaplace(motivos)
        prof=profileLaplaceInc(antprof["occurrences"], oldmot, newmot, len(motivos))#NEW
        #print("profile took",(time.time()-ti))
        #ti=time.time()
        newmot=profileMostProbableKmer(dna[i], k, prof["frequencies"])
        motivos.insert(i, newmot) #perfil más probable en la iésima cadena de Dna
        #print("PMPK took",(time.time()-ti))
        #ti=time.time()
        temp=scoreInc(motivos, prof["occurrences"]) #incremental score (faster)
        s=temp['score']
        #con=temp['consensus']
        #print("score took",(time.time()-ti), "and is ",s)
        antprof=prof
        
        if(s<sb):
            #print(j,") score improvement: ", s)
            bestMotifs=motivos
            sb=s            
    print("Gibbs took", (time.time()-t00))
    print("best motif is ",consensus(bestMotifs), "with score",sb)
    return {"motifs":bestMotifs, "score":sb}

#%%
import time
motifs=["AAGCTTCACCGGCGCAGTCATTCTCATAATCGCCCACGGACTCACATCCTCATTACTATTCTGCCTAGCAAACTCAAACTACGAACGCACTCACAGTCGCATCATAATCCTCTCTCAAGGACTTCAAACTCTACTCCCACTAATAGCTTTTTGATGACTTCTAGCAAGCCTCGCTAACCTCGCCTTACCCCCCACTATTAACCTACTGGGAGAACTCTCTGTGCTAGTAACCACGTTCTCCTGATCAAATATCACTCTCCTACTTACAGGACTCAACATACTAGTCACAGCCCTATACTCCCTCTACATATTTACCACAACACAATGGGGCTCACTCACCCACCACATTAACAACATAAAACCCTCATTCACACGAGAAAACACCCTCATGTTCATACACCTATCCCCCATTCTCCTCCTATCCCTCAACCCCGACATCATTACCGGGTTTTCCTCTTGTAAATATAGTTTAACCAAAACATCAGATTGTGAATCTGACAACAGAGGCTTACGACCCCTTATTTACCGAGAAAGCTCACAAGAACTGCTAACTCATGCCCCCATGTCTAACAACATGGCTTTCTCAACTTTTAAAGGATAACAGCTATCCATTGGTCTTAGGCCCCAAAAATTTTGGTGCAACTCCAAATAAAAGTAATAACCATGCACACTACTATAACCACCCTAACCCTGACTTCCCTAATTCCCCCCATCCTTACCACCCTCGTTAACCCTAACAAAAAAAACTCATACCCCCATTATGTAAAATCCATTGTCGCATCCACCTTTATTATCAGTCTCTTCCCCACAACAATATTCATGTGCCTAGACCAAGAAGTTATTATCTCGAACTGACACTGAGCCACAACCCAAACAACCCAGCTCTCCCTAAGCTT",
"AAGCTTCACCGGCGCAGTCATTCTTATAATCGCCCACGGACTTACATCCTCATTACTATTCTGCCTAGCAAACTCAAATTACGAACGCACCCACAGTCGCATCATAATTCTCTCCCAAGGACTTCAAACTCTACTCCCACTAATAGCCTTTTGATGACTCCTAGCAAGCCTCGCCAACCTCGCCCTACCCCCCACCATTAATCTCCTAGGAGAACTCTCCGTGCTAGTAACCTCATTCTCCTGATCAAATACTACCCTCCTACTCACAGGATTCAACATACTAATTACAGCCCTGTACTCCCTCTACATGTTTACCACAACACAATGAGGCTCACTCACCCACCACATTAATAACATAAAACCCTCATTCACACGAGAAAACACTCTCATATTTATACACCTATCCCCCATCCTCCTCCTATCCCTCAATCCTGATATTATCACTGGATTCACCTCCTGTAAATATAGTTTAACCAAAACATCAGATTGTGAATCTGACAACAGAGGCTCACGACCCCTTATTTACCGAGAAAGCTTATAAGAACTGCTAATTCATATCCCCATGCCTAACAACATGGCTTTCTCAACTTTTAAAGGATAACAGCCATCCGTTGGTCTTAGGCCCCAAAAATTTTGGTGCAACTCCAAATAAAAGTAATAACCATGTATACTACCATAACCACCTTAACCCTAACTCCCTTAATTCTCCCCATCCTCACCACCCTCATTAACCCTAACAAAAAAAACTCATATCCCCATTATGTGAAATCCATTATCGCGTCCACCTTTATCATTAGCCTTTTCCCCACAACAATATTCATATGCCTAGACCAAGAAGCTATTATCTCAAACTGGCACTGAGCAACAACCCAAACAACCCAGCTCTCCCTAAGCTT",
"AAGCTTCACCGGCGCAATTATCCTCATAATCGCCCACGGACTTACATCCTCATTATTATCCTGCCTAGCAAACTCAAATTATGAACGCACCCACAGTCGCATCATAATTCTCTCCCAAGGACTTCAAACTCTACTCCCACTAATAGCCTTTTGATGACTCCTGGCAAGCCTCGCTAACCTCGCCCTACCCCCTACCATTAATCTCCTAGGGGAACTCTCCGTGCTAGTAACCTCATTCTCCTGATCAAATACCACTCTCCTACTCACAGGATTCAACATACTAATCACAGCCCTGTACTCCCTCTACATGTTTACCACAACACAATGAGGCTCACTCACCCACCACATTAATAGCATAAAGCCCTCATTCACACGAGAAAACACTCTCATATTTTTACACCTATCCCCCATCCTCCTTCTATCCCTCAATCCTGATATCATCACTGGATTCACCTCCTGTAAATATAGTTTAACCAAAACATCAGATTGTGAATCTGACAACAGAGGCTCACGACCCCTTATTTACCGAGAAAGCTTATAAGAACTGCTAACTCGTATTCCCATGCCTAACAACATGGCTTTCTCAACTTTTAAAGGATAACAGTTATCCATTGGTCTTAGGCCCCAAAAATTTTGGTGCAACTCCAAATAAAAGTAATAACCATGTATGCTACCATAACCACCTTAGCCCTAACTTCCTTAATTCCCCCCATCCTCGGCGCCCTCATTAACCCTAACAAAAAAAACTCATACCCCCATTACGTGAAATCCATTATCGCATCCACCTTTATCATTAGCCTTTTCCCCACAACAATATTCATATGCCTAGACCAAGAAACTATTATCTCGAACTGACACTGAGCAACAACCCAAACAACCCAACTCTCCCTGAGCTT",
"AAGCTTCACCGGCGCAGTTGTTCTTATAATTGCCCACGGACTTACATCATCATTATTATTCTGCCTAGCAAACTCAAACTACGAACGAACCCACAGCCGCATCATAATTCTCTCTCAAGGACTCCAAACCCTACTCCCACTAATAGCCCTTTGATGACTTCTGGCAAGCCTCGCCAACCTCGCCTTACCCCCCACCATTAACCTACTAGGAGAGCTCTCCGTACTAGTAACCACATTCTCCTGATCAAACACCACCCTTTTACTTACAGGATCTAACATACTAATTACAGCCCTGTACTCCCTTTATATATTTACCACAACACAATGAGGCCCACTCACACACCACATCACCAACATAAAACCCTCATTTACACGAGAAAACATCCTCATATTCATGCACCTATCCCCCATCCTCCTCCTATCCCTCAACCCCGATATTATCACCGGGTTTACCTCCTGTAAATATAGTTTAACCAAAACATCAGATTGTGAATCTGATAACAGAGGCTCACAACCCCTTATTTACCGAGAAAGCTCGTAAGAGCTGCTAACTCATACCCCCGTGCTTAACAACATGGCTTTCTCAACTTTTAAAGGATAACAGCTATCCATTGGTCTTAGGACCCAAAAATTTTGGTGCAACTCCAAATAAAAGTAATAACTATGTACGCTACCATAACCACCTTAGCCCTAACTTCCTTAATTCCCCCTATCCTTACCACCTTCATCAATCCTAACAAAAAAAGCTCATACCCCCATTACGTAAAATCTATCGTCGCATCCACCTTTATCATCAGCCTCTTCCCCACAACAATATTTCTATGCCTAGACCAAGAAGCTATTATCTCAAGCTGACACTGAGCAACAACCCAAACAATTCAACTCTCCCTAAGCTT",
"AAGCTTCACCGGCGCAACCACCCTCATGATTGCCCACGGACTCACATCCTCCCTACTATTCTGCCTAGCAAACTCAAACTACGAACGAACCCACAGCCGCATCATAATCCTCTCTCAAGGCCTTCAAACTCTACTCCCCCTGATAGCCCTCTGATGACTTCTAGCAAGCCTCACTAACCTTGCCCTACCACCCACCATCAACCTACTAGGAGAGCTCTCCGTACTAATAGCCATATTCTCTTGATCTAACATCACCATCCTACTAACAGGACTCAACATATTAATCACAGCCCTATACTCCCTCTACATATTCATCACAACACAACGAGGCACACCCTCACACCACATCAACAACATAAAACCCTCTTCCACACGTGAAAACACCCTCATGCTCATACACCTATTCCCTATCCTCCTCCTATCCCTCAACCCCAGCATCATTGCTGGGCTCACCTACTGTAAATATAGTTTAACCAAAACATTAGATTGTGAATCTAATAATAGGGCCCCACAACCCCTTATTTACCGAGAAAGCTCACAAGAACTGCTAACTCCCAGCCCCATGTATAACAACATGGCTTTCTCGACTTTTAAAGGATAACAGCTATCCCTTGGTCTTAGGACCCAAAAATTTTGGTGCAACTCCAAATAAAAGTAACAGCCATGTTCACCACCATAACCGCCCTCACCTTAACTTCCCTAATCCCCCCCATTACCGCTACCCTCATCAACCCCAACAAAAAGAACTCATACCCCCACTATGTAAAAACGGCTATCGCATCCGCCTTTACTATCAGCCTTATCCCAACAACAATATTTATCTGCCTAGGGCAAGAAACCATCATCACAAACTGATGTTGAACAACCACCCAAACACTGCAACTCTCACTAAGCTT",
"AAGCTTTACAGGTGCAACCGTCCTCATAATCGCCCACGGACTAACCTCTTCCCTGCTATTCTGCCTTGCAAACTCAAACTACGAACGAACTCACAGCCGCATCATAATCCTATCTCGAGGGCTCCAAGCCTTACTCCCACTGATAGCCTTCTGATGACTCGCAGCAAGCCTCGCTAACCTCGCCCTACCCCCCACTATTAACCTCCTAGGTGAACTCTTCGTACTAATGGCCTCCTTCTCCTGGGCAAACACTACTATTACACTCACCGGGCTCAACGTACTAATCACGGCCCTATACTCTCTTTACATATTTATCATAACACAACGAGGCACACTTACACACCACATTAAAAACATAAAACCCTCACTCACACGAGAAAACATATTAATACTTATGCACCTCTTCCCCCTCCTCCTCCTAACCCTCAACCCTAACATCATTACTGGCTTTACTCCCTGTAAACATAGTTTAATCAAAACATTAGATTGTGAATCTAACAATAGAGGCTCGAAACCTCTTGCTTACCGAGAAAGCCCACAAGAACTGCTAACTCACTATCCCATGTATAACAACATGGCTTTCTCAACTTTTAAAGGATAACAGCTATCCATTGGTCTTAGGACCCAAAAATTTTGGTGCAACTCCAAATAAAAGTAATAGCAATGTACACCACCATAGCCATTCTAACGCTAACCTCCCTAATTCCCCCCATTACAGCCACCCTTATTAACCCCAATAAAAAGAACTTATACCCGCACTACGTAAAAATGACCATTGCCTCTACCTTTATAATCAGCCTATTTCCCACAATAATATTCATGTGCACAGACCAAGAAACCATTATTTCAAACTGACACTGAACTGCAACCCAAACGCTAGAACTCTCCCTAAGCTT"]
t0=time.time()
res=gibbsSampler(motifs, 9,100)
print("Gibss took", (time.time()-t0))
#%%
prof=profileLaplace(motifs)

#%%
    


















#%% --------------------------- K-MERS --------------------
# Contar los k-mers que aparecen n o más veces en seq
def countKmers(seq,n,k):
    kmers={}    
    kfinals={}  
    for i in xrange(len(seq)-k+1):
        kmer=seq[i:(i+k)]
        if(kmer in kmers.keys()):	
            kmers[kmer]+=1
        else:			
            kmers[kmer]=1
    for i in kmers.keys():
        if(kmers[i]>=n):
            kfinals[i]=kmers[i]
    return kfinals
    



#%% -------------------- FIND PATTERN ------------------
#Función que busca la secuencia pattern en la secuencia seq
#Retorna un diccionario con la posición como clave y el patrón como valor
def findPattern(pattern,seq):
    pos={}
    for i in xrange(0, len(seq)-len(pattern)+1):
        s=seq[i:(i+len(pattern))]
        if(s==pattern):
            pos[i]=pattern
    return(pos)


#%% ------------------ FIND PATTERN (II) ---------------
# Función que busca la secuencia pattern (o una mutación suya de hasta d SNPs)
# en la secuencia seq
#Retorna un diccionario con la posición como clave y el patrón como valor
def findPattern(pattern,seq, d=3):
    pos={}
    for i in xrange(0, len(seq)-len(pattern)+1):
        s=seq[i:(i+len(pattern))]
        nmismatch=0
        for j in xrange(len(s)):
            if(s[j]!=pattern[j]):
                nmismatch=nmismatch+1
                if(nmismatch>d):
                    break
        if(nmismatch<=d):
            pos[i]=s
    return(pos)
    
    
    #%% -------------------- FIND CONSENSUS -----------------------
# Función similar a findPattern, pero retorna la secuencia consenso en vez de un diccionario
def findConsensus(pattern,seq, d=3):
    matches=findPattern(pattern, seq, 3)
    return consensus(matches.values())


#%% ----------------------- K-MERS (CON MUTACIONES)------------
# Similar a countKmers, retorna para cada k-mer en seq, las veces que aparece
# en seq, bien en su forma normal o con hasta d mutaciones
def countKmersMM(seq,d=3,k=9):
     kmers={}   
     for i in xrange(len(seq)-k+1):
        kmer=seq[i:(i+k)]
        pos=findPattern(kmer,seq,d)
        kmers[kmer]=len(pos)
     return kmers
    


#%% Similar a countKmersMM, pero sólo retorna los k-mers en seq que aparecen
#en más ocasiones (en su forma normal o con hasa d mutaciones)
def mostFrequentKmers(seq,d=3,k=9):
     kmers={}   
     for i in xrange(len(seq)-k+1):
        kmer=seq[i:(i+k)]
        pos=findPattern(kmer,seq,d)
        kmers[kmer]=len(pos)
     kfinals=set()  
     for i in kmers.iterkeys():
        if(kmers[i]==max(kmers.values())):
            kfinals.add(i)
     return kfinals
    




#%% ------------------------- MUTACIONES ---------------------
# Dada una secuencia (word), un número de mutaciones puntuales (num_mismatches)
# y una cadena con los posibles componentes (letters), devuelve una variable 
# generator con todas las posibles mutaciones. 
# OJO: el resultado hay queconvertirlo a list para su porterior uso
#Código extraido de: 
#http://stackoverflow.com/questions/11679855/introducing-mutations-in-a-dna-string-in-python
def mutations(word, num_mismatches, letters="ACGT"):
    import itertools
    for locs in itertools.combinations(xrange(len(word)), num_mismatches):
        this_word = [[char] for char in word]
        for loc in locs:
            orig_char = word[loc]
            this_word[loc] = [l for l in letters if l != orig_char]
        for poss in itertools.product(*this_word):
            yield ''.join(poss)
            
            
#%% ------------------------- MUTACIONES (II) ---------------------
# Como la función anterior, pero ahora retorna una lista con las mutaciones 
#posibles con num_mismathes SNPs *o menos*
def mutationsEqualOrLess(word, num_mismatches, letters="ACGT"):
   matches=set()
   for dd in xrange(num_mismatches,-1,-1): 
       matches.update(list(mutations(word, dd, letters)))
   return matches



#%% ------------------------- MUTACIONES (III) ---------------------
# Retorna como lista todos los posibles k-mers en seq y todas sus mutaciones
#en hasta 2 SNPs
def allMutations(seq, k=9, d=2):
    kmers=set()
    for i in xrange(len(seq)-k+1):
        kmer=seq[i:i+k]
        kmers.update(mutationsEqualOrLess(kmer, d))
    return kmers




#%% ------------------------- K-MERS ESQUIVOS ----------------------
# Cuenta las ocurrencias de cada  k-mer in seq, considerando hast d mutaciones
#puntuales y teniendo en cuenta k-mers que no estén en la secuencia original
def  countKmersMMV2(seq,k=9,d=3, min=3):
     amu=list(allMutations(seq, k, d)) 
     amud=dict(zip(amu, [0]*len(amu))) #we get a dict with all possible muts.
     for key in amud.keys():     #for each possible kmer
        amud[key]=len(findPattern(key,seq,d))
     amud2={}
     for key in amud.keys():
         if(amud[key]>=min):
             amud2[key]=amud[key]
     return amud2

#%% Busca los k-mers más frecuentes en seq, incluyendo k-mers que no estén
#en seq originalmente, con hasta d mutaciones
def mostFrequentKmersV2(seq,k=9,d=3):
     amud=countKmersMMV2(seq,k,d)
     fw=[]
     mvalue=max(amud.values())
     for i in amud.keys():
         if(amud[i]==mvalue):
             fw.append(i)
     return fw



#%% -------------- K-MERS ESQUIVOS Y REVUELTOS --------------------------
# Contamos los k-mers como antes, pero ahora no sólo contando cada k-mer y sus
#hasta d mutaciones, si no también su inverso complementario y sus hasta d muts!
def  countKmersMMRC(seq,k=9,d=3, min=3):
     amu=list(allMutations(seq, k, d)) 
     amud=dict(zip(amu, [0]*len(amu))) #we get a dict with all possible muts.
     for key in amud.keys():     #for each possible kmer
        amud[key]=len(findPattern(key,seq,d))
        amud[key]+=len(findPattern(revcomp(key),seq,d))
     amud2={}
     for key in amud.keys():
         if(amud[key]>=min):
             amud2[key]=amud[key]
     return amud2


#%% Equivalente a frequentWordMM pero ahora con inversos complementarios
def frequentWordMMRC(seq,k=9,d=3):
     amud=countKmersMMRC(seq,k,d,1)
     fw=[]
     mvalue=max(amud.values())
     for i in amud.keys():
         if(amud[i]>=mvalue):
             fw.append(i)
     return fw
     
     
#%%
def consensusProfile(consensus, profile):
    cp=[]
    for i in xrange(len(consensus)):
        cp.append(profile[consensus[i]][i])
    return cp

#%% Entropía del grupo de motivos, a partir de las probabilidades
def entropy(motifs):
    import math
    s=0
    p=profile(motifs)
    for i in xrange(len(p['A'])):
        acc=0
        for k in p.keys():
            if(p[k][i]>0):
                acc+=p[k][i]*math.log(p[k][i], 2)
        s=s-acc
    return s
    
    
#%% Distancia de Hamming de un conjunto de motivos a un patrón
def d(motifs, pattern):
    dtotal=0
    for motif in motifs:
        d=len(pattern)
        for i in xrange(len(motif)-len(pattern)+1):
            dnew=len(pattern)
            chunck=motif[i:i+len(pattern)]
            for j in xrange(len(pattern)):
                if(chunck[j]==pattern[j]):
                    dnew-=1
            if(dnew<d):
                d=dnew
        dtotal+=d
    return dtotal



#%%

def greedyMotifSearch(dna, k):
    t=len(dna)
    bestMotifs= [x[:k] for x in dna]
    for i in xrange(len(dna[0])-k+1):
        mot=[]
        mot.append(dna[0][i:i+k])
        for j in xrange(1,t):
            prof=profileLaplace(mot)
            mot.append(profileMostProbableKmer(dna[j], k, prof)) #perfil más probable en la iésima cadena de Dna
        if(score(mot)<score(bestMotifs)):
            bestMotifs=mot
    return bestMotifs        
    


#%% RANDOMIZED MOTIF SEARCH
#Definimos motifs() como los motivos más probables para un determido perfil
def motifs(perfil, dna):
    return [profileMostProbableKmer(x,len(perfil['A']),perfil) for x in dna]

#%% RANDOMIZED MOTIF SEARCH
#A partir de una semilla aleatoria, va iterando sobre el mismo perfil para
#obtener cada vez motivos con mejor (más bajo) score
def randomizedMotifSearch(dna, k):
    import random
    random.seed()
    motivos=[]
    for i in xrange(len(dna)):
        ri=random.randint(0,len(dna[i])-k) 
        motivos.append(dna[i][ri:ri+k])
    mejores=motivos
    while 1:
        perfil=profileLaplace(motivos)
        motivos=motifs(perfil,dna)
        if(score(motivos)<score(mejores)):
            mejores=motivos
        else:
            return mejores

#Laplace's attempts to improve times and are actually slower
#%%
def profileLaplace1(motifs):
    import re
    import numpy as np
    l=len(motifs[0])
    freqs={'A':np.ones(l), 'C':np.ones(l), 'G':np.ones(l),'T':np.ones(l),'N':np.ones(l)}
    den=(float)(len(motifs)*2)
    for i in xrange(len(motifs)): #para cada fila
        #for k in freqs.keys():
        for k in ['A','C','G','T']:
            freqs[k][[m.start() for m in re.finditer(k, motifs[i])]]+=1
    #for k in freqs.keys():
    for k in ['A','C','G','T']:
        freqs[k]/=den
    return freqs


def profileLaplace2(motifs): #not good
    import numpy as np
    mot=np.matrix([list(x) for x in motifs],dtype="str")
    freqs={'A':[], 'C':[], 'G':[],'T':[]}
    den=(float)(len(motifs)*2)
    for j in xrange(mot[0].size): #para cada columna
        for k in freqs.keys():
            freqs[k].append(("ATTCCG".count(k)+1)/den)
            #freqs[k].append((mot[:,j].tostring().count(k)+1)/den)
    return freqs


#%%
def greedyMotifSearchV0(dna, k):
    t=len(dna)
    bestMotifs= [x[:k] for x in dna]
    #bestProbs=[score(x) for x in bestMotifs]
    for i in xrange(len(dna[0])-k+1):
        mot=[]
        mot.append(dna[0][i:i+k])
        for j in xrange(1,t):
            prof=profile(mot)
            mot.append(profileMostProbableKmer(dna[j], k, prof)) #perfil más probable en la iésima cadena de Dna
        if(score(mot)<score(bestMotifs)):
            print score(mot)
            bestMotifs=mot
    return bestMotifs
    
    
#Consensus with "-" as character
#%% ------------------- CONSENSUS ---------------------
#Retorna la secuencia de consenso de una lista de secuencias
#def consensus(seqs):
#    import operator
#    match=[]
#    for letter in xrange(len(seqs[0])):
#        freq={"A":0,"G":0,"T":0,"C":0, "-":0}
#        for s in seqs:
#            if(s[letter] in freq.keys()):
#                freq[s[letter]]+=1
#            else:
#                freq[s[letter]]=1
#        match.append(max(freq.iteritems(), key=operator.itemgetter(1))[0])
#    return reduce(lambda a,b: '{0}{1}'.format(a,b), match)  #HIC SUNT DRACONES: más info sobre las funciones lambda: http://www.secnetix.de/olli/Python/lambda_functions.hawk
#

    