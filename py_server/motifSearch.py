# -*- coding: utf-8 -*-
"""
SEARCH METHOD FOR GENOME BROWSER
@author: rodri
"""

#%% --------------------------- K-MERS --------------------
# Contar los k-mers que aparecen n o más veces en seq
def countKmers(seq,n,k):
    kmers={}    
    kfinals={}  
    for i in range(len(seq)-k+1):
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
    for i in range(0, len(seq)-len(pattern)+1):
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
    for i in range(0, len(seq)-len(pattern)+1):
        s=seq[i:(i+len(pattern))]
        nmismatch=0
        for j in range(len(s)):
            if(s[j]!=pattern[j]):
                nmismatch=nmismatch+1
                if(nmismatch>d):
                    break
        if(nmismatch<=d):
            pos[i]=s
    return(pos)


#%% ------------------- CONSENSUS ---------------------
#Retorna la secuencia de consenso de una lista de secuencias
def consensus(seqs):
    import operator
    match=[]
    for letter in range(len(seqs[0])):
        freq=dict(A=0,G=0,T=0,C=0)
        for s in seqs:
            if(s[letter] in freq.keys()):
                freq[s[letter]]+=1
            else:
                freq[s[letter]]=1
        match.append(max(freq.iteritems(), key=operator.itemgetter(1))[0])
    return reduce(lambda a,b: '{0}{1}'.format(a,b), match)  #HIC SUNT DRACONES: más info sobre las funciones lambda: http://www.secnetix.de/olli/Python/lambda_functions.hawk




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
     for i in range(len(seq)-k+1):
        kmer=seq[i:(i+k)]
        pos=findPattern(kmer,seq,d)
        kmers[kmer]=len(pos)
     return kmers
    


#%% Similar a countKmersMM, pero sólo retorna los k-mers en seq que aparecen
#en más ocasiones (en su forma normal o con hasa d mutaciones)
def mostFrequentKmers(seq,d=3,k=9):
     kmers={}   
     for i in range(len(seq)-k+1):
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
    for locs in itertools.combinations(range(len(word)), num_mismatches):
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
   for dd in range(num_mismatches,-1,-1): 
       matches.update(list(mutations(word, dd, letters)))
   return matches



#%% ------------------------- MUTACIONES (III) ---------------------
# Retorna como lista todos los posibles k-mers en seq y todas sus mutaciones
#en hasta 2 SNPs
def allMutations(seq, k=9, d=2):
    kmers=set()
    for i in range(len(seq)-k+1):
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
     

# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 10:37:18 2015

@author: rodri
"""

#import sys
#sys.path.append("/Users/rodri/Documents/docencia/bioinformatica/sesionesPractica")
#import sesion2methods
#s2.consensus()=consensus()
#s2.allMutations()=allMutations()


#%% Perfil de los motivos (frecuencias o probabiliades)
def profile(motifs):
    freqs={'A':[], 'C':[], 'G':[],'T':[]}
    for j in range(len(motifs[0])): #para cada columna
        count={'A':0, 'C':0, 'T':0, 'G':0}
        for m in motifs:    #para cada motivo
            count[m[j]]+=1
        for k in count.keys():
            freqs[k].append((float)(count[k])/len(motifs))
    return freqs


#%% Entropía del grupo de motivos, a partir de las probabilidades
def entropy(motifs):
    import math
    s=0
    p=profile(motifs)
    for i in range(len(p['A'])):
        acc=0
        for k in p.keys():
            if(p[k][i]>0):
                acc+=p[k][i]*math.log(p[k][i], 2)
        s=s-acc
    return s

#%% ------------------- CONSENSUS ---------------------
#Retorna la secuencia de consenso de una lista de secuencias
def consensus(seqs):
    import operator
    match=[]
    for letter in range(len(seqs[0])):
        freq=dict(A=0,G=0,T=0,C=0)
        for s in seqs:
            if(s[letter] in freq.keys()):
                freq[s[letter]]+=1
            else:
                freq[s[letter]]=1
        match.append(max(freq.iteritems(), key=operator.itemgetter(1))[0])
    return reduce(lambda a,b: '{0}{1}'.format(a,b), match)  #HIC SUNT DRACONES: más info sobre las funciones lambda: http://www.secnetix.de/olli/Python/lambda_functions.hawk



#%% Puntuación del grupo de motivos, a parrde las probabilidades
def score(motifs):
    s=0
    con=consensus(motifs)
    for m in motifs:
        for j in range(len(m)):
            if(con[j]!=m[j]):
                s+=1
    return s


#%% Distancia de Hamming de un conjunto de motivos a un patrón
def d(motifs, pattern):
    dtotal=0
    for motif in motifs:
        d=len(pattern)
        for i in range(len(motif)-len(pattern)+1):
            dnew=len(pattern)
            chunck=motif[i:i+len(pattern)]
            for j in range(len(pattern)):
                if(chunck[j]==pattern[j]):
                    dnew-=1
            if(dnew<d):
                d=dnew
        dtotal+=d
    return dtotal

#%% MEDIAN STRING
#def medianString(motifs, k):
#    mind=k*len(motifs)
#    bestPattern=""
#    patterns=list(sesion2methods.allMutations("A"*k, k,k))
#    for p in patterns:
#        dp=d(motifs,p)
#        if(dp<mind):
#            bestPattern=p
#            mind=dp
#    return [bestPattern, mind]
    

#%% PROBABILITY
def pr(m,p):
    prod=p[m[0]][0]
    for i in range(1, len(m)):
        prod*=p[m[i]][i]
    return prod


#%% MOST PROBABLE PROFILE
def profileMostProbableKmer(text, k, perfil):
    bestPr=-1
    bestSeq=text[:k]
    for i in range(len(text)-k+1): #para cada posible k-mer en text
        kmer=text[i:(i+k)]
        newPr=pr(kmer, perfil)
        if(newPr>bestPr):
            bestPr=newPr
            bestSeq=kmer
    return bestSeq

#%%
def greedyMotifSearchV0(dna, k):
    t=len(dna)
    bestMotifs= [x[:k] for x in dna]
    #bestProbs=[score(x) for x in bestMotifs]
    for i in range(len(dna[0])-k+1):
        mot=[]
        mot.append(dna[0][i:i+k])
        for j in range(1,t):
            prof=profile(mot)
            mot.append(profileMostProbableKmer(dna[j], k, prof)) #perfil más probable en la iésima cadena de Dna
        if(score(mot)<score(bestMotifs)):
            print score(mot)
            bestMotifs=mot
    return bestMotifs
    

#%% Profile Laplace
def profileLaplace(motifs):
    freqs={'A':[], 'C':[], 'G':[],'T':[]}
    den=len(motifs)*2
    for j in range(len(motifs[0])): #para cada columna
        count={'A':1, 'C':1, 'T':1, 'G':1}
        for m in motifs:    #para cada motivo
            count[m[j]]+=1
        for k in count.keys():
            freqs[k].append((float)(count[k])/den)
    return freqs


#%%
def greedyMotifSearch(dna, k):
    t=len(dna)
    bestMotifs= [x[:k] for x in dna]
    for i in range(len(dna[0])-k+1):
        mot=[]
        mot.append(dna[0][i:i+k])
        for j in range(1,t):
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
    for i in range(len(dna)):
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


#%% GIBBS SAMPLER
#N is the number of iterations we want for the optimization
def gibbsSampler(dna, k, N):
    t=len(dna)
    #Random selection of initial motifs
    import random
    random.seed()
    motivos=[]
    for i in range(len(dna)):
        ri=random.randint(0,len(dna[i])-k) 
        motivos.append(dna[i][ri:ri+k])
    bestMotifs= [x[:k] for x in dna]
    sb=score(bestMotifs)
    #iteration
    for j in range(N):
        i=random.randint(0,t-1)
        del motivos[i]
        prof=profileLaplace(motivos)
        motivos.insert(i, profileMostProbableKmer(dna[i], k, prof)) #perfil más probable en la iésima cadena de Dna
        s=score(motivos)
        if(s<sb):
            bestMotifs=motivos
            sb=s            
    return bestMotifs        



    

