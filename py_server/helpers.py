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
#%%
"""
Parses a text which may contain + and * operations (no parenthesis allowed)
"""   
def convertString(text):
    s=parse(text)
    if(len(s)==1):
        return s[0]
    #solve *
    s2=[]
    for i in range(len(s)-1):
        if(s[i+1]=='*'):
            try:
                int(s[i+2])
                s2.append(s[i]*int(s[i+2]))
            except ValueError:
                s2.append(s[i+2]*int(s[i]))
        elif(s[i]!="*"):
            s2.append(s[i])
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
