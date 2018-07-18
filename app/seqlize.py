# -*- coding: utf-8 -*-
"""
Command line wrapper for BWT sequence searches

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
import os
import time
import string
import pickle #i tried cPickle but is way slower!

import sys

#Our methods
import suffixSearch as ss
import motifSearch as ms
import annotations as ann
import analysis as ana

import helpers
import re

n=3
w=30

if len(sys.argv)==1 or (sys.argv[1]!="process" and sys.argv[1]!="search"):
    print "Usage:\tpython seqlize.py [process|search]" 
    print 
    print " *process*:\tseqlize process [options] inputFile -o processedFile"    
    print "\tinputFile must be a .wig or .bw file"    
    print "\toutputFile is an internal binary indexing you must keep for searches"
    print "\toptions:"
    print "\t\t-n=X number of percentile bins data are being discretized in (default 3)."
    print "\t\t-w=X discretization window (default 30)."
    print "\t\t\tExample: if w=2 and n=3 the wig sequence 1 4 8 9 9 7 6 5 3 0, "
    print "\t\t\tit will be first averaged by w=2 into    2.5 8.5 8.0 5.5 1.5 "
    print "\t\t\tand then assigned to three bins (a,b,c):  a   c   c   b   a"
    print "\t\t\tNote: higher n and lower w increase detail but reduce speed"
    print 
    print " *search*:\tseqlize search [options] inputFile pattern -o searchResults"    
    print "\tinputFile must be a processedFile obtained by seqlize process"
    print "\tsearchResults is the file name were the search will be stored"
    print "\tpattern must be a letter sequence in the range fo the binning detaild by n in seqlize process"
    print "\t\t operators '*' and '+' are allowed"
    print "\t\t Examples: if n=3 and w=30 in seqlize process"
    print "\t\t\t abcba will search for a peak of 150bps"
    print "\t\t\t abcba*3 will search for three consecutive peaks of 150bps"
    print "\t\t\t a*10+abcba will search for a low level 300bps area followed by a 150 bps peak"
    print "\toptions:"
    print "\t\t-d=X number of allowed variations over the pattern (default 0)."
    print "\t\t\tExample: d=1 in pattern search abcba will report abbba or abcca as matches"
    print "\t\t\tNote: d must be in [0,3]. d=2 or 3 significantly reduce time performance"
    print "\t\t-Note: d must be in [0,3]. d=2 or 3 significantly reduce time performance"

if(sys.argv[1]=="process"):
    filenames=[]
    outfile=[]
    
    numBins=3
    windowSize=30

    #These are visualization or client/server parameters that we can leave to default in local console
    maxSize=100000
    stdev=3
    track="None"
    dataName="None"
    interpolation="mean"
    basePath="."
    
    #These are parameters yet to implement, by now defaults
    organism="Saccharomyces cerevisiae"
    
    for i in range(2,len(sys.argv)):
        if(sys.argv[i].startswith("-")):
            option=sys.argv[i][1:].split("=")
            if(option[0]=="n"):
                numBins=int(option[1])
            if(option[0]=="w"):
                windowSize=int(option[1])
        else: 
            break
    while(sys.argv[i].startswith("-")==False):
        filenames.append(sys.argv[i])
        i=i+1
    print sys.argv[i]
    if(sys.argv[i]=="-o"):
        outfile.append(sys.argv[i+1]) #TODO: multiple files
    
    print(outfile)
    print(filenames)
    print(numBins)
    print(windowSize)
               
    ana.batchPreprocessLocal(filenames, outfile, dataName, windowSize, numBins, maxSize, stdev, track, organism, interpolation, basePath)
    print("Data successfuly preprocessed")
    
if(sys.argv[1]=="search"):
    data={}
    d=0
    geo="none"
    for i in range(2,len(sys.argv)):
      if(sys.argv[i].startswith("-")):
        option=sys.argv[i][1:].split("=")
        if(option[0]=="d"):
            d=int(option[1])
            if(d>3):
                print("ERROR: d must be <=3")
                exit
            if(d>1):
                print("WARNING: d>1 might lead to slow searches, be patient")
        if(option[0]=="g"):
            geo=str(option[1])
      else: 
        break

    ifile=sys.argv[i]
    pattern=sys.argv[i+1]
    

    print("Searching",ifile,"for pattern", pattern,"(",d,") restricted to",geo)
    f=open(os.path.join(ifile))
    pdata=pickle.load(f)
    fnames=[]
    print(pdata.keys())
    for k in pdata.keys():
        if(k.find(".")>-1):
            fnames.append(k)
    if(len(fnames)==0):
        print("ERROR: No wig associated to pic data")
        exit
    print fnames
    data=ana.batchCompile(pdata, fnames, fnames[0],"")
    print"DATA LOADED"
    
    ret=ana.searchLocal(data=data, pattern=pattern, d=d, geo=geo, intersect="soft", softMutations="false")
    print"SEARCH FINISHED"
    print ret

   
    
    
    