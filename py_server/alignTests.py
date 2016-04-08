# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 18:26:58 2016

@au#%%thor: rodri
"""

#
#    
# seqs={"100": "AGCCCATGAC", "1345": "GTCCGACTGG"}   
##%%
#from Bio.Alphabet import generic_dna
#from Bio.Seq import Seq
#from Bio.SeqRecord import SeqRecord
#
#seqs={
#    "100": "AGCCCATGAC", 
#    "1345": "GTCCGACTGG"
#  }
#sr=[]
#for k in seqs.keys():
#    sr.append(SeqRecord(Seq(seqs[k], generic_dna), id=k))
#
#from Bio import SeqIO
#output=open("unaligned.fasta", "w")
#SeqIO.write(sr,output, "fasta")
#output.close()
##%%
##TODO: check a way of getting the PATH properly!!
#import subprocess
##import shlex
##cline="t_coffee -output clustalw -infile unaligned.fasta -outfile aligned.aln"
##args=shlex.split(cline)
#p=subprocess.Popen(["/usr/bin/env", "t_coffee","-output "+method+" -infile unaligned.fasta -outfile aligned.aln"], bufsize=-1, cwd=u'/Users/rodri/WebstormProjects/seqview/py_server')
##p=subprocess.Popen(["/usr/bin/env", "/Users/rodri/tcoffee/Version_11.00.8cbe486/bin/t_coffee","-output clustalw -infile unaligned.fasta -outfile aligned.aln"], bufsize=-1, cwd=u'/Users/rodri/WebstormProjects/seqview/py_server')
#print p.returncode
##subprocess.check_output("/Users/rodri/tcoffee/Version_11.00.8cbe486/bin/t_coffee -output clustalw -infile unaligned.fasta -outfile aligned.aln", shell=True, bufsize=-1)
#
##%%
#import re
#lines=open("aligned.aln").readlines()
#aln={}
#for i in range(2,len(seqs)+2):
#    lines[i]=lines[i].replace("\n", "")
#    lines[i]=re.sub("[ ]+", " ", lines[i])
#    lines[i]=lines[i].split(" ")
#    aln[lines[i][0]]=lines[i][1]
#print aln