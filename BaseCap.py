#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: camilla eldridge
"""

''' Trims low qc bases, primers(if found) and output trimmed sequence and it's '.qual file '''

import regex
import os
import sys

#import time
#start_time = time.time()



def reverse_comp(seq):
    """ Return reverse complement """
    result = ""
    seq=seq.lower()
    comp=seq.replace("a", "T").replace("t","A").replace("c","G").replace("g","C").lower()
    result="".join(comp[::-1])
    return result



def find_variable_seq(sequence, primer, error_n):
    """ Return variable string match allowing n errors (del &/or subst"""
    p=str("(" + primer.lower() + ")" + "{e<=" + str(error_n) + "}")
    q=regex.search(p, sequence.lower())
    return q


    
def index_phd(phd):
    """ Index phd file - adds column of base pos """
    indexed = ""
    phd=filter(None, phd.split("\n"))
    
    for x, value in enumerate(phd, 1):
        indexed=indexed + str(x) + " " + str(value) + "\n"
    
    return indexed





class Main(object): 
    
    """   trim low quality bases and primers from .Phd file, output qual and trimmed sequence """
    
    def __init__(self, phd_file, ID, primers, threshold): 
    
        
        self.phd_file = phd_file
        self.ID = ID
        self.primers = primers
        self.threshold = threshold
        
    
    
    def phd_to_sequence(self):
        
        """ Get sequence from phd_file """
        
        self.sequence = ""
        
        with open(self.phd_file) as self.phd1:
            self.phd2=self.phd1.read()
            self.phd=index_phd(self.phd2.split("BEGIN_DNA")[1].split("END_DNA")[0].rstrip())
            self.sequence=self.sequence + "".join(self.phd.split()[1::4]) 
            
        return self.sequence
    
    
    
    def trim_pos_ttuner(self): 
        
      """ Find trim positions from base caller ttuner """
      
      self.trim_pos=""
      
      for line in self.phd2.split("\n"):  
            if "TRIM:" in line:
                self.trim_pos = self.trim_pos + line
                
      self.trim_pos=self.trim_pos.split()[1:3]
      
      self.trim_start = int(self.trim_pos[0])
      self.trim_stop = int(self.trim_pos[1])
      
      
      ''' if both are <0 no trimming needed '''
     
      if ( self.trim_stop < 0 ) and ( self.trim_start < 0 ):
            self.check = 0
      else:
            self.check = 1      
    
      return self.trim_start, self.trim_stop, self.check
          
        
    
    def trim_check(self):
        
        ''' See if trimming is recommended by ttuner 
        we dont know which (if any ) primer will hit so need to keep 
        trim pos as 0 or len(sequence) even if no trimming required '''
        
        
        if self.trim_start < 0:
                self.trim_start = 0
        else:
                self.trim_start = self.trim_start
                
        if self.trim_stop < 0:
            self.trim_stop = int(len(self.sequence))
        else:
            self.trim_stop = self.trim_stop
        
        return self.trim_start, self.trim_stop
        
    
    
    def reverse_primers(self):  
        
        """ reverse complement both primers """
        
        self.primz=[]
        with open(self.primers) as primerz:
            
                primerz=primerz.read().split("\n")[1::2] #list of primers
                
                """ revcomp each primer and make new list """
                self.primz.append(primerz)
                self.primz.append([reverse_comp(prim_seq) for prim_seq in primerz])
        
        self.primz=[a for b in self.primz for a in b]
        
        return self.primz
                          
    
    
    def primer_search(self):
        
        ''' Search for each primer (forward and reverse) - expecting one hit '''
        
        for p in self.primz:
            p=str(p)
            self.search=find_variable_seq(self.sequence,p,3)
            if self.search == None:
                self.primer_pos = None
            else:
                self.primer_pos = self.search
                
                return self.primer_pos, self.search
    


    def locate_primer(self):
        ''' locate which primer hit (assuming orient not known 
        before hand) '''
        
        midpoint=int(len(self.sequence))/2
        
        ''' primer pos: e.g if the position is at the right end
        then then primer hit to end and start is suggested at 0'''
    
    
        if  self.primer_pos.span()[0] > midpoint:
            self.primer_start = 0
            self.primer_stop = int(self.primer_pos.span()[0])
                
        elif self.primer_pos.span()[0] < midpoint:
            self.primer_start=int(self.primer_pos.span()[1])
            self.primer_stop = int(len(self.sequence))
            
        return self.primer_start, self.primer_stop
        

    
    def compare_trim(self):
        
        ''' compare primer and trim positions '''
        ''' we have trim_start and trim_stop, also maybe primer_pos...'''
        
        ''' if trim pos is smaller than primer pos, trim by primer pos etc '''
        
        if self.trim_start < self.primer_start:
            self.start = int(self.primer_start)
        else:
            self.start = self.trim_start
        
        if  self.trim_stop < self.primer_stop:   
            self.stop = self.trim_stop
        else:
            self.stop = int(self.primer_stop)

        return self.start, self.stop
            
    
    
    
    def trim_phd_and_seq(self):
        ''' trim phd and sequence'''
        
        self.trimmed_seq = self.sequence[int(self.start):int(self.stop)]
        
        self.trimmed_phd=self.phd.split("\n")[int(self.start):int(self.stop)]

        return self.trimmed_seq, self.trimmed_phd
    
    

        
    def cap_the_bases(self):
        
        ''' capitalise the low scoring bases
        and get scores for qual
        '''
        
        ''' trimmed and not trimmed have a different structure 
        make sure to remove empty values in list'''
        
        self.trimmed_phd2 = [ele for ele in self.trimmed_phd if ele.strip()]
               
        self.seq_cap = ""
        self.scores = ""
        
        
        for j in self.trimmed_phd2:
            
            self.scores = self.scores + j.split()[2] + " "
            
            if int((j.split()[2])) < int(self.threshold):
                self.seq_cap = self.seq_cap + str(j.split()[1]).upper()
            else:
                self.seq_cap = self.seq_cap + str(j.split()[1]).lower()
                
        self.seq_cap = "".join(self.seq_cap)
        
    
        return self.seq_cap, self.scores

       
       
        
    def write_out(self):
        
        ''' write out trimmed and capped sequence '''
        
        self.qual_file=open(str(self.ID) + ".qual", "w")
        self.qual_file.write(">" + str(self.ID) + "\n" + self.scores)
        self.qual_file.close()
        
        self.seq_file=open(str(self.ID) + ".fasta", "w")
        self.seq_file.write(">" + str(self.ID) + "\n" + self.seq_cap)
        self.seq_file.close()
        
        return self.seq_file, self.qual_file
    
    
    
    def call_methods(self):
        
        ''' to call initial methods ''' 
        
        methods1= ["phd_to_sequence", "trim_pos_ttuner", "reverse_primers" ,"primer_search", "trim_check"]
        for method in methods1:
                getattr(self, method)()
            
        ''' if no primers found 
        miss out compare trim function'''        
        
        if self.primer_pos == None:
            self.start = self.trim_start
            self.stop = self.trim_stop
        else:
            ''' optional methods if primers found'''
            self.locate_primer()
            self.compare_trim()

        
        if self.check == 0:
            self.trimmed_seq = self.sequence
            self.trimmed_phd = self.phd.split("\n")
            
            
        else:
            self.trim_phd_and_seq()

        ''' run final methods '''
        self.cap_the_bases()
        self.write_out()
        
        
        ''' get summary info,print out and save '''
        
        if self.search == None:
            self.search_out = None 
        else:
            self.search = self.primer_pos.group()
        
        ''' if trimmed_seq variable exists, get len of trimmed sequence '''    
        if self.trimmed_seq in locals():
            
            self.if_trim = len(self.sequence)
        else:
            self.if_trim = len(self.trimmed_seq)
        
        
        
        
        self.summary=""
        
        var_list=[self.ID, self.check, self.start, self.stop, len(self.sequence), self.if_trim , self.trim_start, self.trim_stop, self.search]
        
        for variable in var_list:
            self.summary=self.summary + str(variable) + ","
            
        ''' make output in csv format, remove last comma'''
        
        return self.summary[:-1]




''' Input '''
Directory=sys.argv[1] 
Primers= sys.argv[2]
Threshold=sys.argv[3]

header_list=["ID", "trimmed(0/1)" ,"trim_start", "trim_stop", "read_len", "trimmed_read_len", "primer_hit"]

print(",".join(header_list))
           
for filename in os.listdir(Directory):
    
    ''' iterate over all phd.1 files in dir'''
    
    f = os.path.join(Directory, filename)
    ID = str(filename.split(".")[0])
    
    if filename.endswith(".phd.1"):
        Z=Main(f, ID, Primers, Threshold)
        print(Z.call_methods())



#elapsed_time = time.time() - start_time

# test 
#Directory="phd_files/"#sys.argv[1] 
#Primers= "test_primers.txt"#sys.argv[2]
#Threshold=25 #sys.argv[3]

        
        

