# BaseCap 
============== 

`BaseCap.py` is a read QC program for sanger reads, created during my PhD at Kingston University (thesis available soon here: https://eprints.kingston.ac.uk/id/eprint/53583/).  <br /> <br /> <br />



## Background
In both first generation and next generation sequencing methods (NGS) DNA bases are assigned a quality score, known as a Phred score, by base caller programs that analyse the output from sequencing machines. In Sanger sequencing Phred quality scores are assigned based on peak spacing, peak resolution and peak ratios from sequence chromatograms. 

Phred scores, also known as Q scores, represent the probability of an error in the assigned base and are defined by the equation below:

      Q = - 10 log10 Pr (observed allele ≠ true allele)

Phred quality scores above 20 represent a 1% chance a base has been called incorrectly and indicate 99% accuracy in the called base. 



<img src="" width="750" height="600">




Many read filtering pipelines use this threshold to initially remove low scoring reads from a dataset. In addition, to improve sequence quality and increase read-length, amplicons are usually sequenced bi-directionally using a forward and reverse primer, and merged to a final product.




**Problem:**  

- Trawling through 100+ sanger sequence chromatograms takes time.
- It's difficult to identify if SNPs are sequencing artefacts (low scoring bases) or if they represent true variation (alleles) after the QC process.  <br /> <br /> <br />

**Solution:**  

To solve this BaseCap performs the following:  
* Incorporates trimming recommendations from `tracetuner`.  
* Indexes quality scores whilst trimming and identifying primer sequence.
* Determines the best place to trim to minimise sequence loss.  
* Capitalises bases below a phred score threshold in the fasta sequence.  
* Outputs a `.qual` file associated with the trimmed fasta file. <br /> <br /> <br />


**Score indexing:**

<img src="https://github.com/camilla-eldridge/Basecap/blob/main/basecap.jpg" width="750" height="600">


 <br /> <br /> <br />



## Usage  

`BaseCap.py` is run on the commandline and takes 3 arguments: 
1. Directory of phd files.
2. Primer sequence file.
3. Phred score threshold.

             BaseCap.py phd_directory primer_sequence score_threshold
   

The primer sequences should be in a text file in fasta format, as shown in the example below:

    >F
    cggagcatgtaccaaaaatac
    >R
    tcgttttaacttacttcagagt



The following dependancies are required:  
    *`Python3+`  
    *`sys`  
    *`regex`  
    *`os`    

     
<br /> <br /> <br />


## Read QC Workflow  
The pipeline `basecap.sh` uses tracetuner for basecalling sanger reads, BaseCap.py for QC trimming and generation of trimmed .qual files, and finally CAP3 for merging.

Test files provided: 
- `ab1` files.  
- `phd.1` files.  
- primer sequence.  


<br /> <br /> <br />
  
## Example output: 

**Summary table**  
The summary table is in csv format and gives information such as:  

1.Was the read trimmed (0 or 1)?  
2.Trimming positions (start and stop).  
3.The length of the read before and after trimming.  
4.If the primer sequence was found in the read sequence.   

            ID,trimmed(0/1),trim_start,trim_stop,read_len,trimmed_read_len,primer_hit
            SL-CD63-12h_R_069_L18,0,0,168,168,168,0,168,None
            SL-CD63-12e_F_075_E18,1,23,331,357,308,23,342,agtatttttgggacatgactccg
            SL-CD63-12h_F_069_K18,0,0,5,5,5,0,5,None
            SL-CD63-12a_F_052_M16,1,34,355,373,321,34,355,None
            SL-CD63-12e_R_075_F18,1,19,319,320,300,19,319,None
            SL-CD63-12a_R_052_N16,1,30,232,243,202,30,232,None
            SL-CD63-12d_R_077_D18,1,35,323,324,288,35,323,None
            SL-CD63-12d_F_077_C18,1,21,329,354,308,21,343,caagtatttttggtacatgctccg

<br /> <br /> <br />

**Viewing capitalised bases in sequences**  

With low scoring bases capitalised we can now make a more informed decision about our SNPs by viewing sequences them in an alignment viewer. 
If our SNPs represent real variation they will likely be in lower case, if they are a result of sequencing errors we should see them in uppercase.

**(A)** illustrates we have alleleic variation and the low scoring base `G` is more likely to represent real variation.  
In **(B)** we see an SNP `T` within a region of low scoring bases that looks more likely to be a sequencing artefact.  

![alt text](https://github.com/camilla-eldridge/Basecap/blob/main/manual-edit.jpg)

<br /> <br /> <br />


## Citation

If you use this program in your work please cite as:  

Eldridge, C., Majoros, G., Cook, R.T., Kidd, D., Emery, A.M., Lawton, S.P. (2021). BaseCap: A read QC program for Sanger reads [Software]. Unpublished.  


If you use the pipeline in your work, please also cite the authors of these tools:  

**TraceTuner**  
G.A.Denisov, A.B.Arehart and M.D.Curtin (2004). A system and method 
for improving the accuracy of DNA sequencing and error probability estimation
through application of a mathematical model to the analysis of
electropherograms. US Patent 6681186.


**CAP3**  
Huang, X. and Madan, A. (1999) CAP3: A DNA sequence assembly program. Genome Res., 9, 868-877.



<br /> <br /> <br />

<br /> <br /> <br />


Shield: [![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg
