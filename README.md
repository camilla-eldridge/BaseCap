**BaseCap**
==============
BaseCap.py is a read QC program for sanger reads, created during my PhD at Kingston University.

**Problems: **
- Trawling through 100+ sanger sequence chromatograms takes time.
- It's difficult to identify if SNPs are sequencing artefacts (low scoring bases) or if they represent true variation (alleles) after the QC process.

**Solution:**

BaseCap.py does the following:
1. Incorporates trimming recommendations from tracetuner.
2. Indexes quality scores throughout trimming and primer identification.
3. Capitalises bases below a score threshold in the fasta sequence.
4. Outputs a .qual file associated with the trimmed fasta file.

            Usage: BaseCap.py phd_directory primer_sequence score_threshold

**Usual Workflow**
I used tracetuner for basecalling sanger reads, BaseCap.py for QC trimming and generation of trimmed .qual files, then used CAP3 for merging.

test files provided:
- ab1 files
- phd.1 files
- primer sequence

  
Example output:
**Summary csv table**

                ID,trimmed(0/1),trim_start,trim_stop,read_len,trimmed_read_len,primer_hit
                SL-CD63-12h_R_069_L18,0,0,168,168,168,0,168,None
                SL-CD63-12e_F_075_E18,1,23,331,357,308,23,342,agtatttttgggacatgactccg
                SL-CD63-12h_F_069_K18,0,0,5,5,5,0,5,None
                SL-CD63-12a_F_052_M16,1,34,355,373,321,34,355,None
                SL-CD63-12e_R_075_F18,1,19,319,320,300,19,319,None
                SL-CD63-12a_R_052_N16,1,30,232,243,202,30,232,None
                SL-CD63-12d_R_077_D18,1,35,323,324,288,35,323,None
                SL-CD63-12d_F_077_C18,1,21,329,354,308,21,343,caagtatttttggtacatgctccg
**

Then view sequences in an alignment viewer**





<br /> <br /> <br />

<br /> <br /> <br />


Shield: [![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg
