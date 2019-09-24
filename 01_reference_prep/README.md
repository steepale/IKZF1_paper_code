# Custom Alteration of NCBI *Gallus gallus* 5.0 Reference Genome

The Galgal5 reference genome contains a large number of micro-chromosomes and unresolved contigs. Early in our experiment, we chose to focus on only the major chromosomes; therefore, all genomic analyses rely on mapping to an altered version of the NCBI Glagal5 genome. Specifically, we reduced the Galgal5 reference to major contigs (1-28,30-33,W,Z,LGE64, & MT).  

The Galgal5 reference genome was obtained from NCBI on Febuary 22,2016, chromosomes were renamed, and minor contigs were removed.  

The resulting Galgal5 genome was indexed with both BWA (v0.7.12-r1044) and Samtools (v1.3-20-gd49c73b using htslib 1.3-29-g091c89c).  
The sequence dictionary was generated with Picard tools (v1.141).  

Analysis written and performed by [Hongen Xu](https://github.com/hongenxu).  