# fusBreakpoint

A package to simplify the extraction of fusion junction sequences from RNAseq BAM files, using the Biostrings Bioconductor framework. 


## Installation

```
devtools::install_github("jennylsmith/fusBreakpoint")
```

## Usage:

0. You will need the following files as input: 

A. *Sample_ID*
- simply a unique string

B. 

*breakpoint*
a character vector of the breakpoint to provide the coordinates/chromosome to pull the reads from BAM file for the sample in "Sample_ID".
It must be in the format of, "chromosome:basepair_pos|chromosome:basepair_pos" such as "16:4384802|16:88945678". 
It needs to match the same style of the reference genome used to align the BAM file (eg some use "chr16" or "16" in the BAM/SAM). 
You can just check with the command line tool `samtools view my_sample.bam` and look at the first read record.

OR


*fusion_data (optional)*
a data.frame with the colnames, "Sample" (must contain Sample_ID), "Fusion.Category", "breakpoint.TA" (delimiter is a pipe), and 
"Breakpoints.STAR" (delimiter is a pipe). 


C. *file_manifest*
a dataframe with the column names,  "Sample" and "filepath". Must have one bam file per patient ID.

D. *Seqlens*
A named vector with the chromosome lengths of the reference genome used to align the BAM file. 

E. *fusion_theoretical*
a dataframe with the column names, "Name" and "Sequence". Name refers to a unique identifier for each sequence to search (row).

F. *cicero.contig (optional)*
a dataframe with the column names, "Sample"," Fusion.Category", and "contig". 
The dataframe is a very simplified subset of the raw [CICERO](https://github.com/stjude/CICERO) fusion algorithm output. 

G. *TA.contigs (optional)*
a dataframe with the column names, "Sample"," Fusion.Category", and "contig"
The dataframe is a very simplified subset of the raw [TransAbyss](https://github.com/bcgsc/transabyss) fusion algorithm output. 

TO DO: add STAR-fusion contig support as well. 


1. For one sample: 

```
Sample <- "TARGET.20.PAUSXS.09A.01R" 

rnaseq_evidence <- examine_breakpoint_evidence(Sample_ID=Sample,
                                           breakpoint="16:4384802|16:88945678",
                                           file_manifest=bam_files_df,
                                           fusion_theoretical_pdict=pdict_obj,
                                           Seqlens=chromLengths,
                                           cicero.contig=cicero_contigs_df,
                                           TA.contigs=TA_contigs_df)

```


2. for multiple samples: 

The best approach is parallelization since each sample and fusion are treated independently. The example below is part of the vignette "Multiple Samples with RSlurm and Future". I would recommend RSlurm as the sequence search is highly computationally intensive and `FURRR` package can only parallelize so much on a local machine, and it takes anywhere from > 1 hr to > 24hrs to run for a set of ~1,500 RNAsew BAM files (depends on sequencing depth, length of chromosome(s), etc). 


```
outdir <- file.path(SCRATCH,'Fusion_Breakpoints')
sopt <- list(nodes='1', 'cpus-per-task'='2',
             'partition'='campus-new', 'mem'='40G',
             'time' = '48:00:00', 'mail-type'='FAIL,END') 

#define samples as a list and the pdict (pattern dictionary) of sequences to search each BAM file.
samps <- as.list(pull(bam_file_df, Sample))
pdict_obj <- create_custom_pdict(theoretical_seqs_df = theoretical_seqs_df)



jobname <- paste("CBFB_MYH11","polyA_RNAseq_mut20",sep="_")
message(paste0("submitting: ",jobname))

sjob <- rslurm::slurm_map(samps,
                      examine_breakpoint_evidence,
                      breakpoint = "16:4384802|16:88945678",
                      file_manifest = bam_file_df,
                      Seqlens = chromLengths,
                      fusion_theoretical_pdict = pdict_obj,
                      cicero.contig = cicero_contigs_df,
                      TA.contigs = TA_contigs_df,
                      slurm_options=sopt,
                      jobname=jobname,
                      submit = T)

```

[![Travis build status](https://travis-ci.com/jennylsmith/fusBreakpoint.svg?branch=master)](https://travis-ci.com/jennylsmith/fusBreakpoint)
