# Notes on setup  
___

This file aims to give an overview of the tools and files used for the study **Reliable biodiversity metrics from co-occurence based post-clustering curation of amplicon data.**    
All steps/processes for this study can be carried out on the same computer/platform. But, in practise all analyses were carried out on a linux server setup with 64 processors (AMD Opteron(tm) 6380), except R-scripts, which were run on a MacBook Pro (2.6 GHz Intel Core i7, 16 GB 1600 MHz DDR3).
All analyses were carried out in one directory (analyses) and sub-directories of this.

## Bioinformatic tools
### CLI tools were used for this study  

 * VSEARCH v.2.02 (or later) (https://github.com/torognes/vsearch) - used for clustering in the pure VSEARCH pileine as well as the combined DADA2+VSEARCH approach.
 * Cutadapt v 1.10 (https://cutadapt.readthedocs.io/en/stable/)  - used for demultiplexing (assigning reads to samples) using two different scripts, one for VSEARCH, CROP, SWARM, and another for DADA2.
 * Swarm v 2.19 (https://github.com/frederic-mahe/swarm)  - used for single linkage clustering in the SWARM approach.
 * CROP v 1.33 (https://github.com/tingchenlab/CROP)  - used for clustering in the unsupervised bayesian clustering approach.
 * blastn v2.4.0+ (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) - used for assigning taxonomy to sequences, and for making the match lists required for running LULU.  also requires installation of NCBI's taxbd - see https://www.ncbi.nlm.nih.gov/books/NBK279690/?report=reader#CmdLineAppsManual).
 * dbotu3 (https://github.com/swo/dbotu3) - used for distribution based clustering (and in this study also tested as an alternative tool for post-clustering curation).  

### R-packages used for this study  
To replicate the analyses the following packages (and their dependencies) need to be installed.  

 * LULU 0.1.0 - (https://github.com/tobiasgf/lulu) - the algorithm developed for this study. Used to curate the OTU tables produced by the initial clustering methods.
 * taxize v0.8.4 (https://cran.r-project.org/web/packages/taxize/index.html) - used in the custom script for assigning taxonomy to the sequences.
 * rentrez v 1.1.0 - (https://cran.r-project.org/web/packages/rentrez/index.html) - used for searching GenBank.
 * DADA2 v 1.3 (http://benjjneb.github.io/dada2/dada-installation.html) - used for the pure dada2 approach as well as the DADA2 approach with subsequent clustering.
 * stringr 1.2.0 (https://cran.r-project.org/web/packages/stringr/) -  manipulation strings for the analyses.
 * dplyr 0.5.0 - (https://cran.r-project.org/web/packages/dplyr/) -  manipulation of data.frames for analyses.
 * tidyr 0.6.1 - (https://cran.r-project.org/web/packages/tidyr/) -  manipulation of data.frames for analyses.
 * ggplot2 2.2.1 - (https://cran.r-project.org/web/packages/ggplot2/) -  producing the plots for visualization.
 * ggpmisc 0.2.14 - (https://cran.r-project.org/web/packages/ggpmisc/) -  putting regression formulas on the plots.
 * cowplot 0.7.0 - (https://cran.r-project.org/web/packages/cowplot/) - make composite plot (figure1).
 * devtools 1.12.0 - (https://cran.r-project.org/web/packages/devtools/) - used for installing LULU from Github.com.  


### Provided scripts  
Command line scripts  
A number of scripts are provided with this manuscript. Place these in the /bin directory and make them executable with "chmod 755 SCRIPTNAME"" or place the scripts directly in the directory/directories where they should be executed (i.e. the analyses directory). Most of them are simple shell scripts. Their context and use is described in the R markdown (A-N) files documenting the workflow. 

 * Alfa_merge_all_pairs.sh - VSEARCH based script for merging R1/R2 paired end files (raw data).
 * Alfa_demultiplex_universal.sh - demultiplexing the merged paired end fastq files for the VSEARCH, CROP and SWARM approaches.
 * Alfa_demultiplex_for_DADA2.sh -  demultiplexing the raw data (not-merged R1/R2 files!) for the DADA2 approaches.
 * Alfa_concatenate_and_dereplicate_fasta.sh - used to collapse sample triplicates for the VSEARCH, CROP and SWARM approaches.
 * Alfa_concatenate_single_samples_fastq_DADA2.sh- used to collapse sample triplicates for the DADA2 approaches.
 * Alfa_CROP95.sh - the CROP approach at approx 95% clustering level.
 * Alfa_CROP97.sh - the CROP approach at approx 97% clustering level.
 * Alfa_CROP98.sh - the CROP approach at approx 98% clustering level.
 * Alfa_DADA2_vsearch.sh - the subsequent clustering of the initial DADA2 results.
 * Alfa_swarm.sh - the SWARM approach.
 * Alfa_vsearch.sh  - the VSEARCH approach.
 * Alfa_vsearch_dbotu.sh - used to make a 0% clustering table for the dbotu3 algorithm.
 * Alfa_matchlists.sh - the script that makes the match lists required to run the LULU curation.
 * rename.pl  - a perl script used to rename fasta-headers of files.
 * buildOTUtable_simple.sh - used to guide the building of OTU table from the SWARM analyses.
 * OTU_contingency_table_simple.py - python script building the OTU tables from the SWARM analyses.
 * uc2otutab.py  (see http://drive5.com/python/summary.html) - a python script used for making the OTU tables in several pipelines/approaches.

### R-markdown files  
This manuscript includes 9 R markdown files (including this one) documenting the analyses.  

 * A_Preparation_of_sequences.Rmd
 * B_clustering_with_VSEARCH_SWARM_CROP.Rmd
 * C_Processing_with_DADA2.Rmd
 * D_Taxonomic_filtering.Rmd 
 * E_LULU_processing.Rmd
 * F_Calculating_table_statistics.Rmd
 * G_Producing_Figure1.Rmd
 * H_Producing_supplementary_figures_1_6.Rmd
 * I_GenBank_coverage.Rmd
 * J_Dbotu3_onestep_clustering.Rmd
 * K_Benchmarking_dbotu3_curation_and_singleton_culling.Rmd
 * L_Producing_supplementary_figures_7_12.Rmd
 * M_Producing_supplementary_tables.Rmd
 * N_Producing_supplementary_figures_13.Rmd
 * X_setup_notes.Rmd - (this file)  

## Data/Files
A number of files and data provided with this manuscript are necessary for the processing (they need to be placed in the analyses directory):
Files used by Alfa_demultiplex_universal.sh script.
batchfile.list  - file with library information for demultiplexing the R1/R2 fastq paired end files. Using the files with info on primer-tag combinations used for the different samples in the different libraries.

* tags_R1A.list
* tags_R1B.list
* tags_R2A.list
* tags_R2B.list
* tags_R3A.list
* tags_R3B.list

Files used by Alfa_demultiplex_for_DADA2.sh script.
batchfileDADA2.list  - file with library information for demultiplexing the R1/R2 fastq paired end files. Using the files with info on primer-tag combinations used for the different samples in the different libraries. (Different from above, because DADA2 requires unmerged demultiplexing of R1/R2 files).

* DADA2_tags_R1A.list
* DADA2_tags_R1B.list
* DADA2_tags_R2A.list
* DADA2_tags_R2B.list
* DADA2_tags_R3A.list
* DADA2_tags_R3B.list

Table_plants_2014_cleaned.txt - Species/site matrix for the plant survey data.
Raw MiSeq data (not included in supplementary material). (accessible from http://datadryad.org/ accession no XXX).

* ESW-QPKH-ITS4-ITS2-R1A_S1_L001_R1_001.fastq.gz
* ESW-QPKH-ITS4-ITS2-R1A_S1_L001_R2_001.fastq.gz
* ESW-QPKH-ITS4-ITS2-R1B_S2_L001_R1_001.fastq.gz
* ESW-QPKH-ITS4-ITS2-R1B_S2_L001_R2_001.fastq.gz
* ESW-QPKH-ITS4-ITS2-R2A_S3_L001_R1_001.fastq.gz
* ESW-QPKH-ITS4-ITS2-R2A_S3_L001_R2_001.fastq.gz
* ESW-QPKH-ITS4-ITS2-R2B_S4_L001_R1_001.fastq.gz
* ESW-QPKH-ITS4-ITS2-R2B_S4_L001_R2_001.fastq.gz
* ESW-QPKH-ITS4-ITS2-R3A_S5_L001_R1_001.fastq.gz
* ESW-QPKH-ITS4-ITS2-R3A_S5_L001_R2_001.fastq.gz
* ESW-QPKH-ITS4-ITS2-R3B_S6_L001_R1_001.fastq.gz
* ESW-QPKH-ITS4-ITS2-R3B_S6_L001_R2_001.fastq.gz
