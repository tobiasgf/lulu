# LULU
___
A r-package for distribution based post clustering curation of amplicon data.

The purpose of **LULU** is to reduce the number of erroneous OTUs in OTU tables to achieve more realistic biodiversity metrics. By evaluating the co-occurence patterns of OTUs among samples **LULU** identifies OTUs that consistently satisfy some user selected criteria for being errors of more abundant OTUs and merges these. It has been shown that curation with **LULU** consistently result in more realistic diversity metrics. The required input of **LULU** is an OTU table and a corresponding matchlist with all the internal matches of OTUs.   

## Requirements/pre-requisites
___
To be able to run the LULU algorithm the following things are needed:  
1. R (or R-studio) - LULU was developed in R version 3.3.2.  
2. LULU r-package (see below for installation).  
3. An OTU table - produced with any algorithm (see below for format etc).  
4. OTU sequences - a file with representavei sequence of each OTU (see below for format etc.).    
5. Access to a tool for making a match list (see below, e.g. BLASTn or VSEARCH).  

## Practical workflow
___
These are the steps carried out by the user to produce a curated OTU table:  
1. Produce an OTU table and associated file with representative sequences.    
2. Produce a match list with **BLASTn** (or another pairwise sequence dissimilarity algorithm).  
3. Run the **LULU** curation in **R** (by supplying the function with OTU table and match list).  
4. Use the curated OTU table for further analyses.  

#### 1. Produce an OTU table and associated file with representative sequences.    
The initial clustering can be done in **UCLUST, VSEARCH, SWARM, QIIME, DADA2**, etc. The requirements of **LULU** are that the OTU table has samples as columns and OTUs as rows, and that it has unique OTU id's as row names. A file with the representative sequences of each OTU is required to make the match list. The OTU id's has to match those of the OTU table.  

*An example of an OTU table adequate for LULU*  

OTUid|Sample1|Sample2|Sample3|Sample4|...  
--- | --- | --- | --- | ---  | ---
OTU1|11|204|100|299|...   
OTU2|3|2201|100|388|...   
OTU3|0|20|130|10|...   
OTU4|147|0|0|9|...  
...|...|...|...|...

*An example of a fasta file with OTU sequences*  
>OTU1  
AGCGTGGTGSA...  
>OTU2  
GGCGTATGCATGGTA...  
>OTU2  
ATGGTAGGCGTATGC...  
>OTU4  
GCGATGCGAT...  
...  

#### 2. Produce a match list  
The match list can in practice be produced with any tool for pair wise matching of sequences. **BLASTn** is an effective tool for this, and the one that was used in the validation of the **LULU** algorithm. The only requirements is that the match list has three columns with pair wise similarity scores for the OTUs. The first column contains the id of the query OTU, the second column contains the id of the matching OTU, and the third column contains the similarity score (%) of the two OTUs.  

BlastN - a match list can be produced with **blastN** with these commands:  
```
#First produce a blastdatabase with the OTUs
makeblastdb -in OTU_sequences.fasta -parse_seqids -dbtype nucl
# Then blast the OTUs against the database
blastn -db OTU_sequences.fasta -outfmt '6 qseqid sseqid pident' -out match_list.txt -qcov_hsp_perc 80 -perc_identity 84 -query OTU_sequences.fasta
```

Alternatively, a matchlist can also be produced with **VSEARCH** with this command
```
vsearch --usearch_global OTU_sequences.fasta --db OTU_sequences.fasta --self --id .84 --iddef 1 --userout match_list.txt -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10
```

No matter which algorithm is used for internal OTU matching, the resulting match list, should look like this:
```
OTU1  OTU2  0.90
OTU1  OTU3  0.88
OTU3  OTU3  1
OTU3  OTU4  0.82
OTU4  OTU9  0.86
...
```

#### 3. Run the LULU curation in R
With the OTU table and the corresponding match list at hand, we can carry out the curation with **LULU** in **R**.
```
> curated_result <- lulu(otutable_name, matchlist_name)
```
The `lulu` function will return a curated OTU table and some statistics and information on the curation in the object `curated_result`. The curated OTU table can be accessed like this:
```
curated_result$curated_table
```
This table can the be used for further biodiversity analyses, etc.  
For more advanced uses, see the help-file for the function.  
```
>?lulu  
```

## Installation
___
The lulu package can be installed in R (RStudio) using devtools, by typing these commands in R
```
> library(devtools)
> install_github("tobiasgf/lulu")  
```

## The algorithm
___
This is the processing flow employed by the r-function:  
1. Sort the OTUs by decreasing occurrence, secondarily by total read count.  
2. Select one OTU at a time from top to bottom and treat them as a 'potential daughter' (erroneous OTU).  
  - The first OTU will always be accepted as a valid OTU as it is the most widely occurring OTU with the highest read abundance.  
3. Select all OTUs in the match list matching the potential daughter within the dissimilarity threshold selected, and designate these as 'hits'.  
  - The dissimilarity threshold can be set to anything, but if the matchlist was produced with a minimum match %, a lower threshold will not have an effect.  
4. Select OTUs from the 'hits' that have an occurrence at the level of the 'potential daughter' or higher, and designate these as 'potential parents'.  
5. Test all potential parents (one by one, from top to bottom) to see if the distribution of the 'potential daughter' among samples can be explained by co-occurence (satisfying the co-occurence threshold, and the abundance threshold) by the 'potential parent', and if so, flag the 'potential daughter as an error of that 'parent'.  
  - In this step all samples (columns of the OTU table) where the 'potential daughter' occurs are selected also for the 'potential parent'.  
  - If the number of samples where the 'potential parent' has a positive presence is below the minimum_relative_cooccurence (default 95%, meaning that 1 in 20 samples are allowed to have no parent presence), the 'potential parent' is rejected.  
  - If not, the abundance ratio between 'potential daughter' and 'potential parent' is calculated for all samples, and tested against the minimum_ratio.  
  - The minimum ratio  can be set to any number, and the threshold (minimum_ratio_type) can be set to be evaluated as the minimum (min) observed ratio or the average (avg) observed ratio.  
  - If the potential parent satisfies the ratio threshold, the 'potential daughter' is flagged as an error of this OTU.  
6. If no potential parent satisfies the criteria above, the 'potential daughter' is flagged as a valid OTU.  
7. repeat from 2 with the next OTU until the full table has bee processed.  
  - After this step all OTUs will either be flagged as an error of a more abundant OTU, or as a valid OTU.  
8. Process the OTU table from bottom to top, merging flagged errors with their parents to produce a curated OTU table.  
  - By processing the OTU table from the bottom, OTUs can be merged in series, with parents of errors, subsequently being merged with more abundant OTUs.  
  - This allows errors of errors to be merged with their ultimate parents, also it allows for imperfectly assigned 'daughters' in a swarm of errors or biological variants within an abundant and/or genetically variable species to be merged correctly.  
 
 
## Tutorial
___
A step-by-step walk-through with the 97% clustered (VSEARCH) data from the LULU paper.  
The first steps are carried out in Linux/Unix.  
Make a directory for the analyses:    
```
mkdir test_data
cd test_data
```


Download test data from LULU GitHub site (and OTU table and accompanying sequences/centroids)  
```
wget https://raw.githubusercontent.com/tobiasgf/lulu/master/Example_data/centroids_test.txt
wget https://raw.githubusercontent.com/tobiasgf/lulu/master/Example_data/otutable_test.txt
```

Make a match list (using blastN)  
Make blast database
```
makeblastdb -in centroids_test.txt -parse_seqids -dbtype nucl
```
Blast the centoids against themselves
```
blastn -db centroids_test.txt -outfmt '6 qseqid sseqid pident' -out match_list.txt -qcov_hsp_perc 80 -perc_identity 84 -query centroids_test.txt
```

Run the curation with LULU (done in R).  
Read the files
```
otutab <- read.csv("otutable_test.txt",sep='\t',header=TRUE,as.is=TRUE, as.is=TRUE, row.names = 1)
matchlist <- read.table("match_list.txt", header=FALSE,as.is=TRUE, stringsAsFactors=FALSE)
```
Run LULU
``
curated_result <- lulu(otutab, matchlist)
```

