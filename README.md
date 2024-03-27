# LULU
___
A r-package for distribution based post clustering curation of amplicon data.

The purpose of **LULU** is to reduce the number of erroneous OTUs in OTU tables to achieve more realistic biodiversity metrics. By evaluating the co-occurence patterns of OTUs among samples **LULU** identifies OTUs that consistently satisfy some user selected criteria for being errors of more abundant OTUs and merges these. It has been shown that curation with **LULU** consistently result in more realistic diversity metrics. The required input of **LULU** is an OTU table and a corresponding matchlist with all the internal matches of OTUs.  

The method is published here:  
Frøslev, T. G., Kjøller, R., Bruun, H. H., Ejrnæs, R., Brunbjerg, A. K., Pietroni, C., & Hansen, A. J. (2017). Algorithm for post-clustering curation of DNA amplicon data yields reliable biodiversity estimates. Nature Communications, 8(1), 1188.  

https://www.nature.com/articles/s41467-017-01312-x


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

BLASTn - a match list can be produced with **BLASTn** with these commands:  
```
#First produce a blastdatabase with the OTUs
makeblastdb -in OTU_sequences.fasta -parse_seqids -dbtype nucl
# Then blast the OTUs against the database
blastn -db OTU_sequences.fasta -outfmt '6 qseqid sseqid pident' -out match_list.txt -qcov_hsp_perc 80 -perc_identity 84 -query OTU_sequences.fasta
```

**VSEARCH** - Alternatively, a matchlist can also be produced with **VSEARCH** with this command
```
vsearch --usearch_global OTU_sequences.fasta --db OTU_sequences.fasta --self --id .84 --iddef 1 --userout match_list.txt -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10
```

No matter which algorithm is used for internal OTU matching, the resulting match list, should look like this (with match in percent, with or without decimals):
```
OTU1  OTU2  90.000
OTU1  OTU3  88.000
OTU3  OTU3  100.000
OTU3  OTU4  82.123
OTU4  OTU9  86.333
...
```

#### 3. Run the **LULU** curation in **R**
With the OTU table and the corresponding match list at hand, we can carry out the curation with **LULU** in **R**.
```
> curated_result <- lulu(otutable_name, matchlist_name)
```
The `lulu` function will return a curated OTU table and some statistics and information on the curation in the object `curated_result`. The curated OTU table can be accessed like this:
```
curated_result$curated_table
```
This table can then be used for further biodiversity analyses, etc.  
For more advanced uses, see the help-file for the function.  
```
>?lulu  
```

## Installation
___
The lulu package can be installed in **R (RStudio)** using **devtools**, by typing these commands in **R**
```
> library(devtools)
> install_github("adrientaudiere/lulu")  
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
  - If the number of samples where the 'potential parent' has a positive presence is below the `minimum_relative_cooccurence` (default 95%, meaning that 1 in 20 samples are allowed to have no parent presence), the 'potential parent' is rejected.  
  - If not, the abundance ratio between 'potential daughter' and 'potential parent' is calculated for all samples, and tested against the `minimum_ratio`.  
  - The minimum ratio  can be set to any number, and the threshold (`minimum_ratio_type`) can be set to be evaluated as the minimum (`min`) observed ratio or the average (`avg`) observed ratio.  
  - If the potential parent satisfies the ratio threshold, the 'potential daughter' is flagged as an error of this OTU.  
6. If no potential parent satisfies the criteria above, the 'potential daughter' is flagged as a valid OTU.  
7. repeat from 2 with the next OTU until the full table has bee processed.  
  - After this step all OTUs will either be flagged as an error of a more abundant OTU, or as a valid OTU.  
8. Process the OTU table from bottom to top, merging flagged errors with their parents to produce a curated OTU table.  
  - By processing the OTU table from the bottom, OTUs can be merged in series, with parents of errors, subsequently being merged with more abundant OTUs.  
  - This allows errors of errors to be merged with their ultimate parents, also it allows for imperfectly assigned 'daughters' in a swarm of errors or biological variants within an abundant and/or genetically variable species to be merged correctly.  
 
 
## Tutorial
___
A step-by-step walk-through with the 97% clustered (**VSEARCH**) data from the **LULU** paper.  
The first steps are carried out in Linux/Unix.  
Make a directory for the analyses:    
```
mkdir test_data
cd test_data
```


Download test data from this **LULU** GitHub site (and OTU table and accompanying sequences/centroids)  
The OTU table contains 2425 OTUs destributed over 130 sites/samples.
```
wget https://raw.githubusercontent.com/tobiasgf/lulu/master/Example_data/centroids_test.txt
wget https://raw.githubusercontent.com/tobiasgf/lulu/master/Example_data/otutable_test.txt
```

Make a match list (using **BLASTn**)  
Make blast database
```
makeblastdb -in centroids_test.txt -parse_seqids -dbtype nucl
```
Blast the centoids against themselves
```
blastn -db centroids_test.txt -outfmt '6 qseqid sseqid pident' -out match_list.txt -qcov_hsp_perc 80 -perc_identity 84 -query centroids_test.txt
```

Run the curation with **LULU** (done in **R**).  
Read the files
```
otutab <- read.csv("otutable_test.txt",sep='\t',header=TRUE,as.is=TRUE, row.names = 1)
matchlist <- read.table("match_list.txt", header=FALSE,as.is=TRUE, stringsAsFactors=FALSE)
```

Run **LULU**  
```
> curated_result <- lulu(otutab, matchlist)

# ....Which is equivalent of running LULU with default settings for the options minimum_ratio_type, minimum_ratio, minimum_relative_cooccurence  

> curated_result <- lulu(otutab, matchlist, minimum_ratio_type = "min", minimum_ratio = 1, minimum_match = 84, minimum_relative_cooccurence = 0.95)
```

The curated OTU table can now be accessed here:
```
> curated_result$curated_table
```
...and the original table here:
```
curated_result$original_table
```

...and other information on the data.  
Number of OTUs retained:
```
> curated_result$curated_count
[1] 459
```
IDs of OTUs retained (list only first)
```
> head(curated_result$curated_otus)
[1] "001540ee723f903acddaf1c993cefde9c0c43d67" "0039237b99e15bfe7e12f4354d2dfd03b5ae22b0"
[3] "014b752e527009320f56afb0df2caf0591d060ef" "014f57228502f64234316c86f180f555b3151464"
[5] "016bb353b6b859d79a1ef36863cf0850806a3c06" "024f6daf85d8aed1f23c119c31a8acfba5b82a66"
```
(...I use sha1 hash names for my OTUs... hence the strange IDs...)

number of OTUs discarded  
```
> curated_result$discarded_count
[1] 1966
```

IDs of OTUs discarded (list only first)
```
> head(curated_result$discarded_otus)
[1] "7c535c7709639b9ec025858cca671d406966a653" "1ea168de62e8686635707db62629aae301a14b2b"
[3] "0c2c529cbd545bc3675f3433b0160e0cb56c4b2c" "ee7271685168ed084bfcaa5515caf4761012f260"
[5] "2e721d0157683b7e1ab7999fc5c9d22b0a3b4397" "57f40e612102dc39214b2670caca5d8c3f5b7897"
```

Check the computation time
```
> curated_result$runtime
Time difference of 1.053344 mins
```

Check which setting was used for `minimum_match`
```
> curated_result$minimum_match
[1] 84
```
Check which setting was used for `minimum_relative_cooccurence`
```
> curated_result$minimum_relative_cooccurence
[1] 0.95
```
Check how the OTUs were mapped.  
This file includes som basic stats:  
*total* - total read count  
*spread* - the number of samples the OTU is present in  
*parent_id* - ID of OTU with which this OTU was merged (or self)  
*curated* - ("parent" or "merged"), was this OTU kept as a valid OTU (parent) or merged with another  
*rank* - The rank of the OTU in terms of decreasing spread and read count    
```
> head(curated_result$otu_map)
                                          total spread                                parent_id curated rank
ec84eb6504ec23a3fe659c533bf9b3f08f5bd1cb 136715     58 ec84eb6504ec23a3fe659c533bf9b3f08f5bd1cb  parent    6
79a49b866cf4bdc00d11eb1c7b91957ce15a0314 104908     50 79a49b866cf4bdc00d11eb1c7b91957ce15a0314  parent   11
c2f02be9235142d605aaa5170f38d5a9c8a684de  98839     45 c2f02be9235142d605aaa5170f38d5a9c8a684de  parent   13
9b88a08f039c7bfc513e52b4369b4f05857cb1f5 171279     42 9b88a08f039c7bfc513e52b4369b4f05857cb1f5  parent    5
a2e5ad0bd2a99776da541051125b0ad377f7ea6e 634469     41 a2e5ad0bd2a99776da541051125b0ad377f7ea6e  parent    1
aafb7fcf4cfed42eaae4141f2af712b5ca7db7f0 301433     36 aafb7fcf4cfed42eaae4141f2af712b5ca7db7f0  parent    3

# And checking somewat further down the table
> curated_result$otu_map[300:308,]
                                         total spread                                parent_id curated rank
709d050ce8c823a6650c74f085ff034093e3ad42   108      5 ec84eb6504ec23a3fe659c533bf9b3f08f5bd1cb  merged  483
1dbf509b1f6bd8470354e29855459b1c0bf4d033    94      5 0b2e099f3eebf3ef942767f4c190c4ec703bbe30  merged  518
01b7e27549e043e22aebe8c215746fdbbd37a4e4    88      5 ec84eb6504ec23a3fe659c533bf9b3f08f5bd1cb  merged  537
3eda946fc9435377e003bf85089f75ddf7972a7d    87      5 a2e5ad0bd2a99776da541051125b0ad377f7ea6e  merged  539
c623dbeece5ff9df34c56decee695be51d30b5e1    86      5 ec84eb6504ec23a3fe659c533bf9b3f08f5bd1cb  merged  547
212eae5cd8133d47b085ea861a0f6865928a9276    85      5 aafb7fcf4cfed42eaae4141f2af712b5ca7db7f0  merged  550
22740909902c879d2b044a0ac8ac4bbdee2a9bdf    79      5 bd34bf9b277639657f65381c53d7715718a184c7  merged  563
3f17a0b4a4097f5348fa817a1ada92ec3ae7d37e    67      5 aafb7fcf4cfed42eaae4141f2af712b5ca7db7f0  merged  611
9612a49d162af29198945e1b09ddf0616da0288f    65      5 aafb7fcf4cfed42eaae4141f2af712b5ca7db7f0  merged  619
```

The `lulu` function also produces a log file (named something like lulu.log_20171113_162157) which will be placed in the working directory.  
For each OTU processed, the log file contains:  
1) a list of all *hits*, i.e. other OTUs with a sequence similarity above the selected threshold) in the dataset is listed, and   
2) all *potential parents*, i.e. hits with a lower rank number, i.e. higher spread and total read count, and satisfying the selected ratio of read counts, and   
3) *relative co-occurence* of all parents (until a parent satisfying the minimum relative cooccurence (and min avg abundance) thresholds is met, if one such is present), and
4) *min avg abundance* of parents satisfying minimum relative co-occurence, and
5) information whether the OTU was found to have a parent or not ("No parent found!" or "SETTING XXX to be an ERROR of YYY")   

The example below shows parts of the output for an OTU (b168b5b94056f0eef180562cbf6b24bdef011758) that did not have a parent (i.e. was retained as a valid OTU), and another OTU (1ea168de62e8686635707db62629aae301a14b2b) that was found to be an error of another OTU (79a49b866cf4bdc00d11eb1c7b91957ce15a0314).

```
####processing: b168b5b94056f0eef180562cbf6b24bdef011758 #####
<< ... NB: many hits left out ... >>
---hits: 87e8911113447fea3f9cdb75ff86ff3ab89f5add 
---hits: e6accfc1b24e27860d90ef26f2336c24b5dafb75 
---hits: 67b28704e8b6c6489bc43583a4aff4919b7013b3 
---hits: 54c8d513876fa4aa520b9b439d0c2253de80c479 
---hits: a4d4d755b2f698b95d028dccaff3bc58a9488c4f 
---hits: eaaabfb22a7d3beb15ddeb189c00c925a74f4ad9 
---hits: 1fe335b946757d5c28cb3a863cb55dd89db25465
---potential parent: ec84eb6504ec23a3fe659c533bf9b3f08f5bd1cb 
---potential parent: 79a49b866cf4bdc00d11eb1c7b91957ce15a0314 
---potential parent: c2f02be9235142d605aaa5170f38d5a9c8a684de 
---potential parent: 0d87b85358966ea5287480c1e236b06814bd1060 
---potential parent: a666e5d1e6860d8fd9eddb5018020c15124e97a1 
---potential parent: 1ea168de62e8686635707db62629aae301a14b2b
------checking: ec84eb6504ec23a3fe659c533bf9b3f08f5bd1cb
------relative cooccurence: 0.68
------checking: 79a49b866cf4bdc00d11eb1c7b91957ce15a0314
------relative cooccurence: 0.64
------checking: c2f02be9235142d605aaa5170f38d5a9c8a684de
------relative cooccurence: 0.64
------checking: 0d87b85358966ea5287480c1e236b06814bd1060
------relative cooccurence: 0.08
------checking: a666e5d1e6860d8fd9eddb5018020c15124e97a1
------relative cooccurence: 0.52
------checking: 1ea168de62e8686635707db62629aae301a14b2b
------relative cooccurence: 0.36
No parent found!


####processing: 1ea168de62e8686635707db62629aae301a14b2b #####
<< ... NB: many hits left out ... >>
---hits: 21bebcd27a1acb22fe170631aa3d6443a8ab34ba 
---hits: a7d10f70afd816b32dd68b9ae38192b9b73a6404 
---hits: 7801ac4bf8cc9611cb93899aebd0ec6b97c90601 
---hits: f8e6e2ccf691d010b2ddbd6eb29a57334cd029d7 
---hits: 4c609758a47745d53e0783534a2a9bec66771c8b 
---hits: aedf2739938830e87ace58dcaedeaaf5567328d0
---potential parent: ec84eb6504ec23a3fe659c533bf9b3f08f5bd1cb 
---potential parent: 79a49b866cf4bdc00d11eb1c7b91957ce15a0314 
---potential parent: c2f02be9235142d605aaa5170f38d5a9c8a684de 
---potential parent: 0d87b85358966ea5287480c1e236b06814bd1060 
---potential parent: a666e5d1e6860d8fd9eddb5018020c15124e97a1 
---potential parent: b168b5b94056f0eef180562cbf6b24bdef011758
------checking: ec84eb6504ec23a3fe659c533bf9b3f08f5bd1cb
------relative cooccurence: 0.92
------checking: 79a49b866cf4bdc00d11eb1c7b91957ce15a0314
------relative cooccurence: 1 which is sufficient!
------min avg abundance: 27.8888888888889 which is OK!
SETTING 1ea168de62e8686635707db62629aae301a14b2b to be an ERROR of 79a49b866cf4bdc00d11eb1c7b91957ce15a0314
```
