---
  title: "Alignment"
output: html_document
date: "2022-12-18"
---
  
  ```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
#quick update
install.packages("installr")
library(installr)
updateR()
#does not work...


#get link
biomart_link <- 'https://bioconductor.org/packages/2.3/bioc/src/contrib/biomaRt_1.16.0.tar.gz'
install.packages(biomart_link, repo = NULL, type = 'source')

#msa
packageurl <- "https://bioconductor.org/packages/3.16/bioc/src/contrib/Archive/msa/msa_1.30.0.tar.gz"	
install.packages(packageurl, repos = NULL, type = 'source')

```{r}
#TODO: load in libraries
library(biomaRt)
library(dplyr)
library(readxl)
library(stringr)
library(msa)

library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(Homo.sapiens)
```

```{r}
#loading human mart
#reference genome is hg38
mart <- useEnsembl(dataset="hsapiens_gene_ensembl",
                   biomart='ensembl')

#loading cancer excel
dataset <- read_excel('/stor/home/dhp563/Mutation_annotation/Mutation_cancer_type.xlsx')

```



```{r}
#getting list of all genes included
list_all_entrez <- unique(dataset$Entrez_Gene_Id)

output = c('ensembl_gene_id', 'start_position', 'end_position', 
           'ccds', 'cds_start', 'cds_end',  
           'chromosome_name', 'entrezgene_id')

#build a biomart query
#HG38 reference
#Assuming gene name across the old and new one remains the same
#attributes = what you want to get
#filters = the name of the variable used for filter
#values = the stuff labeled by filters

BMtable <- getBM(attributes = output, filters = 'entrezgene_id',
                 values = list_all_entrez, mart = mart)

```

```{r}
#syntax for getting sequence
#options for seq type is 'gene_exon', 'transcript_exon', 'transcript_exon_intron', 
#'gene_exon_intron', 'cdna', 'coding

sequence <- getSequence(chromosome=12, 
                        start=55966781, 
                        end = 55972789,
                        type = 'entrezgene_id',
                        seqType = 'cdna',
                        mart = mart)

getSequence(chromosome=12, 
            start=55966781, 
            end = 55972789,
            type = 'entrezgene_id',
            seqType = 'peptide',
            mart = mart)

```
```{r}
#returning genomic coordinates of coding sequence for entrez_gene_id
get_genomic_coord <- function(entrez_gene_id){
  cds <- cds(Homo.sapiens, columns="TXNAME", filter=list(gene_id=entrez_gene_id))
  cds_grl <- multisplit(cds, cds$TXNAME)
  df_grl <- as.data.frame(cds_grl)
  df_grl %>% group_by(group) %>% 
    summarize(sum_length = sum(width), AA_count = sum(width)/3)
  return(df_grl)
}

```

```{r}
#test dataframe
#how to handle this incase there are multiple versions of cds
#iterate over each version and do alignment ? 
test <- get_genomic_coord(list_all_entrez[1])
#remove chr from chrom string

#potential discrepancies between mart and TxBD?
#this doesn't work, the listed entrez gene id from getSequence doesn't match 
#with the input. 
#misunderstanding of what the output should be? 
#this doesn't work for now
test$seqnames <- substring(test$seqnames, 4)
test <- test %>% mutate(Amino_acid = getSequence(chromosome = seqnames, 
                                                 start = start, 
                                                 end = end,
                                                 type = 'entrezgene_id', 
                                                 seqType = 'cdna', 
                                                 mart = mart))


```

```{r}
#try getting introns and exons using Biomart
#how are the seqs represented? 
#this is still a combo using mutiple libraries

getSequence(id = 142,
            type = 'entrezgene_id', 
            seqType = 'gene_exon', 
            mart = mart)

```

