BiocManager::install("SummarizedExperiment")
BiocManager::install("HDF5Array")
library("TCGAbiolinks")
library("SummarizedExperiment")
#CYU: library()

#Exercise 1.1-2.1 - different because my laptop won't let me download for some reason
setwd("../analysis_data")
sum_exp = HDF5Array::loadHDF5SummarizedExperiment(".")

#exercise 2.2
str(sum_exp)
assays(sum_exp)
rowData(sum_exp)
counts = assays(sum_exp)$"HTSeq - Counts"
head(rowData(sum_exp))
colData(sum_exp)[1:5, 25:29]
metadata( sum_exp )

#exercise 2.3
dim(colData(sum_exp))
dim(rowData(sum_exp))
dim(assays(sum_exp)$"HTSeq - Counts")
# coldata rows are individual patients, which correspond to the cols oof HTSeq - Counts
#similarly, rowdata rows are the genes (in ensemble format), and they correspond to the rows of HTSeq - Counts

#2.4
str(colData(sum_exp))
head(colData(sum_exp))
#rows are individual patients, cols are additional info about each patietn. The rows of colData correspond to the cols of the HTSeq - Counts assay

#2.5
colnames(colData(sum_exp))
#not too sure what you want, so I'm guessing age_at_diagnosis? "find...the age of the patients" is ambiguous - age when? At diagnosis? Or compared to now? Or age when they died if they died? etc. 

#2.6
colData(sum_exp)$age_at_diagnosis[1:10]
#units are days

#2.7
colData(sum_exp)$age_at_diagnosis = colData(sum_exp)$age_at_diagnosis/365 

#2.8
colData(sum_exp)$age_category = ifelse(colData(sum_exp)$age_at_diagnosis < 50, "Young", "Old")

#2.9
head(rowData(sum_exp))
dim(rowData(sum_exp))
#56602 genes
#rows are genes, cols are alternative gene names

#2.10
"RPS20" %in% rowData(sum_exp)$external_gene_name
"IL12RB1" %in% rowData(sum_exp)$external_gene_name

#2.11
assays(sum_exp)$"HTSeq - Counts"[20:25, 30:30]
#rows are genes (ensemble format), cols are individuals

#2.12
geneA_id_mask = (rowData(sum_exp)$external_gene_name == "RPS20") 
sum(geneA_id_mask)
ensembl_geneA = rowData(sum_exp)$ensembl_gene_id[geneA_id_mask]

geneB_id_mask = (rowData(sum_exp)$external_gene_name == "IL12RB1") 
sum(geneB_id_mask)
ensembl_geneB = rowData(sum_exp)$ensembl_gene_id[geneB_id_mask]
#first created mask to find where geneA was, then checked that there was only one row corresponding to geneA, finally extracted geneA by getting corresponding ensemble row


#2.13
#Row

#2.14
min(  assays(sum_exp)$"HTSeq - Counts"[ensembl_geneA, ]  )

max(  assays(sum_exp)$"HTSeq - Counts"[ensembl_geneA, ]  )

summary(assays(sum_exp)$"HTSeq - Counts"[ensembl_geneB, ] )

#2.15
plot(assays(sum_exp)$"HTSeq - Counts"[ensembl_geneA, ],
     assays(sum_exp)$"HTSeq - Counts"[ensembl_geneB, ],
     xlab = "RPS20",
     ylab = "IL12RB1"
     )
#Both genes are clustered towards 0

#2.16
bool_age_na = is.na(colData(sum_exp)$age_category)
num_na = sum(bool_age_na)
num_na

#2/17
age_cat_no_NAs <- colData(sum_exp)$age_category[!bool_age_na]

#2.18
num_na == dim(colData(sum_exp))[1] -length(age_cat_no_NAs)

#2.19
dim(colData(sum_exp))[1]
#521 patients
#They're not the same because there are 4 NAs for patient age categories

#2.20
identical( rownames(colData(sum_exp)),colnames(assays(sum_exp)$"HTSeq - Counts")  )

gene_counts = assays(sum_exp)$"HTSeq - Counts"[ensembl_geneA,(!bool_age_na)]

#2.21
length(gene_counts) == length(age_cat_no_NAs)

#2.22
boxplot(gene_counts~age_cat_no_NAs,
        xlab = "Age",
        ylab = "RPS20 Counts")

#They're pretty similar between young and old patients in the mean and quartiles,
# but the old patients have a lot more outliers

#3.1
#1. Access by calling assays(sum_exp)$"HTSeq - Counts", rows = genes in ensembl, cols = individuals
#2 - Access by calling rowData(assays(sum_exp)$"HTSeq - Counts"). The resulting dataframe stores the genes in ensembl in the rows (same as rows in assays(sum_exp)$"HTSeq - Counts") and othernames for the gene in the cols
#3 - Access by calling colData(assays(sum_exp)$"HTSeq - Counts"). Resulting dataframe stores individuals in the rows (same as cols in assays(assays(sum_exp)$"HTSeq - Counts")) and more info about the patient in the cols
