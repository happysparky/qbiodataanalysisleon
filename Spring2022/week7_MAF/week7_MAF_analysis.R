#1.1
BiocManager::install("maftools")
library("TCGAbiolinks")
library("maftools")
setwd("../analysis_data")

#1.2
clinic <- data.table::fread("coad_clinical_data.csv", data.table=F)
colnames(clinic)[ colnames(clinic) == "bcr_patient_barcode" ] <- "Tumor_Sample_Barcode"

#1.3
length(colnames(clinic))
#colnames(clinic) is of length 78
length(colnames(clinic) == "bcr_patient_barcode")
#colanems(clinic) == "bcr_patient_barcode is also of length 78, and contains boolean data
sum(colnames(clinic) == "bcr_patient_barcode")
#There is 1 true 

#1.4
#my laptop doesn't like GDCquery so I'm just gonna load in the csv data that David gave to me

#1.5
getwd()
list.files()
maf_dataframe = data.table::fread("TCGA.COAD.mutect.03652df4-6090-4f5a-a2ff-ee28a37f9301.DR-10.0.somatic.maf.csv", data.table=F)
clinic <- data.table::fread("coad_clinical_data.csv", data.table=F)
colnames(clinic)[ colnames(clinic) == "bcr_patient_barcode" ] <- "Tumor_Sample_Barcode"
maf_object <- read.maf(maf=maf_dataframe,
                       clinicalData=clinic,
                       isTCGA=TRUE)

#2.1
maf_object
str(maf_object)
maf_object@clinical.data
maf_object@data
colnames(maf_object@clinical.data)
colnames(maf_object@data)

#The two data frames share the Tumor_Sample_Barcode column

#3.1
library("ggplot2")
oncoplot(maf=maf_object,
         top=20)
ggsave("../week7_maf/oncoplot_top_20.png", width=10, height=6)


#3.2
#The APC gene provides instructions for making the APC protein, which plays a critical role in several cellular processes. The APC protein acts as a tumor suppressor, which means that it keeps cells from growing and dividing too fast or in an uncontrolled way. It helps control how often a cell divides, how it attaches to other cells within a tissue, and whether a cell moves within or away from a tissue. This protein also helps ensure that the number of chromosomes in a cell is correct following cell division. The APC protein accomplishes these tasks mainly through association with other proteins, especially those that are involved in cell attachment and signaling

#3.3
clinic <- maf_object@clinical.data
young_patient_ids = clinic$Tumor_Sample_Barcode[(clinic$age_category == "YOUNG")]
young_maf <- subsetMaf(maf=maf_object, tsb=young_patient_ids)

old_maf <- subsetMaf(maf=maf_object, tsb=clinic$Tumor_Sample_Barcode[(clinic$age_category == "OLD")])

#3.4
coOncoplot(m1=young_maf, m2=old_maf, m1Name="Young Patient Oncoplot", m2Name="Old Patient Oncoplot")
#dev.off()


#3.5
#the mutated genes are the same for both groups, which is interesting. I did not quite expect this because I would've thought age has a bigger effect on which genes are mutated.

#3.6
lollipopPlot(maf_object, gene="APC")
ggsave("../week7_maf/apc_lollipop.png")

#3.7
lollipopPlot2(m1=young_maf, m2=old_maf, m1_name="Young Patients", m2_name = "Old Patients", gene="APC")
ggsave("../week7_maf/apc_lollipop_age_comparison.png")

#3.8
#There are 10 samples with no mutation in gene A and no mutation in gene B. 10/40, or 25% of samples with no mutation in gene A have no mutations in gene B.

#3.9
#b=7
#c=2
#d=35
#e=37
#f=42

#3.10
geneA_maf <- subsetMaf(maf=maf_object,genes="TP53")
geneB_maf <- subsetMaf(maf=maf_object,genes="KRAS")

#3.11
str(geneA_maf)
geneA_maf
str(geneB)
geneB_maf
#subsetMaf() takes the subset of the maf_object data structure where the gene is gene A or gene B for each respective variable. 

#No, a sample can have more than one mutation. Looking at geneA_maf, we see that there are 222 counts of mutations, but only 214 samples, so by the pigeonhole principle there has to be a few overlaps. Additionally, when typing str(geneA_maf), we can scroll down to the different kinds of mutations and see that for a sample, there can be multiple mutations

#The number of samples in each data section is not the same - the number of unique samples is the same, but nonunique isn't.  

#3.12
mut_bc_geneA <- geneA_maf@clinical.data$Tumor_Sample_Barcode

mut_bc_geneB <- geneB_maf@clinical.data$Tumor_Sample_Barcode

num_mut_geneA <- length(mut_bc_geneA)
num_mut_geneB <- length(mut_bc_geneB)
#These correspond to the first row third column and third row first column sections of the table

mut_bc_geneAB = intersect(mut_bc_geneA, mut_bc_geneB)
num_mut_geneAB = length(mut_bc_geneAB)

#3.13
num_mut_geneA_only = num_mut_geneA-num_mut_geneAB
num_mut_geneB_only = num_mut_geneB-num_mut_geneAB

#3.14
num_mut_neither <- length(maf_object@clinical.data$Tumor_Sample_Barcode)-num_mut_geneA_only-num_mut_geneB_only+num_mut_geneAB


contig_table <- matrix(c(num_mut_geneAB, 
                         num_mut_geneB_only,
                         num_mut_geneA_only,
                         num_mut_neither), 
                       nrow=2)
contig_table

#3.15
fe_results <- fisher.test(contig_table)
fe_results

#These results mean that it's very likely that these mutations occur together. 