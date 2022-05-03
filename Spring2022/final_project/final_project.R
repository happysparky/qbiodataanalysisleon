#Creating the km-curve 

library(survival)
library(survminer)

setwd("../analysis_data")
#reading in the data
clinic <- read.csv("coad_clinical_data.csv")

#creating survival object
surv_object <- Surv(time = clinic$days_to_death, 
                    event = clinic$death_event)

age_fit <- surv_fit( surv_object ~ clinic$age_category, data = clinic )

#plotting

survplot = ggsurvplot(age_fit, 
                      pval=TRUE,
                      font.x=c(50),
                      font.y=c(50),
                      ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                      xlab="Time #(Days)",
                      legend = "right",
                      )

#making appearance prettier
p = survplot$plot + 
  theme_bw() +  # changes the appearance to be a bit prettier
  theme( # increase font sizes
        axis.text = element_text(size=15),
        legend.title = element_text(size=10),
        legend.text = element_text(size=12),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20))

p
# text = element_text(family = "Times New Roman")
ggsave("../final_project/kmplot_by_age.png", plot = p, width = 10, height = 6)

#creating the bar plots

#get the young and old patients
#remove patients with no stage information
clean_clinic <- clinic[clinic$stage_event_pathologic_stage != "",]

#see different stages available
#use https://www.cancer.net/cancer-types/colorectal-cancer/stages#:~:text=Stage%20IVA%3A%20The%20cancer%20has,has%20spread%20to%20the%20peritoneum.
# as reference
unique(clean_clinic$stage_event_pathologic_stage)

#create new column based on metastatic status
clean_clinic$stage_event_pathologic_stage
clean_clinic$metastatic <- ifelse(
  clean_clinic$stage_event_pathologic_stage == "Stage IIIA" |
    clean_clinic$stage_event_pathologic_stage == "Stage IIIA" |
    clean_clinic$stage_event_pathologic_stage == "Stage IIIB" |
    clean_clinic$stage_event_pathologic_stage == "Stage IIIC" |
    clean_clinic$stage_event_pathologic_stage == "Stage IV" |
    clean_clinic$stage_event_pathologic_stage == "Stage IVA" |
    clean_clinic$stage_event_pathologic_stage == "Stage IVB", 
  "metastatic",
  "nonmetastatic")

young_patients <- clean_clinic[clean_clinic$age_category == "YOUNG",]
old_patients <- clean_clinic[clean_clinic$age_category == "OLD",]
#for each type of patient, count and graph the number of metastatic tumors

young_metastatic <- young_patients$metastatic == "metastatic"
old_metastatic <- old_patients$metastatic == "metastatic"
old_metastatic
metastaticCounts <- c(sum(young_metastatic), 
                        sum(!young_metastatic),
                        sum(old_metastatic),
                        sum(!old_metastatic))

metastaticPercentage <- c(metastaticCounts[1]/(metastaticCounts[1]+metastaticCounts[2]),
                            metastaticCounts[2]/(metastaticCounts[1]+metastaticCounts[2]),
                            metastaticCounts[3]/(metastaticCounts[3]+metastaticCounts[4]),
                            metastaticCounts[4]/(metastaticCounts[3]+metastaticCounts[4]))

par(mar = c(4, 4, 2, 0) + 0.2)
bp_counts = barplot(metastaticCounts, names.arg=c("Metastatic", "Nonmetastatic", "Metastatic", "Nonmetastatic"), cex.names=1.2, xlab="Patients", ylab="Number of Patients", col=c("red", "red", "blue", "blue"))

legend("top", legend=c("Young","Old"), fill=c("red", "blue"))
dev.off()

par(mar = c(4, 4, 2, 0) + 0.2)
bp_percentage = barplot(metastaticPercentage, names.arg=c("Metastatic", "Nonmetastatic", "Metastatic", "Nonmetastatic"), cex.names=1.2, xlab="Patients", ylab="Fraction of Patients", col=c("red", "red", "blue", "blue"))

legend("top", legend=c("Young","Old"), fill=c("red", "blue"))

# ggsave("../final_project/barplot_metastasis_by_age.png", plot = bp, width = 12, height = 9)

#co-oncoplot
library("maftools")

colnames(clinic)[ colnames(clinic) == "bcr_patient_barcode" ] <- "Tumor_Sample_Barcode"

maf_dataframe = data.table::fread("TCGA.COAD.mutect.03652df4-6090-4f5a-a2ff-ee28a37f9301.DR-10.0.somatic.maf.csv", data.table=F)

maf_object <- read.maf(maf=maf_dataframe,
                       clinicalData=clinic,
                       isTCGA=TRUE)

clinic <- maf_object@clinical.data
young_patient_ids = clinic$Tumor_Sample_Barcode[(clinic$age_category == "YOUNG")]
young_maf <- subsetMaf(maf=maf_object, tsb=young_patient_ids)

old_maf <- subsetMaf(maf=maf_object, tsb=clinic$Tumor_Sample_Barcode[(clinic$age_category == "OLD")])

#oncoplot(maf_object, top=10)
oncoplot(young_maf, top=5, titleText="Altered in 49 (100%) of 49 Young Patient Samples")
oncoplot(old_maf, top=5, titleText="Altered in 348 (100%) of 348 Old Patient Samples")

dev.off()

#png("../final_project/cooncoplot1.png", width=2500, height=1000)

#coOncoplot(m1=young_maf, m2=old_maf, m1Name="Young Patient Oncoplot", m2Name="Old Patient Oncoplot")

dev.off()


#lollipop 
lollipopPlot2(m1=young_maf, m2=old_maf, m1_name="Young Patients", m2_name = "Old Patients", gene="TP53")
# ggsave("../final_project/lollipop.png", width=10, height=6)
?lollipopPlot2


library(HDF5Array)
library("SummarizedExperiment")
sum_exp = HDF5Array::loadHDF5SummarizedExperiment(".")

counts = assays(sum_exp)$"HTSeq - Counts"
assays(sum_exp)$"HTSeq - Counts"
gene_id_mask = (rowData(sum_exp)$external_gene_name == "TP53")
ensembl1_gene = rowData(sum_exp)$ensembl_gene_id[gene_id_mask]
ensembl1_gene
colData(sum_exp)

bool_age_na = is.na(colData(sum_exp)$age_at_diagnosis)

age_no_NAs = colData(sum_exp)$age_at_diagnosis[!bool_age_na]
age_no_NAs = age_no_NAs/365
age_category = ifelse(age_no_NAs < 50, "Young", "Old")

gene_counts = assays(sum_exp)$"HTSeq - Counts"[ensembl1_gene, !bool_age_na]

par(mar=c(5,5,2,5))
boxplot(gene_counts ~ age_category,
        xlab = "Age of patient",
        ylab = "TP53 Gene Counts")

