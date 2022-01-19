if (!require(TCGAbiolinks)) BiocManager::install("TCGAbiolinks")

library(TCGAbiolinks)
library(SummarizedExperiment)

#add barcodes argument to query if you want to run on your local machine for smaller files downloaded
barcodes_rnaseq <- c("TCGA-BH-A0DG-01A-21R-A12P-07","TCGA-A2-A0YF-01A-21R-A109-07",
          "TCGA-AN-A04A-01A-21R-A034-07","TCGA-AR-A1AP-01A-11R-A12P-07",
           "TCGA-A2-A0T3-01A-====21R-A115-07", "TCGA-E2-A154-01A-11R-A115-07" )
#barcodes_clinic <- c("TCGA-BH-A0DG","TCGA-A2-A0YF","TCGA-AN-A04A","TCGA-AR-A1AP", "TCGA-A2-A0T3",
#                      "TCGA-E2-A154", "TCGA-AO-A12F", "TCGA-A2-A0YM", "TCGA-BH-A0E0", "TCGA-AN-A0FL")


query <- GDCquery(project = "TCGA-BRCA",
                   data.category = "Transcriptome Profiling",
                   data.type = "Gene Expression Quantification",
                   workflow.type = "HTSeq - Counts",
                   barcode = barcodes_rnaseq)
#GDCdownload(query) #only need this line of code once to download the data
sum_exp <- GDCprepare(query)
# Create a tutorial on SummarizedExperiment

patient_data <- colData(sum_exp)

#get rid of stage iv cancer types
patient_data<-patient_data[!(patient_data$tumor_stage=="stage iv"),]

patient_ages <- patient_data$paper_age_at_initial_pathologic_diagnosis
patient_data$age_category = ifelse(patient_ages < 40, "Young", ifelse(patient_ages >= 60, "Old", "Mid"))

htseq_counts <- assays(sum_exp)$"HTSeq - Counts"

patient_data$GATA3_counts = htseq_counts["ENSG00000107485",]
patient_data$MUC2_counts = htseq_counts["ENSG00000198788",]
patient_data$KMT2C_counts = htseq_counts["ENSG00000055609",]

patient_data$GATA3_counts_log = sapply(htseq_counts["ENSG00000107485",], log10)
patient_data$MUC2_counts_log = sapply(htseq_counts["ENSG00000198788",], log10)
patient_data$KMT2C_counts_log = sapply(htseq_counts["ENSG00000055609",], log10)

# Boxplots by age
boxplot(GATA3_counts_log~age_category, main = "Boxplot of HTSeq - Counts for GATA3 by Age Category", data = patient_data)
boxplot(MUC2_counts_log~age_category, main = "Boxplot of HTSeq - Counts for MUC2 by Age Category", data = patient_data)
boxplot(KMT2C_counts_log~age_category, main = "Boxplot of HTSeq - Counts for KMT2C by Age Category", data = patient_data)
plot(x = patient_data$paper_age_at_initial_pathologic_diagnosis, y = patient_data$GATA3_counts)
plot(x = patient_data$paper_age_at_initial_pathologic_diagnosis, y = patient_data$MUC2_counts)
plot(x = patient_data$paper_age_at_initial_pathologic_diagnosis, y = patient_data$KMT2C_counts)


# Add a new column to colData called "age_category"
# If age_at_initial_pathologic_diagnosis is < 40, define patient as "Young" in new column                               # If age_at_initial_pathologic_diagnosis is >= 60, define patient as "Old" in new column
# Other patients (between 40 and 60), define as "Mid"==
# Choose 3 genes of interest from the paper presentations
# Create 3 different boxplots with age_category on x-axis, counts on the y-axis by repeating the below code for each ge$# remember to give your plot a title and informative axis labels
# png("boxplot_genename.png")
# boxplot(FILL IN HERE)
# *Feel free to add lines here that format your boxplot*
# dev.off()
# Use rsync to copy your create pngfile to local machine for viewing
