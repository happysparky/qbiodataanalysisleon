#1.1

BiocManager::install("DESeq2")
library("DESeq2")
library("SummarizedExperiment")
library("TCGAbiolinks")
library("HDF5Array")

setwd("../analysis_data")
sum_exp = HDF5Array::loadHDF5SummarizedExperiment(".")

#1.2
colnames(colData(sum_exp))
rownames(colData(sum_exp))
patient_NA <- is.na(colData(sum_exp)$age_at_diagnosis)
sum(patient_NA)
#colData has patient information

patient_data <- colData(sum_exp)
counts <- assays(sum_exp)$"HTSeq - Counts"
# counts <- data.frame(counts)
#There are 4 pateitns with NA age

length(patient_data$barcode[!patient_NA])
NA_mask <- patient_data$barcode[!patient_NA]
patient_data <- patient_data[NA_mask, ]
nrow(patient_data)

length(colnames(counts))
length(colnames(counts)[!patient_NA])
counts <- counts[,colnames(counts)[!patient_NA]]
length(colnames(counts))

patient_data$age_category <- ifelse(patient_data$age_at_diagnosis/365 < 50, "young", "old")
patient_data$age_category
patient_data$age_category = factor(patient_data$age_category, levels = c("young", "old"))

#1.3

if (all(rownames(counts) == names(rowRanges(sum_exp)))){
  print("Changing row names!")
  rownames(counts) = rowRanges(sum_exp)$external_gene_name
}
rownames(counts)

counts_row_sums <- rowSums(counts)
length(counts_row_sums)
low_counts_mask <- ifelse(counts_row_sums < 10, FALSE, TRUE)
num_low_count <- length(counts_row_sums)-sum(low_counts_mask)
num_low_count
#There are 4817 low count genes
counts <- counts[low_counts_mask,]
counts

#2.1
dds = DESeqDataSetFromMatrix(countData=counts,
                             colData=patient_data,
                             design=~age_category)
dds_obj = DESeq(dds)
resultsNames(dds_obj)  # see what comparisons got run

# get the young vs. old comparison
results = results(dds_obj, format = "DataFrame", contrast = c("age_category", "young", "old"))

#2.2
head(results)
str(results)

my_df = data.frame(x = c('b', 'd', 'c', 'e', 'a'),
                   y = c(2,4,3,5,1))

order_indices = order(my_df$y)
order_indices
my_df = my_df[order_indices, ]
my_df

#2.3
row_order <- order(results$padj)
results <- results[row_order, ]
head(results, 20)
#I will pick MAP2 - Microtubule-associated protein 2
# This gene is more highly expressed in older patients
# This gene encodes a protein that belongs in the microtubule-associated protein family. The single MAP2 gene produces four major transcripts producing four proteins, MAP2A, MAP2B, MAP2C and MAP2D.

#2.4
log2FoldChange_threshold <- 1
padj_threshold <- 0.05
over_log2foldchange_threshold <- ((results$log2FoldChange > log2FoldChange_threshold) | (results$log2FoldChange < -1*log2FoldChange_threshold)) 
under_p <- results$padj < padj_threshold
under_p[is.na(under_p)] <- FALSE

sum(is.na(over_log2foldchange_threshold))
sum(is.na(under_p))

significant <- results[over_log2foldchange_threshold&under_p,]


#2.5
fc_threshold = 2
p_threshold = 0.05
plot(x=significant$log2FoldChange,
     y=-log(significant$padj,10),
     xlab="young over old log2 fold change",
     ylab="-log10 p-value",
     pch=20)
abline(v=c(-log2(fc_threshold), log2(fc_threshold)), h= c(-log10(p_threshold)), col="green")

fc_threshold = 2  # set a threshold of at least a 2 fold increase (double)
p_threshold = 0.05  # set a threshold of adjusted p-value being <= 0.05


library(ggplot2)

volcano_plot = ggplot(data = data.frame(results), 
                      aes(x = log2FoldChange, y = -log10(padj))) + 
  geom_point(aes(color = ifelse(log2FoldChange < -1 & padj < 0.05, "lower in young",
                                ifelse(log2FoldChange > 1 & padj < 0.05, "higher in young", "NS"))),
             size = 0.5) + 
  theme_minimal() + # make things pretty +
  theme(legend.title = element_blank()) + 
  # next 2 lines draw lines at the thresholds
  geom_vline(xintercept=c(-log2(fc_threshold), log2(fc_threshold)), color="green") + 
  geom_hline(yintercept=-log10(p_threshold), color="green") + 
  scale_color_discrete(type=c("red", "blue", "black")) +
  labs(x = "log2 Fold Change (Young/Old)",
       y = "-log10 Adjusted p-value")


volcano_plot


#2.6
write.csv(x=results, 
          file="../week6_DESeq2/results.csv",
          row.names=FALSE
          )
