library("TCGAbiolinks")
library('survival')
library('survminer')
library("SummarizedExperiment")

getwd()
setwd("../analysis_data")

#getting our clinical data 
clinic <- data.table::fread("coad_clinical_data.csv", data.table=F)


#Getting RNA data
sum_exp = HDF5Array::loadHDF5SummarizedExperiment(".")
counts = assays(sum_exp)$"HTSeq - Counts"

#Extracting the two genes I want to look at
geneA_id_mask = (rowData(sum_exp)$external_gene_name == "TP53") 
sum(geneA_id_mask)
ensembl_geneA = rowData(sum_exp)$ensembl_gene_id[geneA_id_mask]

geneB_id_mask = (rowData(sum_exp)$external_gene_name == "KRAS") 
sum(geneB_id_mask)
ensembl_geneB = rowData(sum_exp)$ensembl_gene_id[geneB_id_mask]

#Getting some stats for expression of TP53
min(  assays(sum_exp)$"HTSeq - Counts"[ensembl_geneA, ]  )
max(  assays(sum_exp)$"HTSeq - Counts"[ensembl_geneA, ]  )
summary(assays(sum_exp)$"HTSeq - Counts"[ensembl_geneA, ] )

#Getting some stats for expression of KRAS
min(  assays(sum_exp)$"HTSeq - Counts"[ensembl_geneB, ]  )
max(  assays(sum_exp)$"HTSeq - Counts"[ensembl_geneB, ]  )
summary(assays(sum_exp)$"HTSeq - Counts"[ensembl_geneB, ] )

#Plotting gene_a and gene_b to see if there's any correlation
plot(assays(sum_exp)$"HTSeq - Counts"[ensembl_geneA, ],
     assays(sum_exp)$"HTSeq - Counts"[ensembl_geneB, ],
     xlab = "TP53",
     ylab = "KRAS",
     main="TP53 vs KRAS Expression in All Patients",
     cex.main=0.75
)

#Getting rows with valid race data
colData(sum_exp)$race
sum(is.na(colData(sum_exp)$race))

no_data_or_na_race <- ((is.na(colData(sum_exp)$race)) | (colData(sum_exp)$race == "not reported"))
sum(no_data_or_na_race)

race_has_data <- colData(sum_exp)$race[!no_data_or_na_race]

#Getting counts with valid race data
identical( rownames(colData(sum_exp)),colnames(assays(sum_exp)$"HTSeq - Counts")  )

geneA_counts = assays(sum_exp)$"HTSeq - Counts"[ensembl_geneA,(!no_data_or_na_race)]
geneB_counts = assays(sum_exp)$"HTSeq - Counts"[ensembl_geneB,(!no_data_or_na_race)]


#These plots were such a pain you would not believe 
par(mar=c(8, 4, 2, 1))

boxplot(geneA_counts~race_has_data,
        xlab = "Race",
        ylab = "TP53 Counts",
        cex.axis=0.8, 
        main="TP53 Counts Across Races",
        cex.main=0.8,
        xaxt = "n")
axis(side=1, labels=FALSE)
unique(race_has_data)
text(x = 1:length(unique(race_has_data)),
     ## Draw labels at the bottom of the chart.
     y = par("usr")[3]-3000,
     ## Use the names from the list.
     labels = unique(race_has_data),
     cex = 0.8,
     srt= 45,
     adj=1,
     xpd=NA)

boxplot(geneB_counts~race_has_data,
        xlab = "Race",
        ylab = "KRAS Counts",
        cex.axis=0.5, 
        main="KRAS Counts Across Races",
        cex.main=0.8,
        xaxt = "n")
axis(side=1, labels=FALSE)

text(x = 1:length(unique(race_has_data)),
     ## Draw labels at the bottom of the chart.
     y = par("usr")[3]-3000,
     ## Use the names from the list.
     labels = unique(race_has_data),
     ## Increase the label size a bit.
     cex = 0.8,
     srt= 45,
     adj=1,
     xpd=NA)



#creating the survival plots

#cleaning up clinic's race list so there's no NA or No Data

sum(is.na(clinic$race_list))
clinic_has_race_mask <- clinic$race_list == "No Data"
clinic_has_race <- clinic[!clinic_has_race_mask, ]

#Since I'm using the dataset we've modified before, the days to death event is replaced with days since last followup if there was no data for the days to death event. 
surv_object <- Surv(time = clinic_has_race$days_to_death, 
                    event = clinic_has_race$death_event)

# We then create a fit object

race_fit <- surv_fit( surv_object ~ clinic_has_race$race_list, data = clinic_has_race )


#the ggtheme and legend arguments are for formatting. 
# Feel free to play around with the margins and legend placement
survplot = ggsurvplot(race_fit, 
                      pval=TRUE, 
                      ggtheme = theme(plot.margin = unit(c(1,2,1,1), "cm")), 
                      legend = "right",
                      title="CRC Survival Curves By Race")

p = survplot$plot + 
  theme_bw() +  # changes the appearance to be a bit prettier
  theme(axis.title = element_text(size=10), # increase font sizes
        axis.text = element_text(size=10),
        legend.title = element_text(size=5),
        legend.text = element_text(size=5))
p

#checking to see why there's no AI/AN data
sum(clinic_has_race$race_list=="AMERICAN INDIAN OR ALASKA NATIVE")
#There's only 1...





