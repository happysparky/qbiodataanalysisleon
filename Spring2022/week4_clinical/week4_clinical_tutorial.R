# this will install packages (if necessary) and load them
if(!require(BiocManager)) install.packages("BiocManager")

# the double colon syntax calls a function from a specific package
# this avoids loading the entire package
# in this case, we need to download TCGAbiolinks from Bioconductor using BiocManager
if(!require(TCGAbiolinks)) BiocManager::install("TCGAbiolinks")

# this just loads a package
library(TCGAbiolinks)

# Can't download stuff so (see david)
clinical_data <- read.csv("clinic.csv")
clinical_data$race_list = as.character(clinic$race_list)
# Exercise 1.1
str(clinical_data)
head(clinical_data)
#str() is compact display of an object (col names, data type, a few examples of ), head() shows first 6 rows. I find days_to_death interesting.

#Exercise 1.2
colnames(clinical_data)
clinical_data$days_to_death

#Exercise 2.1
plot(clinical_data$age_at_initial_pathologic_diagnosis, clinical_data$weight)

# exercise 2.2
unique(clinical_data$race_list)
par(mar=c(10,1,1,1))
boxplot(clinical_data$age_at_initial_pathologic_diagnosis~clinical_data$race_list,
        las = 2,
        cex.axis = 0.5)

#exercise 2.3
num_empty_string <- sum(clinical_data$race_list == "")
num_empty_string
clinical_data$race_list[clinical_data$race_list == ""] <- "No Data"

#exercise 2.4
min(clinical_data$age_at_initial_pathologic_diagnosis)
max(clinical_data$age_at_initial_pathologic_diagnosis)
mean(clinical_data$age_at_initial_pathologic_diagnosis)
median(clinical_data$age_at_initial_pathologic_diagnosis)
summary(clinical_data$age_at_initial_pathologic_diagnosis)

#exercise 2.5
num_young <- sum(clinical_data$age_at_initial_pathologic_diagnosis < 50)
num_old <- sum(clinical_data$age_at_initial_pathologic_diagnosis >= 50)
num_young+num_old == length(clinical_data[,1])


#exercise 2.6
young_patient_ids <- clinical_data$patient_id[clinical_data$age_at_initial_pathologic_diagnosis < 50]
length(young_patient_ids) == num_young

old_patient_ids <- clinical_data$patient_id[clinical_data$age_at_initial_pathologic_diagnosis >= 50]

#exercise 2.7
clinical_data$age_category <- ifelse(clinical_data$age_at_initial_pathologic_diagnosis < 50, "YOUNG", "OLD")

head(clinical_data$age_category)
sum(clinical_data$age_category == "YOUNG") == length(young_patient_ids)
#binary either young or not, so don't have to check the other

#Exercise 2.8
#first row, every column
#rows 2 to 5, every col
#every row, 3rd col
#the comma denotes a multi-dimension array
#leaving a blank before/after a comma means to get every instance of the blank along that one axis

#Exercise 2.9
young_clinic = clinical_data[clinical_data$age_at_initial_pathologic_diagnosis < 50,]

old_clinic = clinical_data[clinical_data$age_at_initial_pathologic_diagnosis >= 50,]
#age_category is a col
#patients correspond to rows

#Exercise 2.10
young_clinic_one_line = clinical_data[clinical_data$age_category == "YOUNG",]

identical(dim(young_clinic), dim(young_clinic_one_line))

# install.packages("survival")
# install.packages("survminer")
library(survival)
library(survminer)

#Exercise 3.1
na_rows_death <- is.na(clinical_data$days_to_death)

clinical_data$days_to_death[na_rows_death] <- 
clinical_data$days_to_last_follow_up[na_rows_death]
#why not use days_to_last_known_alive?

#Exerisie 3.2
clinical_data$vital_status
clinical_data$death_event <- ifelse(clinical_data$vital_status == "Alive",
                                    0,
                                    1)


# We initialize a 'survival' object first, which contains the data we need.
surv_object <- Surv(time = clinical_data$days_to_death, 
                    event = clinical_data$death_event)

# We then create a fit object
race_fit <- surv_fit( surv_object ~ clinical_data$race_list, data = clinical_data )

#the ggtheme and legend arguments are for formatting. 
# Feel free to play around with the margins and legend placement
survplot = ggsurvplot(race_fit, 
                      pval=TRUE, 
                      ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                      legend = "right")

p = survplot$plot + 
  theme_bw() +  # changes the appearance to be a bit prettier
  theme(axis.title = element_text(size=20), # increase font sizes
        axis.text = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))
p

# save the plot as a png
# if you want to present the data, change the strata labels!
ggsave("../week4_clinical/kmplot_by_race.png", plot = p, width = 12, height = 9)

#Exercise 3.3
#key takeaways: clear stratifications based on race
#questions: how can we improve care so that these differences are minimized?
#improvements: plot could use a title, which would improve clarity

# change the file path! make sure it's in your week4 folder
# we set row.names to false since the rows don't contain any relevant info
write.csv(clinical_data, "../week4_clinical/coad_clinical_data.csv", row.names = F)
getwd()
