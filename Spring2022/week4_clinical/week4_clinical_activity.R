clinic <- read.csv("coad_clinical_data.csv")

# Written activity
# 1: A categorical variable is a variable that is split into categories. For example, a column that stored colors might have variables such as "Red", "Green", "Blue". A discrete variable is a variable that are countable in a finite amount of time. For example, number of students in a classroom (you can't have half a student (unless you're REALLY morbid)). A continuous variable is obtained from measuring. For example, height is continuous because you can always add more sigfigs between two heights to get an inifnite, continous, range. 

#2: 
colnames(clinic)
sum(clinic$colon_polyps_present == "")
length(clinic$colon_polyps_present)
# We chose colon_polyps_present. About half is empty, but there's still ~270ish that's present, and that's a pretty sizeable amount.

#3
#Colon polyps are usually found during colonoscopies. Our variable is categorical. 

#4
# https://oce-ovid-com.libproxy2.usc.edu/article/00005792-202112300-00005/HTML
# The differences between fecal microbiota and intestinal fluid microbiota in colon polyps: An observational study
# There are two general types of intestinal microbiota: intestinal cavity (default) and mucosal. Fresh feces and intestinal fluid was sampled from patients with colon polyps and alpha and beta diversities were calculated. Differences in alpha diversity were not statistically significant, but there were statistical differences in beta diveristy. Thus, there are statistical differences in compositionb etween intestinal microbiota and fecal microbiota in colon polyp patients. 


# Article 2: Serrated polyps of the colon and rectum (hyperplastic polyps, sessile serrated adenomas, traditional serrated adenomas, and mixed polyps)—proposal for diagnostic criteria
# Article 2 URL: https://www.proquest.com/docview/750346263/fulltextPDF/42AE1522CDB6419DPQ/1?accountid=14749
# ARticle 2 Summary: Serrated lesions of polyps show characteristic epigenetic alterations not commonly seen in colorectal adenomas and progress to colorectal carcinoma via the so-called serrated pathway. his group of polyps is comprised not only of hyperplastic polyps, but also of sessile serrated adenomas, traditional serrated adenomas and mixed polyps, showing serrated and “classical” adenomatous features. Diagnostic criteria and nomenclature for these lesions are not uniform and, therefore, somewhat confusing


#5
#The second variable we chose is age_at_initial_pathologic_diagnosis. It describes the age at which someone was diagnosed with CRC. It is measured by how old someone is (that is, years after being born).

#6
#Hypothesis 1: increased age at diagnosis correlated with having colon polyps
#hypothesis 2: having colon polyps is correlated with decreased survival
#hypothesis 3: increased age at diagnosis is correlated with decreased survival 

#7 
# a) First, we explored polyp status by age. We found that in general,
# these patients are less likely to have polyps. The ratios of polyps to
# no polyps in both young and old populations were roughly even (12:21 and 80:136,
# respectively.) This suggests our first hypothesis is incorrect - age doesn't
# really affect polyp status. 
# 
# b) Next, we looked how colon polyps were correlated with survival. We found
# that among those that have polyps, the amount of patients that survived vs
# the amount of patients that died were about even. On the other hand,
# for those without polyps, the chance of survival was much higher. These findings
# indicate that if someone survives, it's much more likely they have no polyps, although
# having polyps itself isn't an indicator for death, partially proving our second hypothesis
#
#c) Finally, we looked how age at diagnosis was correlated with survival. Similar 
# to the findings above, the incidents of survival and death were approximately equal
# among the young population. However, among the old population, many more patients
# died. This suggests that if a patient died, they were much more likely to be old. 
# However, if they survived, then they can be either young or old. These results 
# partially prove our third hypothesis. 



#Coding section

#1 - exploring the data
#create a mask of rows with age data
sum(is.na(clinic$age_at_initial_pathologic_diagnosis))
#no empty values for age, so don't have to worry about that

#create mask for patients that are young
is_young <- clinic$age_at_initial_pathologic_diagnosis < 50

#create mask for patients where we have data that polyps are present
has_polyps <- clinic$colon_polyps_present == "YES"
#mask for patients where we have data that polyps aren't present
#can't just use inverse of has_polyps because we have no data cols
no_polyps <- clinic$colon_polyps_present == "NO"
#Create counts of possible polyp and age combinations
polyps_young <- length(clinic$colon_polyps_present[has_polyps & is_young])
polyps_young
no_polyps_young
polyps_old
no_polyps_old
no_polyps_young <- length(clinic$colon_polyps_present[no_polys & is_young])
polyps_old <- length(clinic$colon_polyps_present[has_polyps & !is_young])
no_polyps_old <- length(clinic$colon_polyps_present[no_polyps & !is_young])

#increase margin of plot
par(mar=c(9,4,4,4))
#plot the data
barplot(c(polyps_young, no_polyps_young, polyps_old, no_polyps_old),
        main="Polyp Status By Age",
        names=c('polyps_young', 'no_polyps_young', 'polyps_old', 'no_polyps_old'),
        las=2)
#2 - graphing survival and the first variable (yes polyp or no poly)
# ensure no empty death events
sum(is.na(clinic$death_event))

# there's a death event for every thing, so only need to limit by polyps
#Create counts of possible polyp and survival combinations
polyps_survived <- sum(clinic$death_event == 0 & has_polyps)
polyps_died <- sum(clinic$death_event == 0 & has_polyps)
no_polyps_survived <- sum(clinic$death_event == 0 & no_polyps)
no_polyps_died <- sum(clinic$death_event == 1 & no_polyps)

#plotb the data
barplot(c(polyps_survived, polyps_died, no_polyps_survived, no_polyps_died),
        main="Polyp Status and Survival",
        names=c('polyps_survived', 'polyps_died', 'no_polyps_survived', 'no_polyps_died'),
        las=2)

#3 - graphing survival and the second variable (young or old)
#Create counts of possible age at diagnosis and survival combinations
young_survived <- sum(clinic$death_event == 0 & clinic$age_category == "YOUNG")
young_died <- sum(clinic$death_event == 0 & clinic$age_category == "YOUNG")
old_survived <- sum(clinic$death_event == 1 & clinic$age_category == "OLD")
old_died <- sum(clinic$death_event == 0 & clinic$age_category == "OLD")

#plotb the data
barplot(c(young_survived, young_died, old_survived, old_died),
        main="Age at Diagnosis and Survival",
        names=c('young_survived', 'young_died', 'old_survived', 'old_died'),
        las=2)

