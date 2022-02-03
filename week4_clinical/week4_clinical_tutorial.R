#1 Accessing Data from TCGA Using TCGAbiolinks

# this will install packages (if necessary) and load them
if(!require(BiocManager)) install.packages("BiocManager")

# the double colon syntax calls a function from a specific package
# this avoids loading the entire package
# in this case, we need to download TCGAbiolinks from Bioconductor using BiocManager
if(!require(TCGAbiolinks)) BiocManager::install("TCGAbiolinks")

# this just loads a package
library(TCGAbiolinks)

setwd("/Users/sabrinazhong/desktop/usc")
setwd("spring_2022/qbio490")
setwd("qbio_data_analysis_sabrinaz")
setwd("analysis_data")
getwd()


write.csv(clinic, "/Users/sabrinazhong/Desktop/coad_clinical_data.csv", row.names = F)
clinic_read_in <- read.csv("/Users/sabrinazhong/Desktop/coad_clinical_data.csv")
View(clinic_read_in)

clin_query <- GDCquery(project = "TCGA-COAD", data.category = "Clinical", file.type = "xml")
# Only use this line ONCE! Comment out after you have downloaded the data. 
#GDCdownload(clin_query)
clinic <- GDCprepare_clinic(clin_query, clinical.info = "patient")
# Just adding an underscore between follow and up
names(clinic)[names(clinic)=="days_to_last_followup"] <- "days_to_last_follow_up"

#Exercise 1.1
str(clinic)
head(clinic)

#Exercise 1.2
colnames(clinic)

#Exercise 2.1
plot(clinic$age_at_initial_pathologic_diagnosis, clinic$weight, xlab = "Age at Initial Pathologic Diagnosis", 
     ylab = "Weight")

#Exercise 2.2
unique(clinic$race_list)
#the mar argument in par() sets the plot margins. 
# As the race names may be long, we want to have a large bottom margin. 
# The syntax is par(mar = c(bottom, left, top, right)), where mar is margins
# How does the plot change if you change these margins?
par(mar=c(10,1,1,1))
boxplot(clinic$age_at_initial_pathologic_diagnosis~clinic$race_list, xlab = "Race", 
        ylab = "Age", las = 2, cex.axis = 0.5)

#Exercise 2.3
View(clinic)
clinic$race_list = as.character(clinic$race_list)
nrow(clinic[clinic$race_list == "",])

clinic_read_in$race_list <- ifelse(clinic_read_in$race_list == "", "NO DATA", clinic_read_in$race_list)

#Exercise 2.4
min(clinic_read_in$age_at_initial_pathologic_diagnosis)
max(clinic_read_in$age_at_initial_pathologic_diagnosis)
mean(clinic_read_in$age_at_initial_pathologic_diagnosis)
median(clinic_read_in$age_at_initial_pathologic_diagnosis)

#Exercise 2.5
nrow(clinic_read_in[clinic_read_in$age_at_initial_pathologic_diagnosis < 50, ])
nrow(clinic_read_in[clinic_read_in$age_at_initial_pathologic_diagnosis >= 50, ])
dim(clinic_read_in)

#Exercise 2.6
young_patient_ids = clinic_read_in$patient_id[clinic_read_in$age_at_initial_pathologic_diagnosis < 50,]
old_patient_ids = clinic_read_in$patient_id[clinic_read_in$age_at_initial_pathologic_diagnosis >= 50,]

#Exercise 2.7
clinic_read_in$age_category = ifelse(clinic_read_in$age_at_initial_pathologic_diagnosis < 50, "young", "old")
head(clinic_read_in$age_category)

#Exercise 2.8
clinic[1,1] #This is the top left entry of the dataframe. R has "one-based indexing"
clinic[1,] #This shows everything from the patient in the first row under each column
clinic[2:5,] #This shows the data from row 2 to 5 in each column
clinic[,3] #This shows all the data for the third column of each patient
#Leaving a blank means that you are not specifying the position in that row/column. You simply just print out
#everything from the un specified row or column.

#Exercise 2.9
young_clinic = clinic_read_in[clinic_read_in$age_category == "young",]
head(young_clinic)
old_clinic = clinic_read_in[clinic_read_in$age_category == "old",]
head(old_clinic)

#Exercise 2.10
young_clinic_one_line = clinic_read_in[clinic_read_in$age_at_initial_pathologic_diagnosis < 50,]
identical(dim(young_clinic), dim(young_clinic_one_line))

#3 Kaplan-Meier Curves
install.packages("survival")
library(survival)
install.packages("survminer")
library(survminer)

#Exercise 3.1
View("survival")
View("survminer")

clinic_read_in$days_to_death <- ifelse(is.na(clinic_read_in$days_to_death), 
                                       clinic_read_in$days_to_last_followup, clinic_read_in$days_to_death)

#Exercise 3.2
head(clinic_read_in$vital_status)
clinic_read_in$death_event <- ifelse(clinic_read_in$vital_status == "Alive", 0, 1)
head(clinic_read_in$death_event)

# We initialize a 'survival' object first, which contains the data we need.
surv_object <- Surv(time = clinic_read_in$days_to_death, 
                    event = clinic_read_in$death_event)

# We then create a fit object
race_fit <- surv_fit( surv_object ~ clinic_read_in$race_list, data = clinic )

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
ggsave("qbio_data_analysis_sabrinaz/analysis_data/kmplot_by_race.png", plot = p, width = 12, height = 9)


#Exercise 3.3
#From the plot, it seems that Asian have the highest and a 100% survival rate. Black or African American
#patients have the lowest survival rate. It seems as though the data for Asians does not span over a long
#amount of time as the line for it gets cut shorter than the others. I think this might raise questions as
#to what happens to the sruvival rates for Asians past the time on the graph. A 100% survival does not seem
#very plausible. I think there needs to be more data on Asian patients for a longer time. 

#4 Wrap Up
#Done in the beginning (line 19)

