# Exercise 1.1
BiocManager::install("SummarizedExperiment", force = TRUE)
library(TCGAbiolinks)
library(SummarizedExperiment)
getwd()
setwd("desktop")
setwd("USC")
setwd("spring_2022/qbio490")
setwd("qbio_data_analysis_sabrinaz/analysis_data")

# Exercise 2.1
query <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Transcriptome Profiling", # get the RNA-seq transcriptome
                  data.type = "Gene Expression Quantification", # gets the counts
                  workflow.type = "HTSeq - Counts") # gets the raw counts processed by this method
GDCdownload(query) # only need to download the data once! Comment this out once you have completed it once
sum_exp <- GDCprepare(query)

# Exercise 2.2
assays(sum_exp)$"HTSeq - Counts"[1:5, 1:5]

# Exercise 2.3
dim(colData(sum_exp))
dim(rowData(sum_exp))
dim(assay(sum_exp))

# rows of rowData and assay are the same
# rows of colData and columns of assays are the same

colData(sum_exp)[1:5, 25:29]

# Exercise 2.4
str(colData(sum_exp))
head(colData(sum_exp))

# Exercise 2.5
str(colnames(sum_exp))
head(colnames(sum_exp))

# Exercise 2.6
colData(sum_exp)$age_at_diagnosis[1:10]
# units: days

# Exercise 2.7
colData(sum_exp)$age_at_diagnosis = colData(sum_exp)$age_at_diagnosis / 365
colData(sum_exp)$age_at_diagnosis[1:10]

# Exercise 2.8
colData(sum_exp)$age_category = ifelse(colData(sum_exp)$age_at_diagnosis < 50, "Young", "Old")

# Exercise 2.9
head(rowData(sum_exp))
dim(rowData(sum_exp))

# Exercise 2.10
"MSH2" %in% rowData(sum_exp)$external_gene_name
"MLH1" %in% rowData(sum_exp)$external_gene_name

# Exercise 2.11
assays(sum_exp)$"HTSeq - Counts"[20:25, 30:35]
# column: patient, row: Endsembl gene ID

# Exercise 2.12
which(rowData(sum_exp)$external_gene_name == "MSH2")
which(rowData(sum_exp)$external_gene_name == "MLH1")
# The which function give the position of the element in the row data. 
rowData(sum_exp)[2013,]
rowData(sum_exp)[1338,]

# Exercise 2.13
# Row

# Exercise 2.14
min(assays(sum_exp)$"HTSeq - Counts"[2013,])
max(assays(sum_exp)$"HTSeq - Counts"[2013,])

min(assays(sum_exp)$"HTSeq - Counts"[1338,])
max(assays(sum_exp)$"HTSeq - Counts"[1338,])

summary(assays(sum_exp)$"HTSeq - Counts"[1338,])

# Exercise 2.15
plot(assays(sum_exp)$"HTSeq - Counts"[2013,],
     assays(sum_exp)$"HTSeq - Counts"[1338,],
     xlab = "MSH2",
     ylab = "MLH1"
)
# There seems to be a positive correlation between the two genes. It is more clustered towards the origin of the graph. 

# Exercise 2.16
bool_age_na = is.na(colData(sum_exp)$age_category)
num_na = sum(bool_age_na)
num_na

# Exercise 2.17
age_cat_no_NAs = age_cat_no_NAs[!(is.na(colData(sum_exp)$age_category))]
length(age_cat_no_NAs)

# Exercise 2.18
length(age_cat_no_NAs)

numOfRows = nrow(colData(sum_exp)$age_category)
length = length(age_cat_no_NAs)
numOfRows+length==num_na

# Exercise 2.19
dim(assays(sum_exp)$"HTSeq - Counts")
# 521 Patients. This number matches the age_cat_no_NAs, but i'm not sure why. 

# Exercise 2.20
identical( rownames(colData(sum_exp)), colnames(assays(sum_exp)$"HTSeq - Counts")  )
length(age_cat_no_NAs)
gene_counts = assays(sum_exp)$"HTSeq - Counts"[2013, age_cat_no_NAs]
gene_counts

# Exercise 2.21
# The lengths do match. Gene counts is getting the patient's with an age and their counts
# of the gene I was looking at. So for gene MSH2, patient TCGA-D5-6530-01A-11R-1723-07 has
# 2009 counts of it. 

# Exercise 2.22
boxplot(gene_counts ~ age_cat_no_NAs, xlab = "Gene Counts", ylab = "Age")
# There are a lot of outliers in my boxplot. The median age from diagnosis is around 1300 days.
# The upper quantile seems fairly larger than the lower quantile. 

# Exercise 3.1
# 1) I would access the data frame through the $ sign. It would be something like (sum_exp)$"HTSeq - Counts".
# To access the rows or columns, I would use rowData and colData respectively. 
# 2) I could use the head() function to see the first few rows and columns of the data frame. Then, I can use
# the rownames column to access the rows of the data. The rows of the data usually refer to the gene ID. 
# 3) I could use the colnames function to access the columns of the data frame. The columns store the 
# different patient's information. 


