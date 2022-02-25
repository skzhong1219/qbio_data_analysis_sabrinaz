#Exercise 1.2
#BiocManager::install("maftools")
library(TCGAbiolinks)
library(maftools)

#Exercise 1.2
clinic <- data.table::fread("/Users/sabrinazhong/Desktop/USC/Spring_2022/QBIO490/qbio_data_analysis_sabrinaz/week4_clinical/coad_clinical_data.csv",
                            data.table = F)
colnames(clinic)[ colnames(clinic) == "bcr_patient_barcode" ] <- "Tumor_Sample_Barcode"

#Exercise 1.3
# colnames(clinic) has a length of 76. 
# colnames(clinic) == "bcr_patient_barcode" has a length of 524. It contains elements of type char.
# Theres only 1 true. 

#Exercise 1.4
mutation_query <- GDCquery_Maf(tumor = "COAD", 
                               pipeline = "mutect2",
                               save.csv = TRUE)

maf_object <- read.maf(maf = mutation_query, 
                       clinicalData = clinic, 
                       isTCGA = TRUE)

#Exercise 1.5
getwd()
setwd("/Users/sabrinazhong/Desktop/USC/Spring_2022/QBIO490/qbio_data_analysis_sabrinaz/analysis_data/GDCdata")
list.files()
maf_dataframe = data.table::fread("/Users/sabrinazhong/Desktop/USC/Spring_2022/QBIO490/qbio_data_analysis_sabrinaz/analysis_data/GDCdata/TCGA.COAD.mutect.03652df4-6090-4f5a-a2ff-ee28a37f9301.DR-10.0.somatic.maf.csv",
                                  data.table = F)
clinic <- data.table::fread("/Users/sabrinazhong/Desktop/USC/Spring_2022/QBIO490/qbio_data_analysis_sabrinaz/week4_clinical/coad_clinical_data.csv",
                            data.table = F)

colnames(clinic)[ colnames(clinic) == "bcr_patient_barcode" ] <- "Tumor_Sample_Barcode"
mutation_query <- GDCquery_Maf(tumor = "COAD", 
                               pipeline = "mutect2",
                               save.csv = TRUE)

maf_object <- read.maf(maf = mutation_query, 
                       clinicalData = clinic, 
                       isTCGA = TRUE)

#Exercise 2.1
maf_object
str(maf_object)
maf_object@data
maf_object@clinical.data
# They both have a column for the Tumor Sample Barcode. This makes sense because
# we need to know the information of which patient the tumor sample came from and
# all the information about the tumor's mutation in its gene sequence. 

#Exercise 3.1 
oncoplot(maf = maf_object,
         top = 20) 

library("ggplot2")

ggsave("/Users/sabrinazhong/Desktop/USC/Spring_2022/QBIO490/qbio_data_analysis_sabrinaz/week7_MAF/oncoplot.png")

#Exercise 3.2
# APC is the most mutated gene in my oncoplot. The function of the APC gene is to
# provide a template to create the APC protein. The APC protein is involved in
# suppressing tumors, so if there is a mutation in this gene, cells will start
# growing at an uncontrollable rate. 

#Exercise 3.3
# 1. Write to clinic again
clinic = maf_object@clinical.data

# 2. Create the young_patients_ids vector
clinic$age_category = ifelse(clinic$age_at_initial_pathologic_diagnosis < 50, "young", "old")
young_patients_ids = c(clinic$Tumor_Sample_Barcode[clinic$age_category == "young"])

# 3. Create another 
young_maf = subsetMaf(maf = maf_object,
                      tsb = young_patients_ids)

# 4. Repeat steps 2-3 to create an old_maf! Can you do it in one line?
old_patients_ids = c(clinic$Tumor_Sample_Barcode[clinic$age_category == "old"])
old_maf = subsetMaf(maf = maf_object,
                      tsb = old_patients_ids)

#Exercise 3.4
#install.packages("ggplot2")     # Install ggplot2 package
library("ggplot2") 
coOncoplot(m1 = young_maf, 
           m2 = old_maf, 
           m1Name = "Young Patients", 
           m2Name = "Old Patients")
ggsave("/Users/sabrinazhong/Desktop/USC/Spring_2022/QBIO490/qbio_data_analysis_sabrinaz/week7_MAF/young_vs_old_oncoplot.png")

#Exercise 3.5
#The percentages for the old patients tend to be higher. This means that the older patients
#have more of the APC, TP53, and TTN gene that is mutated. I think this is expected since
#the cell's mechanisms may not be as efficient in older individuals. 

#Exercise 3.6
dev.off()
lollipopPlot(maf_object, gene = "MSH2")

ggsave("/Users/sabrinazhong/Desktop/USC/Spring_2022/QBIO490/qbio_data_analysis_sabrinaz/week7_MAF/lollipopPlot.png")

#Exercise 3.7
lollipopPlot2(m1 = young_maf, 
           m2 = old_maf, 
           m1_name = "Young Patients", 
           m2_name = "Old Patients",
           gene = "MSH2")

ggsave("/Users/sabrinazhong/Desktop/USC/Spring_2022/QBIO490/qbio_data_analysis_sabrinaz/week7_MAF/lollipopPlot2.png")

#Exercise 3.8
#There are 10 patients that do not have a mutation in both genes A and B.

#Exercise 3.9
b = 7
c = 2
d = 35
e = 37
f = 42

#Exercise 3.10
# geneA is TP53
geneA_maf <- subsetMaf(maf = maf_object,
                       genes = "TP53")

# geneB is KRAS
geneB_maf <- subsetMaf(maf = maf_object, 
                       genes = "KRAS")

#Exercise 3.11
# subsetMaf got the data from each specific gene. It gives us how many of that gene 
# went through each type of mutation. 
# No, there are multiple types of mutations that geneA can go through. Looking at
# geneA_maf, it shows each  type of mutation and how many of geneA has that mutation. 
# It is not the same since not every patient has a mutation in geneA.

#Exercise 3.12
# 1. Access the barcodes of the patients with mutations in genes A and B
# bc stands for barcode
mut_bc_geneA = c(geneA_maf@clinical.data$Tumor_Sample_Barcode)
mut_bc_geneB = c(geneB_maf@clinical.data$Tumor_Sample_Barcode)

# 2. Get the lengths of these two vectors
num_mut_geneA = length(mut_bc_geneA)
num_mut_geneB = length(mut_bc_geneB)

# 3. Fill in the intersect here! Then get the nubmer of patients
mut_bc_geneAB = intersect(mut_bc_geneA, mut_bc_geneB)
num_mut_geneAB = length(mut_bc_geneAB)

#Exercise 3.13
num_mut_geneA_only = num_mut_geneA - num_mut_geneAB
num_mut_geneB_only = num_mut_geneB - num_mut_geneAB

#Exercise 3.14
num_neither_mutation = length(maf_object@clinical.data$Tumor_Sample_Barcode) - num_mut_geneAB - num_mut_geneA_only - num_mut_geneB_only

contig_table <- matrix(c(num_mut_geneAB, 
                         num_mut_geneB_only,
                         num_mut_geneA_only,
                         num_neither_mutation), 
                       nrow=2)
contig_table

#Exercise 3.15
fe_results <- fisher.test(contig_table)
fe_results
#The p-value is 0.06543, which is a small number. This means that there is a statistical
#significance between the mutation of gene A and B.


