library(TCGAbiolinks)
library(survival)
library(survminer)
library(SummarizedExperiment)
library(ggplot2)
getwd()
setwd("/Users/sabrinazhong/Desktop/USC/Spring_2022/QBIO490/qbio_data_analysis_sabrinaz/analysis_data")

clin_query <- GDCquery(project = "TCGA-COAD", 
                       data.category = "Clinical", 
                       file.type = "xml")
clinic <- GDCprepare_clinic(clin_query, clinical.info = "patient")
names(clinic)[names(clinic)=="days_to_last_followup"] <- "days_to_last_follow_up"

gene_query <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Transcriptome Profiling", 
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts")
sum_exp <- GDCprepare(gene_query)

# Using boolean indexing to find MSH2's ID
"MSH2" %in% rowData(sum_exp)$external_gene_name
MSH2_id_mask = rowData(sum_exp)$external_gene_name == "MSH2"
sum(MSH2_id_mask)
ensembl_MSH2 = rowData(sum_exp)$external_gene_id[MSH2_id_mask]

"MLH1" %in% rowData(sum_exp)$external_gene_name
MLH1_id_mask = rowData(sum_exp)$external_gene_name == "MLH1"
sum(MLH1_id_mask)
ensembl_MLH1 = rowData(sum_exp)$external_gene_id[MLH1_id_mask]

# Finding the mean, minimum, and maximum counts of MSH2
summary(assays(sum_exp)$"HTSeq - Counts"[MSH2_id_mask,]) 
MSH2_counts = assays(sum_exp)$"HTSeq - Counts"[MSH2_id_mask,]
MSH2_counts_median = median(MSH2_counts)

MLH1_counts = assays(sum_exp)$"HTSeq - Counts"[MLH1_id_mask,]
MLH1_counts_median = median(MLH1_counts)

patient_gender = colData(sum_exp)$gender
gene_clinic = data.frame(colData(sum_exp))

# getting number of patients in study
clinic$gender = as.character(clinic$gender)
clinic$gender[clinic$gender == ""] = "No data"
table(clinic$gender)
# 244 females, 280 males

#scatterplot of gene counts between MSH2 and MLH1 counts
plot(MSH2_counts, MLH1_counts, xlab = "MSH2 Counts", ylab = "MLH1 Counts", main = "MSH2 vs. MLH1 Counts")

boxplot(MSH2_counts ~ patient_gender, 
        xlab = "Gender", 
        ylab = "MSH2 Counts", 
        main = "MSH2 Counts Between Men and Women")

boxplot(MLH1_counts ~ patient_gender, 
        xlab = "Gender", 
        ylab = "MLH1 Counts", 
        main = "MLH1 Counts Between Men and Women")

#Sorting patients into low and high counts of MSH2 and MLH1
gene_clinic$MSH2_counts_category = ifelse(MSH2_counts > MSH2_counts_median, "high", "low")
gene_clinic$MLH1_counts_category = ifelse(MLH1_counts > MLH1_counts_median, "high", "low")


#Creating Kaplain-Meier Plot for Gene Counts of MSH2
#Replacing the NAs to days to last followup
gene_clinic$days_to_death = ifelse(is.na(gene_clinic$days_to_death),
                                   gene_clinic$days_to_last_follow_up,
                                   gene_clinic$days_to_death)

gene_clinic$days_to_death = gene_clinic$days_to_death/365

gene_clinic$death_event = as.integer(colData(sum_exp)$vital_status == "Dead")

gene_surv_object <- Surv(time = gene_clinic$days_to_death,
                    event = gene_clinic$death_event)

MSH2_count_fit <- surv_fit(gene_surv_object ~ gene_clinic$MSH2_counts_category, data = gene_clinic)

MSH2_survplot = ggsurvplot(MSH2_count_fit,
                      pval=TRUE,
                      ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")),
                      legend = "bottom")
MSH2_survplot

#Creating Kaplain-Meier Plot for Gene Counts of MLH1
MLH1_count_fit <- surv_fit(gene_surv_object ~ gene_clinic$MLH1_counts_category, data = gene_clinic)

MLH1_survplot = ggsurvplot(MLH1_count_fit,
                           pval=TRUE,
                           ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")),
                           legend = "bottom")
MLH1_survplot

