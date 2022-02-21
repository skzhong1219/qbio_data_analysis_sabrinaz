
# Exercise 1.1
BiocManager::install("DESeq2")
library(DESeq2)
library(TCGAbiolinks)
library(SummarizedExperiment)

query <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Transcriptome Profiling", # get the RNA-seq transcriptome
                  data.type = "Gene Expression Quantification", # gets the counts
                  workflow.type = "HTSeq - Counts") # gets the raw counts processed by this method
sum_exp <- GDCprepare(query)

# Exercise 1.2
is.na(colData(sum_exp)$age_at_index)

patients_data = colData(sum_exp)
counts = assays(sum_exp)$"HTSeq - Counts"

# counts = data.frame(counts)
# 4 patients with NA age

NApatientdata = is.na(patients_data$age_at_index)

counts = counts[,!NApatientdata]
patients_data = patients_data[!NApatientdata,]

patients_data$age_category = ifelse(patients_data$age_at_index < 50, "young", "old")

patients_data$age_category = factor(patients_data$age_category, levels = c("young", "old"))

# Exercise 1.3
if (all(rownames(counts) == names(rowRanges(sum_exp)))){
  print("Changing row names!")
  rownames(counts) = rowRanges(sum_exp)$external_gene_name
}

counts_row_sums = rowSums(counts)

low_counts_mask = ifelse(counts_row_sums < 10, FALSE, TRUE)
sum(low_counts_mask)

counts = counts[low_counts_mask,]

# Exercise 2.1
dds = DESeqDataSetFromMatrix(countData = counts,
                             colData = patients_data,
                             design = ~age_category)

dds_obj = DESeq(dds)
resultsNames(dds_obj)  # see what comparisons got run

# get the young vs. old comparison
results = results(dds_obj, format = "DataFrame", contrast = c("age_category", "young", "old"))

# Exercise 2.2
str(results)
head(results)
my_df = data.frame(x = c('b', 'd', 'c', 'e', 'a'),
                   y = c(2,4,3,5,1))

order_indices = order(my_df$y)
# we expect c(5, 1, 3, 2, 4) because:
# 1 is index 5
# 2 is index 1
# 3 is index 3
# 4 is index 2
# 5 is index 4
order_indices  # note the order!

my_df = my_df[order_indices, ]
my_df

# Exercise 2.3
row_order = order(results$padj)
results = results[row_order, ]
results
head(results)

# choosing SCARNA5 gene
which(rownames(counts) == "SCARNA5")

# This gene is more highly expressed in old patients. 
# The full name of the gene is Small Cajal Body-Specific RNA 5. They are a type of snoRNA, and their
# primary function involves the biogenesis of snRNPs and guiding the modification of RNA polymerase
# II transcribed spliceosomal RNA U5. 


# Exercise 2.4
log2FoldChange_threshold = 1
padj_threshold = 0.05

log2FoldChange_greater = results$log2FoldChange > log2FoldChange_threshold
padj_lesser = results$padj < padj_threshold
results = log2FoldChange_greater & padj_lesser

# Exercise 2.5
fc_threshold = 2  # set a threshold of at least a 2 fold increase (double)
p_threshold = 0.05  # set a threshold of adjusted p-value being <= 0.05

# fill in your plot code here!
# be sure to relabel the axes!
# note: you can perform the log transformation directly in the plot function
plot(x = results$log2FoldChange,
     y = -log10(results$padj),
     xlab = "Log 2 Fold Change (young/old)", # be sure the specify that it's young over old!
     ylab = "-log10 (adjusted) p-value",
     pch = 20) # smaller solid circles

# these lines put the lines on the plot
# abline() plots straight lines on an R plot.
# v argument is for a vertical line, h argument is for a horizontal line, col argument is color
abline(v=c(-log2(fc_threshold), log2(fc_threshold)), h= c(-log10(p_threshold)), col="green")

# Exercise 2.6
write.csv(x = results,
          file = "/Users/sabrinazhong/Desktop/USC/Spring_2022/QBIO490/qbio_data_analysis_sabrinaz/week6_DESeq2/results.csv",
          row.names = FALSE)
