library(tidyverse)

# Set working directory to where the raw data file is

setwd("C:/Users/pavels/Git/pietzak_nmibc/")


# Import raw data table

a_raw_data <- read.csv("data_mutations_extended.txt", header = TRUE, sep = "\t")


# Extract only columns we need

b_raw_data <- a_raw_data[,c(1,5,6,7,10,11,17,38:44)]


## Count number of patients with a mutation in each gene

# Sort rows by 1) gene name and 2) barcode (patient ID)
# Will show patients with multiple mutations in the same gene in adjacent rows

c_mutations_by_gene <- b_raw_data[order(b_raw_data$Hugo_Symbol, b_raw_data$Tumor_Sample_Barcode),]
c_mutations_by_gene <- as_tibble(c_mutations_by_gene)
# Combine 

c1 <- c_mutations_by_gene %>% group_by(Tumor_Sample_Barcode) %>% count(Hugo_Symbol, Tumor_Sample_Barcode)
c1 <- c_mutations_by_gene %>% count(Hugo_Symbol, .drop=FALSE, sort = TRUE)

# Test if there are duplicates in the raw data
a1_test_for_dups <- a_raw_data[order(a_raw_data$Hugo_Symbol, a_raw_data$Tumor_Sample_Barcode),]
a1a <- duplicated(a1_test_for_dups)
a1a
