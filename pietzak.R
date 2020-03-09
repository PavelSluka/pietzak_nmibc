library(tidyverse)

# Set working directory to where the raw data file is
setwd("C:/Users/pavels/Git/pietzak_nmibc/")


# Import raw data table
a_raw_data <- as_tibble(read.csv("data_mutations_extended.txt", header = TRUE, sep = "\t"))


# Test if there are duplicates in the raw data
a1_test_for_dups <- a_raw_data[order(a_raw_data$Hugo_Symbol, a_raw_data$Tumor_Sample_Barcode),]
a1_test_for_dups <- duplicated(a1_test_for_dups)
sum(a1_test_for_dups) # Should be zero
sum(!a1_test_for_dups) # Should be equal to the number of samples


# Extract only columns we need
b_raw_data <- a_raw_data[,c(1,5,6,7,10:14,17,38:44)]


## PART C - Count number of patients with a mutation in each gene

# Sort rows by 1) gene name and 2) barcode (patient ID)
# Will show patients with multiple mutations in the same gene in adjacent rows
c_mutations_by_gene <- b_raw_data[order(b_raw_data$Hugo_Symbol, b_raw_data$Tumor_Sample_Barcode),]

# Count number of unique patients with a mutation in each gene
c_mutations_by_gene <- c_mutations_by_gene %>% 
  group_by(Hugo_Symbol) %>% 
  summarize(Patients_With_Mutation = n_distinct(Tumor_Sample_Barcode))

# Sort table to show genes with highest number of patients with mutation at top 
c_mutations_by_gene <- arrange(c_mutations_by_gene, desc(Patients_With_Mutation))

# Calculate percentage of patients with mutation for each gene
c_mutations_by_gene <- c_mutations_by_gene %>% mutate(Percent_Patients_With_Mutation = format(round((Patients_With_Mutation / 103 * 100), 1)))

# Plot Highest to lowest % Patents with mutation for each gene
c_mutations_by_gene$Rank <- c(1:nrow(c_mutations_by_gene)) # Add a ranking column
plot(c_mutations_by_gene$rank, c_mutations_by_gene$Percent_Patients_With_Mutation)


## PART D - Count number of patients with each unique mutation

# Sort rows by mutation
d_individual_mutations <- b_raw_data[order(b_raw_data$HGVSc, b_raw_data$Start_Position, b_raw_data$Tumor_Sample_Barcode, b_raw_data$Hugo_Symbol, b_raw_data$Reference_Allele, b_raw_data$Tumor_Seq_Allele2, b_raw_data$HGVSp_Short),]

# Count number of patients with each unique mutation
d_individual_mutations <- d_individual_mutations %>%
  group_by(Hugo_Symbol, HGVSc, Start_Position, Variant_Classification, Variant_Type, Reference_Allele, Tumor_Seq_Allele2, HGVSp_Short ) %>%
  summarize(Patients_With_Mutation = n_distinct(Tumor_Sample_Barcode))

# Sort table to show genes with highest number of patients with mutation at top 
d_individual_mutations <- arrange(d_individual_mutations, desc(Patients_With_Mutation))

# Calculate percentage of patients with each unique mutation
d_individual_mutations <- d_individual_mutations %>% mutate(Percent_Patients_With_Mutation = format(round((Patients_With_Mutation / 103 * 100), 1)))

# Add a ranking column
d_individual_mutations$Rank <- c(1:nrow(d_individual_mutations))

#### Test for duplicates in individual mutation
## Identify columns that identify duplicate mutations (HUGO [Gene], Start site [Mutation location], HGSV [Actual Mutation])
temp <- d_individual_mutations[,c("Hugo_Symbol", "Start_Position", "HGVSc", "Tumor_Seq_Allele2")]
temp1 <- duplicated(temp)
temp <- sum(duplicated(temp))
temp

## Part E - Calculate % Unique Patients

# Take top 1 highest-ranking mutation (d, symbol, )



