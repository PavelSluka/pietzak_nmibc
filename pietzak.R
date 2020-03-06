library(tidyverse)

# Set working directory to where the raw data file is
setwd("C:/Users/pavels/Git/pietzak_nmibc/")


# Import raw data table
a_raw_data <- as_tibble(read.csv("data_mutations_extended.txt", header = TRUE, sep = "\t"))


# Extract only columns we need
b_raw_data <- a_raw_data[,c(1,5,6,7,10,11,17,38:44)]


## Count number of patients with a mutation in each gene

# Sort rows by 1) gene name and 2) barcode (patient ID)
# Will show patients with multiple mutations in the same gene in adjacent rows
c_mutations_by_gene <- b_raw_data[order(b_raw_data$Hugo_Symbol, b_raw_data$Tumor_Sample_Barcode),]

# Count number of unique patients with a mutation in each gene
c1 <- c_mutations_by_gene %>% 
  group_by(Hugo_Symbol) %>% 
  summarize(Patients_With_Mutation = n_distinct(Tumor_Sample_Barcode))

# Sort table to show genes with highest number of patients with mutation at top 
c1 <- arrange(c1, desc(Patients_With_Mutation))

# Calculate percentage of patients with mutation for each gene
c1 <- c1 %>% mutate(Percent_Patients_With_Mutation = format(round((Patients_With_Mutation / 103 * 100), 1)))

# Plot Highest to lowest % Patents with mutation for each gene
c1 <- c1 %>% mutate(rank = c(1:320))
plot(c1$rank, c1$Percent_Patients_With_Mutation)


## Count number of patients with each unique mutation

# Sort rows by mutation
d_individual_mutations <- b_raw_data[order(b_raw_data$HGVSc, b_raw_data$Start_Position, b_raw_data$Tumor_Sample_Barcode, b_raw_data$Hugo_Symbol),]

# Count number of patients with each unique mutation
d1 <- d_individual_mutations %>%
  group_by(Hugo_Symbol, HGVSc, Start_Position) %>%
  summarize(Patients_With_Mutation = n_distinct(Tumor_Sample_Barcode))

# Sort table to show genes with highest number of patients with mutation at top 
d1 <- arrange(d1, desc(Patients_With_Mutation))

# Calculate percentage of patients with each unique mutation
d1 <- d1 %>% mutate(Percent_Patients_With_Mutation = format(round((Patients_With_Mutation / 103 * 100), 1)))



# Test if there are duplicates in the raw data
a1_test_for_dups <- a_raw_data[order(a_raw_data$Hugo_Symbol, a_raw_data$Tumor_Sample_Barcode),]
a1a <- duplicated(a1_test_for_dups)
a1a
