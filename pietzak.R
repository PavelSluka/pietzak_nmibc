library(tidyverse)
library(ggplot2)

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
plot(c_mutations_by_gene$Rank, c_mutations_by_gene$Percent_Patients_With_Mutation)


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
# temp <- d_individual_mutations[,c("Hugo_Symbol", "Start_Position", "HGVSc", "Tumor_Seq_Allele2")]
# temp1 <- duplicated(temp)
# temp <- sum(duplicated(temp))
# temp


## Part E - Calculate % Unique Patients

# Prepare a blank table
e_cumulative_percent_patients_with_unique_mutations <- as.tibble(data.frame(matrix(nrow=0, ncol=length(colnames(d_individual_mutations))+2 )))
colnames(e_cumulative_percent_patients_with_unique_mutations) <- c(colnames(d_individual_mutations), "Cumulative_Patients", "Percent_Cumulative_Patients")

# Define the ranked mutation to use
for (i in 1:100) {
# for (i in 1:1412) {
# for (i in 1:1406) {
# for (i in c(2, 4:10)) {
# for (i in c(1:10, 14, 37, 55, 56, 61)) {
# for (i in 1:11) {
ranked_mutation <- i


# Take top 1 highest-ranking mutation (d, symbol, )
e_mutations_to_analyse <- d_individual_mutations %>% filter(Rank <= ranked_mutation)
# Take only those values that will be compared to the raw data set
e_mutations_to_analyse <- e_mutations_to_analyse[,c("Hugo_Symbol", "HGVSc", "Start_Position", "Variant_Type", "Tumor_Seq_Allele2", "HGVSp_Short")]
# Extract patients with this mutation in the raw data set
e1_patients_with_matching_mutations <- b_raw_data %>% 
  filter(Hugo_Symbol %in% e_mutations_to_analyse$Hugo_Symbol) %>%
  filter(HGVSc %in% e_mutations_to_analyse$HGVSc) %>%
  filter(Start_Position %in% e_mutations_to_analyse$Start_Position) %>%
  filter(Variant_Type %in% e_mutations_to_analyse$Variant_Type) %>%
  filter(Tumor_Seq_Allele2 %in% e_mutations_to_analyse$Tumor_Seq_Allele2) %>%
  filter(HGVSp_Short %in% e_mutations_to_analyse$HGVSp_Short)

# Count number of unique patients
length(unique(e1_patients_with_matching_mutations$Tumor_Sample_Barcode))

# Place count in a table
e_cumulative_percent_patients_with_unique_mutations <- add_row(e_cumulative_percent_patients_with_unique_mutations,
  Hugo_Symbol = d_individual_mutations$Hugo_Symbol[ranked_mutation],
  HGVSc = d_individual_mutations$HGVSc[ranked_mutation],
  Start_Position = d_individual_mutations$Start_Position[ranked_mutation],
  Variant_Classification = d_individual_mutations$Variant_Classification[ranked_mutation],
  Variant_Type = d_individual_mutations$Variant_Type[ranked_mutation],
  Reference_Allele = d_individual_mutations$Reference_Allele[ranked_mutation],
  Tumor_Seq_Allele2 = d_individual_mutations$Tumor_Seq_Allele2[ranked_mutation],
  HGVSp_Short = d_individual_mutations$HGVSp_Short[ranked_mutation],
  Patients_With_Mutation = d_individual_mutations$Patients_With_Mutation[ranked_mutation],
  Percent_Patients_With_Mutation = d_individual_mutations$Percent_Patients_With_Mutation[ranked_mutation],
  Rank = d_individual_mutations$Rank[ranked_mutation],
  Cumulative_Patients = length(unique(e1_patients_with_matching_mutations$Tumor_Sample_Barcode)),
  Percent_Cumulative_Patients = as.numeric(format(round(length(unique(e1_patients_with_matching_mutations$Tumor_Sample_Barcode)) / 103 * 100,1))))
}

# Plot the rankings

ggplot(e_cumulative_percent_patients_with_unique_mutations, aes(x = e_cumulative_percent_patients_with_unique_mutations$Rank, y = e_cumulative_percent_patients_with_unique_mutations$Percent_Cumulative_Patients)) +
    geom_point() +
    geom_line() +
    scale_x_continuous(breaks = seq(0, ranked_mutation, by = ranked_mutation / 10), limits = c(0, ranked_mutation)) +
    scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(0,100)) +
    #scale_y_continuous(breaks = seq(0, 100, by = 10)) +
    labs(x = "Mutation Rank", y = "Cumulative % of Patients")
  
# Part F - Filter data for other focused analyses

# Filtering raw data to focus analysis
b1_raw_data <- b_raw_data %>% filter(Hugo_Symbol != "TERT")
b_raw_data <- b1_raw_data
b2_raw_data_FGFR3 <- b_raw_data %>% filter(Hugo_Symbol == "FGFR3")
b_raw_data <- b2_raw_data_FGFR3
b3_raw_data_TERT_or_FGFR3 <- b_raw_data %>% filter (Hugo_Symbol == "TERT" | Hugo_Symbol == "FGFR3")
b_raw_data <- b3_raw_data_TERT_or_FGFR3