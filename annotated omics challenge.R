## Install/load packages
install.packages(c("httr", "jsonlite"))  # Install required packages
library(httr)  # Load the httr package for HTTP requests
library(jsonlite)  # Load the jsonlite package for JSON manipulation
library(readr)  # Load the readr package for data reading

## API data download for cell frequency
res <- GET("https://www.cmi-pb.org/api/v4/pbmc_cell_frequency")  # Retrieve data from the specified API endpoint
dat <- fromJSON(rawToChar(res$content))  # Parse JSON content into R data structures

## Downloading RDS training data for combined data
harmized_dat <- readRDS(gzcon(url("https://www.cmi-pb.org/downloads/cmipb_challenge_datasets/current/2nd_challenge/processed_datasets/training_dataset/master_harmonized_training_data.RDS")))  # Read RDS data from the specified URL
processed_dat <- readRDS(gzcon(url("https://www.cmi-pb.org/downloads/cmipb_challenge_datasets/current/2nd_challenge/processed_datasets/training_dataset/master_processed_training_data.RDS")))  # Read RDS data from the specified URL

## Testing data
# Read RDS data from the specified URL for prediction
processed_test_dat <- readRDS(gzcon(url("https://www.cmi-pb.org/downloads/cmipb_challenge_datasets/current/2nd_challenge/processed_datasets/prediction_dataset/master_processed_prediction_data.RDS")))

# Extract and manipulate dataframes from processed_test_dat
# Converting matrices to data frames and adding specimen_id as a column
plasma_t <- processed_test_dat$abtiter$processed_similar_to_training
t_plasma_t <- as.data.frame(t(plasma_t))
t_plasma_t$specimen_id <- rownames(t_plasma_t)

gene_t <- processed_test_dat$pbmc_gene_expression$processed_similar_to_training
t_gene_t <- as.data.frame(t(gene_t))
t_gene_t$specimen_id <- rownames(t_gene_t)

freq_t <- processed_test_dat$pbmc_cell_frequency$processed_similar_to_training
t_freq_t <- as.data.frame(t(freq_t))
t_freq_t$specimen_id <- rownames(t_freq_t)

cyto_t <- processed_test_dat$plasma_cytokine_concentrations$processed_similar_to_training
t_cyto_t <- as.data.frame(t(cyto_t))
t_cyto_t$specimen_id <- rownames(t_cyto_t)

# Merging data frames
processed_test_dat2 <- merge(t_plasma_t, t_gene_t, by = "specimen_id", all = TRUE)
processed_test_dat2 <- merge(processed_test_dat2, t_freq_t, by = "specimen_id", all = TRUE)
processed_test_dat2 <- merge(processed_test_dat2, t_cyto_t, by = "specimen_id", all = TRUE)

# Writing data frames to CSV files
write.csv(plasma_t, file = "C:/Users/jans0/Documents/plasma_antibody_levels_test.csv")
write.csv(gene_t, file = "C:/Users/jans0/Documents/pbmc_gene_expression_test.csv")
write.csv(freq_t, file = "C:/Users/jans0/Documents/pbmc_cell_frequency_test.csv")
write.csv(cyto_t, file = "C:/Users/jans0/Documents/plasma_cytokine_concentrations_test.csv")

# Read TSV files from the specified URLs
cell_freq_dat <- readr::read_tsv("https://www.cmi-pb.org/downloads/cmipb_challenge_datasets/current/2nd_challenge/processed_datasets/prediction_dataset/pbmc_cell_frequency_processed_data.tsv")
gene_expr_dat <- readr::read_tsv("https://www.cmi-pb.org/downloads/cmipb_challenge_datasets/current/2nd_challenge/processed_datasets/prediction_dataset/pbmc_gene_expression_processed_data.tsv")
cytokine_dat <- readr::read_tsv("https://www.cmi-pb.org/downloads/cmipb_challenge_datasets/current/2nd_challenge/processed_datasets/prediction_dataset/plasma_cytokine_concentrations_processed_data.tsv")
subject_dat <- readr::read_tsv("https://www.cmi-pb.org/downloads/cmipb_challenge_datasets/current/2nd_challenge/processed_datasets/prediction_dataset/subject_specimen.tsv")

# Downloading TSV files for cell frequency
pred_dat <- readr::read_tsv("https://www.cmi-pb.org/downloads/cmipb_challenge_datasets/current/2nd_challenge/raw_datasets/prediction_data/2022BD_pbmc_cell_frequency.tsv")
training_2020_dat <- readr::read_tsv("https://www.cmi-pb.org/downloads/cmipb_challenge_datasets/current/2nd_challenge/raw_datasets/training_data/2020LD_pbmc_cell_frequency.tsv")
training_2021_dat <- readr::read_tsv("https://www.cmi-pb.org/downloads/cmipb_challenge_datasets/current/2nd_challenge/raw_datasets/training_data/2021LD_pbmc_cell_frequency.tsv")

# Extracting names of columns from the 'harmized_dat' dataframe
plasma <- harmized_dat$plasma_antibody_levels$wide
plasma_antibody_levels <- names(plasma)
gene <- harmized_dat$pbmc_gene_expression_wide$wide
pbmc_gene_expression <- names(gene)
freq <- harmized_dat$pbmc_cell_frequency_wide$wide
pbmc_cell_frequency <- names(freq)
cyto <- harmized_dat$plasma_cytokine_concentrations$wide
plasma_cytokine_concentrations <- names(cyto)

# Writing column names to CSV files
write.csv(plasma_antibody_levels, file = "C:/Users/jans0/Documents/plasma_antibody_levels_names.csv")
write.csv(pbmc_gene_expression, file = "C:/Users/jans0/Documents/pbmc_gene_expression.csv")
write.csv(pbmc_cell_frequency, file = "C:/Users/jans0/Documents/pbmc_cell_frequency.csv")
write.csv(plasma_cytokine_concentrations, file = "C:/Users/jans0/Documents/plasma_cytokine_concentrations_names.csv")

# Further data manipulation and analysis
# Performing imputation
library(missMDA)  # Load the missMDA package for multiple imputation

# Reading RDS files
processed_dat2 <- readRDS("C:/Users/jans0/Downloads/merged_day_0.rds")
proc_test <- readRDS("C:/Users/jans0/Downloads/merged_df.rds")

# Performing imputation with missMDA
# Estimating the number of components for PCA
# Subset the data from processed_dat2 based on specified subject_id values
imp_set <- processed_dat2 |> 
  filter(subject_id %in% c(1,3,4,5,6,9,10,11,13,15,17,18,19,20,21,22,23,24,25,26,27,29,31,32,33,35,36,38,42,43,44,47,48,50,52,53,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96))

# Estimate the number of components for PCA using k-fold cross-validation
imput_est <- missMDA::estim_ncpPCA(imp_set[, c(3:29, 39:8335)],
                                   method.cv = "kfold",
                                   verbose = FALSE)

# Perform multiple imputation using PCA
imp_dat <- missMDA::imputePCA(imp_set[, c(3:29, 39:8335)], ncp = imput_est$ncp)

# Combine imputed data with original data
comp_dat <- cbind.data.frame(imp_set[, c(1:2, 30:38)], imp_dat$completeObs)

# Write the imputed data to a CSV file
write.csv(comp_dat, file = "C:/Users/jans0/Documents/imputed_data_72.csv")

## Test section

# Perform imputation on processed_test_dat2 excluding specified columns
imp_test <- missMDA::imputePCA(processed_test_dat2[, -c(1, 8323:8335)], ncp = 5)

# Combine imputed test data with original data
comp_test2 <- cbind.data.frame(processed_test_dat2[, c(1, 8323:8335)], imp_test$completeObs)

# Extract specific columns from comp_test2 where timepoint is 0
one <- comp_test2[comp_test2$timepoint == 0, c("subject_id", "IgG_PT")]

# Sort the data by IgG_PT values
one[order(one$IgG_PT),]

# Extract specific columns from comp_test2 where timepoint is 0 or -15
two <- comp_test2[comp_test2$timepoint == 0 | comp_test2$timepoint == -15, c("subject_id", "IgG_PT", "timepoint")]
# Spread the data based on timepoint values
two <- two |> 
  spread(timepoint, IgG_PT)
# Rename columns
names(two)[2:3] <- c("before", "after")
# Calculate change
two$change <- two$before - two$after
# Sort the data by change
two[order(two$change),]

# Further data extraction and manipulation
# Extract specific columns from comp_test2 where timepoint is 0
five <- comp_test2[comp_test2$timepoint == 0, c("subject_id", "ENSG00000277632.1")]
# Order the data by ENSG00000277632.1 values
five <- five[order(five$ENSG00000277632.1),]
# Add rank column
five$rank <- 1:21
# Order the data by subject_id
five[order(five$subject_id),]

# Extract specific columns from comp_test2 where timepoint is 0 or -15
six <- comp_test2[comp_test2$timepoint == 0 | comp_test2$timepoint == -15, c("subject_id", "ENSG00000277632.1", "timepoint")]
# Spread the data based on timepoint values
six <- six |> 
  spread(timepoint, ENSG00000277632.1)
# Rename columns
names(six)[2:3] <- c("before", "after")
# Calculate change
six$change <- six$before - six$after
# Order the data by change
six <- six[order(six$change),]
# Add rank column
six$rank <- 1:21
# Order the data by subject_id
six[order(six$subject_id),]

# Write the test data to a CSV file
write.csv(comp_test, file = "C:/Users/jans0/Documents/imputed_test_63.csv")

# Visualize the data using ggplot2
# Convert train_set column to a factor
processed_dat2$train_set <- as.factor(processed_dat2$train_set)
# Create a histogram for Monocytes_1 colored by train_set
ggplot2::ggplot(processed_dat2, aes(x = Monocytes_1, color = train_set)) + geom_histogram()

# Create a histogram for Monocytes_1 colored by train_set using the imputed data
ggplot2::ggplot(comp_dat, aes(x = Monocytes_1, color = as.factor(train_set))) + geom_histogram()
