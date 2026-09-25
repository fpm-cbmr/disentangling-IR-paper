####### validation outlier/ppols #######
library(data.table)
library(limma)
library(PhosR)
#Neat: one Pre sample excluded, sample 57 (column 61 in the report)
#MagNet: pool in H9, outlier in sample H8
#PCA: one Post sample excluded, sample 14 (column 18 in the report)
#Depleted: one Pre sample excluded, sample 9 (column 13 in the report)

#### validation raw file processing ####
Val_Clinical_data <- read.delim("data-raw/Clinical_data_exercise_cohort.txt", dec=",")
Val_Clinical_data <- Val_Clinical_data[!(Val_Clinical_data$New_ID == ""),]
Val_proteomics_design <- read.delim("data-raw/Sample_design_exercise_cohort.txt")
Val_proteomics_design[which(Val_proteomics_design$Before.After == "B"),"Before.After"] <- "Pre"
Val_proteomics_design[which(Val_proteomics_design$Before.After == "A"),"Before.After"] <- "Post"
Val_Clinical_data$ID_Cond <- paste0(Val_Clinical_data$ID,"_",Val_Clinical_data$Condition)
Val_proteomics_design$ID_Cond <- paste0(Val_proteomics_design$Subject.ID,"_",Val_proteomics_design$Before.After)
Val_Clinical_data <- Val_Clinical_data[match(Val_proteomics_design$ID_Cond,Val_Clinical_data$ID_Cond),] ### order clinical data as the order of plasma analysis
# removing one technical outlier sample (row 92); see Methods
Val_Clinical_data <- Val_Clinical_data[-c(92),]
Val_Clinical_data_Olink <- Val_Clinical_data
Val_Clinical_data <- Val_Clinical_data[-c(which(Val_Clinical_data$ID == "PARTICIPANT_1" & Val_Clinical_data$Condition == "Post"),which(Val_Clinical_data$ID == "PARTICIPANT_2" & Val_Clinical_data$Condition == "Pre"),
which(Val_Clinical_data$ID == "PARTICIPANT_3" & Val_Clinical_data$Condition == "Pre")),]

#### Neat ####
Neat_validation_proteome_data <- read.delim("data-raw/Spectronaut_validation_Neat.tsv")
Neat_validation_proteome_data <- Neat_validation_proteome_data[,c(1:95)]
Neat_validation_proteome_data <- Neat_validation_proteome_data[,-c(13,18,61)]

colnames(Neat_validation_proteome_data)[5:ncol(Neat_validation_proteome_data)] <- paste0(Val_Clinical_data$ID,"_",Val_Clinical_data$Condition)

valid_counts <- colSums(!is.na(Neat_validation_proteome_data[, 5:ncol(Neat_validation_proteome_data)]))
valid_counts_df <- data.frame(Sample = names(valid_counts), Valid_Proteins = valid_counts)

ggplot(valid_counts_df, aes(x = Sample, y = Valid_Proteins)) +
    geom_bar(stat = "identity") +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size = 14, color="black"),  # Increase y-axis text size
          axis.title.x = element_text(size = 16), # Increase x-axis title text size
          axis.title.y = element_text(size = 16),
          plot.title = element_text(size = 18, hjust = 0.5)) +
    labs(title = "MagNet validation cohort",
         x = "Sample", y = "Number of Proteins") + geom_hline(yintercept = 1606, color="red")


mean(valid_counts_df$Valid_Proteins) # 907 proteins identified on average for discovery cohort

view(valid_counts_df)

#### MagNet ####
MagNet_validation_proteome_data <- read.delim("data-raw/Spectronaut_validation_MagNet.tsv")
MagNet_validation_proteome_data <- MagNet_validation_proteome_data[,c(1:95)]
MagNet_validation_proteome_data <- MagNet_validation_proteome_data[,-c(13,18,61)]

colnames(MagNet_validation_proteome_data)[5:ncol(MagNet_validation_proteome_data)] <- paste0(Val_Clinical_data$ID,"_",Val_Clinical_data$Condition)

valid_counts <- colSums(!is.na(MagNet_validation_proteome_data[, 5:ncol(MagNet_validation_proteome_data)]))
valid_counts_df <- data.frame(Sample = names(valid_counts), Valid_Proteins = valid_counts)

ggplot(valid_counts_df, aes(x = Sample, y = Valid_Proteins)) +
    geom_bar(stat = "identity") +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size = 14, color="black"),  # Increase y-axis text size
          axis.title.x = element_text(size = 16), # Increase x-axis title text size
          axis.title.y = element_text(size = 16),
          plot.title = element_text(size = 18, hjust = 0.5)) +
    labs(title = "MagNet validation cohort",
         x = "Sample", y = "Number of Proteins") + geom_hline(yintercept = 1606, color="red")


mean(valid_counts_df$Valid_Proteins) # 1606 proteins identified on average for discovery cohort


#### PCA ####
PCA_validation_proteome_data <- read.delim("data-raw/Spectronaut_validation_PCA.tsv")
PCA_validation_proteome_data <- PCA_validation_proteome_data[,c(1:95)]
PCA_validation_proteome_data <- PCA_validation_proteome_data[,-c(13,18,61)]
colnames(PCA_validation_proteome_data)[5:ncol(PCA_validation_proteome_data)] <- paste0(Val_Clinical_data$ID,"_",Val_Clinical_data$Condition)


valid_counts <- colSums(!is.na(PCA_validation_proteome_data[, 5:ncol(PCA_validation_proteome_data)]))
valid_counts_df <- data.frame(Sample = names(valid_counts), Valid_Proteins = valid_counts)

ggplot(valid_counts_df, aes(x = Sample, y = Valid_Proteins)) +
    geom_bar(stat = "identity") +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size = 14, color="black"),  # Increase y-axis text size
          axis.title.x = element_text(size = 16), # Increase x-axis title text size
          axis.title.y = element_text(size = 16),
          plot.title = element_text(size = 18, hjust = 0.5)) +
    labs(title = "MagNet validation cohort",
         x = "Sample", y = "Number of Proteins") + geom_hline(yintercept = 1606, color="red")


mean(valid_counts_df$Valid_Proteins) # 1606 proteins identified on average for discovery cohort


#### Depleted ####
Depleted_validation_proteome_data <- read.delim("data-raw/Spectronaut_validation_Depleted.tsv")
Depleted_validation_proteome_data <- Depleted_validation_proteome_data[,!grepl("pool",colnames(Depleted_validation_proteome_data))]
Depleted_validation_proteome_data <- Depleted_validation_proteome_data[,-96]
Depleted_validation_proteome_data <- Depleted_validation_proteome_data[,-c(13,18,61)]


colnames(Depleted_validation_proteome_data)[5:ncol(Depleted_validation_proteome_data)] <- paste0(Val_Clinical_data$ID,"_",Val_Clinical_data$Condition)



valid_counts <- colSums(!is.na(Depleted_validation_proteome_data[, 5:ncol(Depleted_validation_proteome_data)]))
valid_counts_df <- data.frame(Sample = names(valid_counts), Valid_Proteins = valid_counts)

ggplot(valid_counts_df, aes(x = Sample, y = Valid_Proteins)) +
    geom_bar(stat = "identity") +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size = 14, color="black"),  # Increase y-axis text size
          axis.title.x = element_text(size = 16), # Increase x-axis title text size
          axis.title.y = element_text(size = 16),
          plot.title = element_text(size = 18, hjust = 0.5)) +
    labs(title = "MagNet validation cohort",
         x = "Sample", y = "Number of Proteins") + geom_hline(yintercept = 1606, color="red")


mean(valid_counts_df$Valid_Proteins) # 1606 proteins identified on average for discovery cohort




#### Olink ####
Val_Clinical_data_Olink <- Val_Clinical_data_Olink[-c(which(Val_Clinical_data_Olink$ID == "PARTICIPANT_1"),89:91),] # one participant and four single samples were not measured by Olink

Olink_data <- read.delim("data-raw/Olink_NPX_validation_cohorts.csv", sep=",")

# Keep only protein assay rows
Olink_data_all <- Olink_data[grepl("assay", Olink_data$AssayType), ]
Olink_data_all <- as.data.table(Olink_data_all)
Olink_data_all[, ID := paste(Assay, Block, UniProt, sep = "_")]

# Create unique SampleID to match Olink naming
Val_Clinical_data_Olink$match_col <- paste0(
    Val_Clinical_data_Olink$ID, "_",
    ifelse(Val_Clinical_data_Olink$Condition == "Post", "A", "B")
)

# 4. Filter to samples present in metadata and non-missing
Olink_data_all <- Olink_data_all[SampleID %in% Val_Clinical_data_Olink$match_col & !is.na(SampleID)]
Olink_data_all[, SampleID := factor(SampleID, levels = Val_Clinical_data_Olink$match_col)]

transposed_Olink_data <- dcast(
    Olink_data_all,
    ID ~ SampleID,
    value.var = "PCNormalizedNPX"
)

ordered_samples <- Val_Clinical_data_Olink$match_col[
    Val_Clinical_data_Olink$match_col %in% colnames(transposed_Olink_data)
]

transposed_Olink_data <- transposed_Olink_data[, c("ID", ordered_samples), with = FALSE]

Olink_validation_proteome_data <- as.matrix(transposed_Olink_data[, -1])
rownames(Olink_validation_proteome_data) <- transposed_Olink_data$ID





#### Creating a combined dataset ####
# List of datasets
datasets <- list(Neat_validation_proteome_data, MagNet_validation_proteome_data, PCA_validation_proteome_data, Depleted_validation_proteome_data)
processed_datasets_validation <- list()

# Loop through each dataset and apply transformations
for (i in seq_along(datasets)) {
    dataset <- datasets[[i]]

    # Step 1: Extract raw data (columns 5 to 154)
    raw_data <- as.matrix(dataset[, 5:92])
    rownames(raw_data) <- paste0(dataset[,2],"_",dataset[,3])

    # Step 3: Log2 transformation with NA handling
    log2_data <- apply(raw_data , 2, function(x) ifelse(is.na(x), NA, log2(x)))

    # Step 4: Apply selectOverallPercent to filter based on 20% valid values
    filtered_data <- selectOverallPercent(log2_data, percent = 0.2)

    scaled_data <- medianScaling(filtered_data)

    # Store all stages of processed data in the list
    processed_datasets_validation[[i]] <- list(raw = raw_data, log2 = log2_data, filtered = filtered_data, scaled = scaled_data)
}

# Assign names to each dataset for easy access
names(processed_datasets_validation) <- c("Neat_validation_cohort", "MagNet_validation_cohort", "PCA_validation_cohort", "Depleted_validation_cohort")
processed_datasets_validation$Olink_validation_cohort <- Olink_validation_proteome_data

processed_datasets_validation$Neat_validation_cohort$complete_corrected <- medianScaling(selectOverallPercent(processed_datasets_validation$Neat_validation_cohort$log2,1))
processed_datasets_validation$MagNet_validation_cohort$complete_corrected <- medianScaling(selectOverallPercent(processed_datasets_validation$MagNet_validation_cohort$log2,1))
processed_datasets_validation$PCA_validation_cohort$complete_corrected <- medianScaling(selectOverallPercent(processed_datasets_validation$PCA_validation_cohort$log2,1))
processed_datasets_validation$Depleted_validation_cohort$complete_corrected <- medianScaling(selectOverallPercent(processed_datasets_validation$Depleted_validation_cohort$log2,1))


# Save the processed data into "data" folder.
saveRDS(processed_datasets_validation, file = "data/processed_datasets_validation.rds")
saveRDS(Val_Clinical_data, file = "data/Val_Clinical_data.rds")
saveRDS(Val_Clinical_data_Olink, file = "data/Val_Clinical_data_Olink.rds")
