###############################################################################
# NOTE ON DATA AVAILABILITY
# The input file(s) this script reads contain individual-level proteomic and/or
# clinical data and cannot be shared publicly.
#
###############################################################################

library(tidyverse)
library(ggplot2)
library(PhosR)

#### SAMPLE ID ####
Sample_ID <- read.delim("data-raw/Sample_design_discovery_cohort.txt")
Sample_ID$Gender[grepl("Male",Sample_ID$Group)] <- "Male"
Sample_ID$Gender[grepl("Female",Sample_ID$Group)] <- "Female"
Sample_ID$disease[grepl("T2D",Sample_ID$Group)] <- "T2D"
Sample_ID$disease[grepl("NGT",Sample_ID$Group)] <- "NGT"
Sample_ID$Group <- factor(Sample_ID$Group)
Sample_ID$Subject_ID <- sub("IRS","IRS1",Sample_ID$Subject_ID)
#Sample_ID$New_subject_ID <- sub("Pre", "1", sub("Post", "2", paste0(Sample_ID$Subject_ID, "-", Sample_ID$Clamp)))

#overwriting Clamp ID as all samples were run "Pre" "Post"
Sample_ID$Clamp <- rep(c("Pre","Post"),nrow(Sample_ID)/2)

Sample_ID <- Sample_ID %>%
    mutate(
        New_subject_ID = ifelse(Clamp == "Pre", paste0(Subject_ID, "-1"),
                                paste0(Subject_ID, "-2"))
    )


#filtering sample ID matrix
Sample_ID_filt <- Sample_ID[-c(8,28,127,148),]


#### Clinical data ####
Clinical_data <- read.delim("data-raw/Clinical_data_discovery_cohort.txt", dec=",")
Our_clin_data <- Clinical_data[Clinical_data$Subject.ID %in% Sample_ID$Subject_ID,]
colnames(Our_clin_data)[1] <- "Subject_ID"
Our_clin_data <- Our_clin_data[order(match(Our_clin_data[,1],Sample_ID[,5])),]
duplicated_Our_clin_data <- Our_clin_data %>%
    uncount(weights = 2, .remove = FALSE)

duplicated_Our_clin_data_filt <- duplicated_Our_clin_data[-c(8,28,127,148),]

#### MagNet ####
MagNet_plasma_data <- read.delim("data-raw/Spectronaut_discovery_MagNet.tsv")
MagNet_plasma_data$PG.ProteinGroups <- word(MagNet_plasma_data$PG.ProteinGroups,1,sep=";")
MagNet_plasma_data$PG.Genes <- word(MagNet_plasma_data$PG.Genes,1,sep=";")

# Correcting number labeling
colnames(MagNet_plasma_data) <- sapply(colnames(MagNet_plasma_data), function(name) {
  number <- gsub(".*_S(\\d+)\\.raw.*", "\\1", name)
  padded_number <- sprintf("%03d", as.numeric(number))
  new_name <- gsub("_S\\d+\\.raw", paste0("_S", padded_number, ".raw"), name)
  return(new_name)
})
first_four_cols <- MagNet_plasma_data[, 1:4]

# Extract and sort the remaining columns that need to be sorted
remaining_cols <- MagNet_plasma_data[, 5:ncol(MagNet_plasma_data)]
sorted_remaining_cols <- remaining_cols[, order(colnames(remaining_cols))]

# Combine the first three columns with the sorted remaining columns
MagNet_plasma_data <- cbind(first_four_cols, sorted_remaining_cols)
colnames(MagNet_plasma_data)

sample_order <- MagNet_plasma_data[,-c(1:4)]
sample_order <- sample_order[,order(word(colnames(sample_order),10,sep="_"))]
MagNet_plasma_data_ordered <- cbind(MagNet_plasma_data[,1:4],sample_order)
colnames(MagNet_plasma_data_ordered)[5:ncol(MagNet_plasma_data_ordered)] <- word(word(colnames(MagNet_plasma_data_ordered)[5:ncol(MagNet_plasma_data_ordered)],10,sep="_"),1,sep=".raw")

# split data in discovery and validation cohorts. Run 1:155 are discovery cohort. Rest is validation cohort
#MagNet_discovery_cohort <- MagNet_plasma_data_ordered[,c(1:159)]

MagNet_discovery_cohort_filt <- MagNet_plasma_data_ordered[,-c(12,32,131,152,159)]

colnames(MagNet_discovery_cohort_filt)[5:ncol(MagNet_discovery_cohort_filt)] <- Sample_ID_filt$New_subject_ID

#### Neat PLASMA ####
Neat_plasma_data <- read.delim("data-raw/Spectronaut_discovery_Neat.tsv")
Neat_plasma_data$PG.ProteinGroups <- word(Neat_plasma_data$PG.ProteinGroups,1,sep=";")
Neat_plasma_data$PG.Genes <- word(Neat_plasma_data$PG.Genes,1,sep=";")

# split data in discovery and validation cohorts. Run 1:155 are discovery cohort. Rest is validation cohort. S28 and and S155 are QC pools
#Neat_discovery_cohort <- Neat_plasma_data[,c(1:159)]
Neat_discovery_cohort <- Neat_plasma_data

colnames(Neat_discovery_cohort)[5:159] <- paste0("S",1:155)

#Neat_validation_cohort <- Neat_plasma_data[,c(1:4,160:ncol(Neat_plasma_data))]
#colnames(Neat_validation_cohort)[5:ncol(Neat_validation_cohort)] <- word(word(colnames(Neat_validation_cohort[5:ncol(Neat_validation_cohort)]),10,sep="_"),1,sep=".raw")

old_colnames <- colnames(Neat_discovery_cohort)[5:159]
new_colnames <- sprintf("S%03d", as.numeric(gsub("S", "", old_colnames)))
colnames(Neat_discovery_cohort)[5:159] <- new_colnames

Neat_discovery_cohort_filt <- Neat_discovery_cohort[,-c(12,32,131,152,159)]

colnames(Neat_discovery_cohort_filt)[5:ncol(Neat_discovery_cohort_filt)] <- Sample_ID_filt$New_subject_ID



#### PCA PLASMA ####
PCA_plasma_data <- read.delim("data-raw/Spectronaut_discovery_PCA.tsv")
PCA_plasma_data$PG.ProteinGroups <- word(PCA_plasma_data$PG.ProteinGroups,1,sep=";")
PCA_plasma_data$PG.Genes <- word(PCA_plasma_data$PG.Genes,1,sep=";")

head(PCA_plasma_data$PG.Genes)

# split data in discovery and validation cohorts. Run 1:155 are discovery cohort. Rest is validation cohort. S95 and S155 are QC pools. Note the first 94 samples are from 47 individuals
#plate two we have 61 samples where the last one was a QC pool
#PCA_discovery_cohort <- PCA_plasma_data[,c(1:160)]
PCA_discovery_cohort <- PCA_plasma_data

colnames(PCA_discovery_cohort)[5:ncol(PCA_discovery_cohort)] <- paste0("S",word(word(colnames(PCA_discovery_cohort)[5:ncol(PCA_discovery_cohort)],11,sep="_"),1,sep=".raw"))

PCA_discovery_cohort_filt <- PCA_discovery_cohort[,-c(12,32,99,131,152,160)]

colnames(PCA_discovery_cohort_filt)[5:ncol(PCA_discovery_cohort_filt)] <- Sample_ID_filt$New_subject_ID


#### Depleted depleted PLASMA ####
Depleted_plasma_data <- read.delim("data-raw/Spectronaut_discovery_Depleted.tsv")
Depleted_discovery_cohort <- Depleted_plasma_data
Depleted_discovery_cohort$PG.ProteinGroups <- word(Depleted_discovery_cohort$PG.ProteinGroups,1,sep=";")
Depleted_discovery_cohort$PG.Genes <- word(Depleted_discovery_cohort$PG.Genes,1,sep=";")


Depleted_discovery_cohort_filt <- Depleted_discovery_cohort [,-c(5,6,14,34,79:82,137,158,165:166)]
colnames(Depleted_discovery_cohort_filt)[5:ncol(Depleted_discovery_cohort_filt)] <- word(word(colnames(Depleted_discovery_cohort_filt)[5:ncol(Depleted_discovery_cohort_filt)],11,sep="_"),1,sep=".raw")
colnames(Depleted_discovery_cohort_filt)[5:ncol(Depleted_discovery_cohort_filt)] <- Sample_ID_filt$New_subject_ID

Depleted_plasma_data$PG.ProteinGroups
#### Olink ####
library(data.table)
Olink_data <- read.delim("data-raw/Olink_NPX_discovery_cohort.csv", sep=",")
# one sample was mislabelled and is relabelled here (identifier withheld)
Olink_data[which(Olink_data$SampleID == "PARTICIPANT_D1"),]$SampleID <- "PARTICIPANT_D1_RELABELLED"
Olink_data_all <- Olink_data[which(grepl("assay",Olink_data$AssayType)),]
Olink_data_all <- Olink_data_all[Olink_data_all$SampleType == "SAMPLE",]
Olink_data_all$Gene_ID <- paste0(Olink_data_all$Assay,"_",Olink_data_all$Block)

setDT(Olink_data_all)
Olink_data_all[, ID := paste(Assay, Block, UniProt, sep = "_")]
transposed_Olink_data <- dcast(Olink_data_all, ID ~ SampleID, value.var = "PCNormalizedNPX")

#transposed_Olink_data <- dcast(Olink_data_all, Gene_ID ~ SampleID, value.var = "PCNormalizedNPX")
transposed_Olink_data <- as.data.frame(transposed_Olink_data)

ordered_subjects <- as.character(Sample_ID$New_subject_ID)
ordered_subjects <- intersect(ordered_subjects, colnames(transposed_Olink_data)) # extracting intersection

# Reorder columns of transposed_data, keeping 'Gene_ID' as the first column
#transposed_Olink_data <- transposed_Olink_data[, c("Gene_ID", ordered_subjects)]
transposed_Olink_data <- transposed_Olink_data[, c("ID", ordered_subjects)]

Olink_discovery_cohort_filt <- transposed_Olink_data[,-c(9,29,128,149)]
rownames(Olink_discovery_cohort_filt) <- Olink_discovery_cohort_filt$ID


#### Creating a combined dataset ####
# List of datasets
datasets <- list(MagNet_discovery_cohort_filt, Neat_discovery_cohort_filt, PCA_discovery_cohort_filt, Depleted_discovery_cohort_filt)
processed_datasets <- list()

# Loop through each dataset and apply transformations
for (i in seq_along(datasets)) {
    dataset <- datasets[[i]]

    # Step 1: Extract raw data (columns 5 to 154)
    raw_data <- as.matrix(dataset[, 5:154])
    rownames(raw_data) <- paste0(dataset[,2],"_",dataset[,3])

    # Step 3: Log2 transformation with NA handling
    log2_data <- apply(raw_data , 2, function(x) ifelse(is.na(x), NA, log2(x)))

    # Step 4: Apply selectOverallPercent to filter based on 20% valid values
    filtered_data <- selectOverallPercent(log2_data, percent = 0.2)

    scaled_data <- medianScaling(filtered_data)

    # Store all stages of processed data in the list
    processed_datasets[[i]] <- list(raw = raw_data, log2 = log2_data, filtered = filtered_data, scaled = scaled_data)
}

# Assign names to each dataset for easy access
names(processed_datasets) <- c("MagNet_discovery_cohort", "Neat_discovery_cohort", "PCA_discovery_cohort", "Depleted_discovery_cohort")
processed_datasets$Olink_discovery_cohort <- Olink_discovery_cohort_filt

# The plate info used to corrrect batch effects
Plate_info <- list(MagNet = factor(c(rep(1,94),rep(2,56)), levels=c(1,2)), Neat = factor(c(rep(1,94),rep(2,56)), levels=c(1,2)), PCA = factor(c(rep(1,92),rep(2,58)), levels=c(1,2)))

processed_datasets$MagNet_discovery_cohort$corrected <- medianScaling(removeBatchEffect(processed_datasets$MagNet_discovery_cohort$filtered,Plate_info$MagNet))
processed_datasets$Neat_discovery_cohort$corrected <- medianScaling(removeBatchEffect(processed_datasets$Neat_discovery_cohort$filtered,Plate_info$Neat))
processed_datasets$PCA_discovery_cohort$corrected <- medianScaling(removeBatchEffect(processed_datasets$PCA_discovery_cohort$filtered,Plate_info$PCA))
processed_datasets$Depleted_discovery_cohort$corrected <- processed_datasets$Depleted_discovery_cohort$scaled

processed_datasets$MagNet_discovery_cohort$complete_corrected <- medianScaling(removeBatchEffect(selectOverallPercent(processed_datasets$MagNet_discovery_cohort$log2,1),Plate_info$MagNet))
processed_datasets$Neat_discovery_cohort$complete_corrected <- medianScaling(removeBatchEffect(selectOverallPercent(processed_datasets$Neat_discovery_cohort$log2,1),Plate_info$Neat))
processed_datasets$PCA_discovery_cohort$complete_corrected <- medianScaling(removeBatchEffect(selectOverallPercent(processed_datasets$PCA_discovery_cohort$log2,1),Plate_info$PCA))
processed_datasets$Depleted_discovery_cohort$complete_corrected <- medianScaling(selectOverallPercent(processed_datasets$Depleted_discovery_cohort$log2,1))


rownames(processed_datasets$PCA_discovery_cohort$raw)[829]
rownames(processed_datasets$Depleted_discovery_cohort$raw)[829]

rownames(processed_datasets$Depleted_discovery_cohort$complete_corrected)
rownames(PCA_discovery_cohort_filt)

# Save the processed data into "data" folder.
saveRDS(processed_datasets, file = "data/processed_datasets.rds")
saveRDS(Sample_ID_filt, file = "data/Sample_ID_filt.rds")
saveRDS(duplicated_Our_clin_data_filt, file = "data/Clinical_data_filt")
