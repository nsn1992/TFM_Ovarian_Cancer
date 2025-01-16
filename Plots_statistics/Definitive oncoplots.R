################# ONCOPLOTS #####################################3

# Load necessary libraries
library(maftools)  # For mutational analysis and visualization
library(data.table)  # For data manipulation
library(readxl)  # To read Excel files
library(tibble)  # For data frames and tibbles
library(dplyr)  # For data manipulation (filter, mutate, etc.)
library(tidyr)  # For data reshaping

# Load the mutation file with the corresponding filter to obtain driver genes
mutations <- read_excel('/home/vant/TFM/High_Filter_Davidaltreads40VAF0.1.xlsx')

# Remove samples with artifacts and mixed case and rename duplicates
mutations <- mutations[!(mutations$Sample %in% c('OVE44', 'OVE16','RVB3')), ]
mutations$Sample <- gsub('LP16_old', 'LP16', mutations$Sample)
mutations$Sample <- gsub('LP17_old', 'LP17', mutations$Sample)
mutations$Sample <- gsub('MDA51T1', 'MDA51_T1_CCC', mutations$Sample)
mutations$Sample <- gsub('MDA51T2', 'MDA51_T2_EC', mutations$Sample)

# Rename columns to match maftools format
colnames(mutations)[colnames(mutations) == "Sample"] <- "Tumor_Sample_Barcode"
colnames(mutations)[colnames(mutations) == "Gene"] <- "Hugo_Symbol"
colnames(mutations)[colnames(mutations) == "Function"] <- "Variant_Classification"
colnames(mutations)[colnames(mutations) == "Pos"] <- "Chromosome"
colnames(mutations)[colnames(mutations) == "Ref"] <- "Reference_Allele"
colnames(mutations)[colnames(mutations) == "Alt"] <- "Tumor_Seq_Allele2"

# Type of varinant in my mutations file
unique(mutations$Variant_Classification)

# Replace values in Variant_Classification to match specific mutation categories recognized by maftools
mutations$Variant_Classification <- gsub("^frameshift_deletion$", "Frame_Shift_Del", mutations$Variant_Classification, ignore.case = TRUE)
mutations$Variant_Classification <- gsub("^missense$", "Missense_Mutation", mutations$Variant_Classification, ignore.case = TRUE)
mutations$Variant_Classification <- gsub("^stopgain$", "Nonsense_Mutation", mutations$Variant_Classification, ignore.case = TRUE)
mutations$Variant_Classification <- gsub("^frameshift_insertion$", "Frame_Shift_Ins", mutations$Variant_Classification, ignore.case = TRUE)
mutations$Variant_Classification <- gsub("^splicing$", "Splice_Site", mutations$Variant_Classification, ignore.case = TRUE)
mutations$Variant_Classification <- gsub("^frameshift_substitution$", "Frame_Shift_Del", mutations$Variant_Classification, ignore.case = TRUE) ##All 'frameshift_subtitution' mutations are deletions
mutations$Variant_Classification <- gsub("^nonframeshift_deletion$", "In_Frame_Del", mutations$Variant_Classification, ignore.case = TRUE)
#mutations$Variant_Classification <- gsub("^startloss$", "Translation_Start_Site", mutations$Variant_Classification, ignore.case = TRUE)

## Check changes
unique(mutations$Variant_Classification)

# Deduce the type of variants (SNP, INS, DEL)
mutations$Variant_Type <- NA
mutations$Variant_Type[mutations$Reference_Allele != mutations$Tumor_Seq_Allele2 &
                         nchar(mutations$Reference_Allele) == 1 &
                         nchar(mutations$Tumor_Seq_Allele2) == 1] <- "SNP"
mutations$Variant_Type[nchar(mutations$Tumor_Seq_Allele2) > nchar(mutations$Reference_Allele)] <- "INS"
mutations$Variant_Type[nchar(mutations$Reference_Allele) > nchar(mutations$Tumor_Seq_Allele2)] <- "DEL"

# Separate chromosome and position information
mutations <- mutations %>%
  separate(Chromosome, into = c("Chromosome", "Position"), sep = ":", remove = TRUE)
mutations$Position <- as.numeric(mutations$Position)

# Function to calculate the final position of a mutation
calculate_end_position <- function(ref, alt, pos) {
  ref_length <- nchar(ref)
  alt_length <- nchar(alt)
  if (alt_length < ref_length) {
    # For deletions
    end_position <- pos + ref_length - 1
  } else {
    # For insertions and SNP
    end_position <- pos + alt_length - 1
  }
  return(end_position)
}

mutations <- mutations %>%
  rowwise() %>%
  mutate(
    Start_Position = Position,
    End_Position = calculate_end_position(Reference_Allele, Tumor_Seq_Allele2, Position)
  ) %>%
  ungroup()


# Load the histology data with additional information
histology_data <- read_excel('/home/vant/TFM/Final_db/Samples_alltypedata_annotated.xlsx')

# Eliminate rows in df_variables with no necessary information for the analysis
histology_data <- histology_data[-c(1, 2), ]

# Change the misclassified samples value in MMR_final_status column (according to molecular classifier)
histology_data <- histology_data %>%
  mutate(MMR_final_status = case_when(
    ID_CNIO %in% c('LP20', 'LP26', 'MDA15') ~ 'MMRp',
    ID_CNIO %in% c('MDA26') ~ 'MMRd',
    TRUE ~ MMR_final_status
  ))
# Identify samples with MMR_final_status == "MMRd"
mmrd_samples <- histology_data %>%
  filter(MMR_final_status == "MMRd") %>%
  pull(ID_CNIO)  # Get the names of the MMRd samples


## EOC SAMPLES
# Filter mutations to match histology data and add the HISTOLOGY column
mutations_filtered <- mutations %>%
  filter(Tumor_Sample_Barcode %in% histology_data$ID_CNIO) %>%
  left_join(histology_data %>% dplyr::select(ID_CNIO, HISTOLOGY), by = c("Tumor_Sample_Barcode" = "ID_CNIO")) %>%
  filter(HISTOLOGY == 0)


# Create a MAF object
maf_object <- tryCatch({
  maf <- read.maf(maf = mutations_filtered)
  cat("MAF object created successfully!\n")
  maf
}, error = function(e) {
  message("Error creating the MAF object: ", e$message)
  NULL
})

# Add a new 'MMR_status' column to the MAF object
if (!is.null(maf_object)) {
  # Add MMR_status column indicating whether the sample is MMRd or not
  maf_object@clinical.data$MMR_status <- ifelse(
    maf_object@clinical.data$Tumor_Sample_Barcode %in% mmrd_samples, 
    "MMRd", 
    "MMRp"
  )
  
  # Define colors for annotations
  annotation_colors <- list(
    MMR_status = c("MMRd" = "red", "MMRp" = "lightgreen")
  )
  
  # Generate the oncoplot
  png(filename = "Oncoplot_EOC_top8.png", width = 3500, height = 1500, res = 300)
  
  oncoplot(
    maf = maf_object,
    top = 8,  # Adjust the number of genes to show
    showTumorSampleBarcodes = TRUE,
    drawRowBar = FALSE,
    drawColBar = FALSE,
    draw_titv = FALSE,
    SampleNamefontSize = 0.8,
    barcode_mar = 12,  # Adjust the margin for the tumor sample barcodes
    sepwd_genes = 0.15,
    clinicalFeatures = "MMR_status",  # Show the MMR_status annotation
    annotationColor = annotation_colors,
    sortByAnnotation = TRUE  # Sort the samples based on the annotation
  )
  dev.off()  # Save the plot to file
  cat("Oncoplot generated successfully with MMRd samples highlighted in red and annotation bar adjusted.\n")
} else {
  cat("Error: MAF object is NULL.\n")
}
# Summary plots
png(filename = "plotmafSummary_EOC.png", width = 3500, height = 2000, res = 300)
plotmafSummary(maf = maf_object, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()
# Transversions-transitions plots results 
png(filename = "plottitv_EOC.png", width = 3500, height = 2000, res = 300)
maf_object_H1.titv = titv(maf = maf_object, plot = FALSE, useSyn = TRUE)
plotTiTv(res = maf_object_H1.titv)
dev.off()
# Somatic interactions analysis (Co-occurrence, exclusivity)
png(filename = "Somaticinteractions_EOC.png", width = 3500, height = 2000, res = 300)
somaticInteractions(maf = maf_object, top = 20, pvalue = c(0.05, 0.1))
dev.off()
# Enriched pathways analysis
png(filename = "Pathways_EOC.png", width = 3500, height = 2000, res = 300)
OncogenicPathways(maf = maf_object)
dev.off()


#### CCOC SAMPLES ##############

# Filter mutations for samples where HISTOLOGY == 1 (proficient MMR)
mutations_filtered_H1 <- mutations %>%
  filter(Tumor_Sample_Barcode %in% histology_data$ID_CNIO) %>%
  left_join(histology_data %>% dplyr::select(ID_CNIO, HISTOLOGY), by = c("Tumor_Sample_Barcode" = "ID_CNIO")) %>%
  filter(HISTOLOGY == 1)

# Create a MAF object for these filtered samples
maf_object_H1 <- tryCatch({
  maf <- read.maf(maf = mutations_filtered_H1)
  cat("MAF object for HISTOLOGY == 1 created successfully!\n")
  maf
}, error = function(e) {
  message("Error creating the MAF object: ", e$message)
  NULL
})

# Add a new 'MMR_status' column to the MAF object for these samples
if (!is.null(maf_object_H1)) {
  # Add MMR_status column indicating whether the sample is MMRd or not
  maf_object_H1@clinical.data$MMR_status <- ifelse(
    maf_object_H1@clinical.data$Tumor_Sample_Barcode %in% mmrd_samples, 
    "MMRd", 
    "MMRp"
  )
  
  # Define colors for annotations
  annotation_colors_H1 <- list(
    MMR_status = c("MMRd" = "red", "MMRp" = "lightgreen")
  )
  
  # Generate the oncoplot for HISTOLOGY == 1
  png(filename = "Oncoplot_CCOC_top8.png", width = 3500, height = 1500, res = 300)
  
  oncoplot(
    maf = maf_object_H1,
    top = 8,  # Adjust the number of genes to show
    showTumorSampleBarcodes = TRUE,
    drawRowBar = FALSE,
    drawColBar = FALSE,
    draw_titv = FALSE,
    SampleNamefontSize = 0.8,
    barcode_mar = 12,  # Adjust the margin for the tumor sample barcodes
    sepwd_genes = 0.15,
    clinicalFeatures = "MMR_status",  # Show the MMR_status annotation
    annotationColor = annotation_colors_H1,
    sortByAnnotation = TRUE  # Sort the samples based on the annotation
  )
  
  dev.off()  # Save the plot to file
  cat("Oncoplot CCOC adjusted.\n")
} else {
  cat("Error: MAF object for HISTOLOGY == 1 is NULL.\n")
}
# Summary plots
png(filename = "plotmafSummary_CCOC.png", width = 3500, height = 2000, res = 300)
plotmafSummary(maf = maf_object_H1, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()
# Transversions-transitions plots results 
png(filename = "plottitv_CCOC.png", width = 3500, height = 2000, res = 300)
maf_object_H1.titv = titv(maf = maf_object_H1, plot = FALSE, useSyn = TRUE)
plotTiTv(res = maf_object_H1.titv)
dev.off()
# Somatic interactions analysis (Co-occurrence, exclusivity)
png(filename = "Somaticinteractions_CCOC.png", width = 3500, height = 2000, res = 300)
somaticInteractions(maf = maf_object_H1, top = 20, pvalue = c(0.05, 0.1))
dev.off()
# Enriched pathways analysis
png(filename = "Pathways_CCOC.png", width = 3500, height = 2000, res = 300)
OncogenicPathways(maf = maf_object_H1)
dev.off()

