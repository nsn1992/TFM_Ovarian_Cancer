
#### Heatmaps/Oncoplot for MMRd cases ##############

# Load libraries for the analysis
library(readxl)
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(patchwork)
library (tibble)

# Load mutations data obtained after specific filtering
df_mutations <- read_excel('/home/vant/TFM/High_Filter_Davidaltreads40VAF0.1.xlsx')

# Check structure of the dataset
head(df_mutations) # Check initial rows
colnames(df_mutations)  # Check column names
unique(df_mutations$Sample) # Unique samples

# Eliminate samples with artefacts in variant calling (OVE16 and OVE44) and mixed case (RVB3) in the mutations file
df_mutations <- df_mutations %>% 
  filter(!Sample %in% c('OVE44', 'OVE16', 'RVB3'))

# Rename specific samples in df_mutations
df_mutations <- df_mutations %>%
  mutate(Sample = recode(Sample,
                         'LP16_old' = 'LP16',
                         'LP17_old' = 'LP17',
                         'MDA51T1' = 'MDA51_T1_CCC',
                         'MDA51T2' = 'MDA51_T2_EC'))

# Load file with additional information 
df_variables <- read_excel('/home/vant/TFM/Final_db/Samples_alltypedata_annotated.xlsx')

# Eliminate rows in df_variables with no necessary information for the analysis
df_variables <- df_variables[-c(1, 2), ]

# Adjust misclassified sample values in the MMR_final_status column based on molecular classifiers
df_variables <- df_variables %>%
  mutate(MMR_final_status = case_when(
    ID_CNIO %in% c('LP20', 'LP26', 'MDA15') ~ 'MMRp',
    ID_CNIO %in% c('MDA26') ~ 'MMRd',
    TRUE ~ MMR_final_status
  ))

# Eliminate the same samples in the variables file
df_variables <- df_variables %>% 
  filter(!ID_CNIO %in% c('OVE44', 'OVE16', 'RVB3'))

# Select and clean columns of interest
columns_of_interest <- c('ID_CNIO', 'HISTOLOGY', 'MMR_final_status')
df_variables_filtered <- df_variables %>%
  select(all_of(columns_of_interest)) %>%
  mutate(ID_CNIO = trimws(toupper(ID_CNIO))) 

# Clean values in df_mutations
df_mutations <- df_mutations %>%
  mutate(Sample = trimws(toupper(Sample)))

# Verify common samples between df_mutations and df_variables_filtered
common_samples <- intersect(df_mutations$Sample, df_variables_filtered$ID_CNIO)
cat("Número de muestras comunes entre df_mutations y df_variables:", length(common_samples), "\n")

# Melt DataFrames based on the column with common values (Sample and ID_CNIO). Keep only common samples
df_combined <- inner_join(df_mutations, df_variables_filtered, by = c("Sample" = "ID_CNIO"))

# Print the length of the merged dataframe
cat("Length of df_combined is:", nrow(df_combined), "\n")

# Check missing values in the merged dataframe
head(df_combined)
colSums(is.na(df_combined))  

# Filter samples with 'MMR_final_status' equal to 'MMRd'
df_mmrd <- df_combined %>%
  filter(MMR_final_status == "MMRd")
# Imprimir la longitud del DataFrame 
cat("Mutations in MMR genes", nrow(df_mmrd), "\n")

# Filter to keep only genes of interest in MMRd samples
genes_of_interest <- c("MLH1", "MSH2", "MSH6", "PMS2")
df_mmrd_genes <- df_mmrd %>%
  filter(Gene %in% genes_of_interest)

# Create a new column to eliminate duplicates
df_mmrd_genes <- df_mmrd_genes %>%
  mutate(Sample_Gene = paste(Sample, Gene, sep = "_")) %>%
  distinct(Sample_Gene, .keep_all = TRUE)

# Verify the content of the new 'Sample_Gene' column
head(df_mmrd_genes$Sample_Gene)

# Print samples and combinations Sample-Gene in the MMRd samples
df_mmrd_genes$Sample
df_mmrd_genes$Sample_Gene

# Calculate the total number of samples with mutations in MMRd genes
total_samples_mmrd_mut <- n_distinct(df_mmrd_genes$Sample)  
print(paste("The number of MMRd samples mutated in MMR genes is:",total_samples_mmrd_mut))

# Calculate the total number of MMRd samples with and without mutations
total_rows_mmrd <- df_variables %>% filter (MMR_final_status=='MMRd')
total_samples_mmrd <- nrow(total_rows_mmrd)
cat('Number of MMRd samples is:', total_samples_mmrd)

# Filter out synonymous and intronic mutations
df_mmrd_genes <- df_mmrd_genes %>%
  filter(!Function %in% c('synonymous', 'intronic')) 

# Calculate the percentage of samples mutated in the genes of interest
gene_mutation_percent <- df_mmrd_genes %>%
  group_by(Gene) %>%
  summarize(Mutated_Samples = n_distinct(Sample_Gene)) %>%
  mutate(Percent = (Mutated_Samples / total_samples_mmrd) * 100) %>%
  arrange(desc(Percent))  # Order genes by mutation frequency

# See the results
print(gene_mutation_percent)

#  Assign an order to the genes based on mutation frequency
df_mmrd_genes <- df_mmrd_genes %>%
  mutate(Gene = factor(Gene, levels = rev(gene_mutation_percent$Gene)))  # Invert the order

# Create an order for the samples based on mutations frequency
sample_order <- df_mmrd_genes %>%
  arrange(Gene) %>%
  group_by(Sample) %>%
  summarize(First_Gene_Mutated = first(Gene)) %>%  
  arrange(First_Gene_Mutated, Sample) %>%         
  pull(Sample)

# Create heatmap using df_mmrd_genes
heatmap_data <- df_mmrd_genes %>%
  filter(Gene %in% genes_of_interest) %>%
  mutate(
    Sample = factor(Sample, levels = sample_order),  # Set the samples order
    Gene = factor(Gene, levels = rev(gene_mutation_percent$Gene))  # Order of the genes
  )

# I change the order of the 'Sample' factor to group samples on the left according to the most frequently mutated gene
heatmap_data <- heatmap_data %>%
  mutate(Sample = factor(Sample, levels = rev(levels(Sample))))  # Invert level orders

# Update histology values for heatmap
histology_values <- df_combined %>%
  select(Sample, HISTOLOGY) %>%
  mutate(
    Sample = factor(Sample, levels = sample_order),  # Order of the samples
    HISTOLOGY = recode(HISTOLOGY, "0" = "EOC", "1" = "CCOC")  # Map HISTOLOGY values
  )


# Update order in HISTOLOGY 
histology_values <- histology_values %>%
  mutate(Sample = factor(Sample, levels = rev(levels(Sample)))) 

# Check common samples
common_samples <- intersect(heatmap_data$Sample, histology_values$Sample)
print(common_samples)

# Filter both dataframes to keep only common samples
heatmap_data <- heatmap_data %>%
  filter(Sample %in% common_samples)

histology_values <- histology_values %>%
  filter(Sample %in% common_samples)

# Create heatmap for mutations
heatmap_mutations <- ggplot(heatmap_data, aes(x = Sample, y = Gene, fill = Function)) +
  geom_tile(color = "white",width=1,height=0.8) +
  scale_fill_manual(
    values = c(
      "frameshift_deletion" = "green",
      "frameshift_substitution" = "darkgreen",
      "frameshift_insertion" = "lightgreen",
      "nonframeshift_deletion" = "pink",
      "nonframeshift_insertion" = "purple",
      "nonframeshift_substitution" = "brown",
      "missense"="blue",
      "startloss"="yellow",
      "splicing" = "red",
      "stopgain" = "orange",
      "no_mutation" = "white"
    ),
    na.value = "gray"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2,size = 9,face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold",vjust = 0),
    legend.position = "right",
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5),
    panel.spacing.y = unit(0, "lines"),
    plot.margin = margin(0, 0, 0, 0)

  ) +
  labs(
    title = "Mutations in MMR genes in MMRd samples",
    x =NULL,
    y = "",
    fill = "Mutation type"
  ) +
  scale_y_discrete(
    labels = function(labels) {
      paste0(labels, " (", round(gene_mutation_percent$Percent[match(labels, gene_mutation_percent$Gene)], 0), "%)")
    }
  ) +
  coord_cartesian(clip = "off")

# Create heatmap for histology
heatmap_histology <- ggplot(histology_values,  aes(x = Sample, y = 1, fill = as.factor(HISTOLOGY))) + 
  geom_tile(color = "white", width = 1) +
  scale_fill_manual(values = c("EOC" = "lightblue", "CCOC" = "lightgreen")) +
  theme_minimal() + 
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid =  element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  labs(fill = "HISTOLOGY") +
  coord_cartesian(ylim = c(0, 1))

# Unify both heatmaps ambos using `patchwork`
final_plot <- heatmap_mutations / heatmap_histology +
  plot_layout(heights = c(2, 0.4)) 
  
# Show the combined heatmap
print(final_plot)

# Save the figure
ggsave(
  "Heatmap_MMRd_DRIVERS_HIGH_FILTER.png",
  plot = final_plot,
  width = 8,
  height = 6,
  dpi = 300,
  #bg = "white"
)
