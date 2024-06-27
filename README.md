# Infinium Methylation Data Analysis

This project focuses on analyzing Infinium methylation data using R. We will preprocess the raw data, perform normalization, visualize the results, and interpret the findings. The project involves several steps, including loading data, preprocessing, visualization, and statistical analysis.

## Table of Contents
1. [Introduction](#introduction)
2. [Data Preparation](#data-preparation)
3. [Preprocessing and Normalization](#preprocessing-and-normalization)
4. [Visualization and Analysis](#visualization-and-analysis)
5. [Conclusions](#conclusions)

## Introduction

Methylation analysis is crucial for understanding epigenetic modifications that affect gene expression without changing the DNA sequence. This project uses Infinium methylation arrays to measure DNA methylation at specific CpG sites. The workflow includes preprocessing raw data, normalization, visualization, and statistical analysis to identify patterns and differences between groups (e.g., wild type (WT) and mutant (MUT)).

## Data Preparation

First, we load the necessary libraries and data files. The raw data files include methylation intensity values for different samples.

```r
# Load necessary libraries
library(minfi)
library(knitr)

# Load data
beta_values <- getBeta(MSet.raw)
M_values <- getM(MSet.raw)
pheno_data <- read.csv("C:/Users/rober/OneDrive/Desktop/Uni/Bioinformatics/DNA-RNA/report/input_data/Samplesheet_report_2024.csv", header=TRUE, stringsAsFactors=TRUE)
```

### Sample Group Definition

We define sample groups based on the phenotype data. The WT and MUT groups are created using the sample sheet.

```r
# Define WT and MUT sample groups
WT_samples <- apply(pheno_data[pheno_data$Group == "WT", c("Sentrix_ID", "Sentrix_Position")], 1, function(row) paste(row[1], row[2], sep = "_"))
MUT_samples <- apply(pheno_data[pheno_data$Group == "MUT", c("Sentrix_ID", "Sentrix_Position")], 1, function(row) paste(row[1], row[2], sep = "_"))
```

## Preprocessing and Normalization

### Raw Data Visualization

We create a function to plot the density of beta values and M values for the WT and MUT groups.

```r
# Define a function to make plots
create_plot <- function(data1, data2, func, title) {
  aggregated_data1 <- apply(data1, 1, func)
  aggregated_data2 <- apply(data2, 1, func)
  density_data1 <- density(aggregated_data1, na.rm=TRUE)
  density_data2 <- density(aggregated_data2, na.rm=TRUE)
  plot(density_data1, col="blue", main=title)
  lines(density_data2, col="red")
}

# Plot beta value means
create_plot(beta_values[, WT_samples], beta_values[, MUT_samples], mean, "Distribution of Beta-value Means")
legend("topright", legend = c("WT", "MUT"), col = c("blue", "red"), lty = 1)

# Plot M value means
create_plot(M_values[, WT_samples], M_values[, MUT_samples], mean, "Distribution of M-value Means")
legend("topright", legend = c("WT", "MUT"), col = c("blue", "red"), lty = 1)
```

### Normalization

We preprocess the raw data using the Noob normalization method, which corrects for background fluorescence.

```r
# Normalize data using preprocessNoob
norm_result <- preprocessNoob(RGset)
beta_norm <- getBeta(norm_result) 
```

### Visualization of Normalized Data

We create plots to compare the distribution of normalized beta values between Infinium I and Infinium II probes.

```r
# Raw and normalized beta values divided by Infinium design
df_design_I <- Illumina450Manifest_clean[Illumina450Manifest_clean$Infinium_Design_Type == "I",]
df_design_I <- droplevels(df_design_I)
df_design_II <- Illumina450Manifest_clean[Illumina450Manifest_clean$Infinium_Design_Type == "II",]
df_design_II <- droplevels(df_design_II)

beta_I <- beta_values[rownames(beta_values) %in% df_design_I$IlmnID,]
beta_II <- beta_values[rownames(beta_values) %in% df_design_II$IlmnID,]
beta_norm_I <- beta_norm[rownames(beta_norm) %in% df_design_I$IlmnID,]
beta_norm_II <- beta_norm[rownames(beta_norm) %in% df_design_II$IlmnID,]

group_colors <- ifelse(colnames(beta_values) %in% WT_samples, "blue", "red")
par(mfrow = c(2, 3))

create_plot(beta_I, beta_II, mean, "Beta-value Mean Distribution")
legend("topright", legend = c("Infinium I", "Infinium II"), col = c("blue", "red"), lty = 1)

create_plot(beta_I, beta_II, sd, "Beta-value StD Distribution")
legend("topright", legend = c("Infinium I", "Infinium II"), col = c("blue", "red"), lty = 1)

boxplot(beta_values, col = group_colors, las = 2, names = pheno_data[,"Sentrix_Position"], main="Boxplots of Beta Values")

create_plot(beta_norm_I, beta_norm_II, mean, "Normalized Beta-value Mean Distribution")
legend("topright", legend = c("Infinium I", "Infinium II"), col = c("blue", "red"), lty = 1)

create_plot(beta_norm_I, beta_norm_II, sd, "Normalized Beta-value StD Distribution")
legend("topright", legend = c("Infinium I", "Infinium II"), col = c("blue", "red"), lty = 1)

boxplot(beta_norm, col = group_colors, las = 2, names = pheno_data[,"Sentrix_Position"], main="Boxplots of Normalized Beta Values")
```

### Save Plots to PDF

We save a plot to a PDF file for documentation.

```r
# Save plot to a PDF
pdf("negative_control_plot.pdf")
controlStripPlot(RGset, controls = 'NEGATIVE')
dev.off()
```

### Failed Positions Analysis

We create a dataframe to count the number of failed positions for each sample.

```r
# Create a dataframe for failed positions
failed_summary <- summary(failed)
failed_counts <- sapply(failed_summary, function(x) sum(x == TRUE))
sample_names <- gsub(".*_", "", names(failed_summary))
failed_positions_df <- data.frame(
  Sample = sample_names,
  n_Failed_positions = failed_counts
)

# Display the dataframe aesthetically
kable(failed_positions_df)
```

## Visualization and Analysis

### PCA Analysis

Principal Component Analysis (PCA) is performed to reduce the dimensionality of the data and visualize the differences between WT and MUT groups.

```r
# Perform PCA
pca_result <- prcomp(t(beta_values), scale. = TRUE)
plot(pca_result$x[, 1], pca_result$x[, 2], col = group_colors, main = "PCA of Beta Values", xlab = "PC1", ylab = "PC2")
legend("topright", legend = c("WT", "MUT"), col = c("blue", "red"), pch = 1)
```

**Interpretation:** PCA helps in visualizing the separation between WT and MUT groups based on their methylation profiles. Clear clustering indicates distinct methylation patterns between the groups.

### Boxplots and Density Plots

Boxplots and density plots are used to visualize the distribution of beta values and M values.

```r
# Boxplots and density plots
boxplot(beta_values, col = group_colors, las = 2, names = pheno_data[,"Sentrix_Position"], main="Boxplots of Beta Values")
plot(density(rowMeans(beta_values[, WT_samples], na.rm = TRUE)), col = "blue", main = "Density of Beta Values")
lines(density(rowMeans(beta_values[, MUT_samples], na.rm = TRUE)), col = "red")
legend("topright", legend = c("WT", "MUT"), col = c("blue", "red"), lwd = 2)
```

**Interpretation:** Boxplots show the range and distribution of beta values for each sample, while density plots compare the overall distribution between WT and MUT groups.

## Conclusions

This project provides a comprehensive analysis of Infinium methylation data, including preprocessing, normalization, visualization, and statistical analysis. Key steps include loading data, defining sample groups, normalizing data, and creating informative plots. The analysis reveals differences in methylation patterns between WT and MUT groups, which can be further explored for biological insights.

### Additional Resources
- [minfi package documentation](https://bioconductor.org/packages/release/bioc/html/minfi.html)
- [Infinium Methylation Assay overview](https://www.illumina.com/techniques/microarrays/microarray-kits/infinium-methylation-assay.html)

By following this guide, you can replicate the analysis and gain insights into methylation patterns in your data.
