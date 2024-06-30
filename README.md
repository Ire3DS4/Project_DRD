# Infinium Methylation Data Analysis

This project for the exam of DNA/RNA Dynamics focuses on analyzing Infinium methylation data using R.

P.S: in the README file there are going to be code snippets as well as interpretations of the results and various plots; for the complete code please look at 'R code'.

## Table of Contents
1. [Introduction](#introduction)
2. [Data Preparation](#data-preparation)
3. [Preprocessing and Normalization](#preprocessing-and-normalization)
4. [Visualization and Analysis](#visualization-and-analysis)
5. [Conclusions](#conclusions)

## Introduction

The following pipeline is designed to analyze DNA methylation data from the Illumina HumanMethylation450 BeadChip platform, a popular tool for genome-wide methylation studies. The Infinium assay provides insights into DNA methylation patterns, which are crucial for understanding various biological processes and disease mechanisms.

This repository offers a comprehensive, user-friendly pipeline for processing, analyzing, and interpreting methylation data from the Infinium platform.

### Workflow Overview
Quality Control: Assess data quality to identify potential issues or biases.
Preprocessing: Apply normalization techniques to reduce technical variations and batch effects.
Differential Methylation Analysis: Identify differentially methylated regions (DMRs) or CpG sites associated with specific conditions or phenotypes.
Visualization: Create informative plots to visualize DNA methylation patterns and results.

## Data Preparation

### 1. Load Raw Data

Begin by loading the raw data files, which contain intensity values for both methylated and unmethylated probes across samples. This step ensures that the data is correctly formatted for subsequent analyses.

```r
# Load necessary libraries
library(minfi)
library(knitr)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

# Load raw data
baseDir <- ('Input')
targets <- read.metharray.sheet(baseDir)
RGset <- read.metharray.exp(targets = targets)
```

### 2. Create R/G Dataframes

Separate the raw intensity data into red and green channels. This distinction helps in identifying the chemistry of the probes and is essential for analyzing the fluorescence signals accurately.

```r
# Create R/G dataframes
Red <- data.frame(getRed(RGset))
dim(Red)
```
```r
## [1] 622399      8
```

```r
Green <- data.frame(getGreen(RGset))
dim(Green)
```
```r
## [1] 622399      8
```

### 3. Check Probe Info by Address

To determine the probe type, refer to the manifest file from Illumina, which categorizes probes into Type I or Type II. This classification is important because Type II probes use only one color channel, whereas Type I probes use both.

```r
# Check probe info by address
address <- "39802405"

if (address %in% rownames(Red) & address %in% rownames(Green)) 
    Red_fluor <- Red[address, ]
    Green_fluor <- Green[address, ]

load("Illumina450Manifest_clean.RData")
probe_type = Illumina450Manifest_clean[Illumina450Manifest_clean$AddressA_ID==address, 'Infinium_Design_Type']
```
```r
## [1] II
## Levels: I II
```

```r
# Create and fill the fluorescence_data dataframe

fluorescence_data <- data.frame(
  Sample = colnames(Red),
  Red_fluor = Red_fluor,
  Green_fluor = Green_fluor,
  Type = probe_type,
  Color = 'NA'
)
fluorescence_data$Sample <- sapply(strsplit(rownames(fluorescence_data), "_"), `[`, 2)
rownames(fluorescence_data) <- NULL
```

```r
|Sample | Red_fluor| Green_fluor|Type |Color |
|:------|---------:|-----------:|:----|:-----|
|R01C01 |      4254|        8361|II   |NA    |
|R02C01 |      4584|       10343|II   |NA    |
|R03C01 |      4201|        9859|II   |NA    |
|R04C01 |      3627|        8552|II   |NA    |
|R02C02 |      5669|        1003|II   |NA    |
|R03C02 |      7689|        1041|II   |NA    |
|R04C02 |      5954|        6336|II   |NA    |
|R05C02 |      5989|         761|II   |NA    |
```
Notice that the assigned probes are of Infinium II design, hence no color channel needs to be specified since type II probes use a single bead type for both methylated and unmethylated states and measure the intensities of the two states using the same color channel. As a result, Type II probes do not require separate color information for the red or green channels because they are not differentiated by color.

### 4. Create the Object MSet.raw

The `MSet.raw` object contains methylated and unmethylated signal intensities, facilitating further analysis of methylation levels.

```r
# Create MSet.raw object
MSet.raw <- preprocessRaw(RGset)
```


## Preprocessing and Normalization

### 5. Quality Check

QCplot visualizes the median intensities of methylation and unmethylation signals. High median values in both distributions indicate good data quality. However, the plot has limitations, such as not accounting for background signals or potential sample preparation failures.

##### 5.1 QCplot

The plot shows clustering of samples with high median values, indicating good quality data. Low median values suggest poor quality, which can affect downstream analyses.

```r
# Quality check
qc <- qcReport(RGset, pdf = "QCReport.pdf")
```
![QC plot](plots/QCplot.png)

##### 5.2 Negative control intensity check

Negative control probes estimate system background intensity. Normally, these values range from 100 to 1000 units. High values may indicate degraded DNA, affecting data quality. Using controlStripPlot(), we can check the intensity levels, with higher values indicating potential issues.

```r
controlStripPlot(RGset, controls="NEGATIVE")
```
![Control strip plot](plots/ControlStripPlot.png)

##### 5.3 Failed positions
```r
# Calculate detection pValues

failed_positions <- colSums(failed)

# Summarize failed positions per sample
failed_positions_summary <- data.frame(
  Sample = colnames(failed),
  n_Failed_positions = failed_positions
)
failed_positions_summary$Sample <- sapply(strsplit(rownames(failed_positions_summary), "_"), `[`, 2)
rownames(failed_positions_summary) <- NULL
```

```r
|Sample | n_Failed_positions|
|:------|------------------:|
|R01C01 |                 45|
|R02C01 |                 26|
|R03C01 |                 28|
|R04C01 |                 32|
|R02C02 |                190|
|R03C02 |                130|
|R04C02 |                 17|
|R05C02 |                406|
```

### 6. Beta and M Values

Beta values represent the proportion of methylation at each CpG site, while M values are logit-transformed beta values.

```r
# Calculate Beta and M values
beta_values <- getBeta(MSet.raw)
M_values <- getM(MSet.raw)
```
![WT beta and M plots](plots/WT_beta_m.png)

### 7. Functional Normalization

Normalization adjusts for technical variation between samples.

```r
# Perform functional normalization
norm_result <- preprocessNoob(RGset)
beta_norm <- getBeta(norm_result)
```

```r
# Filter Infinium I and II data
dfI <- Illumina450Manifest_clean %>% filter(Infinium_Design_Type == "I") %>% droplevels()
dfII <- Illumina450Manifest_clean %>% filter(Infinium_Design_Type == "II") %>% droplevels()

# Get beta values for Infinium I and II
beta <- getBeta(MSet.raw)
beta_I <- beta[rownames(beta) %in% dfI$IlmnID, ]
beta_II <- beta[rownames(beta) %in% dfII$IlmnID, ]

# Remove rows with all NAs
beta_I <- beta_I[rowSums(is.na(beta_I)) != ncol(beta_I), ]
beta_II <- beta_II[rowSums(is.na(beta_II)) != ncol(beta_II), ]

# Calculate densities for raw means and sds
mean_beta_I <- rowMeans(beta_I, na.rm = TRUE)
mean_beta_I <- na.omit(mean_beta_I)
density_mean_beta_I <- density(mean_beta_I)

mean_beta_II <- rowMeans(beta_II, na.rm = TRUE)
mean_beta_II <- na.omit(mean_beta_II)
density_mean_beta_II <- density(mean_beta_II)

sd_beta_I <- apply(beta_I, 1, sd, na.rm = TRUE)
sd_beta_I <- na.omit(sd_beta_I)
density_sd_of_beta_I <- density(sd_beta_I)

sd_beta_II <- apply(beta_II, 1, sd, na.rm = TRUE)
sd_beta_II <- na.omit(sd_beta_II)
density_sd_of_beta_II <- density(sd_beta_II)

# Apply SWAN normalization
RGSet_SWAN <- preprocessSWAN(RGset)

# Get beta values for Infinium I and II after SWAN normalization
beta_SWAN <- getBeta(RGSet_SWAN)
beta_I_SWAN <- beta_SWAN[rownames(beta_SWAN) %in% dfI$IlmnID, ]
beta_II_SWAN <- beta_SWAN[rownames(beta_SWAN) %in% dfII$IlmnID, ]

# Remove rows with all NAs
beta_I_SWAN <- beta_I_SWAN[rowSums(is.na(beta_I_SWAN)) != ncol(beta_I_SWAN), ]
beta_II_SWAN <- beta_II_SWAN[rowSums(is.na(beta_II_SWAN)) != ncol(beta_II_SWAN), ]

# Calculate densities for SWAN means and sds
mean_beta_I_SWAN <- rowMeans(beta_I_SWAN, na.rm = TRUE)
mean_beta_I_SWAN <- na.omit(mean_beta_I_SWAN)
density_mean_beta_I_SWAN <- density(mean_beta_I_SWAN)

mean_beta_II_SWAN <- rowMeans(beta_II_SWAN, na.rm = TRUE)
mean_beta_II_SWAN <- na.omit(mean_beta_II_SWAN)
density_mean_beta_II_SWAN <- density(mean_beta_II_SWAN)

sd_beta_I_SWAN <- apply(beta_I_SWAN, 1, sd, na.rm = TRUE)
sd_beta_I_SWAN <- na.omit(sd_beta_I_SWAN)
density_sd_of_beta_I_SWAN <- density(sd_beta_I_SWAN)

sd_beta_II_SWAN <- apply(beta_II_SWAN, 1, sd, na.rm = TRUE)
sd_beta_II_SWAN <- na.omit(sd_beta_II_SWAN)
density_sd_of_beta_II_SWAN <- density(sd_beta_II_SWAN)
```
![Density plots](plots/density_plots.png)

### 8. PCA

Principal Component Analysis (PCA) reduces data dimensionality and helps visualize sample clustering.

```r
# Perform PCA
par(mfrow = c(1, 2))

# Define the palette for the first plot
palette(c("#FF0000", "#9966FF"))

# Plot the first PCA plot
plot(pca_results$x[, 1], pca_results$x[, 2], cex = 2, pch = 2, col = targets$Group,
     xlab = "PC1", ylab = "PC2", xlim = c(-1000, 1000), ylim = c(-1000, 1000))
text(pca_results$x[, 1], pca_results$x[, 2], labels = rownames(pca_results$x), cex = 0.5, pos = 1)
legend("bottomright", legend = levels(targets$Group), col = 1:nlevels(targets$Group), pch = 2)

# Define the palette for the second plot
palette(c("#FF9999","#0066CC"))

# Set shapes for sexes
gender <- c("M" =-0x2640L, "F" = -0x2642L)

# Plot the second PCA plot
plot(pca_results$x[, 1], pca_results$x[, 2], cex = 1.5, pch = gender[samplesheet$Sex], col = samplesheet$Sex,
     xlab = "PC1", ylab = "PC2", xlim = c(-1000, 1000), ylim = c(-1000, 1000))
text(pca_results$x[, 1], pca_results$x[, 2], labels = rownames(pca_results$x), cex = 0.5, pos = 1)
legend("bottomright", legend=levels(samplesheet$Sex),  col = c(1:nlevels(samplesheet$Sex)), pch = gender)
```
![PCA](plots/PCA.png)

## Visualization and Analysis

### 9. Differential Methylation Analysis

We identify differentially methylated positions (DMPs) between WT and MUT groups.

```r
# Define the Mann-Whitney test function
mann_whitney_function <- function(x) {
  mw_test <- wilcox.test(x ~ samplesheet$Group)
  return(mw_test$p.value)
}

# Apply the Mann-Whitney test function to each row
pValues_mw <- apply(beta_SWAN, 1, mann_whitney_function)

# Create a data frame with beta_SWAN values and p-values
final_mw <- data.frame(beta_SWAN, pValues_mw)
final_mw <- final_mw[order(final_mw$pValues_mw),]

kable(head(final_mw))
```

```r
|           |    R01C01|    R02C01|    R03C01|    R04C01|    R02C02|    R03C02|    R04C02|    R05C02| pValues_mw|
|:----------|---------:|---------:|---------:|---------:|---------:|---------:|---------:|---------:|----------:|
|cg00213748 | 0.2800583| 0.3392271| 0.2065237| 0.2220701| 0.4043062| 0.5243037| 0.1556169| 0.5729493|  0.0285714|
|cg03695421 | 0.6014599| 0.5250000| 0.5522746| 0.6032083| 0.3615561| 0.3838308| 0.5667055| 0.4803304|  0.0285714|
|cg03750315 | 0.0373791| 0.0835118| 0.0783275| 0.0479203| 0.4081164| 0.5822222| 0.0828804| 0.6153310|  0.0285714|
|cg04462340 | 0.5786193| 0.6238414| 0.6119760| 0.5876453| 0.6516410| 0.7722205| 0.5696010| 0.7369912|  0.0285714|
|cg05865243 | 0.9346293| 0.8883589| 0.9297575| 0.9287217| 0.8526497| 0.7607313| 0.9129813| 0.4261603|  0.0285714|
|cg07939587 | 0.9090038| 0.8314741| 0.9360740| 0.9045395| 0.8650881| 0.7309030| 0.9311156| 0.5170557|  0.0285714|
```

### 10. Multiple Test Correction

Adjust p-values for multiple testing to control the false discovery rate (FDR).

```r
# Multiple test correction
corrected_pValues_Bonf <- p.adjust(final_mw$pValues_mw,"bonferroni")
corrected_pValues_BH <- p.adjust(final_mw$pValues_mw,"BH")
final_mw_corr <- data.frame(sum(final_mw <= 0.05), sum(corrected_pValues_Bonf <= 0.05), sum(corrected_pValues_BH <= 0.05))
kable(final_mw_corr)
```

```r
|                                 | before_corr| after_Bonf| after_BH|
|:--------------------------------|-----------:|----------:|--------:|
|Differentially methylated probes |       53343|          0|        0|
```

### 11. Volcano and Manhattan Plots

Visualize DMPs with volcano and Manhattan plots.

```r
# WT group mean
WT_group <- final_mw_corr[,targets$Group=="WT"]
WT_group_mean <- apply(WT_group, 1, mean)

# MUT group mean
MUT_group <- final_mw_corr[,targets$Group=="MUT"]
MUT_group_mean <- apply(MUT_group, 1, mean)

# Compute delta
delta <- WT_group_mean - MUT_group_mean

#create a dataframe
toVolcPlot <- data.frame(delta, -log10(final_mw_corr$pValues_mw))
kable(head(toVolcPlot))
```

```r
|           |      delta| X.log10.final_mw_corr.pValues_mw.|
|:----------|----------:|---------------------------------:|
|cg00213748 | -0.1046692|                          1.544068|
|cg03695421 |  0.1565745|                          1.544068|
|cg03750315 | -0.1773085|                          1.544068|
|cg04462340 | -0.0461887|                          1.544068|
|cg05865243 |  0.1515448|                          1.544068|
|cg07939587 |  0.1439945|                          1.544068|
```
![Volcano Plot](plots/VolcanoPlot.png)

### 12. Heatmap

Heatmaps visualize methylation patterns across samples.

```r
# Heatmap
library(pheatmap)
significant_dmp <- dmp[dmp$adj.P.Val < 0.05, ]
pheatmap(beta_norm[rownames(significant_dmp), ])
```

## Conclusions

This project provides a comprehensive analysis of Infinium methylation data, including preprocessing, normalization, visualization, and statistical analysis. The key steps are loading data, defining sample groups, normalizing data, creating informative plots, and performing differential methylation analysis. The analysis reveals differences in methylation patterns between WT and MUT groups, which can be further explored for biological insights.

### Additional Resources
- [minfi package documentation](https://bioconductor.org/packages/release/bioc/html/minfi.html)
- [Infinium Methylation Assay overview](https://www.illumina.com/techniques/microarrays/microarray-kits/infinium-methylation-assay.html)

By following this guide, you can replicate the analysis and gain insights into methylation patterns in your data.
