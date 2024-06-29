# 0
library(minfi)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(knitr)
library(dplyr)
library(factoextra) ###

# 1
baseDir <- ('Input')
targets <- read.metharray.sheet(baseDir)
# RGset <- read.metharray.exp(targets = targets)
# save(RGset, file='RGset.RData')
load("RGset.RData")
samplesheet <- read.csv("Input/Samplesheet.csv",header=T, stringsAsFactors=T)

# 2
Red <- getRed(RGset)
Green <- getGreen(RGset)

print(dim(Red))      # [1] 622399      8
print(dim(Green))    # [1] 622399      8

# 3
address <- "39802405"

if (address %in% rownames(Red) & address %in% rownames(Green)) 
    Red_fluor <- Red[address, ]
    Green_fluor <- Green[address, ]

load("Illumina450Manifest_clean.RData")
probe_type = Illumina450Manifest_clean[Illumina450Manifest_clean$AddressA_ID==address, 'Infinium_Design_Type']

fluorescence_data <- data.frame(
  Sample = colnames(Red),
  Red_fluor = Red_fluor,
  Green_fluor = Green_fluor,
  Type = probe_type,
  Color = 'NA'
)
fluorescence_data$Sample <- sapply(strsplit(rownames(fluorescence_data), "_"), `[`, 2)
rownames(fluorescence_data) <- NULL

kable(fluorescence_data)

# 4
# MSet.raw <- preprocessRaw(RGset)
# save(MSet.raw,file='MSet_raw.RData')
load("MSet_raw.RData")

# 5
qc <- getQC(MSet.raw)
# pdf('QCplot.pdf')
plotQC(qc)
dev.off()

# pdf(file = 'controlStripPlot_NEGATIVE.pdf')
controlStripPlot(RGset, controls = 'NEGATIVE')
dev.off()

# detP <- detectionP(RGset) 
# save(detP, file = 'detP.RData')
load("detP.RData")

threshold <- 0.05
failed <- detP > threshold
summary(failed) 

failed_positions <- colSums(failed)
failed_positions_summary <- data.frame(
  Sample = colnames(failed),
  n_Failed_positions = colSums(failed),
  means_of_columns = colMeans(failed)*100
)
failed_positions_summary$Sample <- sapply(strsplit(rownames(failed_positions_summary), "_"), `[`, 2)
rownames(failed_positions_summary) <- NULL

kable(failed_positions_summary)

# 6
beta_values <- getBeta(MSet.raw)
M_values <- getM(MSet.raw)

wt_indices <- which(samplesheet$Group == "WT")
mut_indices <- which(samplesheet$Group == "MUT")

mean_beta_wt <- rowMeans(beta_values[, wt_indices], na.rm = TRUE)
mean_beta_mut <- rowMeans(beta_values[, mut_indices], na.rm = TRUE)

mean_beta_wt <- mean_beta_wt[!is.na(mean_beta_wt)]
mean_beta_mut <- mean_beta_mut[!is.na(mean_beta_mut)]

M_values_wt <- as.vector(M_values[, wt_indices])
M_values_mut <- as.vector(M_values[, mut_indices])

M_values_wt <- M_values_wt[!is.na(M_values_wt)]
M_values_mut <- M_values_mut[!is.na(M_values_mut)]

par(mfrow = c(1, 2))

plot(density(mean_beta_wt), main = "WT Beta Values", col = "#FF9900")
lines(density(mean_beta_mut), col = "#9966FF")
legend("topright", legend = c("WT", "MUT"), col = c("#FF9900", "#9966FF"), lwd = 3, cex = 0.8)

plot(density(M_values_wt), xlim = c(-8, 8), ylim = c(0, 0.2), main = "WT M Values", col = "#FF9900")
lines(density(M_values_mut), col = "#9966FF")
legend("topright", legend = c("WT", "MUT"), col = c("#FF9900", "#9966FF"), lwd = 3, cex = 0.8)

par(mfrow = c(1, 1))

# 7
dfI <- Illumina450Manifest_clean %>% filter(Infinium_Design_Type == "I") %>% droplevels()
dfII <- Illumina450Manifest_clean %>% filter(Infinium_Design_Type == "II") %>% droplevels()

beta <- getBeta(MSet.raw)
beta_I <- beta[rownames(beta) %in% dfI$IlmnID, ]
beta_II <- beta[rownames(beta) %in% dfII$IlmnID, ]

beta_I <- beta_I[rowSums(is.na(beta_I)) != ncol(beta_I), ]
beta_II <- beta_II[rowSums(is.na(beta_II)) != ncol(beta_II), ]

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

RGSet_SWAN <- preprocessSWAN(RGset)

beta_SWAN <- getBeta(RGSet_SWAN)
beta_I_SWAN <- beta_SWAN[rownames(beta_SWAN) %in% dfI$IlmnID, ]
beta_II_SWAN <- beta_SWAN[rownames(beta_SWAN) %in% dfII$IlmnID, ]

beta_I_SWAN <- beta_I_SWAN[rowSums(is.na(beta_I_SWAN)) != ncol(beta_I_SWAN), ]
beta_II_SWAN <- beta_II_SWAN[rowSums(is.na(beta_II_SWAN)) != ncol(beta_II_SWAN), ]

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

colnames(beta) <- sub("^[^_]+_", "", colnames(beta))
colnames(beta_SWAN) <- sub("^[^_]+_", "", colnames(beta_SWAN))

# png("density_plots.png", width = 1200, height = 800)
par(mfrow = c(2, 3), xpd = TRUE)
targets$Group <- as.factor(targets$Group)

palette(c("#D6604D", "#4393C3"))

plot(density_mean_beta_I, col = "#FF6600", main = "Raw Beta Density", lwd = 2.5, xlab = 'Beta')
lines(density_mean_beta_II, col = "#0033FF", lwd = 2.5)
legend('topright', legend = c("Type I", "Type II"), fill = c("#FF6600", "#0033FF"), cex = 1)

plot(density_sd_of_beta_I, col = "#FF6600", main = "Raw SD Density", lwd = 2.5, xlab = 'Beta SD')
lines(density_sd_of_beta_II, col = "#0033FF", lwd = 2.5)
legend('topright', legend = c("Type I", "Type II"), fill = c("#FF6600", "#0033FF"), cex = 1)

boxplot(beta, col = targets$Group, names = colnames(beta), main = "Raw Beta by Sample", las = 2)

plot(density_mean_beta_I_SWAN, col = "#FF6600", main = "SWAN Beta Density", lwd = 2.5, xlab = 'Beta')
lines(density_mean_beta_II_SWAN, col = "#0033FF", lwd = 2.5)
legend('topright', legend = c("Type I", "Type II"), fill = c("#FF6600", "#0033FF"), cex = 1)

plot(density_sd_of_beta_I_SWAN, col = "#FF6600", main = "SWAN SD Density", lwd = 2.5, xlab = 'Beta SD')
lines(density_sd_of_beta_II_SWAN, col = "#0033FF", lwd = 2.5)
legend('topright', legend = c("Type I", "Type II"), fill = c("#FF6600", "#0033FF"), cex = 1)

boxplot(beta_SWAN, col = targets$Group, names = colnames(beta_SWAN), main = "SWAN Beta by Sample", las = 2)

par(mfrow = c(1, 1))
dev.off()

pca_results <- prcomp(t(beta_SWAN), scale. = TRUE)

# png("PCA.png", width = 1200, height = 600)

par(mfrow = c(1, 2))

palette(c("#FF0000", "#9966FF"))

plot(pca_results$x[, 1], pca_results$x[, 2], cex = 2, pch = 2, col = targets$Group,
     xlab = "PC1", ylab = "PC2", xlim = c(-1000, 1000), ylim = c(-1000, 1000))
text(pca_results$x[, 1], pca_results$x[, 2], labels = rownames(pca_results$x), cex = 0.5, pos = 1)
legend("bottomright", legend = levels(targets$Group), col = 1:nlevels(targets$Group), pch = 2)

palette(c("#FF9999","#0066CC"))

gender <- c("M" =-0x2640L, "F" = -0x2642L)

plot(pca_results$x[, 1], pca_results$x[, 2], cex = 1.5, pch = gender[samplesheet$Sex], col = samplesheet$Sex,
     xlab = "PC1", ylab = "PC2", xlim = c(-1000, 1000), ylim = c(-1000, 1000))
text(pca_results$x[, 1], pca_results$x[, 2], labels = rownames(pca_results$x), cex = 0.5, pos = 1)
legend("bottomright", legend=levels(samplesheet$Sex),  col = c(1:nlevels(samplesheet$Sex)), pch = gender)

dev.off()

fviz_eig(pca_results, addlabels = T,xlab='PC number',ylab='% of variance', barfill = "#0063A6", barcolor = "black") ###

# 9
mann_whitney_function <- function(x) {
  mw_test <- wilcox.test(x ~ samplesheet$Group)
  return(mw_test$p.value)
}

pValues_mw <- apply(beta_SWAN, 1, mann_whitney_function)

final_mw <- data.frame(beta_SWAN, pValues_mw)
final_mw <- final_mw[order(final_mw$pValues_mw),]

kable(head(final_mw))

# save(final_mw, file = "final_mw.RData")
load("final_mw.RData")

hist_res <- hist(final_mw$pValues_mw, main = "P-value distribution (Mann-Whitney test)", xlab = 'p-val', col = "#A9A9A9", border = "#A9A9A9")
abline(v = 0.05, col = "#EF233C", lwd = 2)

# 10
corrected_pValues_Bonf <- p.adjust(final_mw$pValues_mw,"bonferroni")
corrected_pValues_BH <- p.adjust(final_mw$pValues_mw,"BH")
final_mw_corr <- data.frame(sum(final_mw <= 0.05), sum(corrected_pValues_Bonf <= 0.05), sum(corrected_pValues_BH <= 0.05))

kable(final_mw_corr)

par(mfrow=c(1,1))
boxplot(final_mw_corr [,9:11], col = c("seashell2", "#C7FF38", "red"), names=NA) ###

# 11
WT_group <- final_mw_corr[,targets$Group=="WT"]
WT_group_mean <- apply(WT_group, 1, mean)

MUT_group <- final_mw_corr[,targets$Group=="MUT"]
MUT_group_mean <- apply(MUT_group, 1, mean)

delta <- WT_group_mean - MUT_group_mean

toVolcPlot <- data.frame(delta, -log10(final_mw_corr$pValues_mw))

kable(head(toVolcPlot))

# png(filename = "VolcanoPlot.png", width = 800, height = 600)
par(mfrow=c(1,1))

HighLight <- toVolcPlot[abs(toVolcPlot[,1])>0.1 & toVolcPlot[,2]>(-log10(0.05)),]
plot(toVolcPlot[,1],toVolcPlot[,2], pch=19, cex=0.6 ,xlab="delta",ylab="-log10(Pvalue)", col="black")
abline(a=-log10(0.05),b=0, col="#fb8b24")
points(HighLight[,1], HighLight[,2], pch=17, cex=0.6, col="#CC00FF")

dev.off()













