# 0
library(minfi)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(knitr)
library(dplyr)
library(qqman)  
library(gplots)
library(ggplot2)
library(viridis)

# 1
baseDir <- ('Input')
targets <- read.metharray.sheet(baseDir)
RGset <- read.metharray.exp(targets = targets)
save(RGset, file='RGset.RData')                    # note that I could not upload it on GitHubs since even the zipped file is >25MB
# load("RGset.RData")
samplesheet <- read.csv("Input/Samplesheet_report_2024.csv",header=T, stringsAsFactors=T)

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
# png("QCplot.png", width = 1200, height = 800)
plotQC(qc)
dev.off()

# png("controlStripPlot_NEGATIVE.png", width = 1200, height = 800)
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
beta <- getBeta(MSet.raw)
M <- getM(MSet.raw)

wt_indices <- which(samplesheet$Group == "WT")
mut_indices <- which(samplesheet$Group == "MUT")

beta_wt <- rowMeans(beta[, wt_indices], na.rm = TRUE)
beta_mut <- rowMeans(beta[, mut_indices], na.rm = TRUE)

beta_wt <- beta_wt[!is.na(beta_wt)]
beta_mut <- beta_mut[!is.na(beta_mut)]

M_values_wt <- as.vector(M[, wt_indices])
M_values_mut <- as.vector(M[, mut_indices])

M_values_wt <- M_values_wt[!is.na(M_values_wt)]
M_values_mut <- M_values_mut[!is.na(M_values_mut)]

# png("BetaMvalues.png", width = 1200, height = 800)
par(mfrow = c(1, 2))

plot(density(beta_wt), main = "Beta Values", col = "#FF9900")
lines(density(beta_mut), col = "#9966FF")
legend("topright", legend = c("WT", "MUT"), col = c("#FF9900", "#9966FF"), lwd = 3, cex = 0.8)

plot(density(M_values_wt), xlim = c(-8, 8), ylim = c(0, 0.2), main = "M Values", col = "#FF9900")
lines(density(M_values_mut), col = "#9966FF")
legend("topright", legend = c("WT", "MUT"), col = c("#FF9900", "#9966FF"), lwd = 3, cex = 0.8)

dev.off()

# 7
dfI <- Illumina450Manifest_clean %>% filter(Infinium_Design_Type == "I") %>% droplevels()
dfII <- Illumina450Manifest_clean %>% filter(Infinium_Design_Type == "II") %>% droplevels()

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

# 8
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

var_explained <- pca_results$sdev^2 / sum(pca_results$sdev^2) * 100

eig_df <- data.frame(
  PC = seq_along(var_explained),
  Variance = var_explained
)

ggplot(eig_df, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "#00CCCC", color = "black") +
  geom_text(aes(label = round(Variance, 1)), vjust = -0.5) +
  labs(x = "PCs", y = "% of variance") +
  theme_minimal()

# ggsave("PCA_variance.png", width = 8, height = 6)

# 9
mann_whitney_function <- function(x) {
  mw_test <- wilcox.test(x ~ samplesheet$Group)
  return(mw_test$p.value)
}

pValues_mw <- apply(beta_SWAN, 1, mann_whitney_function)
final_mw <- data.frame(beta_SWAN, pValues_mw)
final_mw <- final_mw[order(final_mw$pValues_mw),]

save(final_mw, file = "final_mw.RData")                       # note that I could not upload it on GitHubs since even the zipped file is >25MB
# load("final_mw.RData")

kable(head(final_mw))

# png(filename = "pvalue_distribution.png", width = 800, height = 600)

hist_res <- hist(final_mw$pValues_mw, main = "P-value distribution (Mann-Whitney test)", xlab = 'p-val', col = "#A9A9A9", border = "#A9A9A9")
abline(v = 0.05, col = "#EF233C", lwd = 2)

dev.off()

# 10
corr_pValues_Bonf <- p.adjust(final_mw$pValues_mw, "bonferroni")
corr_pValues_BH <- p.adjust(final_mw$pValues_mw, "BH")
final_mw_corr <- data.frame(final_mw, corr_pValues_Bonf, corr_pValues_BH)

before_corr <- nrow(final_mw[final_mw$pValues_mw <= 0.05,])
after_Bonf <- nrow(final_mw[final_mw$corr_pValues_Bonf <= 0.05,])
after_BH <- nrow(final_mw_corr[final_mw_corr$corr_pValues_BH <= 0.05,])

diff_meth_probes <- data.frame(before_corr, after_Bonf, after_BH)
rownames(diff_meth_probes) <- c("Differentially methylated probes")

kable(diff_meth_probes)

# png("boxplot_diff_meth_probes.png", width = 800, height = 600)
par(mfrow=c(1,1))
boxplot(final_mw_corr [,9:11], col = c("darkslategray1", "cornsilk", "deepskyblue"), names=NA) 
legend("bottomright", legend=c("raw", "Bonferroni", "BH"),col=c("darkslategray1", "cornsilk", "deepskyblue"), lty=1:1, cex=0.8, xpd=TRUE, pch=15)
dev.off()

# 11
WT_group <- final_mw_corr[, targets$Group == "WT"]
MUT_group <- final_mw_corr[, targets$Group == "MUT"]
mean_WT_group <- apply(WT_group, 1, mean)
mean_MUT_group <- apply(MUT_group, 1, mean)
delta <- mean_WT_group - mean_MUT_group

toVolcPlot <- data.frame(delta, -log10(final_mw_corr$pValues_mw))
kable(head(toVolcPlot))

# png(filename = "VolcanoPlot.png", width = 800, height = 600)
par(mfrow=c(1,1))

HighLight <- toVolcPlot[abs(toVolcPlot[,1])>0.1 & toVolcPlot[,2]>(-log10(0.05)),]
plot(toVolcPlot[,1],toVolcPlot[,2], pch=19, cex=0.6 ,xlab="delta",ylab="-log10(Pvalue)", col="black")
abline(a=-log10(0.05),b=0, col="#FB8B24")
points(HighLight[,1], HighLight[,2], pch=17, cex=0.6, col="#CC00FF")

dev.off()

final_mw_corr_df <- data.frame(IlmnID = rownames(final_mw_corr), final_mw_corr)
final_mw_corr_ann <- merge(final_mw_corr_df, Illumina450Manifest_clean, by = "IlmnID")

input_Manhattan <- data.frame(id = final_mw_corr_ann$IlmnID, 
                              chr = final_mw_corr_ann$CHR, 
                              map = final_mw_corr_ann$MAPINFO, 
                              pval = final_mw_corr_ann$pValues_mw)

input_Manhattan$chr <- as.character(input_Manhattan$chr)
input_Manhattan$chr[input_Manhattan$chr == "X"] <- "23"
input_Manhattan$chr[input_Manhattan$chr == "Y"] <- "24"
input_Manhattan$chr <- as.numeric(input_Manhattan$chr)

input_Manhattan <- input_Manhattan[is.finite(input_Manhattan$pval) & input_Manhattan$pval > 0 & input_Manhattan$pval <= 1, ]

kable(head(input_Manhattan))

# png(filename = "ManhattanPlot.png", width = 800, height = 600)
par(mfrow=c(1,1))
manhattan(input_Manhattan, snp = "id", chr = "chr", bp = "map", p = "pval", suggestiveline = -log10(0.05), col = rainbow(24), annotateTop = F)
dev.off()

# 12
input_heatmap = as.matrix(final_mw[1:100, 1:8])

wt <- targets[targets$Group == "WT", "Basename"]
wt <- sub("^[^_]+_", "", wt)
group_color = c()
i = 1

for (name in colnames(beta)) {
  if (name %in% wt) { group_color[i] = "#FF6600" }
  else { group_color[i] = "#CC00CC" }
  i = i + 1
}

create_heatmap <- function(input_heatmap, group_color, targets, filename, main_title, linkage_method) {
  png(filename = filename, width = 800, height = 600)
  heatmap.2(
    input_heatmap, col = plasma(100), Rowv = TRUE, Colv = TRUE,
    hclustfun = function(x) hclust(x, method = linkage_method),
    dendrogram = "both", key = TRUE, ColSideColors = group_color,
    density.info = "none", trace = "none", scale = "none",
    symm = FALSE, main = main_title,
    key.xlab = 'beta-val', key.title = NA, keysize = 1, labRow = NA
  )
  legend("topright", legend = levels(targets$Group),
         col = c('#CC00CC', '#FF6600'), pch = 19, cex = 0.7)
  dev.off()
}

# note that the three .png files are going to be saved
create_heatmap(input_heatmap, group_color, targets, "CompleteLinkage.png", "Complete linkage", "complete")
create_heatmap(input_heatmap, group_color, targets, "SingleLinkage.png", "Single linkage", "single")
create_heatmap(input_heatmap, group_color, targets, "AverageLinkage.png", "Average linkage", "average")
