# =========================== #
#
# Maerl Whole Genome Re-sequencing Project 2024
#
# SNP Genetic Differentiation
#
# Species:
# Phymatolithon calcareum
# Lithothamnion corallioides
#
# =========================== #

# In RStudio set working directory to the path where this R script is located
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
set.seed(123)

# Load packages
library(vcfR)
library(adegenet)
library(stringr)
library(dplyr)
library(tidyr)
library(tibble)
library(poppr)
library(hierfstat)
library(mmod)
library(pegas)
library(RColorBrewer)
library(scales)
library(patchwork)


# ----------------- #
# Phymatolithon calcareum ####
# ----------------- #

# Read in Pcalcareum SNPs
pcal_genind <- vcfR2genind(read.vcfR("outputs/pcalcareum_SNPs_Britain.vcf.gz"))
pcal_genind

# Merge the following samples into the same group (population) for Fst analysis
# • St Mawes coarse maerl from 2011 and 2022
# • St Mawes maerl from 2015 and 2022
# • All samples from St Austell area (St Austell Bay and Little Gribbin)
pcal_genind$pop <- as.factor(sub("\\_.*", "", indNames(pcal_genind)))
pcal_genind$pop <- as.factor(str_replace(pcal_genind$pop, "Maw22C|Maw11C", "MawC"))
pcal_genind$pop <- as.factor(str_replace(pcal_genind$pop, "Maw15|Maw22", "Maw"))
pcal_genind$pop <- as.factor(str_replace(pcal_genind$pop, "Gri|AusII", "Aus"))
summary(pcal_genind$pop)

# Remove single sample from Swanage
pcal_genind <- popsub(pcal_genind, exclude = c("Swa"))
summary(pcal_genind$pop)

# Compute pairwise Fst
# pcal_Fst <- hierfstat::genet.dist(pcal_genind, method = "WC84")
# save(pcal_Fst, file = "outputs/Pcalcareum_Fst.RData")
load("outputs/Pcalcareum_Fst.RData")
pcal_Fst |> round(2)


# ----------------- #
# Lithothamnion corallioides ####
# ----------------- #

# Read in Lcorallioides SNPs
lcor_genind <- vcfR2genind(read.vcfR("./outputs/lcorallioides_SNPs.vcf.gz"))

# Merge the following samples into the same group (population)
# • St Mawes maerl from 2015 and 2022
# • All samples from St Austell area (St Austell Bay and Little Gribbin)
lcor_genind$pop <- as.factor(sub("\\_.*", "", indNames(lcor_genind)))
lcor_genind$pop <- as.factor(str_replace(lcor_genind$pop, "Maw15|Maw22", "Maw"))
lcor_genind$pop <- as.factor(str_replace(lcor_genind$pop, "Aus[A-Z]+|Gri", "Aus"))
# lcor_genind$pop <- as.factor(str_replace(lcor_genind$pop, "Wey|Swa", "WeySwa"))
# lcor_genind$pop <- as.factor(str_replace(lcor_genind$pop, "Mil[1-2]*", "Mil"))
summary(lcor_genind$pop)

# Remove sites with fewer than three individuals
lcor_genind <- popsub(lcor_genind, exclude = "Tud")
summary(lcor_genind$pop)

# Compute pairwise Fst
# lcor_Fst <- hierfstat::genet.dist(lcor_genind, method = "WC84")
# save(lcor_Fst, file = "outputs/Lcorallioides_Fst.RData")
load("outputs/Lcorallioides_Fst.RData")
lcor_Fst |> round(3)


# ----------------- #
# Supporting Figure ####
# ----------------- #

# Function to plot Fst heatmap
plot_fst_heatmap <- function(mat, lab_order, digits = 2) {
  
  # Change order of rows and cols
  fst.mat = as.matrix(round(mat, digits))
  fst.mat1 = fst.mat[lab_order, ]
  fst.mat2 = fst.mat1[, lab_order]
  
  # Create a data.frame
  ind = which(upper.tri(fst.mat2), arr.ind = TRUE)
  fst.df = data.frame(Site1 = dimnames(fst.mat2)[[2]][ind[,2]],
                      Site2 = dimnames(fst.mat2)[[1]][ind[,1]],
                      Fst = fst.mat2[ ind ])
  
  # Keep the order of the levels in the data.frame for plotting 
  fst.df$Site1 = factor(fst.df$Site1, levels = unique(fst.df$Site1))
  fst.df$Site2 = factor(fst.df$Site2, levels = unique(fst.df$Site2))
  
  # Convert minus values to zero
  fst.df$Fst[fst.df$Fst < 0] = 0
  
  # Fst italic label
  fst.label = expression(italic("F")[ST])
  
  # Extract middle Fst value for gradient argument
  mid = max(fst.df$Fst) / 2
  
  # Plot heatmap
  ggplot(data = fst.df, aes(x = Site1, y = Site2, fill = Fst))+
    geom_tile(colour = "black")+
    geom_text(aes(label = Fst), color="black", size = 3)+
    scale_fill_gradient2(low = "blue", mid = "pink", high = "red",
                         midpoint = mid, name = fst.label,
                         limits = c(0, max(fst.df$Fst)))+
    scale_x_discrete(expand = c(0,0))+
    scale_y_discrete(expand = c(0,0), position = "right")+
    theme(axis.text = element_text(colour = "black", size = 10, face = "bold"),
          axis.title = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          legend.position = "right",
          legend.title = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 10),
          plot.title = element_text(hjust = 0.5, face = "bold.italic"),
    )
}

# Pcalareum
(plt1 <- plot_fst_heatmap(pcal_Fst, popNames(pcal_genind))+ggtitle("Phymatolithon calcareum"))
  
# Lcorallioides
(plt2 <- plot_fst_heatmap(lcor_Fst, popNames(lcor_genind))+ggtitle("Lithothamnion corallioides"))

# Composer
(fig <- plt1 + plt2)

# Export
ggsave("supporting_figures/Figure_SX_FST.png", dpi = 300)
