# =========================== #
#
# Maerl Whole Genome Resequencing Project 2024
#
# Metagenomics: Classify Unmapped Reads
#
# Figure 2
#
# =========================== #

# In RStudio set working directory to the path where this R script is located
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load packages
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(performance)
library(see)

# Read in Kraken2 report output file
kraken <- read_csv("data/allsamples_kraken2.csv")
kraken

# Add column for site
kraken$Site <- substr(kraken$Prefix, 1, 3)
kraken

# Change order of columns
kraken <- select(kraken, Prefix, Site, Rank_code, everything()) 

#--------------#
# Domain ####
#--------------#

# Extract all records for Domain, convert to wide format and (re)calculate %
kraken_D <- kraken |> 
  filter(Rank_code == "D") |> 
  pivot_wider(id_cols = 1:3, names_from = "Scientific_name", values_from = "Number_of_fragments") |> 
  select(!Eukaryota) |> 
  mutate(Bacteria_perc = Bacteria / (Bacteria + Archaea) * 100) |> 
  mutate(Archaea_perc = Archaea / (Bacteria + Archaea) * 100)
kraken_D

# Convert to format for plotting
kraken_D_ind <- kraken_D |> 
  pivot_longer(cols = 6:7, names_to = "Domain", values_to = "Percentage")
kraken_D_ind

# Plot all samples
ggplot(data = kraken_D_ind)+
  geom_bar(
    mapping = aes(x = Prefix, y = Percentage, fill = Domain),
    stat = "identity", position="dodge"
  )+
  scale_y_continuous(expand = c(0,0))+
  theme(
    axis.text.x = element_text(angle = 90)
  )

# Convert to mean % Bacteria and Archaea per site
kraken_D_site <- kraken_D |> 
  group_by(Site) |> 
  summarise(
    Bacteria_mean = mean(Bacteria_perc),
    # Bacteria_sd = sd(Bacteria_perc),
    Archaea_mean = mean(Archaea_perc),
    # Archaea_sd = sd(Archaea_perc)
  ) |> 
  pivot_longer(cols = ends_with("mean"), names_to = "Domain", values_to = "Percentage_mean")
kraken_D_site

# Plot site means
ggplot(data = kraken_D_site)+
  geom_bar(
    mapping = aes(x = Site, y = Percentage_mean, fill = Domain),
    stat = "identity", position="dodge"
  )+
  scale_y_continuous(expand = c(0,0), limits = c(0, 100))

#--------------#
# Groups: D1/D2/D3/P ####
#--------------#

# Extract all records for D1, convert to wide format and add total number of fragments column
kraken_grps <- kraken |> 
  filter(Rank_code == "D1" | Rank_code == "D2" | Rank_code == "D3" | Rank_code == "P") |> 
  filter(Scientific_name != "environmental samples") |> 
  pivot_wider(id_cols = 1:2, names_from = "Scientific_name", values_from = "Number_of_fragments")
kraken_grps

# Select groups of interest
groups <- c("Nitrososphaerota", "Pseudomonadota", "Bacteroidota", "Nitrospirota")
kraken_grps <- select(kraken_grps, Prefix, Site, all_of(groups))
kraken_grps

# Add Bacteria and Archaea total to data frame
kraken_grps$Total <- kraken_D$Bacteria + kraken_D$Archaea
kraken_grps

# Calculate percentages for each group
kraken_grps_per <- kraken_grps |> 
  mutate(Nitrososphaerota_per = Nitrososphaerota / Total * 100) |> 
  mutate(Pseudomonadota_per = Pseudomonadota / Total * 100) |> 
  mutate(Bacteroidota_per = Bacteroidota / Total * 100) |> 
  mutate(Nitrospirota_per = Nitrospirota / Total * 100) |> 
  mutate(Other_per = (Total - (Nitrososphaerota + Pseudomonadota + Bacteroidota + Nitrospirota)) / Total * 100) |> 
  select(Prefix, Site, contains("per"))
kraken_grps_per

# Convert to format for plotting individuals
kraken_grps_per_ind <- kraken_grps_per |> 
  pivot_longer(cols = 3:7, names_to = "Domain", values_to = "Percentage")
kraken_grps_per_ind

# Plot all samples
ggplot(data = kraken_grps_per_ind)+
  geom_bar(
    mapping = aes(x = Prefix, y = Percentage, fill = Domain),
    stat = "identity", position="dodge"
  )+
  scale_y_continuous(expand = c(0,0))+
  theme(
    axis.text.x = element_text(angle = 90)
  )

# Convert to format for plotting site means
kraken_grps_per_site <- kraken_grps_per_ind |> 
  group_by(Site, Domain) |> 
  summarise(Percentage_mean = mean(Percentage), Percentage_sd = sd(Percentage))
kraken_grps_per_site

# Reorder domain groups
kraken_grps_per_site$Domain <- factor(
   x = kraken_grps_per_site$Domain,
   levels = c("Pseudomonadota_per","Bacteroidota_per","Nitrospirota_per","Nitrososphaerota_per","Other_per"),
   labels = c("Pseudomonadota","Bacteroidota","Nitrospirota","Nitrososphaerota","Other")
)

# Reorder maerl bed groups
kraken_grps_per_site$Site <- factor(
  x = kraken_grps_per_site$Site,
  levels = c("Maw","Hel","Man","Bem"),
  labels = c("St Mawes Bank", "Helford River", "The Manacles", "Bembridge")
)

# Parameters
segment_x1 <- 0.55
segment_x2 <- 3.45
segment_y <- -6
text_x <- (segment_x1+segment_x2)/2
text_y <- segment_y - 2.5
linetype <- 2
linecol <- "grey50"
textcol <- "grey50"

# Plot site means
ggplot(data = kraken_grps_per_site, aes(x = Domain, y = Percentage_mean, fill = Site))+
  geom_bar(stat = "identity", position = position_dodge(), colour = "black", linewidth = 0.1)+
  geom_errorbar(
    mapping = aes(ymin = pmax(0, Percentage_mean-Percentage_sd), ymax = Percentage_mean+Percentage_sd),
    width = 0.0, position = position_dodge(0.90),
  )+
  coord_cartesian(expand = FALSE, ylim = c(0, 70), clip = "off")+
  scale_y_continuous(breaks = seq(0,100,10))+
  annotate("segment", x = segment_x1, xend = segment_x2, y = segment_y, yend = segment_y, linetype = linetype, colour = linecol)+
  annotate("text", x = text_x, y = text_y, label = "Bacteria", colour = textcol)+
  annotate("segment", x = segment_x2+0.1, xend = segment_x2+1.01, y = segment_y, yend = segment_y, linetype = linetype, colour = linecol)+
  annotate("text", x = (segment_x2+0.1+segment_x2+1.01)/2, y = text_y, label = "Archaea", colour = textcol)+
  annotate("segment", x = segment_x2+1.11, xend = segment_x2+2.01, y = segment_y, yend = segment_y, linetype = linetype, colour = linecol)+
  annotate("text", x = (segment_x2+1.11+segment_x2+2.01)/2, y = text_y, label = "Both", colour = textcol)+
  scale_fill_manual(values = c("deeppink", "pink", "red4", "grey60"))+
  ylab("Microbes Detected Per Site (mean %)\n")+
  theme(
    plot.margin = margin(t = 10, l = 10, r = 15, b = 50),
    panel.background = element_rect(fill = "white"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(colour = "grey90"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 12, face = "bold"),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
  )

ggsave("Figure_02.png", width = 10, height = 7, units = "in", dpi = 600)
ggsave("Figure_02.pdf", width = 10, height = 7, units = "in")

#--------------#
# Modelling ####
#--------------#

# Statistics
kraken_grps_per_ind |> filter(Domain == "Pseudomonadota_per") |> pull(Percentage) |> range()
kraken_grps_per_ind |> filter(Domain == "Bacteroidota_per") |> pull(Percentage) |> range()
kraken_grps_per_ind |> filter(Domain == "Nitrospirota_per") |> pull(Percentage) |> range()
kraken_grps_per_ind |> filter(Domain == "Nitrososphaerota_per") |> pull(Percentage) |> range()

# ANOVA: Bacterioidota model percentage as a function of site (factor)
model1 <- lm(Percentage ~ Site, data = subset(kraken_grps_per_ind, Domain == "Bacteroidota_per"))
performance::check_model(model1)
summary(model1)

# Subset data for Nitrospirota model
nitro1 <- subset(kraken_grps_per_ind, Domain == "Nitrospirota_per")
nitro1$Site <- factor(nitro1$Site, levels = c("Maw","Man","Bem"))

# ANOVA: Nitrospirota model percentage as a function of site (factor)
model2 <- lm(Percentage ~ Site, data = nitro1)
performance::check_model(model2)
summary(model2)

# Subset data for Nitrospirota model
nitro2 <- subset(kraken_grps_per_ind, Domain == "Nitrososphaerota_per")
nitro2$Site <- factor(nitro1$Site, levels = c("Maw","Man","Bem"))

# ANOVA: Nitrospirota model percentage as a function of site (factor)
model3 <- lm(Percentage ~ Site, data = nitro2)
performance::check_model(model3)
summary(model3)
