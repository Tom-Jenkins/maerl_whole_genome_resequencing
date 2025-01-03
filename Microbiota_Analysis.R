# =========================== #
#
# Maerl Whole Genome Re-sequencing Project 2024
#
# Metagenomics: Classify Unmapped Reads
#
# =========================== #

# In RStudio set working directory to the path where this R script is located
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
set.seed(123)

# Load packages
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(performance)
library(see)
library(car)

# Read in Kraken2 report output file
kraken <- read_csv("data/kraken_results/allsamples_kraken2.csv")
kraken

# Add column for site
kraken$Site <- substr(kraken$Prefix, 1, 3)

# Filter out sites not for analysis
kraken <- kraken |> filter(Site != "Gri", Site != "Her", Site != "Wey", Site != "Biz")

# Change order of columns
kraken <- dplyr::select(kraken, Prefix, Site, Rank_code, everything()) 

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
    Bacteria = mean(Bacteria_perc),
    # Bacteria_sd = sd(Bacteria_perc),
    Archaea = mean(Archaea_perc),
    # Archaea_sd = sd(Archaea_perc)
  ) |> 
  pivot_longer(cols = 2:3, names_to = "Domain", values_to = "Percentage_mean")
kraken_D_site

# Stats
kraken_D_site |> filter(Domain == "Bacteria") |> arrange(Percentage_mean)
kraken_D_site |> filter(Domain == "Archaea") |> arrange(Percentage_mean)

# Plot site means
ggplot(data = kraken_D_site)+
  geom_bar(
    mapping = aes(x = Site, y = Percentage_mean, fill = Domain),
    stat = "identity", position="dodge"
  )+
  scale_y_continuous(expand = c(0,0), limits = c(0, 100))+
  ylab("Microbiota per site average (%)")+
  theme(
    axis.title.x = element_blank()
  )

# Export to supporting information
ggsave("supporting_figures/microbiota_bacteria_versus_archaea.png")

#--------------#
# Main analysis ####
#--------------#

# Select groups of interest (Class level or above)
groups <- c("Pseudomonadota", "Alphaproteobacteria", "Gammaproteobacteria",
            "Actinomycetota", "Cyanobacteriota", "Flavobacteriia",
            "Nitrososphaerota", "Nitrospirota")

# Extract all records for groups of interest
kraken_grps <- kraken |> 
  filter(Rank_code == "D1" | Rank_code == "D2" | Rank_code == "D3" | Rank_code == "P" | Rank_code == "C") |> 
  filter(Scientific_name != "environmental samples") |>
  pivot_wider(id_cols = 1:2, names_from = "Scientific_name", values_from = "Number_of_fragments") |> 
  select(Prefix, Site, all_of(groups))
kraken_grps

# Add column for total Bacteria and Archaea (to calculate proportions next)
kraken_grps$Total <- kraken_D$Bacteria + kraken_D$Archaea

# Add proportions of each group by dividing by the total
kraken_grps_per <- mutate(kraken_grps, across(all_of(groups), ~ . / Total * 100))

# Convert to format for plotting individuals
kraken_grps_per_ind <- kraken_grps_per |> 
  select(-Total) |> 
  pivot_longer(cols = 3:(ncol(kraken_grps_per)-1), names_to = "Taxa", values_to = "Percentage")
kraken_grps_per_ind

# Plot all samples
ggplot(data = kraken_grps_per_ind)+
  geom_bar(
    mapping = aes(x = Prefix, y = Percentage, fill = Taxa),
    stat = "identity", position="dodge"
  )+
  scale_y_continuous(expand = c(0,0))+
  theme(
    axis.text.x = element_text(angle = 90)
  )

# Change Mil samples to Mil1 or Mil2 depending on individual name
# kraken_grps_per_ind$Site[which(str_detect(kraken_grps_per_ind$Prefix, "Mil1"))] <- "Mil1"
# kraken_grps_per_ind$Site[which(str_detect(kraken_grps_per_ind$Prefix, "Mil2"))] <- "Mil2"

# Convert to format for plotting site means
kraken_grps_per_site <- kraken_grps_per_ind |> 
  group_by(Site, Taxa) |> 
  summarise(Percentage_mean = mean(Percentage), Percentage_sd = sd(Percentage))
kraken_grps_per_site

#--------------#
# Visualisation ####
#--------------#

# Parameters
segment_x1 <- 0.40
segment_x2 <- 6.46
segment_y <- -2.0
text_x <- (segment_x1+segment_x2)/2
text_y <- segment_y - 1.10
linetype <- 2
linecol <- "grey50"
textcol <- "grey50"
colours <- c("deeppink","#FFD300","#F8DE7E","#FFF100","red3","red4")

# Filter data frame
kraken_vis <- filter(kraken_grps_per_site, Taxa != "Pseudomonadota")

# Reorder domain groups
kraken_vis$Taxa <- factor(
   x = kraken_vis$Taxa,
   levels = c("Alphaproteobacteria","Gammaproteobacteria", "Actinomycetota","Cyanobacteriota", "Flavobacteriia", "Nitrospirota","Nitrososphaerota"),
   labels = c("Alphaproteobacteria","Gammaproteobacteria", "Actinobacteria", "Cyanobacteria", "Flavobacteria", "Nitrospirota","Nitrososphaerota")
)

# Reorder maerl bed groups to represent classification
kraken_vis$Site <- factor(
  x = kraken_vis$Site,
  levels = c("Maw","Hel","Mil","Man","Swa","Bem"),
  labels = c("St Mawes (A1)", "Helford River (A3)", "Milford Haven (B1)",
             "The Manacles (B3)", "Swanage (B3)", "Bembridge (C2)")
)

# Plot site means
ggplot(data = kraken_vis, aes(x = Taxa, y = Percentage_mean, fill = Site))+
  geom_bar(stat = "identity", position = position_dodge(), colour = "black", linewidth = 0.1)+
  geom_errorbar(
    mapping = aes(ymin = pmax(0, Percentage_mean-Percentage_sd), ymax = Percentage_mean+Percentage_sd),
    width = 0.0, position = position_dodge(0.90), linewidth = 0.2
  )+
  coord_cartesian(expand = FALSE, ylim = c(0, 30), clip = "off")+
  scale_y_continuous(breaks = seq(0,30,10))+
  annotate("segment", x = segment_x1, xend = segment_x2, y = segment_y, yend = segment_y, linetype = linetype, colour = linecol, linewidth = 0.2)+
  annotate("text", x = text_x, y = text_y, label = "Bacteria", colour = textcol, size = 3.2)+
  annotate("segment", x = segment_x2+0.1, xend = segment_x2+1.01, y = segment_y, yend = segment_y, linetype = linetype, colour = linecol, linewidth = 0.2)+
  annotate("text", x = (segment_x2+0.1+segment_x2+1.01)/2, y = text_y, label = "Archaea", colour = textcol, size = 3.2)+
  # annotate("segment", x = segment_x2+1.11, xend = segment_x2+2.01, y = segment_y, yend = segment_y, linetype = linetype, colour = linecol)+
  # annotate("text", x = (segment_x2+1.11+segment_x2+2.01)/2, y = text_y, label = "Both", colour = textcol)+
  scale_fill_manual(values = colours)+
  ylab("Microbiota per site average (%)\n")+
  theme(
    plot.margin = margin(t = 10, l = 10, r = 15, b = 40),
    panel.background = element_rect(fill = "white"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(colour = "grey90"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 9, face = "bold"),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
  )
ggsave("figures/Figure_07.png", width = 10, height = 7, units = "in", dpi = 600)
ggsave("figures/Figure_07.pdf", width = 10, height = 7, units = "in")


#--------------#
# Modelling ####
#--------------#

# Type II ANOVA
# https://x.com/donaldmcknight2/status/1853508145286414381
?car::Anova

# Statistics
kraken_grps_per_ind |> filter(Taxa == "Pseudomonadota") |> pull(Percentage) |> range()
kraken_grps_per_ind |> filter(Taxa == "Alphaproteobacteria") |> pull(Percentage) |> range()
kraken_grps_per_ind |> filter(Taxa == "Gammaproteobacteria") |> pull(Percentage) |> range()
kraken_grps_per_ind |> filter(Taxa == "Actinomycetota") |> pull(Percentage) |> range()
kraken_grps_per_ind |> filter(Taxa == "Cyanobacteriota") |> pull(Percentage) |> range()
kraken_grps_per_ind |> filter(Taxa == "Flavobacteriia") |> pull(Percentage) |> range()
kraken_grps_per_ind |> filter(Taxa == "Nitrospirota") |> pull(Percentage) |> range()
kraken_grps_per_ind |> filter(Taxa == "Nitrososphaerota") |> pull(Percentage) |> range()

# ANOVA: Pseudomonadota model percentage as a function of site (factor)
model1 <- aov(Percentage ~ Site, data = subset(kraken_grps_per_ind, Taxa == "Pseudomonadota"))
performance::check_model(model1)
summary(model1)
TukeyHSD(model1, conf.level=.95)
plot(TukeyHSD(model1, conf.level=.95), las = 2)

# ANOVA: Flavobacteriaceae model percentage as a function of site (factor)
model2 <- aov(Percentage ~ Site, data = subset(kraken_grps_per_ind, Taxa == "Flavobacteriia"))
performance::check_model(model2)
summary(model2)
TukeyHSD(model2)
plot(TukeyHSD(model2, conf.level=.95), las = 2)

# Subset data for Nitrospirota model
nitro1 <- subset(kraken_grps_per_ind, Taxa == "Nitrospirota")
# nitro1$Site <- factor(nitro1$Site, levels = c("Maw","Man","Bem"))

# ANOVA: Nitrospirota model percentage as a function of site (factor)
model3 <- aov(Percentage ~ Site, data = nitro1)
performance::check_model(model3)
summary(model3)
TukeyHSD(model3)
plot(TukeyHSD(model3, conf.level=.95), las = 2)

# Subset data for Nitrospirota model
nitro2 <- subset(kraken_grps_per_ind, Taxa == "Nitrososphaerota")
# nitro2$Site <- factor(nitro1$Site, levels = c("Maw","Man","Bem"))

# ANOVA: Nitrospirota model percentage as a function of site (factor)
model4 <- aov(Percentage ~ Site, data = nitro2)
performance::check_model(model4)
summary(model4)
TukeyHSD(model4)
plot(TukeyHSD(model4, conf.level=.95), las = 2)

