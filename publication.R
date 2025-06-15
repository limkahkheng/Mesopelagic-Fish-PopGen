#Title: Mitochondrial DNA markers reveal panmixia and new sister lineage in Red Sea mesopelagic fish
#Author: Kah Kheng LIM
#Date: 29 May 2025

## Load required packages
library(tidyverse)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(rstatix)
library(car)
library(patchwork)

## Load data
allfish_data<-read.csv("/Users/limk/Desktop/Ongoing manuscript/Popgen/R script for length-weight relationship/Fish_LW_data_noBP.csv")

#################################################
################## Figure 2 #####################
#################################################
## Scatter plots with LOESS smoothing
raw_data <- ggplot(data = allfish_data, aes(x = SL, y = Weight)) +
  geom_point() +
  geom_smooth(method = "loess", level = 0.95) +
  facet_wrap(~Species, scales = "free") +
  theme_minimal() +  # Apply a minimal theme
  geom_vline(xintercept = 0, linetype = "solid", color = "black") +  # Vertical line at x = 0
  geom_hline(yintercept = 0, linetype = "solid", color = "black")   # Horizontal line at y = 0
raw_data

## Prepare species list
Species.list<-unique(allfish_data$Species)
Species.list<-Species.list[order(Species.list)]

## Initialize output dataframe for LWR parameters
param.output<-matrix(nrow=2,ncol=12)
rownames(param.output)<-Species.list
colnames(param.output)<-c("N","Size.min","Size.max",
                          "a","a.lower","a.upper",
                          "b","b.lower","b.upper",
                          "r2","Weight.Min","Weight.Max")

## Loop through species to fit LWR models
for (i in 1:length(Species.list)){
  data.sp <- allfish_data %>% filter(Species == Species.list[i]) %>%
    select(Species, SL, Weight) %>% droplevels()
  
  # Descriptive data
  count.sp<-nrow(data.sp)
  min.Size.sp<-min(data.sp$SL)
  max.Size.sp<-max(data.sp$SL)
  # Save the data
  param.output[i,1]<-count.sp
  param.output[i,2]<-min.Size.sp
  param.output[i,3]<-max.Size.sp
  
  # Descriptive data
  count.sp<-nrow(data.sp)
  min.Weight.sp<-min(data.sp$Weight)
  max.Weight.sp<-max(data.sp$Weight)
  # Save the data
  param.output[i,11]<-min.Weight.sp
  param.output[i,12]<-max.Weight.sp
  
  # Model
  model<-nls(Weight~I(exp(1)^(a+b*log(SL))),
             start=list(a=0,b=1),
             data=data.sp)
  
  #simplified code for model
  model<-nls(Weight ~ exp(a+b*log(SL)),
             start=list(a=0,b=1),
             data=data.sp)
  
  # Extract the parameters
  a.sp<-exp(coef(model)[1])
  a.lower.sp<-exp(confint(model)[1,1])
  a.upper.sp<-exp(confint(model)[1,2])
  
  b.sp<-coef(model)[2]
  b.lower.sp<-confint(model)[2,1]
  b.upper.sp<-confint(model)[2,2]
  
  # Calculate the r2
  RSS.sp<-deviance(model)
  
  Added.Sq.total<-sum(data.sp$Weight^2,na.rm=T)
  CF.total<-((sum(data.sp$Weight,na.rm=T))^2)/nrow(data.sp)
  SS.total<-Added.Sq.total-CF.total
  
  # So the SS model is:
  SS.sp<-SS.total-RSS.sp
  r2.sp<-SS.sp/SS.total
  
  # Save the parameters:
  
  param.output[i,4]<-a.sp
  param.output[i,5]<-a.lower.sp
  param.output[i,6]<-a.upper.sp
  param.output[i,7]<-b.sp
  param.output[i,8]<-b.lower.sp
  param.output[i,9]<-b.upper.sp
  param.output[i,10]<-r2.sp
}

param.output<-as.data.frame(param.output)
param.output$Species<-rownames(param.output)

# Print the histogram
plot_histogram <- function(data, species_name, fill_color) {
  ggplot(data, aes(x = SL)) +
    geom_histogram(fill = fill_color, color = "black", binwidth = 1.5) +
    labs(
      x = "Standard Length (mm)", y = "Frequency",
      title = paste("Histogram of SL -", species_name)
    ) +
    theme_minimal()
}

bp_hist <- plot_histogram(subset(allfish_data, Species == "Benthosema pterotum"), "Benthosema", "pink")
vm_hist <- plot_histogram(subset(allfish_data, Species == "Vinciguerria mabahiss"), "Vinciguerria", "skyblue")

print(bp_hist)
print(vm_hist)

#################################################
######## Morphometric stats starts here ########
#################################################
# Factor levels for consistency
allfish_data$Region <- factor(allfish_data$Region, levels = c("North", "NorthCentral", "Central", "SouthCentral", "South"))
allfish_data$Species <- factor(allfish_data$Species, levels = c("Vinciguerria mabahiss", "Benthosema pterotum"))

##Define color for each species
species_colors <- c("Vinciguerria mabahiss" = "#3B7DB5", "Benthosema pterotum" = "#FF69B4")

#################################################
########### SL comparison by region #############
#################################################
## ANOVA by Region for each species
anova_SL <- allfish_data %>%
  group_by(Species) %>%
  anova_test(SL ~ Region) %>%
  get_anova_table()
print(anova_SL)

## Post-hoc Tukey (Bonferroni adjusted)
tukey_SL <- allfish_data %>%
  group_by(Species) %>%
  emmeans_test(SL ~ Region, p.adjust.method = "bonferroni") %>%
  add_xy_position(x = "Region", dodge = 0.8)

## Boxplot for SL with significance labels
bxp_SL <- ggboxplot(
  allfish_data, x = "Region", y = "SL",
  color = "Species", palette = species_colors,  # apply colour
  ylab = "Standard length (mm)", xlab = "Region",
  add = "jitter", facet.by = "Species"
)
a<-bxp_SL + 
  stat_pvalue_manual(tukey_SL, label = "p.adj.signif", tip.length = 0.01) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#################################################
######### Weight comparison by region ###########
#################################################
## ANOVA by Weight for each species
anova_W <- allfish_data %>%
  group_by(Species) %>%
  anova_test(Weight ~ Region) %>%
  get_anova_table()
print(anova_W)

## Post-hoc Tukey (Bonferroni adjusted)
tukey_W <- allfish_data %>%
  group_by(Species) %>%
  emmeans_test(Weight ~ Region, p.adjust.method = "bonferroni") %>%
  add_xy_position(x = "Region", dodge = 0.8)

## Boxplot for Weight with significance labels
bxp_W <- ggboxplot(
  allfish_data, x = "Region", y = "Weight",
  color = "Species", palette = species_colors,  # apply colour
  ylab = "Wet weight (g)", xlab = "Region",
  add = "jitter", facet.by = "Species"
)
b<-bxp_W +  
  stat_pvalue_manual(tukey_W, label = "p.adj.signif", tip.length = 0.01) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#################################################
########### SL comparison by species ############
#################################################
## Summary stats
allfish_data %>%
  group_by(Species) %>%
  get_summary_stats(SL, type ="mean_sd")

## Boxplot for SL by species
bxp_SL <- ggboxplot(
  allfish_data, x = "Species", y = "SL", 
  color = "Species", palette = species_colors,  # apply colour
  ylab = "Standard length (mm)", xlab = "Species", add = "jitter"
)

## Levene's Test for equal variance
leveneTest(SL ~ Species, data = allfish_data)

## Welch's t-test
stat.test_SL <- allfish_data %>%
  t_test(SL ~ Species) %>%
  add_significance() %>%
  add_xy_position(x = "Species")

## Boxplot for SL by species
c<-bxp_SL + 
  stat_pvalue_manual(stat.test_SL, tip.length = 0) +
  labs(subtitle = get_test_label(stat.test_SL, detailed = TRUE))

#################################################
######### Weight comparison by species ##########
#################################################
## Summary stats
allfish_data %>%
  group_by(Species) %>%
  get_summary_stats(Weight, type ="mean_sd")

## Boxplot for Weight by species
bxp_W <- ggboxplot(
  allfish_data, x = "Species", y = "Weight", 
  color = "Species", palette = species_colors,  # apply colour
  ylab = "Wet weight (g)", xlab = "Species", add = "jitter"
)

## Levene's Test for equal variance
leveneTest(Weight ~ Species, data = allfish_data)

## Welch's t-test
stat.test_W <- allfish_data %>%
  t_test(Weight ~ Species) %>%
  add_significance() %>%
  add_xy_position(x = "Species")

## Boxplot for Weight by species
d<-bxp_W + 
  stat_pvalue_manual(stat.test_W, tip.length = 0) +
  labs(subtitle = get_test_label(stat.test_W, detailed = TRUE))

#################################################
########## Morphometric stats end here #########
#################################################

#################################################
################# Appendix S3 ###################
#################################################
gridExtra::grid.arrange(c,d,ncol=2)
gridExtra::grid.arrange(a,b,ncol=2)

