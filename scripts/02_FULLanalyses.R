################################################################################
#################  Chapter II Ecological successions ###########################
##################          BETAPART ANALYSES        ###########################
##################       Julie Morgane Guenat        ###########################
#################              March 2026             ###########################             
################################################################################

setwd("C:/Users/jguenat/OneDrive - Université de Lausanne/PhD_Analyses/Chapter_2_ecological_successions")

load("../Graph_settings.RData")

################################################################################
# 1. Packages ####

library(tidyverse)
library(glmmTMB)
library(vegan)
library(cowplot)
library(betapart)
library(ggpubr)
library(DHARMa)      
library(performance)
library(effects)
library(sjPlot)
library(ggeffects)
library(patchwork)
library(corrplot)
library(lubridate)
library(car)
library(emmeans)
library(broom.mixed) 
library(flextable)   
library(officer)
library(ggdist)

################################################################################
# 2. Load the datasets ####
################################################################################
## 2.1 Metadata ####
metadata <- read.csv("data/samples_with_environmental_data.csv", header=T, sep= ",")
#selecting the data I need: 
metadata <- metadata %>% dplyr::filter(Sampling_type=="Monthly" & filter_pore_size == "0.45")

## 2.2 Family datasets ####
taxo_data <- read.csv("../Cleaning_data/final_datasets/1.Euka_DNA_Taxo_verified_pres_abs.csv", header = T, sep = ",")

#Selecting aquatic family only
aquatic <- taxo_data %>% 
  dplyr::filter(aquatic == 1)

family_all <- aquatic %>%
  dplyr::select(c(5, 8:ncol(.))) %>%  # Select columns dynamically
  group_by(family) %>% 
  summarize(across(where(is.numeric), sum, na.rm = TRUE)) %>%  # Summarize numeric columns
  mutate(across(where(is.numeric), ~ ifelse(. > 0, 1, 0)))  # Apply mutation to numeric columns

family_all<- as.data.frame(family_all)
family_all<- family_all[-278,] #rm NA
row.names(family_all)<-family_all$family
family_all<- family_all[,-(1:3)] #rm family_name, aquatic, terrestrial
family_all<- t(family_all)
family_all<-as.data.frame(family_all)

family_all <- family_all %>%
  rownames_to_column("sample_ID") %>%
  select(!"Percidae")


## 2.3 adding the metadata ####
# 1st was add the correct names
samp <- read.csv("../Cleaning_data/final_datasets/sample_mapping_names.csv", header = T, sep = ";")
samples_pla <- read.csv("../Cleaning_data/sample_names.csv", sep= ";", header=T)

samp_unique <- samp[!duplicated(samp$Simplified_Name), ]
family_all <- merge(samp_unique, family_all, by.x = "Simplified_Name", by.y = "sample_ID")


family_all <- merge(family_all, samples_pla[, c("X", "Samples_names")], 
                    by.x = "Sample_ID", by.y = "X", 
                    all.x = TRUE)

family_metadat <- merge(metadata, family_all, by.x = "Samples_ID", by.y = "Samples_names")


## 2.4 Column predators and nutrient ####
family_metadat_prednutr <- family_metadat %>%
  mutate(
    predator = case_when(
      Treatment %in% c("NP", "P") ~ 1,
      Treatment %in% c("none", "N") ~ 0
    ),
    nutrients = case_when(
      Treatment %in% c("NP", "N") ~ 1,
      Treatment %in% c("none", "P") ~ 0
    )
  )

family_metadat_prednutr$predator <- factor(family_metadat_prednutr$predator, levels = c(0, 1))  
family_metadat_prednutr$nutrients <- factor(family_metadat_prednutr$nutrients, levels = c(0, 1))

#write.csv(family_metadat_prednutr, "data/Family_metadata_pred_nutr.csv", row.names = F)


rm(aquatic, taxo_data, family_metadat_prednutr, family_all, family_metadat, metadata, samp, samples_pla, samp_unique)

################################################
# 3. BETA.MULTI analysis ####
################################################
#beta.multi provides an overall insight of whether the changes in sp. composition
#is due to turnover (sp. replacement --> beta.SIM) or to nestedness 
#(sp. gain or loss --> beta.SNE). 
#And the overall beta diversity (Sorensen, beta.SOR)

## 3.1 Single Lake analysis ####
# this part is only for me to test the method and is not meant to be run!
# To be able to analyse for all the lakes, I need first to understand how it works 
# for one single lake: 
# For one Lake, I've chose L5 since there is value for each of the month: 
l5_metadata <- metadata %>%
  filter(Lake_ID == "L5") %>%
  arrange(collection_date)

# Extract just the L5 community data
l5_sample_ids <- l5_metadata$sample_id
l5_comm_data <- community_data[l5_sample_ids, ]
dim(l5_comm_data)

#Check the number of species: 
species_counts <- rowSums(l5_comm_data)
print(species_counts)

#Perform the computation of the different component of beta diversity: 
beta_multi <- beta.multi(l5_comm_data, index.family = "sorensen")


## 3.2 beta.multi ALL Lakes analysis ####
## 3.3 Data preparation ####

data<-read.csv("data/Family_metadata_pred_nutr.csv", header = T, sep = ",")

family_cols <- 36:311
family_names <- names(data)[family_cols]

data <- data %>% 
  mutate(
    sample_id = Samples_ID,
    lake_id = Lake_ID,
    collection_date = as.Date(Collection_date, format = "%d.%m.%Y")
  ) %>%
  arrange(lake_id, collection_date)


metadata <- data %>%
  select(-all_of(family_cols)) 

is.Date(metadata$collection_date)

# Create community matrix (samples x families)
community_data <- data %>%
  select(Samples_ID, all_of(family_cols)) %>%
  column_to_rownames("Samples_ID")


# Create data frame to store lake-level summary with beta diversity metrics
beta_multi_dataset <- metadata %>%
  select(Lake_ID, predator, nutrients) %>%
  distinct() %>%
  arrange(Lake_ID)

# Add columns for beta diversity components
beta_multi_dataset$Beta_SOR <- NA
beta_multi_dataset$Beta_SIM <- NA  
beta_multi_dataset$Beta_SNE <- NA


# Loop through each lake to calculate beta diversity metrics
for(i in 1:nrow(beta_multi_dataset)) {
  lake <- beta_multi_dataset$Lake_ID[i]
  
  # Get data for this lake, sorted by date
  lake_metadata <- metadata %>%
    filter(Lake_ID == lake) %>%
    arrange(collection_date)
  
  # Extract community data for this lake
  lake_sample_ids <- lake_metadata$sample_id
  lake_comm_data <- community_data[lake_sample_ids, ]
  
  # Calculate beta.multi
  beta_multi_result <- beta.multi(lake_comm_data, index.family = "sorensen")
  
  # Add results to summary dataframe
  beta_multi_dataset$Beta_SOR[i] <- beta_multi_result$beta.SOR
  beta_multi_dataset$Beta_SIM[i] <- beta_multi_result$beta.SIM
  beta_multi_dataset$Beta_SNE[i] <- beta_multi_result$beta.SNE
}

rm(beta_multi_result, data, lake_comm_data, lake_metadata, 
   family_cols, family_names, i, lake, lake_sample_ids)

## 3.4 Statistics on the beta.multi ####

# We set the treatment as factor and ensure the levels are in the good order
beta_multi_dataset$predator<-factor(beta_multi_dataset$predator, levels = c("0", "1"))
beta_multi_dataset$nutrients<-factor(beta_multi_dataset$nutrients, levels = c("0", "1"))

# we add the index == being the nestedess or the turnover as a variable
# to be able to have a value of their contribution. 

beta_multi_index<-rbind(beta_multi_dataset, beta_multi_dataset)
beta_multi_index$index <- gl(2, 12, label=c("n", "t"))
beta_multi_index$value <- c(beta_multi_dataset$Beta_SNE, beta_multi_dataset$Beta_SIM)

### 3.4.2 INDEX and Treatment effects ####

m.0 <- glmmTMB(value ~index * predator * nutrients, 
               data = beta_multi_index, 
               family = beta_family(link = "logit"))

simulateResiduals(m.0, plot = T)

summary(m.0)
Anova(m.0)


m.1 <- glmmTMB(value ~index * predator * nutrients + (1|Lake_ID), 
               data = beta_multi_index, 
               family = beta_family(link = "logit"))

simulateResiduals(m.1, plot = T)

anova(m.0, m.1) 
# the random factor is not significant. Unlike for the beta.pair analysis
# we do not need to take into account the repeatability of measures on the 
# same lake. So we do not take this into account here. 

rm(m.1)

### 3.4.1 Turnovers: Simpson index ####
m.SIM.0 <- glmmTMB(Beta_SIM ~ predator * nutrients + (1 | Lake_ID), 
                   family = beta_family(link = "logit"),
                   data = beta_multi_dataset)

summary(m.SIM.0)
simulateResiduals(m.SIM.0, plot = T)


m.SIM.1 <- glmmTMB(Beta_SIM ~ predator * nutrients, 
                   family = beta_family(link = "logit"),
                   data = beta_multi_dataset)

#Checking if the lake_ID is necessary as random effect:
anova(m.SIM.0, m.SIM.1)
# no significant, removing it improve slightly the model AIC. 
simulateResiduals(m.SIM.1, plot = T)
summary(m.SIM.1)
Anova(m.SIM.1)

#Why to prefer Anova (particularly Type II or Type III) for factorial designs because:
# It tests the significance of terms as a whole rather than individual coefficients
# It's more appropriate for testing effects in models with interactions
# It follows the principle of marginality (testing higher-order terms before lower-order terms)
# It aligns with the traditional experimental design approach

# here our results suggest that removal of top-predator does not impact the turnover, 
# neither does the nutrient enrichment, BUT the interaction does! 
# this means that the effect of pred rm depends on the level of the nutrient enrichment.
# We observe a higher turnovers in lakes with pred rm and nutri add. 

rm(m.SIM.0)

### 3.4.2 Nestedness: beta.SNE ####

m.SNE.0 <- glmmTMB(Beta_SNE ~ predator * nutrients + (1 | Lake_ID), 
                   family = beta_family(link = "logit"),
                   data = beta_multi_dataset)
summary(m.SNE.0)
simulateResiduals(m.SNE.0, plot = T)


m.SNE.1 <- glmmTMB(Beta_SNE ~ predator * nutrients, 
                   family = beta_family(link= "logit"),
                   data = beta_multi_dataset)
summary(m.SNE.1)
simulateResiduals(m.SNE.1, plot = T)
Anova(m.SNE.1)

anova(m.SNE.0, m.SNE.1)

rm(m.SNE.0)

### 3.4.3 Overall: Sorenson ####

m.SOR.0 <- glmmTMB(Beta_SOR ~ predator * nutrients + (1 | Lake_ID), 
                   family = beta_family(link = "logit"),
                   data = beta_multi_dataset)
summary(m.SOR.0)
simulateResiduals(m.SOR.0, plot = T)

m.SOR.1 <- glmmTMB(Beta_SOR ~ predator * nutrients, 
                   family = beta_family(link = "logit"),
                   data = beta_multi_dataset)
summary(m.SOR.1)
simulateResiduals(m.SOR.1, plot = T)
anova(m.SOR.0, m.SOR.1)

Anova(m.SOR.1)

# nothing significant! 
rm(m.SOR.0)

################################################
# 4. BETA.MULTI Table and Visualization results ####
################################################
## 4.1 Tables ####
### 4.1.1 TABLE S4: Create summary statistics ####

m0_tidy <- tidy(m.0, conf.int = TRUE)
m0_tidy$model <- "Combined (Nestedness + Turnover)"

mSIM_tidy <- tidy(m.SIM.1, conf.int = TRUE)
mSIM_tidy$model <- "Turnover Component"

mSNE_tidy <- tidy(m.SNE.1, conf.int = TRUE)
mSNE_tidy$model <- "Nestedness Component"

mSOR_tidy <- tidy(m.SOR.1, conf.int = TRUE)
mSOR_tidy$model <- "Total Beta"

# Combine into a single table
model_results <- bind_rows(m0_tidy, mSIM_tidy, mSNE_tidy, mSOR_tidy)

rm(m0_tidy, mSIM_tidy, mSNE_tidy, mSOR_tidy)

# Rename columns for clarity and select relevant ones
model_results <- model_results %>%
  select(model, term, estimate, std.error, statistic, p.value, conf.low, conf.high) %>%
  rename(
    "Component" = model,
    "Term" = term,
    "Estimate" = estimate,
    "Std. Error" = std.error,
    "z value" = statistic,
    "p-value" = p.value,
    "Lower CI" = conf.low,
    "Upper CI" = conf.high
  )

# Add significance indicators
model_results$Significance <- ifelse(model_results$`p-value` < 0.001, '***',
                                     ifelse(model_results$`p-value` < 0.01, '**',
                                            ifelse(model_results$`p-value` < 0.05, '*',
                                                   ifelse(model_results$`p-value` < 0.1, '.', ''))))

# Combine p-value with significance
model_results$FormattedPvalue <- paste0(
  sprintf("%.3f", model_results$`p-value`), 
  " ", 
  model_results$Significance
)

# First select and rename columns in the data frame
model_results_formatted <- model_results %>%
  select(-Significance, -`p-value`) %>%  # Remove separate columns
  rename(`p-value` = FormattedPvalue)

# Then create the flextable with the formatted data frame
ft_models <- flextable(model_results_formatted) %>%
  colformat_double(digits = 3) %>%
  bold(i = model_results$`p-value` < 0.05, j = "p-value") %>%
  theme_vanilla() %>%
  add_header_lines("Beta Diversity Components Model Results") %>%
  autofit()

# Save to Word
doc <- read_docx()
doc <- body_add_flextable(doc, ft_models)
print(doc, target = "tables/TableS4.beta.multi_model_results.docx")

rm(ft_models, model_results, model_results_formatted)

### 4.1.2 TABLE S5: Create ANOVA results ####
# Run Anova on each model
m0_anova <- Anova(m.0)
mSIM_anova <- Anova(m.SIM.1)
mSNE_anova <- Anova(m.SNE.1)
mSOR_anova <- Anova(m.SOR.1)

# Convert Anova results to data frames
m0_anova_df <- as.data.frame(m0_anova)
m0_anova_df$Effect <- rownames(m0_anova_df)
rownames(m0_anova_df) <- NULL
m0_anova_df$Component <- "Combined (Nestedness + Turnover)"

mSIM_anova_df <- as.data.frame(mSIM_anova)
mSIM_anova_df$Effect <- rownames(mSIM_anova_df)
rownames(mSIM_anova_df) <- NULL
mSIM_anova_df$Component <- "Turnover Component"

mSNE_anova_df <- as.data.frame(mSNE_anova)
mSNE_anova_df$Effect <- rownames(mSNE_anova_df)
rownames(mSNE_anova_df) <- NULL
mSNE_anova_df$Component <- "Nestedness Component"

mSOR_anova_df <- as.data.frame(mSOR_anova)
mSOR_anova_df$Effect <- rownames(mSOR_anova_df)
rownames(mSOR_anova_df) <- NULL
mSOR_anova_df$Component <- "Total Beta"


# Combine into a single ANOVA results table
anova_results <- bind_rows(m0_anova_df, mSIM_anova_df, mSNE_anova_df, mSOR_anova_df)

rm(m0_anova, mSIM_anova, mSNE_anova, m0_anova_df, mSIM_anova_df, mSNE_anova_df, mSOR_anova_df)

# Rename columns for clarity
names(anova_results)[names(anova_results) == "Chisq"] <- "Chi-Square"
names(anova_results)[names(anova_results) == "Pr(>Chisq)"] <- "p-value"

# Reorder columns
anova_results <- anova_results[, c("Component", "Effect", "Chi-Square", "Df", "p-value")]

# Add significance indicators
anova_results$Significance <- ifelse(anova_results$`p-value` < 0.001, '***',
                                     ifelse(anova_results$`p-value` < 0.01, '**',
                                            ifelse(anova_results$`p-value` < 0.05, '*',
                                                   ifelse(anova_results$`p-value` < 0.1, '.', ''))))

# Combine p-value with significance
anova_results$FormattedPvalue <- paste0(
  sprintf("%.3f", anova_results$`p-value`), 
  " ", 
  anova_results$Significance
)

# First select and rename columns in the data frame
anova_results_formatted <- anova_results %>%
  select(-Significance, -`p-value`) %>%  # Remove separate columns
  rename(`p-value` = FormattedPvalue)

# Create a flextable
ft_anova <- flextable(anova_results_formatted) %>%
  colformat_double(digits = 3) %>%
  bold(i = anova_results$`p-value` < 0.05, j = "p-value") %>%
  theme_vanilla() %>%
  add_header_lines("ANOVA Results for Beta Diversity Components") %>%
  autofit()

# Save to Word
doc <- read_docx()
doc <- body_add_flextable(doc, ft_anova)
#print(doc, target = "tables/TableS5.beta.multi_anova_results.docx")

rm(doc, ft_anova, anova_results, anova_results_formatted)

## 4.2 Graphical representation ####
### 4.2.1 Interaction plots ####
# TURNOVERS
emm_sim <- emmeans(m.SIM.1, ~ predator * nutrients)
emm_sim_df <- as.data.frame(emm_sim)

p_sim_emm <- ggplot(emm_sim_df, aes(x = predator, y = emmean, group = nutrients, color = nutrients)) +
  geom_point(size = 2.5) +
  geom_line(linewidth = 1.1) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.1, linewidth = 1.1) +
  labs(x = "", y = "Estimated Marginal Means \n(logit)",
       tag = "B",
       color = "", 
       title = "Turnover") +
  scale_x_discrete(labels = predator_labels) +
  scale_color_manual(values = nutrient_colors, labels = nutrient_labels) +
  theme_jg_right +
  theme(legend.position = c(1, 1.2), 
        legend.justification = c("right", "top"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 11), 
        axis.title.x = element_text(size = 12, margin = margin(t = 15, b = 5)), # Add top margin
        axis.title.y = element_text(size = 12, margin = margin(r = 15, l = 5)), # Add right margin
        axis.text.y = element_text(size = 11),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11),
        title = element_text(size = 12),
        plot.margin = margin(3, 3, 3, 3))

p_sim_emm

# NESTEDNESS 
emm_sne <- emmeans(m.SNE.1, ~ predator * nutrients)
emm_sne_df <- as.data.frame(emm_sne)

p_sne_emm <- ggplot(emm_sne_df, aes(x = predator, y = emmean, group = nutrients, color = nutrients)) +
  geom_point(size = 2.5) +
  geom_line(linewidth = 1.1) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.1, linewidth = 1.1) +
  labs(x = "", y = "",
       tag = "",
       color = "", 
       title = "Nestedness") +
  scale_x_discrete(labels = predator_labels) +
  scale_color_manual(values = nutrient_colors, labels = nutrient_labels) +
  theme_jg_right +
  theme(legend.position = c(1.0, 1.2), 
        legend.justification = c("right", "top"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 11), 
        axis.title.x = element_text(size = 12, margin = margin(t = 15, b = 5)), # Add top margin
        axis.title.y = element_text(size = 12, margin = margin(r = 15, l = 5)), # Add right margin
        axis.text.y = element_text(size = 11),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11),
        title = element_text(size = 12), 
        plot.margin = margin(3, 3, 3, 3), )

p_sne_emm

#### 4.2.1.1 FIGURE S5: interaction plot tot beta #### 
# Total beta diversity  
emm_sor <- emmeans(m.SOR.1, ~ predator * nutrients)
emm_sor_df <- as.data.frame(emm_sor)

p_sor_emm <- ggplot(emm_sor_df, aes(x = predator, y = emmean, group = nutrients, color = nutrients)) +
  geom_point(size = 2) +
  geom_line(linewidth = 1) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.1, linewidth = 1) +
  labs(x = "", y = "Estimated Marginal Means \n(logit)",
       tag = "",
       color = "", 
       title = "Total beta Diversity",
       subtitle = "") +
  scale_x_discrete(labels = predator_labels) +
  scale_color_manual(values = nutrient_colors, labels = nutrient_labels) +
  theme_classic() + theme_jg_right +
  theme(legend.text = element_text(size = 12),      
        legend.position = c(1.0, 1.20), 
        legend.justification = c("right", "top")
  )

p_sor_emm

#ggsave("figure/S5.FigS5_Beta.multi_totDiversity.pdf", p_sor_emm, width = 6, height = 6, dpi = 300)
#ggsave("figure/S5.FigS5_Beta.multi_totDiversity.png", p_sor_emm, width = 6, height = 6, dpi = 300)


### 4.2.2 Stacked bar chart for beta div comp ####
# Create the summary_data data frame from beta_multi_dataset
summary_data <- beta_multi_dataset %>%
  group_by(predator, nutrients) %>%
  summarize(
    mean_SIM = mean(Beta_SIM),
    mean_SNE = mean(Beta_SNE),
    mean_SOR = mean(Beta_SOR),
    se_SIM = sd(Beta_SIM)/sqrt(n()),
    se_SNE = sd(Beta_SNE)/sqrt(n()),
    se_SOR = sd(Beta_SOR)/sqrt(n()),
    .groups = 'drop'
  )

# Create treatment combinations and restructure data for plotting
summary_long <- summary_data %>%
  # Create a treatment factor with the desired order
  mutate(
    treatment = case_when(
      predator == "1" & nutrients == "0" ~ "Pred+/Nutr-",
      predator == "0" & nutrients == "0" ~ "Pred-/Nutr-",
      predator == "1" & nutrients == "1" ~ "Pred+/Nutr+",
      predator == "0" & nutrients == "1" ~ "Pred-/Nutr+"
    ),
    # Create ordered factor
    treatment = factor(treatment, levels = c(
      "Pred-/Nutr-",
      "Pred-/Nutr+",
      "Pred+/Nutr-",
      "Pred+/Nutr+"
    ))
  ) %>%
  # Pivot to long format for component means
  select(treatment, mean_SIM, mean_SNE, mean_SOR, se_SIM, se_SNE, se_SOR) %>%
  pivot_longer(
    cols = starts_with("mean_"),
    names_to = "component",
    values_to = "value"
  ) %>%
  # Pivot to long format for component standard errors
  pivot_longer(
    cols = starts_with("se_"),
    names_to = "se_component",
    values_to = "se"
  ) %>%
  # Match means with their corresponding standard errors
  filter(substring(component, 6) == substring(se_component, 4)) %>%
  # Clean up component names
  mutate(
    component = case_when(
      component == "mean_SIM" ~ "Turnover",
      component == "mean_SNE" ~ "Nestedness",
      component == "mean_SOR" ~ "Total"
    ),
    component = factor(component, levels = c("Total", "Turnover", "Nestedness"))
  ) %>%
  # Remove the se_component column
  select(-se_component)

# Filter only Turnover and Nestedness for stacking, and Total for error bars
components_data <- summary_long %>%
  filter(component %in% c("Turnover", "Nestedness")) %>%
  mutate(component = factor(component, levels = c("Turnover", "Nestedness")))

total_data <- summary_long %>%
  filter(component == "Total")

# Create improved stacked bar chart
barplot_beta_comp <- ggplot() +
  # Base layer: stacked bars for Turnover and Nestedness
  geom_bar(
    data = components_data, 
    aes(x = treatment, y = value, fill = component),
    stat = "identity", 
    position = "stack", 
    width = 0.5
  ) +
  # Add error bars for total beta diversity
  geom_errorbar(
    data = total_data,
    aes(y = value, ymin = value - se, ymax = value + se, x = treatment),
    width = 0.2, 
    color = "black"
  ) +
  # Custom colors for components
  scale_fill_manual(values = component_colors[c("Turnover", "Nestedness")]) +
  # Labels
  labs(tag = "A",
       x = "",
       y = "Beta Diversity Value",
       fill = "Component"
  ) +
  # Theme customization
  theme_jg_bottom+ 
  theme(axis.text.x = element_text( size = 11), 
        legend.margin = margin(t = -15), 
        axis.title.x = element_text(size = 12, margin = margin(t = 15, b = 5)), # Add top margin
        axis.title.y = element_text(size = 12, margin = margin(r = 15, l = 5)), # Add right margin
        axis.text.y = element_text(size = 11),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11),
        title = element_text(size = 12),
        plot.margin = margin(3, 3, 3, 3))


barplot_beta_comp


### 4.2.3 FIGURE 1: Combine the plots with external panel labels ####

# Now combine the plots without any additional annotation
Fig.1 <- (barplot_beta_comp / (p_sim_emm + p_sne_emm)) +
  plot_layout(heights = c(1, 1.2))


# View the combined figure
Fig.1

# Save high-quality figure
ggsave("figure/R1.Fig1_Beta.multi_Diversity.pdf", Fig.1, width = 8, height = 7, dpi = 600)
ggsave("figure/R1.Fig1_Beta.multi_Diversity.png", Fig.1, width = 8, height = 7, dpi = 600)

rm(beta_multi_dataset, beta_multi_index, components_data, emm_sim_df, emm_sim, 
   emm_sne_df, emm_sne, m.0, m.SIM.1, m.SNE.1, summary_data, summary_long, total_data,
   p_sim_emm, p_sne_emm, barplot_beta_comp)


################################################
# 5. BETA.PAIR dataset preparation ####
################################################
# beta.pair shows how diversity component changes between specific time point. 
# gives dist matrices with pairwise comparison. 
# used to examine specific temporal transition or identify when the biggest 
# community change occurs. 

## 5.1 Single Lake analysis ####
# To be able to analyse for all the lakes, I need first to understand how it works 
# for one single lake: 
# For one Lake, I've chose L5 since there is value for each of the month: 
l5_metadata <- metadata %>%
  filter(Lake_ID == "L5") %>%
  arrange(collection_date)

# Extract just the L5 community data
l5_sample_ids <- l5_metadata$sample_id
l5_comm_data <- community_data[l5_sample_ids, ]
dim(l5_comm_data)

#Check the number of species: 
species_counts <- rowSums(l5_comm_data)
print(species_counts)

#Perform the computation of the different component of beta diversity: 
beta_pair <- beta.pair(l5_comm_data, index.family = "sorensen")


## 5.2 beta.pair computation ALL Lakes analysis ####

# Create data frame to store lake-level summary with beta diversity metrics
beta_pair_dataset <- metadata %>%
  select(Lake_ID, predator, nutrients) %>%
  distinct() %>%
  arrange(Lake_ID)

# we need to work with list here since the output of beta.pair is 
# distance matrices. 

beta_pair_results <- list()

# Calculate the beta.pair: 
# First loop: Just calculate beta.pair for each lake and store the results

for(i in 1:nrow(beta_pair_dataset)) {
  lake <- beta_pair_dataset$Lake_ID[i]
  
  # Get data for this lake, sorted by date
  lake_metadata <- metadata %>%
    filter(Lake_ID == lake) %>%
    arrange(collection_date)
  
  # Extract community data for this lake
  lake_sample_ids <- lake_metadata$sample_id
  lake_comm_data <- community_data[lake_sample_ids, ]
  
  # Try to calculate beta.pair within a tryCatch to handle errors
  tryCatch({
    # Calculate beta.pair
    beta_pair_result <- beta.pair(lake_comm_data, index.family = "sorensen")
    
    # Store in the list
    beta_pair_results[[lake]] <- list(
      metadata = lake_metadata,
      beta_pair = beta_pair_result,
      treatment = data.frame(
        predator = unique(lake_metadata$predator),
        nutrients = unique(lake_metadata$nutrients)
      )
    )
  })
}


## Extract pairwise distances for analysis (only if we have results)

# Create dataframes to store all pairwise comparisons
all_pairs_sor <- data.frame()
all_pairs_sim <- data.frame()
all_pairs_sne <- data.frame()


# Second loop: Process each lake with beta.pair results
for(lake in names(beta_pair_results)) {
  lake_data <- beta_pair_results[[lake]]
  
  # Get treatment info
  predator_status <- lake_data$treatment$predator
  nutrient_status <- lake_data$treatment$nutrients
  
  # Get dates for each sample
  dates <- lake_data$metadata$collection_date
  
  # Extract distance matrices
  sor_mat <- as.matrix(lake_data$beta_pair$beta.sor)
  sim_mat <- as.matrix(lake_data$beta_pair$beta.sim)
  sne_mat <- as.matrix(lake_data$beta_pair$beta.sne)
  
  # For each pair of samples - careful with the indices
  for(row_idx in 1:(nrow(sor_mat)-1)) {
    for(col_idx in (row_idx+1):ncol(sor_mat)) {
      
      # Calculate time difference between samples in days
      time_diff <- as.numeric(difftime(dates[col_idx], dates[row_idx], units = "days"))
      
      # Add to dataframes - using explicit indices from matrix
      pair_data <- data.frame(
        Lake_ID = lake,
        Sample1 = as.character(rownames(sor_mat)[row_idx]),
        Sample2 = as.character(colnames(sor_mat)[col_idx]),
        Date1 = dates[row_idx],
        Date2 = dates[col_idx],
        TimeDiff_days = time_diff,
        predator = predator_status,
        nutrients = nutrient_status
      )
      
      # Add specific distance metrics - FIXED: using row_idx and col_idx instead of i,j
      sor_row <- pair_data
      sor_row$distance <- sor_mat[row_idx, col_idx]
      sor_row$component <- "Total (SOR)"
      
      sim_row <- pair_data
      sim_row$distance <- sim_mat[row_idx, col_idx]
      sim_row$component <- "Turnover (SIM)"
      
      sne_row <- pair_data
      sne_row$distance <- sne_mat[row_idx, col_idx]
      sne_row$component <- "Nestedness (SNE)"
      
      # Append to the respective dataframes
      all_pairs_sor <- rbind(all_pairs_sor, sor_row)
      all_pairs_sim <- rbind(all_pairs_sim, sim_row)
      all_pairs_sne <- rbind(all_pairs_sne, sne_row)
    }
  }
}



# Combine all components into one dataframe for easier analysis

all_pairs_combined <- rbind(all_pairs_sor, all_pairs_sim, all_pairs_sne)

#write.csv(all_pairs_combined, "data/beta.pair_dataset.csv", row.names = F)

rm(all_pairs_combined, all_pairs_sim, all_pairs_sne, all_pairs_sor, 
   beta_pair_dataset, beta_pair_result, beta_pair_results, community_data, Fig.1, 
   lake_comm_data, lake_metadata, pair_data, sim_mat, sne_mat, sim_row, sne_row, 
   sor_row, sor_mat, lake_data, metadata, col_idx, dates, i, lake, lake_sample_ids,
   nutrient_status, predator_status, row_idx, time_diff)

## 5.3 Visual exploration  without summarizing info ####
# data loading: 
beta_data<-read.table("data/beta.pair_dataset.csv", header= T, sep = ",")

beta_data$predator<-factor(beta_data$predator, levels = c("0", "1"))
beta_data$nutrients<-factor(beta_data$nutrients, levels = c("0", "1"))

beta_data <- beta_data %>%
  mutate(treatment = factor(paste(
    ifelse(predator == "1", "Pred+", "Pred-"),
    ifelse(nutrients == "0", "Nutr-", "Nutr+"),
    sep = "/"))) %>%
  mutate(
    Date1 = as.Date(Date1), 
    Date2 = as.Date(Date2), 
    Month1 = format(as.Date(Date1), "%b"),
    Month1 = factor(Month1, levels = month.abb)
  ) %>% 
  mutate(component = factor(component, 
                            levels= c("Total (SOR)", "Turnover (SIM)", "Nestedness (SNE)")))

str(beta_data)

# Basic temporal trends by treatment and component
ggplot(beta_data, aes(x = TimeDiff_days, y = distance, color = treatment)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", se = TRUE) +
  facet_wrap(~ component, scales = "free_y") +
  labs(x = "Time difference (days)", y = "Beta diversity", color = "Treatment") +
  theme_minimal() +
  scale_color_brewer(palette = "Set1") +
  theme(legend.position = "bottom")

#This graph represent all the beta-pair comparisons, with all the delta comparison
#I don't think it is really information yet for our analyses. 

## 5.4 Stacked area chart for turnover vs. nestedness component ####

# First create a dataset with the correct time bins
time_bins <- seq(0, max(beta_data$TimeDiff_days), by = 30)
beta_data$time_bin <- cut(beta_data$TimeDiff_days, 
                          breaks = time_bins, 
                          labels = time_bins[-length(time_bins)] + 15)

# Get component-specific data
component_summary <- beta_data %>%
  filter(component %in% c("Turnover (SIM)", "Nestedness (SNE)")) %>%
  group_by(treatment, time_bin, component) %>%
  summarize(mean_distance = mean(distance, na.rm = TRUE), .groups = "drop") %>%
  # Convert to wide format for stacking
  pivot_wider(id_cols = c(treatment, time_bin),
              names_from = component,
              values_from = mean_distance) %>%
  # Rename for easier handling
  rename(mean_sim = `Turnover (SIM)`, mean_sne = `Nestedness (SNE)`) %>%
  # Long format for ggplot
  pivot_longer(cols = c(mean_sim, mean_sne),
               names_to = "component",
               values_to = "value") %>%
  # Make sure time_bin is numeric for plotting
  mutate(time_bin = as.numeric(as.character(time_bin)),
         # Create nicer labels for the components
         component = factor(component, 
                            levels = c("mean_sim", "mean_sne"),
                            labels = c("Turnover", "Nestedness")))


# Create the stacked area chart
ggplot(component_summary, aes(x = time_bin, y = value, fill = component)) +
  geom_area() +
  facet_wrap(~ treatment) +
  labs(x = "Time (days)", y = "Beta diversity component", fill = "Component") +
  scale_fill_manual(values = component_colors) +
  theme_minimal() +
  theme(legend.position = "bottom")


#This graphs also shows all the comparisons between all the pairs of the dataset
#I put it aside for now. 
rm(time_bins, component_summary)

## 5.5 Annotate delta t + create dataset ####
# To understand the patterns of the succession in Lakes, I would like to 
# analyse the consecutive time points rather than the whole pairwise
# comparison. 

beta_data <- read.table("data/beta.pair_dataset.csv", header = T, sep = ",")
beta_data$predator <- factor(beta_data$predator, levels = c("0", "1"))
beta_data$nutrients <- factor(beta_data$nutrients, levels = c("0", "1"))
beta_data$Date1 <- as.Date(beta_data$Date1)
beta_data$Date2 <- as.Date(beta_data$Date2)
beta_data$month1 <- as.numeric(format(beta_data$Date1, "%m"))

# Identify all sampling dates for each lake
unique_dates <- beta_data %>%
  select(Lake_ID, Date1, Date2) %>%
  pivot_longer(cols = c(Date1, Date2), names_to = "date_type", values_to = "date") %>%
  select(-date_type) %>%
  distinct() %>%
  arrange(Lake_ID, date)

# Create time points for each lake
lake_timepoints <- unique_dates %>%
  group_by(Lake_ID) %>%
  mutate(timepoint = row_number()) %>%
  select(Lake_ID, date, timepoint)

# Join timepoints back to the original data
beta_data_with_timepoints <- beta_data %>%
  # Join timepoint for Date1
  left_join(lake_timepoints, by = c("Lake_ID", "Date1" = "date")) %>%
  rename(timepoint1 = timepoint) %>%
  # Join timepoint for Date2
  left_join(lake_timepoints, by = c("Lake_ID", "Date2" = "date")) %>%
  rename(timepoint2 = timepoint)

# Calculate delta category based on difference in timepoints and months
beta_data_with_delta <- beta_data_with_timepoints %>%
  mutate(
    # Ensure timepoint2 > timepoint1
    timePointsApart = timepoint2 - timepoint1,
    # Calculate months between dates
    monthsApart = (year(Date2) - year(Date1)) * 12 + (month(Date2) - month(Date1)),
    # Only keep consecutive timepoints (where difference is 1)
    isConsecutive = timePointsApart == 1,
    # Create delta category based on actual months, for consecutive timepoints
    deltaCategory = ifelse(isConsecutive, paste0("delta_", monthsApart), NA)
  ) %>%
  # Filter to keep only consecutive timepoints
  filter(isConsecutive)


# Add manual correction for April23 sampling (which was done on May 3rd)
beta_data_with_delta <- beta_data_with_delta %>%
  mutate(
    # Correct monthsApart for April23 sampling date (03.05.2023)
    monthsApart = case_when(
      # If Date2 is May 3rd, 2023 and monthsApart is 2, change to 1
      Date2 == as.Date("2023-05-03") & monthsApart == 2 ~ 1,
      # Otherwise keep original value
      TRUE ~ monthsApart
    ),
    # Update deltaCategory to reflect the corrected monthsApart
    deltaCategory = paste0("delta_", monthsApart)
  )

#write.csv(beta_data_with_delta, "data/BetaPair_deltat_consecutive.csv", row.names = FALSE)

# Clean up environment
rm(beta_data, beta_data_with_timepoints, lake_timepoints, unique_dates)

################################################
# 6. BETA.PAIR dataset and raw visualization ####
################################################
## 6.1 Load datasets ####

beta <- read.csv("data/BetaPair_deltat_consecutive.csv", header = T, sep=",")

beta$Lake_ID<-as.factor(beta$Lake_ID)
beta$Date1<-as.Date(beta$Date1)
beta$Date2<- as.Date(beta$Date2)
beta$component<-factor(beta$component, 
                       levels = c("Total (SOR)", "Turnover (SIM)", "Nestedness (SNE)"),
                       labels = c("Tot", "Turn", "Nest"))

beta$predator<-factor(beta$predator, levels = c("0", "1"))
beta$nutrients<-factor(beta$nutrients, levels = c("0", "1"))
beta$deltaCategory<-as.factor(beta$deltaCategory)
beta$timepoint1<-as.factor(beta$timepoint1)

### 6.1.1 Delta 1 dataset preparation  ####

delta1 <- beta %>% 
  filter(monthsApart == 1) %>%
  mutate(days_since_start = as.numeric(Date1 - min(Date1)))

## 6.2 Visualization delta 1 data ####
### 6.2.1 Total diversity (sor) ####
tot_d1<-delta1 %>% 
  filter(component == "Tot") %>% 
  mutate(treatment = factor(paste(
    ifelse(predator == "1", "Pred+", "Pred-"),
    ifelse(nutrients == "0", "Nutr-", "Nutr+"),
    sep = "/")))

tot_d1$treatment <- factor(tot_d1$treatment, 
                           levels = c("Pred-/Nutr-", "Pred+/Nutr-", 
                                      "Pred-/Nutr+", "Pred+/Nutr+"))


# Create the plot with customized elements
raw_tot <- ggplot(tot_d1, aes(x = Date1, y = distance, color = treatment, 
                              fill = treatment)) + 
  geom_point(size = 2.5, alpha = 0.7) + 
  geom_smooth(method = loess, span = 0.7, alpha = 0.15, linewidth = 1.2) + 
  scale_color_manual(values = custom_colors) +  
  scale_fill_manual(values = custom_colors) +
  labs(title = "Total Beta Diversity",
       tag = "A",
       x = "",
       y = " Total Beta Diversity  \nvalue",
       color = "Treatment combination"
  ) +
  scale_x_date(
    date_breaks = "1 month",
    labels = function(x) {
      month_num <- as.numeric(format(x, "%m"))
      year <- format(x, "%y")  # Use %y for 2-digit year
      paste0(english_months[month_num], year)  # Format like "Jan'22"
    },
    expand = c(0.01, 0.01)
  ) +
  guides(fill = "none") +
  theme_jg_right +
  theme(panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "right")

# Display the plot
raw_tot


### 6.2.2 Turnover (sim)  ####

turn_d1<-delta1 %>% 
  filter(component == "Turn") %>% 
  mutate(treatment = factor(paste(
    ifelse(predator == "1", "Pred+", "Pred-"),
    ifelse(nutrients == "0", "Nutr-", "Nutr+"),
    sep = "/")))

turn_d1$treatment <- factor(turn_d1$treatment, 
                            levels = c("Pred-/Nutr-", "Pred+/Nutr-", 
                                       "Pred-/Nutr+", "Pred+/Nutr+"))


# Create the plot with customized elements
raw_turn <- ggplot(turn_d1, aes(x = Date1, y = distance, color = treatment, 
                                fill = treatment)) + 
  geom_point(size = 2.5, alpha = 0.7) + 
  geom_smooth(method = loess, span = 0.7, alpha = 0.15, linewidth = 1.2) + 
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  labs(title = "Turnover",
       tag = "B",
       x = "",
       y = "Turnover \nvalue"
  ) +
  scale_x_date(
    date_breaks = "1 month",
    labels = function(x) {
      month_num <- as.numeric(format(x, "%m"))
      year <- format(x, "%y")  # Use %y for 2-digit year
      paste0(english_months[month_num], year)  # Format like "Jan'22"
    },
    expand = c(0.01, 0.01)
  ) +
  theme_jg_right +
  theme(legend.position = "none",
        panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))

# Display the plot
raw_turn


### 6.2.3 Nestedness(sne) ####

nest_d1<-delta1 %>% 
  filter(component == "Nest") %>% 
  mutate(treatment = factor(paste(
    ifelse(predator == "1", "Pred+", "Pred-"),
    ifelse(nutrients == "0", "Nutr-", "Nutr+"),
    sep = "/")))

nest_d1$treatment <- factor(nest_d1$treatment, 
                            levels = c("Pred-/Nutr-", "Pred+/Nutr-", 
                                       "Pred-/Nutr+", "Pred+/Nutr+"))


# Create the plot with customized elements
raw_nest <- ggplot(nest_d1, aes(x = Date1, y = distance, color = treatment, 
                                fill = treatment)) + 
  geom_point(size = 2.5, alpha = 0.7) + 
  geom_smooth(method = loess, span = 0.7, alpha = 0.15, linewidth = 1.2) + 
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  labs(title = "Nestedness",
       tag = "C",
       x = "",
       y = "Nestedness \nvalue"
  ) +
  scale_x_date(
    date_breaks = "1 month",
    labels = function(x) {
      month_num <- as.numeric(format(x, "%m"))
      year <- format(x, "%y")  # Use %y for 2-digit year
      paste0(english_months[month_num], year)  # Format like "Jan'22"
    },
    expand = c(0.01, 0.01)
  ) +
  theme_jg_right + 
  theme(legend.position = "none",
        panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))

# Display the plot
raw_nest


### 6.2.4 FIGURE S6: Assembling raw data ####

# Then assemble your plot
FigS6 <- cowplot::plot_grid(
  raw_tot,
  cowplot::plot_grid(raw_turn, raw_nest, ncol = 2, align = "h"),
  ncol = 1,
  rel_heights = c(1, 1, 0.2)
)

# Check the result
print(FigS6)

# Save using ggsave
#ggsave("figure/S6.FigS6_Beta.pair_delta1.png", FigS6, width = 10, height = 12, dpi = 300)
#ggsave("figure/S6.FigS6_Beta.pair_delta1.pdf", FigS6, width = 10, height = 12, dpi = 300)

rm(raw_nest, raw_tot, raw_turn, FigS6)

### 6.2.5 Heat map of beta div total (sor)  NOT ADDED IN THE PAPER ####
# This provides a general vision of the dataset (tot_d1)
# heatmap dataset: 
heatmap_data_tot <- tot_d1 %>%
  select(Lake_ID, timepoint1, timepoint2, distance, treatment, days_since_start)

# Create heatmap
heat_tot <- ggplot(heatmap_data_tot, aes(x = factor(days_since_start), y = factor(treatment), fill = distance)) +
  geom_tile() +
  scale_fill_gradient2(low = "white", mid = "#2A9D8F", high = "#264653", 
                       midpoint = median(heatmap_data_tot$distance),
                       name = "") +
  labs(x = "Days since start", y = "",
       tag="A",
       title = "Total beta-Diversity") +
  theme_jg_right 

heat_tot

### 6.2.6 Heat map of turnover (sim) ####

# Join with original data
heatmap_data_turn <- turn_d1 %>%
  select(Lake_ID, timepoint1, timepoint2, distance, treatment, days_since_start) 

# Create heatmap
heat_turn <- ggplot(heatmap_data_turn, aes(x = factor(days_since_start), y = factor(treatment), fill = distance)) +
  geom_tile() +
  scale_fill_gradient2(low = "white", mid = "#F4A261", high = "#E76F51", 
                       midpoint = median(heatmap_data_turn$distance),
                       name = "") +
  labs(x = "Days since start", y = "",
       tag ="B",
       title = "Turnover") +
  theme_jg_right


heat_turn


### 6.2.7 Heat map of Nestedness (sne) ####

# Join with original data
heatmap_data_nest <- nest_d1 %>%
  select(Lake_ID, timepoint1, timepoint2, distance, treatment, days_since_start)

# Create heatmap
heat_nest <- ggplot(heatmap_data_nest, aes(x = factor(days_since_start), y = factor(treatment), fill = distance)) +
  geom_tile() +
  scale_fill_gradient2(low = "white", mid = "#C77DFF", high = "#7B2CBF", 
                       midpoint = median(heatmap_data_nest$distance),
                       name = "") +
  labs(x = "Days since start", y = "",
       tag ="C",
       title = "Nestedness") +
  theme_jg_right

heat_nest

### 6.2.8 NOT INCLUDED: combined heatmaps ####
FigS3 <- heat_tot / heat_turn / heat_nest

# Display the combined plot
FigS3

#ggsave("figure/S3.FigS3_Beta.pair_delta1_heatmaps.png", FigS3, width = 10, height = 12, dpi = 300)
#ggsave("figure/S3.FigS3_Beta.pair_delta1_heatmaps.pdf", FigS3, width = 10, height = 12, dpi = 300)

rm(beta, FigS2, heat_nest, heat_tot, heat_turn, heatmap_data_nest, 
   heatmap_data_tot, heatmap_data_turn)

################################################
# 6b. Environmental variable dataset and raw visualization ####
################################################
## 6b.1 load dataset ####
env <- read.csv("data/samples_with_environmental_data.csv", header = T, sep = ",")

env <- env %>% filter(filter_pore_size == "0.45",
                      startsWith(Samples_ID, "L")) %>%
  distinct(Lake_ID, Collection_month, .keep_all = TRUE) %>%
  mutate(date_env_data = case_when(
    # December: -01 in wrong date → correct day was 4th, 5th, or 6th January
    Collection_month == "December" & date_env_data == "2023-04-01" ~ "2023-01-04",
    Collection_month == "December" & date_env_data == "2023-05-01" ~ "2023-01-05",
    Collection_month == "December" & date_env_data == "2023-06-01" ~ "2023-01-06",
    # March: all lakes have same wrong date
    Collection_month == "March"  ~ "2023-03-06",
    # April23: two groups
    Collection_month == "April23" & date_env_data == "2023-02-05" ~ "2023-05-02",
    Collection_month == "April23" & date_env_data == "2023-03-05" ~ "2023-05-03",
    TRUE ~ date_env_data
  ))

env$Collection_date <- as.Date(env$Collection_date, "%d.%m.%Y")
env$date_env_data   <- as.Date(env$date_env_data)

str(env)

## 6b.2 FIGURE S2: Visualization of env distribution ####

env_var_raw <- c("temperature", "do_percent", "do_mg_l", "ph",
                 "conductivity", "turbidity", "chlorophyll_rfu", "chlorophyll_ug_l")

env_summary <- env %>%
  select(all_of(env_var_raw)) %>%
  summary()

print(env_summary)

#par(mfrow = c(3, 3))
#for (v in env_var_raw) {
#  hist(env[[v]], main = v, xlab = v, col = "steelblue", breaks = 20)
#}
#par(mfrow = c(1, 1))

# we can see that there is a outlier in the conductivity!

FigS2 <- env %>%
  select(Lake_ID, Collection_month, all_of(env_var_raw)) %>%
  pivot_longer(cols = all_of(env_var_raw),
               names_to  = "variable",
               values_to = "value") %>%
  mutate(variable = factor(variable,
                           levels = env_var_raw,
                           labels = c("Temperature", "DO (%)", "DO (mg/L)",
                                      "pH", "Conductivity", "Turbidity",
                                      "Chlorophyll (RFU)", "Chlorophyll (µg/L)"))) %>%
  ggplot(aes(x = value)) +
  geom_histogram(fill = "#2166ac", colour = "white", bins = 20, alpha = 0.85) +
  facet_wrap(~ variable, scales = "free", ncol = 3) +
  labs(x = "Value", y = "Count",
       title = "Raw distributions of environmental variables") +
  theme_bw() +
  theme(strip.text = element_text(size = 11),
        panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5))

FigS2

#ggsave("figure/S2.FigS2_raw_distributions.png", FigS2, width = 10, height = 9, dpi = 300)
#ggsave("figure/S2.FigS2_raw_distributions.pdf", FigS2, width = 10, height = 9, dpi = 300)

### 6b.2.1 Flagging conductivity outlier ####

conductivity_mean <- mean(env$conductivity, na.rm = TRUE)
conductivity_sd   <- sd(env$conductivity,   na.rm = TRUE)

env <- env %>%
  mutate(conductivity = if_else(
    conductivity > conductivity_mean + 3 * conductivity_sd,
    NA_real_,
    conductivity
  ))

### 6b.2.2 Flag June 2022 as truly missing ####
# June 2022 measurements failed entirely — set to NA so they are interpolated
env <- env %>%
  mutate(across(all_of(env_var_raw),
                ~ if_else(Collection_month == "June", NA_real_, .)))

## 6b.3 Interpolate to exact eDNA sampling date ####
# For each variable, use all available sonde measurements (anchored at their
# actual measurement date) to interpolate the value at the exact eDNA
# collection date. This uses all available information regardless of the
# gap between sonde and eDNA dates.

env_interpolation <- env %>%
  group_by(Lake_ID) %>%
  mutate(
    # Convert dates to numeric for interpolation
    date_env_numeric        = as.numeric(date_env_data),
    date_collection_numeric = as.numeric(Collection_date),
    
    # Interpolate each variable at the exact eDNA collection date
    # x = actual sonde measurement dates (NAs excluded automatically)
    # xout = eDNA collection dates we want values for
    # rule = 2: carry boundary values forward/backward for edge months
    temperature_interp = approx(
      x = date_env_numeric[!is.na(temperature)],
      y = temperature[!is.na(temperature)],
      xout = date_collection_numeric, rule = 2)$y,
    
    do_percent_interp = approx(
      x = date_env_numeric[!is.na(do_percent)],
      y = do_percent[!is.na(do_percent)],
      xout = date_collection_numeric, rule = 2)$y,
    
    do_mg_l_interp = approx(
      x = date_env_numeric[!is.na(do_mg_l)],
      y = do_mg_l[!is.na(do_mg_l)],
      xout = date_collection_numeric, rule = 2)$y,
    
    ph_interp = approx(
      x = date_env_numeric[!is.na(ph)],
      y = ph[!is.na(ph)],
      xout = date_collection_numeric, rule = 2)$y,
    
    conductivity_interp = approx(
      x = date_env_numeric[!is.na(conductivity)],
      y = conductivity[!is.na(conductivity)],
      xout = date_collection_numeric, rule = 2)$y,
    
    turbidity_interp = approx(
      x = date_env_numeric[!is.na(turbidity)],
      y = turbidity[!is.na(turbidity)],
      xout = date_collection_numeric, rule = 2)$y,
    
    chlorophyll_rfu_interp = approx(
      x = date_env_numeric[!is.na(chlorophyll_rfu)],
      y = chlorophyll_rfu[!is.na(chlorophyll_rfu)],
      xout = date_collection_numeric, rule = 2)$y,
    
    chlorophyll_ug_l_interp = approx(
      x = date_env_numeric[!is.na(chlorophyll_ug_l)],
      y = chlorophyll_ug_l[!is.na(chlorophyll_ug_l)],
      xout = date_collection_numeric, rule = 2)$y
  ) %>%
  ungroup()


### 6b.3.1 FIGURE S3: Visualization raw vs interpolated values ####

env_long <- env_interpolation %>%
  select(Lake_ID, Collection_month,
         temperature, do_percent, do_mg_l, ph,
         conductivity, turbidity, chlorophyll_rfu, chlorophyll_ug_l,
         temperature_interp, do_percent_interp, do_mg_l_interp, ph_interp,
         conductivity_interp, turbidity_interp, chlorophyll_rfu_interp, chlorophyll_ug_l_interp) %>%
  pivot_longer(
    cols      = c(temperature:chlorophyll_ug_l),
    names_to  = "variable",
    values_to = "raw"
  ) %>%
  pivot_longer(
    cols      = c(temperature_interp:chlorophyll_ug_l_interp),
    names_to  = "variable_interp",
    values_to = "interpolated"
  ) %>%
  filter(paste0(variable, "_interp") == variable_interp) %>%
  select(-variable_interp)

FigS3 <- env_long %>%
  pivot_longer(cols = c(raw, interpolated),
               names_to  = "type",
               values_to = "value") %>%
  mutate(type = factor(type, levels = c("raw", "interpolated"))) %>%
  ggplot(aes(x = value, fill = type, colour = type)) +
  geom_density(alpha = 0.4) +
  facet_wrap(~ variable, scales = "free") +
  scale_fill_manual(values   = c("raw" = "#2166ac", "interpolated" = "#d6604d"),
                    labels   = c("raw" = "Raw (NA flagged)", "interpolated" = "Interpolated")) +
  scale_colour_manual(values = c("raw" = "#2166ac", "interpolated" = "#d6604d"),
                      labels = c("raw" = "Raw (NA flagged)", "interpolated" = "Interpolated")) +
  labs(x = "Value", y = "Density", fill = "Type", colour = "Type",
       title = "Distribution: raw vs interpolated environmental variables") +
  theme_bw() +
  theme(legend.position = "bottom")

FigS3

#ggsave("figure/S3.FigS3_env_interpolation_vs_raw.png", FigS3, width = 10, height = 12, dpi = 300)
#ggsave("figure/S3.FigS3_env_interpolation_vs_raw.pdf", FigS3, width = 10, height = 12, dpi = 300)


### 6b.3.2 FIGURE S4: Correlation check among variables ####

cor_check <- env_interpolation %>%
  select(temperature_interp, do_percent_interp, do_mg_l_interp,
         ph_interp, conductivity_interp, turbidity_interp,
         chlorophyll_rfu_interp, chlorophyll_ug_l_interp) %>%
  drop_na()

corrplot(cor(cor_check, method = "pearson"),
         method = "color", type = "upper", addCoef.col = "black",
         number.cex = 0.7, tl.cex = 0.8, tl.col = "black",
         col = COL2("RdBu", 200), diag = FALSE,
         title = "Correlation among raw interpolated env variables", mar = c(0, 0, 2, 0))


# Variable selection: T, DO% (less correlated than mg/L), pH, conductivity,
# turbidity, chlorophyll ug/L (conceptually clearer).
# do_mg_l and chlorophyll_rfu dropped.

png("figure/S4.FigS4_env_correlation.png", width = 9, height = 9, 
    units = "in", res = 300)
corrplot(cor(cor_check, method = "pearson"),
         method = "color", type = "upper", addCoef.col = "black",
         number.cex = 0.7, tl.cex = 0.8, tl.col = "black",
         col = COL2("RdBu", 200), diag = FALSE,
         title = "Correlation among raw interpolated env variables", mar = c(0, 0, 2, 0))
dev.off()

pdf("figure/S4.FigS4_env_correlation.pdf", width = 9, height = 9)
corrplot(cor(cor_check, method = "pearson"),
         method = "color", type = "upper", addCoef.col = "black",
         number.cex = 0.7, tl.cex = 0.8, tl.col = "black",
         col = COL2("RdBu", 200), diag = FALSE,
         title = "Correlation among raw interpolated env variables", mar = c(0, 0, 2, 0))
dev.off()

## 6b.4 Computation of delta — log-ratio ####
### 6b.4.1 Build table for delta computation ####

env_raw <- env_interpolation %>%
  mutate(time_index = c("April22"=1,"May"=2,"June"=3,"July"=4,"August"=5,
                        "September"=6,"October"=7,"November"=8,"December"=9,
                        "January"=10,"February"=11,"March"=12,"April23"=13
  )[Collection_month]) %>%
  select(Lake_ID, time_index,
         temperature_interp, do_percent_interp, ph_interp,
         conductivity_interp, turbidity_interp, chlorophyll_ug_l_interp)

env_t1 <- env_raw %>%
  rename_with(~ paste0(.x, "_t1"), -c(Lake_ID, time_index)) %>%
  rename(timepoint1 = time_index)

env_t2 <- env_raw %>%
  rename_with(~ paste0(.x, "_t2"), -c(Lake_ID, time_index)) %>%
  rename(timepoint2 = time_index)

env_paired <- env_t1 %>%
  inner_join(env_t2, by = "Lake_ID") %>%
  filter(timepoint2 == timepoint1 + 1)

### 6b.4.2 Log-ratio: delta = log(raw_t2 / raw_t1) ####
# No z-scoring: the log scale already standardises the values

env_delta_log <- env_paired %>%
  mutate(
    delta_temperature  = log(temperature_interp_t2   / temperature_interp_t1),
    delta_do_percent   = log(do_percent_interp_t2     / do_percent_interp_t1),
    delta_ph           = log(ph_interp_t2             / ph_interp_t1),
    delta_conductivity = log(conductivity_interp_t2   / conductivity_interp_t1),
    delta_turbidity    = log(turbidity_interp_t2      / turbidity_interp_t1),
    delta_chlorophyll  = log(chlorophyll_ug_l_interp_t2 / chlorophyll_ug_l_interp_t1)
  ) %>%
  select(Lake_ID, timepoint1, timepoint2, starts_with("delta_"))

# Check for NA/Inf introduced by log of non-positive values
env_delta_log %>% summarise(across(starts_with("delta_"), ~ sum(is.na(.) | is.infinite(.))))

### 6b.4.3 Join env delta with beta ####

beta_env_log <- delta1 %>%
  mutate(timepoint1 = as.integer(as.character(timepoint1)),
         timepoint2 = as.integer(as.character(timepoint2))) %>%
  inner_join(env_delta_log, by = c("Lake_ID", "timepoint1", "timepoint2")) %>%
  mutate(
    days_since_start = as.numeric(Date1 - min(Date1)),
    midpoint_date    = Date1 + (Date2 - Date1) / 2,
    midmonth_numeric = as.numeric(format(midpoint_date, "%m")),
    sin_year         = sin(2 * pi * midmonth_numeric / 12),
    cos_year         = cos(2 * pi * midmonth_numeric / 12),
    sin_year2        = sin(4 * pi * midmonth_numeric / 12),
    cos_year2        = cos(4 * pi * midmonth_numeric / 12)
  )

rm(env, env_interpolation, env_long, env_raw, env_t1, env_t2, env_paired, env_delta_log, 
   FigS3, FigS2)

### 6b.4.4 Create component subsets ####
# tot_d1, turn_d1, nest_d1 are created here from beta_env_log so that
# the environmental delta columns are available for the residual analyses
# in sections 7, 8 and 9.

tot_d1 <- beta_env_log %>%
  filter(component == "Tot") %>%
  mutate(timepoint1 = as.factor(timepoint1),
         timepoint2 = as.factor(timepoint2),
         treatment  = factor(paste(
           ifelse(predator == "1", "Pred+", "Pred-"),
           ifelse(nutrients == "0", "Nutr-", "Nutr+"),
           sep = "/"),
           levels = c("Pred-/Nutr-", "Pred+/Nutr-", "Pred-/Nutr+", "Pred+/Nutr+")))

turn_d1 <- beta_env_log %>%
  filter(component == "Turn") %>%
  mutate(timepoint1 = as.factor(timepoint1),
         timepoint2 = as.factor(timepoint2),
         treatment  = factor(paste(
           ifelse(predator == "1", "Pred+", "Pred-"),
           ifelse(nutrients == "0", "Nutr-", "Nutr+"),
           sep = "/"),
           levels = c("Pred-/Nutr-", "Pred+/Nutr-", "Pred-/Nutr+", "Pred+/Nutr+")))

nest_d1 <- beta_env_log %>%
  filter(component == "Nest") %>%
  mutate(timepoint1 = as.factor(timepoint1),
         timepoint2 = as.factor(timepoint2),
         treatment  = factor(paste(
           ifelse(predator == "1", "Pred+", "Pred-"),
           ifelse(nutrients == "0", "Nutr-", "Nutr+"),
           sep = "/"),
           levels = c("Pred-/Nutr-", "Pred+/Nutr-", "Pred-/Nutr+", "Pred+/Nutr+")))

### 6b 4.5 effect of treatments on the delta env ####
# Temperature 
t1 <- glmmTMB(delta_temperature ~ predator * nutrients + 
                (1|Lake_ID), 
              data = tot_d1)
simulateResiduals(t1, plot = T)
summary(t1)
Anova(t1)

# DO 
d1 <- glmmTMB(delta_do_percent ~ predator * nutrients + 
                (1|Lake_ID), 
              data = tot_d1)
simulateResiduals(d1, plot = T)
summary(d1)
Anova(d1)

# pH 
p1 <- glmmTMB(delta_ph ~ predator * nutrients + 
                (1|Lake_ID), 
              data = tot_d1)
simulateResiduals(p1, plot = T)
summary(p1)
Anova(p1)

# Conductivity
c1 <- glmmTMB(delta_conductivity ~ predator * nutrients + 
                (1|Lake_ID), 
              data = tot_d1)
simulateResiduals(c1, plot = T)
summary(c1)
Anova(c1)


# Turbidity 
tu1 <- glmmTMB(delta_turbidity ~ predator * nutrients + 
                (1|Lake_ID), 
              data = tot_d1)
simulateResiduals(tu1, plot = T)
summary(tu1)
Anova(tu1)

# Chlorophyll
cl1 <- glmmTMB(delta_chlorophyll ~ predator * nutrients + 
                (1|Lake_ID), 
              data = tot_d1)
simulateResiduals(cl1, plot = T)
summary(cl1)
Anova(cl1)

#Since none
#### 6b 4.5.1 TABLE S2: glmm summary ####
t1_tidy  <- tidy(t1,  conf.int = TRUE)
t1_tidy$model  <- "Temperature"

d1_tidy  <- tidy(d1,  conf.int = TRUE)
d1_tidy$model  <- "Dissolved oxygen (%)"

p1_tidy  <- tidy(p1,  conf.int = TRUE)
p1_tidy$model  <- "pH"

c1_tidy  <- tidy(c1,  conf.int = TRUE)
c1_tidy$model  <- "Conductivity"

tu1_tidy <- tidy(tu1, conf.int = TRUE)
tu1_tidy$model <- "Turbidity"

cl1_tidy <- tidy(cl1, conf.int = TRUE)
cl1_tidy$model <- "Chlorophyll (µg/L)"

# Combine into a single table
model_results_env <- bind_rows(t1_tidy, d1_tidy, p1_tidy, c1_tidy, tu1_tidy, cl1_tidy)

rm(t1_tidy, d1_tidy, p1_tidy, c1_tidy, tu1_tidy, cl1_tidy)

# Rename columns for clarity and select relevant ones
model_results_env <- model_results_env %>%
  select(model, term, estimate, std.error, statistic, p.value, conf.low, conf.high) %>%
  rename(
    "Component"   = model,
    "Term"        = term,
    "Estimate"    = estimate,
    "Std. Error"  = std.error,
    "z value"     = statistic,
    "p-value"     = p.value,
    "Lower CI"    = conf.low,
    "Upper CI"    = conf.high
  )

# Add significance indicators
model_results_env$Significance <- ifelse(model_results_env$`p-value` < 0.001, '***',
                                         ifelse(model_results_env$`p-value` < 0.01,  '**',
                                                ifelse(model_results_env$`p-value` < 0.05,  '*',
                                                       ifelse(model_results_env$`p-value` < 0.1,   '.', ''))))

# Combine p-value with significance
model_results_env$FormattedPvalue <- paste0(
  sprintf("%.3f", model_results_env$`p-value`),
  " ",
  model_results_env$Significance
)

# Select and rename columns
model_results_env_formatted <- model_results_env %>%
  select(-Significance, -`p-value`) %>%
  rename(`p-value` = FormattedPvalue)

# Create flextable
ft_models_env <- flextable(model_results_env_formatted) %>%
  colformat_double(digits = 3) %>%
  bold(i = model_results_env$`p-value` < 0.05, j = "p-value") %>%
  theme_vanilla() %>%
  add_header_lines("GLMM Summary: Effect of Treatments on Monthly Changes in Environmental Variables (δenv)") %>%
  autofit()

# Save to Word
doc <- read_docx()
doc <- body_add_flextable(doc, ft_models_env)
#print(doc, target = "tables/TableS2.delta_env_glmm_results.docx")

rm(ft_models_env, model_results_env, model_results_env_formatted)

#### 6b 4.5.2 TABLES S3: anova results ####
t1_anova  <- Anova(t1)
d1_anova  <- Anova(d1)
p1_anova  <- Anova(p1)
c1_anova  <- Anova(c1)
tu1_anova <- Anova(tu1)
cl1_anova <- Anova(cl1)

# Convert to data frames
t1_anova_df        <- as.data.frame(t1_anova)
t1_anova_df$Effect <- rownames(t1_anova_df)
rownames(t1_anova_df) <- NULL
t1_anova_df$Component <- "Temperature"

d1_anova_df        <- as.data.frame(d1_anova)
d1_anova_df$Effect <- rownames(d1_anova_df)
rownames(d1_anova_df) <- NULL
d1_anova_df$Component <- "Dissolved oxygen (%)"

p1_anova_df        <- as.data.frame(p1_anova)
p1_anova_df$Effect <- rownames(p1_anova_df)
rownames(p1_anova_df) <- NULL
p1_anova_df$Component <- "pH"

c1_anova_df        <- as.data.frame(c1_anova)
c1_anova_df$Effect <- rownames(c1_anova_df)
rownames(c1_anova_df) <- NULL
c1_anova_df$Component <- "Conductivity"

tu1_anova_df        <- as.data.frame(tu1_anova)
tu1_anova_df$Effect <- rownames(tu1_anova_df)
rownames(tu1_anova_df) <- NULL
tu1_anova_df$Component <- "Turbidity"

cl1_anova_df        <- as.data.frame(cl1_anova)
cl1_anova_df$Effect <- rownames(cl1_anova_df)
rownames(cl1_anova_df) <- NULL
cl1_anova_df$Component <- "Chlorophyll (µg/L)"

# Combine
anova_results_env <- bind_rows(t1_anova_df, d1_anova_df, p1_anova_df,
                               c1_anova_df, tu1_anova_df, cl1_anova_df)

rm(t1_anova, d1_anova, p1_anova, c1_anova, tu1_anova, cl1_anova,
   t1_anova_df, d1_anova_df, p1_anova_df, c1_anova_df, tu1_anova_df, cl1_anova_df)

# Rename columns
names(anova_results_env)[names(anova_results_env) == "Chisq"]       <- "Chi-Square"
names(anova_results_env)[names(anova_results_env) == "Pr(>Chisq)"]  <- "p-value"

# Reorder columns
anova_results_env <- anova_results_env[, c("Component", "Effect", "Chi-Square", "Df", "p-value")]

# Add significance indicators
anova_results_env$Significance <- ifelse(anova_results_env$`p-value` < 0.001, '***',
                                         ifelse(anova_results_env$`p-value` < 0.01,  '**',
                                                ifelse(anova_results_env$`p-value` < 0.05,  '*',
                                                       ifelse(anova_results_env$`p-value` < 0.1,   '.', ''))))

# Combine p-value with significance
anova_results_env$FormattedPvalue <- paste0(
  sprintf("%.3f", anova_results_env$`p-value`),
  " ",
  anova_results_env$Significance
)

# Select and rename columns
anova_results_env_formatted <- anova_results_env %>%
  select(-Significance, -`p-value`) %>%
  rename(`p-value` = FormattedPvalue)

# Create flextable
ft_anova_env <- flextable(anova_results_env_formatted) %>%
  colformat_double(digits = 3) %>%
  bold(i = anova_results_env$`p-value` < 0.05, j = "p-value") %>%
  theme_vanilla() %>%
  add_header_lines("ANOVA Results: Effect of Treatments on Monthly Changes in Environmental Variables (δenv)") %>%
  autofit()

# Save to Word
doc <- read_docx()
doc <- body_add_flextable(doc, ft_anova_env)
#print(doc, target = "tables/TableS3.delta_env_anova_results.docx")

rm(doc, ft_anova_env, anova_results_env, anova_results_env_formatted)

## 6b.5 Visualization of environmental delta variables ####
# Mirrors the beta diversity visualizations in section 6.2:
# temporal trend colored by treatment + heatmap, faceted by variable.

### 6b.5.1 Prepare long format dataset ####

# Use tot_d1 since env deltas are identical across all three components
env_delta_long <- tot_d1 %>%
  select(Lake_ID, Date1, days_since_start, treatment,
         delta_temperature, delta_do_percent, delta_ph,
         delta_conductivity, delta_turbidity, delta_chlorophyll) %>%
  pivot_longer(
    cols      = starts_with("delta_"),
    names_to  = "variable",
    values_to = "value"
  ) %>%
  mutate(variable = factor(variable,
                           levels = c("delta_temperature", "delta_do_percent", "delta_ph",
                                      "delta_conductivity", "delta_turbidity", "delta_chlorophyll"),
                           labels = c("Temperature", "DO%", "pH",
                                      "Conductivity", "Turbidity", "Chlorophyll")))

### 6b.6.2 FIGURE 7 -Temporal trend plot — faceted by variable ####

FigS7 <- ggplot(env_delta_long, aes(x = Date1, y = value,
                                            color = treatment, fill = treatment)) +
  geom_point(size = 2, alpha = 0.6) +
  geom_smooth(method = loess, span = 0.7, alpha = 0.15, linewidth = 1.1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  facet_wrap(~ variable, scales = "free_y", ncol = 2) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  scale_x_date(
    date_breaks = "1 month",
    labels = function(x) {
      month_num <- as.numeric(format(x, "%m"))
      year <- format(x, "%y")
      paste0(english_months[month_num], year)
    },
    expand = c(0.01, 0.01)
  ) +
  guides(fill = "none") +
  labs(x = "",
       y = "Log-ratio delta (log t2/t1)",
       color = "Treatment combination",
       tag = "") +
  theme_jg_right +
  theme(panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 11),
        legend.position = "right")

FigS7

#ggsave("figure/S7.FigS6_deltaEnv_throughtime.png", FigS7, width = 10, height = 12, dpi = 300)
#ggsave("figure/S7.FigS6_deltaEnv_throughtime.pdf", FigS7, width = 10, height = 12, dpi = 300)

### 6b.6.3 FIGURE S7: Heatmap — faceted by variable ####

FigS7 <- ggplot(env_delta_long,
                            aes(x = factor(Date1),
                                y = factor(treatment),
                                fill = value)) +
  geom_tile(color = "white", linewidth = 0.3) +
  facet_wrap(~ variable, scales = "free", ncol = 2) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#d6604d",
                       midpoint = 0,
                       name = "Log-ratio\ndelta") +
  labs(x = "Days since start",
       y = "",
       tag = "B") +
  theme_jg_right +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 11),
        legend.position = "right")

FigS7


#ggsave("figure/S7_FigS7_heatmaps_env_throughtime.png", FigS7, width = 10, height = 14, dpi = 300)
#ggsave("figure/S7_FigS7_heatmaps_env_throughtime.pdf", FigS7, width = 10, height = 14, dpi = 300)

rm(cor_check, env_delta_long, FigS7)


###########################################################
# 7. Beta.pair Analyses glmmTMB  TOTAL DIVERSITY (sor) ####
###########################################################

range(tot_d1$distance)

## 7.1 models through the year ####
### 7.1.1 Analysis treatment only ####
# I selected the first model to analyse the diversity through the year, since it 
# allows for reducing the d.f. and lower AIC, while accounting for the multiple
# sampling of the same lake. 

tot_m1 <- glmmTMB(distance ~ predator * nutrients * as.numeric(days_since_start) 
                  + (1|Lake_ID), 
                  family = beta_family(link= "logit"), 
                  data = tot_d1)

simulateResiduals(tot_m1, plot = T)
summary(tot_m1)
Anova(tot_m1)

tot_m2 <- glmmTMB(distance ~ predator * nutrients * as.numeric(days_since_start), 
                  family = beta_family(link= "logit"), 
                  data = tot_d1)

simulateResiduals(tot_m2, plot = T)
anova(tot_m1, tot_m2)
# interaction not significant however I will still include Lake_ID to take into 
# account that we are smapling the sample lake along the year. 

summary(tot_m2)# AIC lower 
Anova(tot_m2)

rm(tot_m2)

### 7.1.2 Taking into account seasonality using sine cosine ####
# according to Stolwijk et al. 1999 we can take account for seasonality as confunding 
# factor in  multivariate analyses. 
# this is kind of circular analyses which was proposed to be genuinely more used 
# in ecological studies as proposed by Chen et al. 2015

# Modified model with seasonal adjustment - replace your tot_m1 with this
tot_m1_seasonal <- glmmTMB(distance ~ predator * nutrients * as.numeric(days_since_start) + 
                             (1| Lake_ID) + 
                             sin_year + cos_year, 
                           family = beta_family(link= "logit"), 
                           data = tot_d1)

simulateResiduals(tot_m1_seasonal, plot = T)
summary(tot_m1_seasonal)

# Compare models to see if seasonality is important
anova(tot_m1, tot_m1_seasonal)

# season seems not having any impact on the total beta_diversity, 
#no need to include it. 

rm(tot_m1_seasonal)

### 7.2.2 Residual analysis environmental variable ####

tot_d1$res_tot <- residuals(tot_m1)

res_tot <- glmmTMB(res_tot ~ delta_temperature + delta_do_percent + delta_ph +
                delta_conductivity + delta_turbidity + delta_chlorophyll + (1|Lake_ID),
              data = tot_d1)
simulateResiduals(res_tot, plot = T)
summary(res_tot)
Anova(res_tot)

## 7.2 models per 3 months ####
### phase 1: April2022 - June2022 ####
phase1_tot <- tot_d1 %>% 
  filter(Date1 %in% as.Date(c("2022-04-21", "2022-05-19", "2022-06-23")))


p1_tot_m1<- glmmTMB(distance ~ predator * nutrients * as.factor(days_since_start)+
                      (1|Lake_ID), 
                    family = beta_family(link = "logit"), 
                    data = phase1_tot)

simulateResiduals(p1_tot_m1, plot = T)
summary(p1_tot_m1)
Anova(p1_tot_m1)

#### Post-hoc phase 1 ####

p1_tot_emm <- emmeans(p1_tot_m1, ~ predator * nutrients * days_since_start)

# Compare phases within each predator x nutrient combination
pairs_p1_tot <- p1_tot_emm %>% 
  contrast(interaction = c("pairwise"), by = c("predator", "nutrients"), 
           adjust = "tukey")
summary(pairs_p1_tot)


### phase 2: July2022 - September2022 ####
phase2_tot <- tot_d1 %>% 
  filter(Date1 %in% as.Date(c("2022-07-21", "2022-08-18", "2022-09-22")))


p2_tot_m1<- glmmTMB(distance ~ predator * nutrients * as.factor(days_since_start) +
                      (1|Lake_ID), 
                    family = beta_family(link = "logit"), 
                    data = phase2_tot)

simulateResiduals(p2_tot_m1, plot = T)
summary(p2_tot_m1)
Anova(p2_tot_m1)

#### Post-hoc phase 2 ####

p2_tot_emm <- emmeans(p2_tot_m1, ~ predator * nutrients * days_since_start)

# Compare phases within each predator x nutrient combination
pairs_p2_tot <- p2_tot_emm %>% 
  contrast(interaction = c("pairwise"), by = c("predator", "nutrients"), 
           adjust = "tukey")
summary(pairs_p2_tot)

### phase 3: October2022 - December2022 ####
phase3_tot <- tot_d1 %>% 
  filter(Date1 %in% as.Date(c("2022-10-21", "2022-11-18", "2022-12-15")))


p3_tot_m1<- glmmTMB(distance ~ predator * nutrients * as.factor(days_since_start)+
                      (1|Lake_ID), 
                    family = beta_family(link = "logit"), 
                    data = phase3_tot)

simulateResiduals(p3_tot_m1, plot = T)
summary(p3_tot_m1)
Anova(p3_tot_m1)

#### Post-hoc phase 3 ####

p3_tot_emm <- emmeans(p3_tot_m1, ~ predator * nutrients * days_since_start)

# Compare phases within each predator x nutrient combination
pairs_p3_tot <- p3_tot_emm %>% 
  contrast(interaction = c("pairwise"), by = c("predator", "nutrients"), 
           adjust = "tukey")
summary(pairs_p3_tot)

### phase 4: January2023 - March2023 ####
phase4_tot <- tot_d1 %>% 
  filter(Date1 %in% as.Date(c("2023-01-25", "2023-02-21", "2023-03-28")))


p4_tot_m1<- glmmTMB(distance ~ predator * nutrients * as.factor(days_since_start) +
                      (1|Lake_ID), 
                    family = beta_family(link = "logit"), 
                    data = phase4_tot)

simulateResiduals(p4_tot_m1, plot = T)
summary(p4_tot_m1)
Anova(p4_tot_m1)

#### Post-hoc phase 4 ####

p4_tot_emm <- emmeans(p4_tot_m1, ~ predator * nutrients * days_since_start)

# Compare phases within each predator x nutrient combination
pairs_p4_tot <- p4_tot_emm %>% 
  contrast(interaction = c("pairwise"), by = c("predator", "nutrients"), 
           adjust = "tukey")
summary(pairs_p4_tot)


## 7.3 testing whether the treatment effect change through phases ####
# First, combine all your phase data into one dataset
all_phases_data_tot <- bind_rows(
  mutate(phase1_tot, Phase = "Phase1"),
  mutate(phase2_tot, Phase = "Phase2"),
  mutate(phase3_tot, Phase = "Phase3"),
  mutate(phase4_tot, Phase = "Phase4")
)

# Make Phase a factor with the correct order
all_phases_data_tot$Phase <- factor(all_phases_data_tot$Phase, 
                                    levels = c("Phase1", "Phase2", "Phase3", "Phase4"))

# Create unified model with Phase interactions
unified_model_tot <- glmmTMB(distance ~ predator * nutrients * Phase + (1|Lake_ID), 
                             family = beta_family(link = "logit"), 
                             data = all_phases_data_tot)

# Check model
simulateResiduals(unified_model_tot, plot = T)
summary(unified_model_tot)
Anova(unified_model_tot)

# phases is in interaction with treatment which means that the effect of treatment 
# on beta diversity indeed depends on the "phase" we are looking. 
### 7.3.1 Post-hoc comparisons ####
# Get predicted values from the model
emm_tot <- emmeans(unified_model_tot, ~ predator * nutrients * Phase)
emm_tot_df <- as.data.frame(emm_tot)

# Compare phases within each predator x nutrient combination
pairs_phase_tot <- emm_tot %>% 
  contrast(interaction = c("pairwise"), by = c("predator", "nutrients"), 
           adjust = "tukey")
summary(pairs_phase_tot)


# Compare predator effects within each phase x nutrient combination
pairs_pred_tot <- emm_tot %>% 
  contrast(interaction = c("pairwise"), by = c("Phase", "nutrients"), 
           adjust = "tukey")
summary(pairs_pred_tot)

# Compare nutrient effects within each phase x predator combination
pairs_nutr <- emm_tot %>% 
  contrast(interaction = c("pairwise"), by = c("Phase", "predator"), 
           adjust = "tukey")
summary(pairs_nutr)

# Compare nutrient effects within each phase x predator combination
pairs_nutr_tot <- emm_tot %>% 
  contrast(interaction = c("pairwise"), by = c("Phase", "predator"), 
           adjust = "tukey")
summary(pairs_nutr_tot)

## 7.4 Rate of change ####
# How quickly is community composition changing over time in different treatments?
# This addresses whether predators or nutrients affect the speed of community 
# turnover, not just the magnitude.

### 7.4.1 Compute rate of change ####
# Calculate rate of change as : beta-div / time diff 
tot_d1 <- tot_d1 %>%
  # Ensure data is ordered correctly
  arrange(Lake_ID, as.factor(days_since_start)) %>%
  # Calculate rate of change as distance/time
  mutate(rate_change = distance / TimeDiff_days)

# Check for extreme values or zeros in TimeDiff_days that might cause problems
summary(tot_d1$TimeDiff_days)
summary(tot_d1$rate_change)

### 7.4.2 Models rate of change overall ####

rate_tot_m1 <- glmmTMB(rate_change ~ predator * nutrients + 
                         (1|Lake_ID), 
                       data = tot_d1,
                       family = beta_family(link = "logit"))

simulateResiduals(rate_tot_m1, plot = T)
summary(rate_tot_m1)

rate_tot_m2 <- glmmTMB(rate_change ~ predator * nutrients + sin_year + cos_year +
                         (1|Lake_ID), 
                       data = tot_d1,
                       family = beta_family(link = "logit"))

simulateResiduals(rate_tot_m2, plot = T)
summary(rate_tot_m2)

anova(rate_tot_m1, rate_tot_m2)

#treatment do not impact the rate of change, However, seasonality does. 

### 7.4.3 Visualization of rate of change ####
# Create sequence of dates for the middle of each month
monthly_dates <- seq(from = as.Date("2022-04-15"), 
                     to = as.Date("2023-04-15"), 
                     by = "month")

# Create prediction data frame with all treatment combinations
monthly_pred_data <- expand.grid(
  Date = monthly_dates,
  predator = c("1", "0"),
  nutrients = c("0", "1")
)

# Add necessary columns for prediction
monthly_pred_data <- monthly_pred_data %>%
  mutate(
    # Month numeric for seasonal components
    month_numeric = as.numeric(format(Date, "%m")),
    # Calculate sine and cosine components matching your model
    sin_year = sin(2 * pi * month_numeric / 12),
    cos_year = cos(2 * pi * month_numeric / 12),
    # Add Lake_ID as NA to use population-level predictions
    Lake_ID = NA
  )

# Create the prediction of the rate of change using the output model 
monthly_pred_data$predicted_rate <- predict(rate_tot_m2, 
                                            newdata = monthly_pred_data, 
                                            type = "response", 
                                            allow.new.levels = TRUE)

# Add treatment labels matching your format
monthly_pred_data <- monthly_pred_data %>%
  mutate(treatment = factor(paste(
    ifelse(predator == "1", "Pred+", "Pred-"),
    ifelse(nutrients == "0", "Nutr-", "Nutr+"),
    sep = "/")))

# Reorder factor levels to match your existing plots
monthly_pred_data$treatment <- factor(monthly_pred_data$treatment, 
                                      levels = c("Pred+/Nutr-", "Pred-/Nutr-", 
                                                 "Pred+/Nutr+", "Pred-/Nutr+"))

#Now we can create a plot where the point represent the real data computed and
#lines represent the predictions. 

actual_dates_plot <- ggplot() +
  # First add the model predictions as smooth lines
  geom_smooth(data = monthly_pred_data, 
              aes(x = Date, y = predicted_rate, color = treatment, group = treatment),
              method = "gam", formula = y ~ s(x, k = 8), 
              se = FALSE, linewidth = 1.2) +
  # Then add the actual data points from tot_d1
  geom_point(data = tot_d1, 
             aes(x = midpoint_date, y = rate_change, color = treatment),
             size = 3, alpha = 0.7) +
  scale_color_manual(values = custom_colors) +
  labs(title = "Rate of Change in Beta Diversity Over Time",
       tag= "D",
       x = "",
       y = "Rate of Change \n(Beta Diversity / Day)",
       color = "Treatment combination") +
  scale_x_date(
    date_breaks = "1 month",
    labels = function(x) {
      month_num <- as.numeric(format(x, "%m"))
      year <- format(x, "%y")
      paste0(english_months[month_num], year)
    },
    limits = c(as.Date("2022-04-01"), as.Date("2023-04-30")),
    expand = c(0.05, 0.05)
  ) +
  theme_jg_right +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

actual_dates_plot

###########################################################
# 8. Beta.pair Analyses glmmTMB  TURNOVER  (sim) ####
###########################################################

range(turn_d1$distance)

# since we have values for turnover equal to 0 we need to add 0.00001 to be able
# to perform the analysis using beta family. 

turn_d1$distance <- turn_d1$distance + 0.00001
range(turn_d1$distance)

## 8.1 models through the year ####
### 8.1.1 Analysis treatment only ####
turn_m1 <- glmmTMB(distance ~ predator * nutrients * as.numeric(days_since_start) 
                   + (1|Lake_ID), 
                   family = beta_family(link= "logit"), 
                   data = turn_d1)

simulateResiduals(turn_m1, plot = T)
summary(turn_m1)
Anova(turn_m1)

turn_m2 <- glmmTMB(distance ~ predator * nutrients * as.numeric(days_since_start), 
                   family = beta_family(link= "logit"), 
                   data = turn_d1)

simulateResiduals(turn_m2, plot = T)
anova(turn_m1, turn_m2)
# interaction significant 
summary(turn_m2)# AIC lower 
Anova(turn_m2)

rm(turn_m2)

### 8.1.2 Taking into account seasonality using sine cosine ####

simulateResiduals(turn_m1, plot = T)
summary(turn_m1)


turn_m1_seasonal <- glmmTMB(distance ~ predator * nutrients * as.numeric(days_since_start) + 
                              sin_year + cos_year + (1|Lake_ID), 
                            family = beta_family(link= "logit"), 
                            data = turn_d1)

simulateResiduals(turn_m1_seasonal, plot = T)
summary(turn_m1_seasonal)
Anova(turn_m1_seasonal)

# Compare models to see if seasonality is important
anova(turn_m1, turn_m1_seasonal)

# season has an impact since cosine is significant. 
rm(turn_m1)

### 8.1.3 Residual analysis environmental variable ####

turn_d1$res_turn <- residuals(turn_m1_seasonal)

res_turn <- glmmTMB(res_turn ~ delta_temperature + delta_do_percent + delta_ph +
                delta_conductivity + delta_turbidity + delta_chlorophyll + (1|Lake_ID),
              data = turn_d1)
summary(res_turn)
Anova(res_turn)

## 8.2 models per 3 months ####
### phase 1: April2022 - June2022 ####
phase1_turn <- turn_d1 %>% 
  filter(Date1 %in% as.Date(c("2022-04-21", "2022-05-19", "2022-06-23")))


p1_turn_m1<- glmmTMB(distance ~ predator * nutrients * as.factor(days_since_start) +
                       (1|Lake_ID), 
                     family = beta_family(link = "logit"), 
                     data = phase1_turn)

simulateResiduals(p1_turn_m1, plot = T)
summary(p1_turn_m1)
Anova(p1_turn_m1)


#### Post-hoc phase 1 ####

p1_turn_emm <- emmeans(p1_turn_m1, ~ predator * nutrients * days_since_start)

# Compare phases within each predator x nutrient combination
pairs_p1_turn <- p1_turn_emm %>% 
  contrast(interaction = c("pairwise"), by = c("predator", "nutrients"), 
           adjust = "tukey")
summary(pairs_p1_turn)


### phase 2: July2022 - September2022 ####
phase2_turn <- turn_d1 %>% 
  filter(Date1 %in% as.Date(c("2022-07-21", "2022-08-18", "2022-09-22")))


p2_turn_m1<- glmmTMB(distance ~ predator * nutrients * as.factor(days_since_start) +
                       (1|Lake_ID), 
                     family = beta_family(link = "logit"), 
                     data = phase2_turn)

simulateResiduals(p2_turn_m1, plot = T)
summary(p2_turn_m1)
Anova(p2_turn_m1)

#### Post-hoc phase 2 ####

p2_turn_emm <- emmeans(p2_turn_m1, ~ predator * nutrients * days_since_start)

# Compare phases within each predator x nutrient combination
pairs_p2_turn <- p2_turn_emm %>% 
  contrast(interaction = c("pairwise"), by = c("predator", "nutrients"), 
           adjust = "tukey")
summary(pairs_p2_turn)


### phase 3: October2022 - December2022 ####

phase3_turn <- turn_d1 %>% 
  filter(Date1 %in% as.Date(c("2022-10-21", "2022-11-18", "2022-12-15")))


p3_turn_m1<- glmmTMB(distance ~ predator * nutrients * as.factor(days_since_start) +
                       (1|Lake_ID), 
                     family = beta_family(link = "logit"), 
                     data = phase3_turn)

simulateResiduals(p3_turn_m1, plot = T)
summary(p3_turn_m1)
Anova(p3_turn_m1)

#### Post-hoc phase 3 ####

p3_turn_emm <- emmeans(p3_turn_m1, ~ predator * nutrients * days_since_start)

# Compare phases within each predator x nutrient combination
pairs_p3_turn <- p3_turn_emm %>% 
  contrast(interaction = c("pairwise"), by = c("predator", "nutrients"), 
           adjust = "tukey")
summary(pairs_p3_turn)


### phase 4: January2023 - March2023 ####
phase4_turn <- turn_d1 %>% 
  filter(Date1 %in% as.Date(c("2023-01-25", "2023-02-21", "2023-03-28")))


p4_turn_m1<- glmmTMB(distance ~ predator * nutrients * as.factor(days_since_start) +
                       (1|Lake_ID), 
                     family = beta_family(link = "logit"), 
                     data = phase4_turn)

simulateResiduals(p4_turn_m1, plot = T)
summary(p4_turn_m1)
Anova(p4_turn_m1)

#### Post-hoc phase 4 ####

p4_turn_emm <- emmeans(p4_turn_m1, ~ predator * nutrients * days_since_start)

# Compare phases within each predator x nutrient combination
pairs_p4_turn <- p4_turn_emm %>% 
  contrast(interaction = c("pairwise"), by = c("predator", "nutrients"), 
           adjust = "tukey")
summary(pairs_p4_turn)


## 8.3 testing whether the treatment effect change through phases ####
# First, combine all your phase data into one dataset
all_phases_turn <- bind_rows(
  mutate(phase1_turn, Phase = "Phase1"),
  mutate(phase2_turn, Phase = "Phase2"),
  mutate(phase3_turn, Phase = "Phase3"),
  mutate(phase4_turn, Phase = "Phase4")
)

# Make Phase a factor with the correct order
all_phases_turn$Phase <- factor(all_phases_turn$Phase, 
                                levels = c("Phase1", "Phase2", "Phase3", "Phase4"))

# Create unified model with Phase interactions
unified_model_turn <- glmmTMB(distance ~ predator * nutrients * Phase + (1|Lake_ID), 
                              family = beta_family(link = "logit"), 
                              data = all_phases_turn)

# Check model
simulateResiduals(unified_model_turn, plot = T)
summary(unified_model_turn)
Anova(unified_model_turn)

# phases is in interaction with treatment which means that the effect of treatment 
# on beta diversity indeed depends on the "phase" we are looking. 

### 8.3.1 Post-hoc comparisons ####
# Get predicted values from the model
emm_turn <- emmeans(unified_model_turn, ~ predator * nutrients * Phase)
emm_turn_df <- as.data.frame(emm_turn)

# Compare phases within each predator x nutrient combination
pairs_phase_turn <- emm_turn %>% 
  contrast(interaction = c("pairwise"), by = c("predator", "nutrients"), 
           adjust = "tukey")
summary(pairs_phase_turn)


# Compare predator effects within each phase x nutrient combination
pairs_pred_turn <- emm_turn %>% 
  contrast(interaction = c("pairwise"), by = c("Phase", "nutrients"), 
           adjust = "tukey")
summary(pairs_pred_turn)

# Compare nutrient effects within each phase x predator combination
pairs_nutr_turn <- emm_turn %>% 
  contrast(interaction = c("pairwise"), by = c("Phase", "predator"), 
           adjust = "tukey")
summary(pairs_nutr_turn)


###########################################################
# 9. Beta.pair Analyses glmmTMB  NESTEDNESS (sne) ####
###########################################################

range(nest_d1$distance)

# since we have values for nestedness equal to 0 we need to add 0.00001 to be able
# to perform the analysis using beta family. 

nest_d1$distance <- nest_d1$distance + 0.00001
range(nest_d1$distance)

## 9.1 models through the year ####
### 9.1.1 Analysis treatment only ####
nest_m1 <- glmmTMB(distance ~ predator * nutrients * as.numeric(days_since_start) 
                   + (1|Lake_ID), 
                   family = beta_family(link= "logit"), 
                   data = nest_d1)

simulateResiduals(nest_m1, plot = T)
summary(nest_m1)
Anova(nest_m1)

nest_m2 <- glmmTMB(distance ~ predator * nutrients * as.numeric(days_since_start), 
                   family = beta_family(link= "logit"), 
                   data = nest_d1)

simulateResiduals(nest_m2, plot = T)
anova(nest_m1, nest_m2)
# interaction not significant however I will still include Lake_ID to take into 
# account that we are smapling the sample lake along the year. 
summary(nest_m2)# AIC lower 
Anova(nest_m2)

rm(nest_m2)

### 9.1.2 Taking into account seasonality using sine cosine ####

simulateResiduals(nest_m1, plot = T)
summary(nest_m1)


nest_m1_seasonal <- glmmTMB(distance ~ predator * nutrients * as.numeric(days_since_start) + 
                              sin_year + cos_year + (1|Lake_ID), 
                            family = beta_family(link= "logit"), 
                            data = nest_d1)

simulateResiduals(nest_m1_seasonal, plot = T)
summary(nest_m1_seasonal)
Anova(nest_m1_seasonal)

# Compare models to see if seasonality is important
anova(nest_m1, nest_m1_seasonal)

# season has an impact since cosine is significant. 
rm(nest_m1)

### 9.1.3 Residual analysis environmental variable ####

nest_d1$res_nest <- residuals(nest_m1_seasonal)

res_nest <- glmmTMB(res_nest ~ delta_temperature + delta_do_percent + delta_ph +
                 delta_conductivity + delta_turbidity + delta_chlorophyll +(1|Lake_ID),
               data = nest_d1)
simulateResiduals(res_nest, plot=T)
summary(res_nest)
Anova(res_nest)

## 9.2 models per 3 months ####
### phase 1: April2022 - June2022 ####
phase1_nest <- nest_d1 %>% 
  filter(Date1 %in% as.Date(c("2022-04-21", "2022-05-19", "2022-06-23")))


p1_nest_m1<- glmmTMB(distance ~ predator * nutrients * as.factor(days_since_start) +
                       (1|Lake_ID), 
                     family = beta_family(link = "logit"), 
                     data = phase1_nest)

simulateResiduals(p1_nest_m1, plot = T)
summary(p1_nest_m1)
Anova(p1_nest_m1)

#### Post-hoc phase 1 ####

p1_nest_emm <- emmeans(p1_nest_m1, ~ predator * nutrients * days_since_start)

# Compare phases within each predator x nutrient combination
pairs_p1_nest <- p1_nest_emm %>% 
  contrast(interaction = c("pairwise"), by = c("predator", "nutrients"), 
           adjust = "tukey")
summary(pairs_p1_nest)


### phase 2: July2022 - September2022 ####
phase2_nest <- nest_d1 %>% 
  filter(Date1 %in% as.Date(c("2022-07-21", "2022-08-18", "2022-09-22")))


p2_nest_m1<- glmmTMB(distance ~ predator * nutrients * as.factor(days_since_start) +
                       (1|Lake_ID), 
                     family = beta_family(link = "logit"), 
                     data = phase2_nest)

simulateResiduals(p2_nest_m1, plot = T)
summary(p2_nest_m1)
Anova(p2_nest_m1)


#### Post-hoc phase 2 ####

p2_nest_emm <- emmeans(p2_nest_m1, ~ predator * nutrients * days_since_start)

# Compare phases within each predator x nutrient combination
pairs_p2_nest <- p2_nest_emm %>% 
  contrast(interaction = c("pairwise"), by = c("predator", "nutrients"), 
           adjust = "tukey")
summary(pairs_p2_nest)


### phase 3: October2022 - December2022 ####
phase3_nest <- nest_d1 %>% 
  filter(Date1 %in% as.Date(c("2022-10-21", "2022-11-18", "2022-12-15")))


p3_nest_m1<- glmmTMB(distance ~ predator * nutrients * as.factor(days_since_start) +
                       (1|Lake_ID), 
                     family = beta_family(link = "logit"), 
                     data = phase3_nest)

simulateResiduals(p3_nest_m1, plot = T)
summary(p3_nest_m1)
Anova(p3_nest_m1)

#### Post-hoc phase 3 ####

p3_nest_emm <- emmeans(p3_nest_m1, ~ predator * nutrients * days_since_start)

# Compare phases within each predator x nutrient combination
pairs_p3_nest <- p3_nest_emm %>% 
  contrast(interaction = c("pairwise"), by = c("predator", "nutrients"), 
           adjust = "tukey")
summary(pairs_p3_nest)


### phase 4: January2023 - March2023 ####
phase4_nest <- nest_d1 %>% 
  filter(Date1 %in% as.Date(c("2023-01-25", "2023-02-21", "2023-03-28")))


p4_nest_m1<- glmmTMB(distance ~ predator * nutrients * as.factor(days_since_start) +
                       (1|Lake_ID), 
                     family = beta_family(link = "logit"), 
                     data = phase4_nest)

simulateResiduals(p4_nest_m1, plot = T)
summary(p4_nest_m1)
Anova(p4_nest_m1)

#### Post-hoc phase 4 ####

p4_nest_emm <- emmeans(p4_nest_m1, ~ predator * nutrients * days_since_start)

# Compare phases within each predator x nutrient combination
pairs_p4_nest <- p4_nest_emm %>% 
  contrast(interaction = c("pairwise"), by = c("predator", "nutrients"), 
           adjust = "tukey")
summary(pairs_p4_nest)


## 9.3 testing whether the treatment effect change through phases ####
# First, combine all your phase data into one dataset
all_phases_nest <- bind_rows(
  mutate(phase1_nest, Phase = "Phase1"),
  mutate(phase2_nest, Phase = "Phase2"),
  mutate(phase3_nest, Phase = "Phase3"),
  mutate(phase4_nest, Phase = "Phase4")
)

# Make Phase a factor with the correct order
all_phases_nest$Phase <- factor(all_phases_nest$Phase, 
                                levels = c("Phase1", "Phase2", "Phase3", "Phase4"))

# Create unified model with Phase interactions
unified_model_nest <- glmmTMB(distance ~ predator * nutrients * Phase + (1|Lake_ID), 
                              family = beta_family(link = "logit"), 
                              data = all_phases_nest)

# Check model
simulateResiduals(unified_model_nest, plot = T)
summary(unified_model_nest)
Anova(unified_model_nest)

# phases is in interaction with treatment which means that the effect of treatment 
# on beta diversity indeed depends on the "phase" we are looking. 


### 9.3.1 Post-hoc comparisons ####
# Get predicted values from the model
emm_nest <- emmeans(unified_model_nest, ~ predator * nutrients * Phase)
emm_nest_df <- as.data.frame(emm_nest)

# Compare phases within each predator x nutrient combination
pairs_phase_nest <- emm_nest %>% 
  contrast(interaction = c("pairwise"), by = c("predator", "nutrients"), 
           adjust = "tukey")
summary(pairs_phase_nest)


# Compare predator effects within each phase x nutrient combination
pairs_pred_nest <- emm_nest %>% 
  contrast(interaction = c("pairwise"), by = c("Phase", "nutrients"), 
           adjust = "tukey")
summary(pairs_pred_nest)

# Compare nutrient effects within each phase x predator combination
pairs_nutr_nest <- emm_nest %>% 
  contrast(interaction = c("pairwise"), by = c("Phase", "predator"), 
           adjust = "tukey")
summary(pairs_nutr_nest)

# Compare nutrient effects within each phase x predator combination
pairs_nutr_nest <- emm_nest %>% 
  contrast(interaction = c("pairwise"), by = c("Phase", "predator"), 
           adjust = "tukey")
summary(pairs_nutr_nest)


###########################################################
# 10. TABLES: Beta.pair Analyses Through the year result report ####
###########################################################
## 10.1 Model selected ####

# Regarding the total diversity, I selected the model tot_m1
# distance ~ predator * nutrients * as.numeric(timepoint1) + (1 |Lake_ID)
# Since the seasonality was not significant 
summary(tot_m1)
Anova(tot_m1)

# Regarding the Turnover I've selected the model turn_m1_seasonal, since the 
# season seems to impact the trunover
summary(turn_m1_seasonal)
Anova(turn_m1_seasonal)

# Regarding the Nestedness I've also selected the model nest_m1_seasonal, 
# since season seems to impact nestedness component. 
summary(nest_m1_seasonal)
Anova(nest_m1_seasonal)

## 10.2 TABLE S6: Report summary output ####
# Combine into a single table
# Extract model summary statistics 
tot_m1_tidy <- tidy(tot_m1, conf.int = TRUE)  # No seasonal effect for total
turn_m1_tidy <- tidy(turn_m1_seasonal, conf.int = TRUE)  # Use seasonal model for turnover
nest_m1_tidy <- tidy(nest_m1_seasonal, conf.int = TRUE)  # Use seasonal model for nestedness

# Combine into a single table
model_results <- bind_rows(
  mutate(tot_m1_tidy, model = "Total Beta Diversity"),
  mutate(turn_m1_tidy, model = "Turnover Component (with seasonality)"),
  mutate(nest_m1_tidy, model = "Nestedness Component (with seasonality)")
)

# Rename columns for clarity and select relevant ones
model_results <- model_results %>%
  select(model, term, estimate, std.error, statistic, p.value, conf.low, conf.high) %>%
  rename(
    "Component" = model,
    "Term" = term,
    "Estimate" = estimate,
    "Std. Error" = std.error,
    "z value" = statistic,
    "p-value" = p.value,
    "Lower CI" = conf.low,
    "Upper CI" = conf.high
  )

# Add significance indicators
model_results$Significance <- ifelse(model_results$`p-value` < 0.001, '***',
                                     ifelse(model_results$`p-value` < 0.01, '**',
                                            ifelse(model_results$`p-value` < 0.05, '*',
                                                   ifelse(model_results$`p-value` < 0.1, '.', ''))))

# Combine p-value with significance
model_results$FormattedPvalue <- paste0(
  sprintf("%.3f", model_results$`p-value`), 
  " ", 
  model_results$Significance
)

# First select and rename columns in the data frame
model_results_formatted <- model_results %>%
  select(-Significance, -`p-value`) %>%  # Remove separate columns
  rename(`p-value` = FormattedPvalue)

# Then create the flextable with the formatted data frame
ft_models <- flextable(model_results_formatted) %>%
  colformat_double(digits = 3) %>%
  bold(i = model_results$`p-value` < 0.05, j = "p-value") %>%
  theme_vanilla() %>%
  add_header_lines("Beta Diversity Components Model Results") %>%
  autofit()


# Save to Word
doc <- read_docx()
doc <- body_add_flextable(doc, ft_models)
#print(doc, target = "tables/TableS6.beta.pair_throughyear_model_results.docx")

## 10.3 TABLE S7: Report Anova output ####
# Anova on each model to get the test statistics
tot_anova <- Anova(tot_m1)
turn_seasonal_anova <- Anova(turn_m1_seasonal)
nest_seasonal_anova <- Anova(nest_m1_seasonal)

# Convert Anova results to data frames
tot_anova_df <- as.data.frame(tot_anova)
tot_anova_df$Effect <- rownames(tot_anova_df)
rownames(tot_anova_df) <- NULL

turn_seasonal_anova_df <- as.data.frame(turn_seasonal_anova)
turn_seasonal_anova_df$Effect <- rownames(turn_seasonal_anova_df)
rownames(turn_seasonal_anova_df) <- NULL

nest_seasonal_anova_df <- as.data.frame(nest_seasonal_anova)
nest_seasonal_anova_df$Effect <- rownames(nest_seasonal_anova_df)
rownames(nest_seasonal_anova_df) <- NULL

# Create a combined ANOVA results table
anova_results <- bind_rows(
  mutate(tot_anova_df, Component = "Total Beta Diversity"),
  mutate(turn_seasonal_anova_df, Component = "Turnover Component (seasonal)"),
  mutate(nest_seasonal_anova_df, Component = "Nestedness Component (seasonal)")
)

# Rename columns for clarity
names(anova_results)[names(anova_results) == "Chisq"] <- "Chi-Square"
names(anova_results)[names(anova_results) == "Pr(>Chisq)"] <- "p-value"

# Reorder columns
anova_results <- anova_results[, c("Component", "Effect", "Chi-Square", "Df", "p-value")]

# Add significance indicators
anova_results$Significance <- ifelse(anova_results$`p-value` < 0.001, '***',
                                     ifelse(anova_results$`p-value` < 0.01, '**',
                                            ifelse(anova_results$`p-value` < 0.05, '*',
                                                   ifelse(anova_results$`p-value` < 0.1, '.', ''))))

# Combine p-value with significance
anova_results$FormattedPvalue <- paste0(
  sprintf("%.3f", anova_results$`p-value`), 
  " ", 
  anova_results$Significance
)

# First select and rename columns in the data frame
anova_results_formatted <- anova_results %>%
  select(-Significance, -`p-value`) %>%  # Remove separate columns
  rename(`p-value` = FormattedPvalue)


# Create a flextable
ft_anova <- flextable(anova_results_formatted) %>%
  colformat_double(digits = 3) %>%
  bold(i = anova_results$`p-value` < 0.05, j = "p-value") %>%
  theme_vanilla() %>%
  add_header_lines("ANOVA Results for Beta Diversity Components") %>%
  autofit()

# Save to Word
doc <- read_docx()
doc <- body_add_flextable(doc, ft_anova)
#print(doc, target = "tables/TableS7.beta_pair_through_year_anova_results.docx")

###########################################################
# 10b. TABLES: Residual analyses environmental variables ####
###########################################################
# R² helper function for Gaussian GLMMs
r2_glmm <- function(model) {
  vc        <- VarCorr(model)
  # glmmTMB returns a nested list: $cond$Lake_ID etc.
  var_rand  <- sum(sapply(unlist(vc, recursive = FALSE), 
                          function(x) attr(x, "stddev")^2))
  var_fixed <- var(predict(model, re.form = NA))
  var_resid <- sigma(model)^2
  var_total <- var_fixed + var_rand + var_resid
  data.frame(R2m = round(var_fixed / var_total, 3),
             R2c = round((var_fixed + var_rand) / var_total, 3))
}

## 10b.1 TABLE S8: Residual GLMM summary (with R²) ####

res_tot_tidy  <- tidy(res_tot,  effects = "fixed", conf.int = TRUE) %>%
  mutate(model = "Total Beta Diversity")
res_turn_tidy <- tidy(res_turn, effects = "fixed", conf.int = TRUE) %>%
  mutate(model = "Turnover Component")
res_nest_tidy <- tidy(res_nest, effects = "fixed", conf.int = TRUE) %>%
  mutate(model = "Nestedness Component")

res_results <- bind_rows(res_tot_tidy, res_turn_tidy, res_nest_tidy)

rm(res_tot_tidy, res_turn_tidy, res_nest_tidy)

# Compute R² for each residual model
r2_res_df <- data.frame(
  model          = c("Total Beta Diversity", "Turnover Component", "Nestedness Component"),
  R2_marginal    = c(r2_glmm(res_tot)$R2m,  r2_glmm(res_turn)$R2m,  r2_glmm(res_nest)$R2m),
  R2_conditional = c(r2_glmm(res_tot)$R2c,  r2_glmm(res_turn)$R2c,  r2_glmm(res_nest)$R2c)
)

res_results <- res_results %>%
  select(model, term, estimate, std.error, statistic, p.value, conf.low, conf.high) %>%
  rename(
    "Component"  = model,
    "Term"       = term,
    "Estimate"   = estimate,
    "Std. Error" = std.error,
    "z value"    = statistic,
    "p-value"    = p.value,
    "Lower CI"   = conf.low,
    "Upper CI"   = conf.high
  ) %>%
  left_join(r2_res_df, by = c("Component" = "model")) %>%
  group_by(Component) %>%
  mutate(
    R2_marginal    = ifelse(row_number() == 1, R2_marginal,    NA),
    R2_conditional = ifelse(row_number() == 1, R2_conditional, NA)
  ) %>%
  ungroup() %>%
  rename("R²m" = R2_marginal, "R²c" = R2_conditional)

res_results$Significance <- ifelse(res_results$`p-value` < 0.001, '***',
                                   ifelse(res_results$`p-value` < 0.01,  '**',
                                          ifelse(res_results$`p-value` < 0.05,  '*',
                                                 ifelse(res_results$`p-value` < 0.1,   '.', ''))))

res_results$FormattedPvalue <- paste0(
  sprintf("%.3f", res_results$`p-value`), " ", res_results$Significance)

res_results_formatted <- res_results %>%
  select(-Significance, -`p-value`) %>%
  rename(`p-value` = FormattedPvalue)

ft_res <- flextable(res_results_formatted) %>%
  colformat_double(digits = 3) %>%
  bold(i = res_results$`p-value` < 0.05, j = "p-value") %>%
  theme_vanilla() %>%
  add_header_lines("Residual Analysis: Effect of Environmental Variables on Beta Diversity Residuals") %>%
  autofit()

doc <- read_docx()
doc <- body_add_flextable(doc, ft_res)
#print(doc, target = "tables/TableS8_res_env_glmm_summary.docx")

rm(ft_res, res_results, res_results_formatted, r2_res_df, doc)


## 10b.2 TABLE S9: Residual ANOVA ####

res_tot_anova  <- Anova(res_tot)
res_turn_anova <- Anova(res_turn)
res_nest_anova <- Anova(res_nest)

res_tot_anova_df          <- as.data.frame(res_tot_anova)
res_tot_anova_df$Effect   <- rownames(res_tot_anova_df)
rownames(res_tot_anova_df) <- NULL
res_tot_anova_df$Component <- "Total Beta Diversity"

res_turn_anova_df          <- as.data.frame(res_turn_anova)
res_turn_anova_df$Effect   <- rownames(res_turn_anova_df)
rownames(res_turn_anova_df) <- NULL
res_turn_anova_df$Component <- "Turnover Component"

res_nest_anova_df          <- as.data.frame(res_nest_anova)
res_nest_anova_df$Effect   <- rownames(res_nest_anova_df)
rownames(res_nest_anova_df) <- NULL
res_nest_anova_df$Component <- "Nestedness Component"

res_anova_results <- bind_rows(res_tot_anova_df, res_turn_anova_df, res_nest_anova_df)

rm(res_tot_anova, res_turn_anova, res_nest_anova,
   res_tot_anova_df, res_turn_anova_df, res_nest_anova_df)

names(res_anova_results)[names(res_anova_results) == "Chisq"]      <- "Chi-Square"
names(res_anova_results)[names(res_anova_results) == "Pr(>Chisq)"] <- "p-value"

res_anova_results <- res_anova_results[, c("Component", "Effect", "Chi-Square", "Df", "p-value")]

res_anova_results$Significance <- ifelse(res_anova_results$`p-value` < 0.001, '***',
                                         ifelse(res_anova_results$`p-value` < 0.01,  '**',
                                                ifelse(res_anova_results$`p-value` < 0.05,  '*',
                                                       ifelse(res_anova_results$`p-value` < 0.1,   '.', ''))))

res_anova_results$FormattedPvalue <- paste0(
  sprintf("%.3f", res_anova_results$`p-value`), " ", res_anova_results$Significance)

res_anova_results_formatted <- res_anova_results %>%
  select(-Significance, -`p-value`) %>%
  rename(`p-value` = FormattedPvalue)

ft_res_anova <- flextable(res_anova_results_formatted) %>%
  colformat_double(digits = 3) %>%
  bold(i = res_anova_results$`p-value` < 0.05, j = "p-value") %>%
  theme_vanilla() %>%
  add_header_lines("ANOVA Results: Effect of Environmental Variables on Beta Diversity Residuals") %>%
  autofit()

doc <- read_docx()
doc <- body_add_flextable(doc, ft_res_anova)
#print(doc, target = "tables/TableS9_res_env_anova.docx")

rm(ft_res_anova, res_anova_results, res_anova_results_formatted, doc)

## 10b.2 FIGURE S8: Coefficient plot residual analyses ####

res_tot_tidy  <- tidy(res_tot,  effects = "fixed", conf.int = TRUE) %>% mutate(Component = "Total Beta Diversity")
res_turn_tidy <- tidy(res_turn, effects = "fixed", conf.int = TRUE) %>% mutate(Component = "Turnover")
res_nest_tidy <- tidy(res_nest, effects = "fixed", conf.int = TRUE) %>% mutate(Component = "Nestedness")

res_plot_data <- bind_rows(res_tot_tidy, res_turn_tidy, res_nest_tidy) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    term = factor(term,
                  levels = c("delta_temperature", "delta_do_percent", "delta_ph",
                             "delta_conductivity", "delta_turbidity", "delta_chlorophyll"),
                  labels = c("Temperature", "DO%", "pH",
                             "Conductivity", "Turbidity", "Chlorophyll")),
    Component = factor(Component,
                       levels = c("Total Beta Diversity", "Turnover", "Nestedness")),
    significant = p.value < 0.05
  )

FigS8 <- ggplot(res_plot_data, aes(x = estimate, y = term, color = significant)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2, linewidth = 0.8) +
  facet_wrap(~ Component, ncol = 3) +
  scale_color_manual(values = c("TRUE" = "#d6604d", "FALSE" = "grey40"),
                     labels = c("TRUE" = "p < 0.05", "FALSE" = "n.s."),
                     name = "") +
  labs(x = "Regression coefficient (± 95% CI)",
       y = "",
       tag = "") +
  theme_jg_right +
  theme(panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5),
        strip.text = element_text(size = 11),
        legend.position = "bottom")

FigS8

#ggsave("figure/S8.FigS8_res_env_coefficient_plot.png", FigS8, width = 10, height = 5, dpi = 300)
#ggsave("figure/S8.FigS8_res_env_coefficient_plot.pdf", FigS8, width = 10, height = 5, dpi = 300)

rm(res_tot_tidy, res_turn_tidy, res_nest_tidy, res_plot_data, FigS8)


###########################################################
# 11. TABLES: Beta.pair Analyses per phase result report ####
###########################################################
## 11.1 Functions to extract model outputs ####

# First, extract the ANOVA results for all models
extract_anova_results <- function(model, component, phase) {
  # Extract ANOVA results
  anova_result <- car::Anova(model)
  
  # Convert to data frame
  anova_df <- as.data.frame(anova_result)
  anova_df$Effect <- rownames(anova_df)
  rownames(anova_df) <- NULL
  
  # Add component and phase info
  anova_df$Component <- component
  anova_df$Phase <- phase
  
  return(anova_df)
}

# Extract model summary results
extract_model_summary <- function(model, component, phase) {
  # Try to extract using tidy() - this handles potential errors
  tryCatch({
    model_summary <- tidy(model)
    model_summary$Component <- component
    model_summary$Phase <- phase
    return(model_summary)
  }, error = function(e) {
    warning(paste("Could not extract summary for", component, phase, ":", e$message))
    return(NULL)
  })
}

# Extract summary for each model, including all parameters
extract_model_results <- function(model, component, phase) {
  # Get model summary
  model_tidy <- extract_model_summary(model, component, phase)
  
  # Return all parameters if summary extraction was successful
  if (!is.null(model_tidy)) {
    results <- model_tidy %>%
      select(Component, Phase, term, estimate, std.error, statistic, p.value)
    
    return(results)
  } else {
    return(NULL)
  }
}

## 11.2 TABLE S10: Summary tables for all components BETWEEN phases ####
# Combine all model results
all_model_results <- bind_rows(
  # Total Beta Diversity
  extract_model_results(unified_model_tot, "Total Beta Diversity", "Between Phases"),
  
  # Turnover Component
  extract_model_results(unified_model_turn, "Turnover", "Between Phases"),
  
  # Nestedness Component
  extract_model_results(unified_model_nest, "Nestedness", "Between Phases")
)

# Add significance indicators
all_model_results$Significance <- ifelse(all_model_results$p.value < 0.001, '***',
                                         ifelse(all_model_results$p.value < 0.01, '**',
                                                ifelse(all_model_results$p.value < 0.05, '*',
                                                       ifelse(all_model_results$p.value < 0.1, '.', ''))))

# Combine p-value with significance
all_model_results$FormattedPvalue <- paste0(
  sprintf("%.3f", all_model_results$p.value), 
  " ", 
  all_model_results$Significance
)

# Format the full model results table
model_table_formatted <- all_model_results %>%
  select(Component, Phase, term, estimate, std.error, statistic, FormattedPvalue) %>%
  rename(
    "Term" = term,
    "Estimate" = estimate,
    "Std. Error" = std.error,
    "z value" = statistic,
    "p-value" = FormattedPvalue
  ) %>%
  # Reorder phases and components for better readability
  mutate(
    Phase = factor(Phase, 
                   levels = c("Phase 1 (Apr-Jun)", "Phase 2 (Jul-Sep)", 
                              "Phase 3 (Oct-Dec)", "Phase 4 (Jan-Mar)",
                              "Unified (All Phases)")),
    Component = factor(Component, 
                       levels = c("Total Beta Diversity", "Turnover", "Nestedness"))
  ) %>%
  arrange(Component, Phase, Term)

# Create a flextable for the full model results
ft_model_results <- flextable(model_table_formatted) %>%
  colformat_double(digits = 3) %>%
  bold(i = grepl("\\*", model_table_formatted$`p-value`), j = "p-value") %>%
  theme_vanilla() %>%
  add_header_lines("Full Model Results for Beta Diversity Components between Phase") %>%
  autofit()


# Save the full model results table
doc_model_results <- read_docx()
doc_model_results <- body_add_flextable(doc_model_results, ft_model_results)
#print(doc_model_results, target = "tables/TableS10_beta.pair_BETWEENphase_full_model_results.docx")


## 11.3 TABLE S11: Create ANOVA tables all components between phases ####

# Extract ANOVA results for each model
all_anova_results <- bind_rows(
  # Total Beta Diversity
  extract_anova_results(unified_model_tot, "Total Beta Diversity", "Between Phases"),
  
  # Turnover Component
  extract_anova_results(unified_model_turn, "Turnover", "Between Phases"),
  
  # Nestedness Component
  extract_anova_results(unified_model_nest, "Nestedness", "Between Phases")
)

# Rename columns for clarity
names(all_anova_results)[names(all_anova_results) == "Chisq"] <- "Chi-Square"
names(all_anova_results)[names(all_anova_results) == "Pr(>Chisq)"] <- "p-value"

# Add significance indicators
all_anova_results$Significance <- ifelse(all_anova_results$`p-value` < 0.001, '***',
                                         ifelse(all_anova_results$`p-value` < 0.01, '**',
                                                ifelse(all_anova_results$`p-value` < 0.05, '*',
                                                       ifelse(all_anova_results$`p-value` < 0.1, '.', ''))))

# Combine p-value with significance
all_anova_results$FormattedPvalue <- paste0(
  sprintf("%.3f", all_anova_results$`p-value`), 
  " ", 
  all_anova_results$Significance
)

# Format the ANOVA table
anova_table_formatted <- all_anova_results %>%
  select(Component, Phase, Effect, `Chi-Square`, Df, FormattedPvalue) %>%
  rename(`p-value` = FormattedPvalue) %>%
  # Reorder phases and components for better readability
  mutate(
    Phase = factor(Phase),
    Component = factor(Component, 
                       levels = c("Total Beta Diversity", "Turnover", "Nestedness"))
  ) %>%
  arrange(Component, Phase, Effect)

# Create a flextable for the ANOVA results
ft_anova <- flextable(anova_table_formatted) %>%
  colformat_double(digits = 3) %>%
  bold(i = grepl("\\*", anova_table_formatted$`p-value`), j = "p-value") %>%
  theme_vanilla() %>%
  add_header_lines("ANOVA Results for Beta Diversity Components between Phase") %>%
  autofit()


# Save the ANOVA table
doc_anova <- read_docx()
doc_anova <- body_add_flextable(doc_anova, ft_anova)
print(doc_anova, target = "tables/TableS11_beta.pair_BETWEENphase_anova_results.docx")


## 11.5 TABLE S12: Create combined summary tables for all components within phases  ####
# Combine all model results
all_model_results <- bind_rows(
  # Total Beta Diversity
  extract_model_results(p1_tot_m1, "Total Beta Diversity", "Phase 1 (Apr-Jun)"),
  extract_model_results(p2_tot_m1, "Total Beta Diversity", "Phase 2 (Jul-Sep)"),
  extract_model_results(p3_tot_m1, "Total Beta Diversity", "Phase 3 (Oct-Dec)"),
  extract_model_results(p4_tot_m1, "Total Beta Diversity", "Phase 4 (Jan-Mar)"),
  
  # Turnover Component
  extract_model_results(p1_turn_m1, "Turnover", "Phase 1 (Apr-Jun)"),
  extract_model_results(p2_turn_m1, "Turnover", "Phase 2 (Jul-Sep)"),
  extract_model_results(p3_turn_m1, "Turnover", "Phase 3 (Oct-Dec)"),
  extract_model_results(p4_turn_m1, "Turnover", "Phase 4 (Jan-Mar)"),
  
  # Nestedness Component
  extract_model_results(p1_nest_m1, "Nestedness", "Phase 1 (Apr-Jun)"),
  extract_model_results(p2_nest_m1, "Nestedness", "Phase 2 (Jul-Sep)"),
  extract_model_results(p3_nest_m1, "Nestedness", "Phase 3 (Oct-Dec)"),
  extract_model_results(p4_nest_m1, "Nestedness", "Phase 4 (Jan-Mar)")
)

# Add significance indicators
all_model_results$Significance <- ifelse(all_model_results$p.value < 0.001, '***',
                                         ifelse(all_model_results$p.value < 0.01, '**',
                                                ifelse(all_model_results$p.value < 0.05, '*',
                                                       ifelse(all_model_results$p.value < 0.1, '.', ''))))

# Combine p-value with significance
all_model_results$FormattedPvalue <- paste0(
  sprintf("%.3f", all_model_results$p.value), 
  " ", 
  all_model_results$Significance
)

# Format the full model results table
model_table_formatted <- all_model_results %>%
  select(Component, Phase, term, estimate, std.error, statistic, FormattedPvalue) %>%
  rename(
    "Term" = term,
    "Estimate" = estimate,
    "Std. Error" = std.error,
    "z value" = statistic,
    "p-value" = FormattedPvalue
  ) %>%
  # Reorder phases and components for better readability
  mutate(
    Phase = factor(Phase, 
                   levels = c("Phase 1 (Apr-Jun)", "Phase 2 (Jul-Sep)", 
                              "Phase 3 (Oct-Dec)", "Phase 4 (Jan-Mar)",
                              "Unified (All Phases)")),
    Component = factor(Component, 
                       levels = c("Total Beta Diversity", "Turnover", "Nestedness"))
  ) %>%
  arrange(Component, Phase, Term)

# Create a flextable for the full model results
ft_model_results <- flextable(model_table_formatted) %>%
  colformat_double(digits = 3) %>%
  bold(i = grepl("\\*", model_table_formatted$`p-value`), j = "p-value") %>%
  theme_vanilla() %>%
  add_header_lines("Full Model Results for Beta Diversity Components by Phase") %>%
  autofit()


# Save the full model results table
doc_model_results <- read_docx()
doc_model_results <- body_add_flextable(doc_model_results, ft_model_results)
#print(doc_model_results, target = "tables/TableS12_beta.pair_phase_full_model_results.docx")

## 11.6 TABLE S13: Create combined ANOVA tables for all components within phases ####

# Extract ANOVA results for each model
all_anova_results <- bind_rows(
  # Total Beta Diversity
  extract_anova_results(p1_tot_m1, "Total Beta Diversity", "Phase 1 (Apr-Jun)"),
  extract_anova_results(p2_tot_m1, "Total Beta Diversity", "Phase 2 (Jul-Sep)"),
  extract_anova_results(p3_tot_m1, "Total Beta Diversity", "Phase 3 (Oct-Dec)"),
  extract_anova_results(p4_tot_m1, "Total Beta Diversity", "Phase 4 (Jan-Mar)"),
  
  # Turnover Component
  extract_anova_results(p1_turn_m1, "Turnover", "Phase 1 (Apr-Jun)"),
  extract_anova_results(p2_turn_m1, "Turnover", "Phase 2 (Jul-Sep)"),
  extract_anova_results(p3_turn_m1, "Turnover", "Phase 3 (Oct-Dec)"),
  extract_anova_results(p4_turn_m1, "Turnover", "Phase 4 (Jan-Mar)"),
  
  # Nestedness Component
  extract_anova_results(p1_nest_m1, "Nestedness", "Phase 1 (Apr-Jun)"),
  extract_anova_results(p2_nest_m1, "Nestedness", "Phase 2 (Jul-Sep)"),
  extract_anova_results(p3_nest_m1, "Nestedness", "Phase 3 (Oct-Dec)"),
  extract_anova_results(p4_nest_m1, "Nestedness", "Phase 4 (Jan-Mar)")
)

# Rename columns for clarity
names(all_anova_results)[names(all_anova_results) == "Chisq"] <- "Chi-Square"
names(all_anova_results)[names(all_anova_results) == "Pr(>Chisq)"] <- "p-value"

# Add significance indicators
all_anova_results$Significance <- ifelse(all_anova_results$`p-value` < 0.001, '***',
                                         ifelse(all_anova_results$`p-value` < 0.01, '**',
                                                ifelse(all_anova_results$`p-value` < 0.05, '*',
                                                       ifelse(all_anova_results$`p-value` < 0.1, '.', ''))))

# Combine p-value with significance
all_anova_results$FormattedPvalue <- paste0(
  sprintf("%.3f", all_anova_results$`p-value`), 
  " ", 
  all_anova_results$Significance
)

# Format the ANOVA table
anova_table_formatted <- all_anova_results %>%
  select(Component, Phase, Effect, `Chi-Square`, Df, FormattedPvalue) %>%
  rename(`p-value` = FormattedPvalue) %>%
  # Reorder phases and components for better readability
  mutate(
    Phase = factor(Phase, 
                   levels = c("Phase 1 (Apr-Jun)", "Phase 2 (Jul-Sep)", 
                              "Phase 3 (Oct-Dec)", "Phase 4 (Jan-Mar)",
                              "Unified (All Phases)")),
    Component = factor(Component, 
                       levels = c("Total Beta Diversity", "Turnover", "Nestedness"))
  ) %>%
  arrange(Component, Phase, Effect)

# Create a flextable for the ANOVA results
ft_anova <- flextable(anova_table_formatted) %>%
  colformat_double(digits = 3) %>%
  bold(i = grepl("\\*", anova_table_formatted$`p-value`), j = "p-value") %>%
  theme_vanilla() %>%
  add_header_lines("ANOVA Results for Beta Diversity Components by Phase") %>%
  autofit()


# Save the ANOVA table
doc_anova <- read_docx()
doc_anova <- body_add_flextable(doc_anova, ft_anova)
print(doc_anova, target = "tables/TableS13_beta.pair_phase_anova_results.docx")

## 11.7 TABLE S14: Post_hoc within treatment within phases ####
create_phase_specific_posthoc_table <- function() {
  # Function to format p-values with significance indicators
  add_significance <- function(p_value) {
    significance <- ifelse(p_value < 0.001, '***',
                           ifelse(p_value < 0.01, '**',
                                  ifelse(p_value < 0.05, '*',
                                         ifelse(p_value < 0.1, '.', ''))))
    return(paste0(sprintf("%.4f", p_value), " ", significance))
  }
  
  # Function to convert day numbers to month labels
  day_to_month <- function(day) {
    case_when(
      day == 0 ~ "Apr/May",
      day == 28 ~ "May/Jun",
      day == 63 ~ "Jun/Jul",
      day == 91 ~ "Jul/Aug",
      day == 119 ~ "Aug/Sep",
      day == 154 ~ "Sep/Oct",
      day == 183 ~ "Oct/Nov",
      day == 211 ~ "Nov/Dec",
      day == 238 ~ "Dec/Jan",
      day == 279 ~ "Jan/Feb",
      day == 306 ~ "Feb/Mar",
      day == 341 ~ "Mar/Apr",
      TRUE ~ as.character(day)
    )
  }
  
  # Function to convert contrast strings (e.g., "28 - 0") to month labels
  convert_contrast_to_months <- function(contrast_str) {
    # Extract the day numbers from the contrast string
    days <- as.numeric(unlist(strsplit(gsub(" ", "", contrast_str), "-")))
    
    if(length(days) == 2) {
      # Convert days to month labels
      month1 <- day_to_month(days[1])
      month2 <- day_to_month(days[2])
      
      # Return the formatted month comparison
      return(paste(month1, "-", month2))
    } else {
      # If parsing fails, return the original string
      return(contrast_str)
    }
  }
  
  # Process post-hoc results for each phase and component
  process_phase_pairs <- function(pairs_df, phase, component) {
    # Create a copy to avoid modifying the original
    result_df <- pairs_df
    
    # Add phase and component information
    result_df$Phase <- phase
    result_df$Component <- component
    result_df$FormattedPvalue <- sapply(result_df$p.value, add_significance)
    
    # Rename columns to standardize
    names(result_df)[1] <- "Contrast"  # This will rename days_since_start_pairwise to Contrast
    
    # Convert day contrasts to month labels
    result_df$Contrast <- sapply(result_df$Contrast, convert_contrast_to_months)
    
    # Create Treatment group column - will be updated for actual data
    result_df$Treatment <- paste("Predator:", result_df$predator, 
                                 "| Nutrients:", result_df$nutrients)
    
    return(result_df)
  }
  
  # Convert pairs objects to data frames
  # Phase 1
  pairs_p1_tot_df <- as.data.frame(pairs_p1_tot)
  pairs_p1_turn_df <- as.data.frame(pairs_p1_turn)
  pairs_p1_nest_df <- as.data.frame(pairs_p1_nest)
  
  # Phase 2
  pairs_p2_tot_df <- as.data.frame(pairs_p2_tot)
  pairs_p2_turn_df <- as.data.frame(pairs_p2_turn)
  pairs_p2_nest_df <- as.data.frame(pairs_p2_nest)
  
  # Phase 3
  pairs_p3_tot_df <- as.data.frame(pairs_p3_tot)
  pairs_p3_turn_df <- as.data.frame(pairs_p3_turn)
  pairs_p3_nest_df <- as.data.frame(pairs_p3_nest)
  
  # Phase 4
  pairs_p4_tot_df <- as.data.frame(pairs_p4_tot)
  pairs_p4_turn_df <- as.data.frame(pairs_p4_turn)
  pairs_p4_nest_df <- as.data.frame(pairs_p4_nest)
  
  # Combine all post-hoc results
  all_posthoc_results <- bind_rows(
    # Phase 1 results
    process_phase_pairs(pairs_p1_tot_df, "Phase 1", "Total Beta Diversity"),
    process_phase_pairs(pairs_p1_turn_df, "Phase 1", "Turnover"),
    process_phase_pairs(pairs_p1_nest_df, "Phase 1", "Nestedness"),
    
    # Phase 2 results
    process_phase_pairs(pairs_p2_tot_df, "Phase 2", "Total Beta Diversity"),
    process_phase_pairs(pairs_p2_turn_df, "Phase 2", "Turnover"),
    process_phase_pairs(pairs_p2_nest_df, "Phase 2", "Nestedness"),
    
    # Phase 3 results
    process_phase_pairs(pairs_p3_tot_df, "Phase 3", "Total Beta Diversity"),
    process_phase_pairs(pairs_p3_turn_df, "Phase 3", "Turnover"),
    process_phase_pairs(pairs_p3_nest_df, "Phase 3", "Nestedness"),
    
    # Phase 4 results
    process_phase_pairs(pairs_p4_tot_df, "Phase 4", "Total Beta Diversity"),
    process_phase_pairs(pairs_p4_turn_df, "Phase 4", "Turnover"),
    process_phase_pairs(pairs_p4_nest_df, "Phase 4", "Nestedness")
  )
  
  # Create a custom ordering for components
  component_order <- c("Total Beta Diversity", "Turnover", "Nestedness")
  
  # Create a custom ordering for treatments
  treatment_order <- c(
    "Predator: 0 | Nutrients: 0",
    "Predator: 1 | Nutrients: 0",
    "Predator: 0 | Nutrients: 1",
    "Predator: 1 | Nutrients: 1"
  )
  
  # Create a function to order contrasts within phases
  order_contrasts <- function(contrasts, phase) {
    if(phase == "Phase 1") {
      contrast_order <- c("Apr/May - May/Jun", "Apr/May - Jun/Jul", "May/Jun - Jun/Jul")
    } else if(phase == "Phase 2") {
      contrast_order <- c("Jul/Aug - Aug/Sep", "Jul/Aug - Sep/Oct", "Aug/Sep - Sep/Oct")
    } else if(phase == "Phase 3") {
      contrast_order <- c("Oct/Nov - Nov/Dec", "Oct/Nov - Dec/Jan", "Nov/Dec - Dec/Jan")
    } else if(phase == "Phase 4") {
      contrast_order <- c("Jan/Feb - Feb/Mar", "Jan/Feb - Mar/Apr", "Feb/Mar - Mar/Apr")
    } else {
      return(contrasts)  # Return as is if phase is unknown
    }
    
    # Create a factor with the desired order
    factor(contrasts, levels = contrast_order)
  }
  
  # Format table for output with the desired ordering
  posthoc_table <- all_posthoc_results %>%
    select(Component, Phase, Treatment, Contrast, estimate, SE, FormattedPvalue) %>%
    rename(
      "Estimate" = estimate,
      "Std. Error" = SE,
      "p-value" = FormattedPvalue
    ) %>%
    # Apply custom ordering
    mutate(
      Component = factor(Component, levels = component_order),
      Treatment = factor(Treatment, levels = treatment_order),
      ContrastOrder = pmap(list(Contrast, Phase), function(c, p) order_contrasts(c, p))
    ) %>%
    # Sort by the new factors
    arrange(Component, Phase, Treatment, ContrastOrder) %>%
    # Remove the helper column used for sorting
    select(-ContrastOrder)
  
  # Create a flextable
  ft_posthoc <- flextable(posthoc_table) %>%
    colformat_double(j = c("Estimate", "Std. Error"), digits = 3) %>%
    bold(i = grepl("\\*", posthoc_table$`p-value`), j = "p-value") %>%
    theme_vanilla() %>%
    add_header_lines("Post-hoc Comparisons of Days Within Each Phase") %>%
    merge_v(j = c("Component", "Phase", "Treatment")) %>%
    valign(j = c("Component", "Phase", "Treatment"), valign = "top") %>%
    autofit()
  
  # Save the table
  doc_posthoc <- read_docx()
  doc_posthoc <- body_add_flextable(doc_posthoc, ft_posthoc)
  print(doc_posthoc, target = "tables/TableS14_phase_specific_posthoc_results.docx")
  
  return(ft_posthoc)
}

# Call the function to create the table
phase_posthoc_table <- create_phase_specific_posthoc_table()

##########################################################
# 12. FIGURES (2,3,4):Combining all figures into a comprehensive figure ####
###########################################################
## 12.1 FIGURE 2: Beta Diversity Time Series  ####
### 12.1.1 Create prediction dataframes for each model ####

# Total Beta Diversity (no seasonal effect)
pred_tot_df <- expand.grid(
  days_since_start = c("0", "28", "63", "91", "119", "154","183",
                       "211", "238", "279", "306", "341","377"),  # Name this column properly
  predator = c(0, 1),
  nutrients = c(0, 1),
  Lake_ID = unique(tot_d1$Lake_ID)[1]  # Use first lake for reference
)

# Add treatment labels
pred_tot_df$treatment <- paste0(
  ifelse(pred_tot_df$predator == 1, "Pred+", "Pred-"), "/",
  ifelse(pred_tot_df$nutrients == 1, "Nutr+", "Nutr-")
)

# Generate predictions
pred_tot_df$predicted <- predict(tot_m1, newdata = pred_tot_df, type = "response")


# Turnover Component (with seasonality)
#here I have to add manually the days from start: 
#levels(as.factor(turn_d1$days_since_start))
#[1] "0"   "28"  "63"  "91"  "119" "154" "183" "211" "238" "279" "306" "341"
# + "377"
pred_turn_df <- expand.grid(
  days_since_start = c("0", "28", "63", "91", "119", "154","183",
                       "211", "238", "279", "306", "341","377"),
  predator = c(0, 1),
  nutrients = c(0, 1),
  Lake_ID = unique(turn_d1$Lake_ID)[1]
)

# Add seasonal components based on actual months (4-16 represents Apr 22 - Apr 23)
# April is month 4, so we use 4:16 (April year 1 to April year 2)
pred_turn_df$midmonth_numeric <- 4:16
pred_turn_df$sin_year <- sin(2 * pi * ((pred_turn_df$midmonth_numeric - 1) %% 12 + 1) / 12)
pred_turn_df$cos_year <- cos(2 * pi * ((pred_turn_df$midmonth_numeric - 1) %% 12 + 1) / 12)

# Add treatment labels
pred_turn_df$treatment <- paste0(
  ifelse(pred_turn_df$predator == 1, "Pred+", "Pred-"), "/",
  ifelse(pred_turn_df$nutrients == 1, "Nutr+", "Nutr-")
)

# Generate predictions from seasonal model
pred_turn_df$predicted <- predict(turn_m1_seasonal, newdata = pred_turn_df, type = "response")


# Nestedness Component (with seasonality)
pred_nest_df <- expand.grid(
  days_since_start = c("0", "28", "63", "91", "119", "154","183",
                       "211", "238", "279", "306", "341","377"),
  predator = c(0, 1),
  nutrients = c(0, 1),
  Lake_ID = unique(nest_d1$Lake_ID)[1]
)

# Add seasonal components based on actual months
pred_nest_df$midmonth_numeric <- 4:16
pred_nest_df$sin_year <- sin(2 * pi * ((pred_nest_df$midmonth_numeric - 1) %% 12 + 1) / 12)
pred_nest_df$cos_year <- cos(2 * pi * ((pred_nest_df$midmonth_numeric - 1) %% 12 + 1) / 12)

# Add treatment labels
pred_nest_df$treatment <- paste0(
  ifelse(pred_nest_df$predator == 1, "Pred+", "Pred-"), "/",
  ifelse(pred_nest_df$nutrients == 1, "Nutr+", "Nutr-")
)

# Generate predictions from seasonal model
pred_nest_df$predicted <- predict(nest_m1_seasonal, newdata = pred_nest_df, type = "response")

### 12.1.2 Graphs ####

month_labels=c("Apr22", "May22", "Jun22", "Jul22", "Aug22", "Sep22",
               "Oct22", "Nov22", "Dec22", "Jan23", "Feb23", "March23",
               "Apr23")

# Total Beta Diversity Plot
pred_tot_df$treatment <- factor(pred_tot_df$treatment, 
                                levels = c("Pred-/Nutr-", "Pred+/Nutr-", 
                                           "Pred-/Nutr+", "Pred+/Nutr+"))

tot_plot <- ggplot(pred_tot_df, aes(x = as.numeric(days_since_start), y = predicted, 
                                    color = treatment, group = treatment)) +
  geom_line(linewidth = 1, alpha = 0.8) +
  geom_point(size = 2) +
  labs(x = "", y = "Model Estimate",  
       title = "Total beta diversity") +
  scale_color_manual(values = custom_colors, name = "Treatment") +
  scale_x_continuous(breaks = seq(1, 13), labels = month_labels[seq(1, 13, by = 1)]) +
  theme_jg_right +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1), 
        plot.margin = margin(5, 5, 5, 5))

tot_plot

# Turnover Plot
pred_turn_df$treatment <- factor(pred_turn_df$treatment, 
                                 levels = c("Pred-/Nutr-", "Pred+/Nutr-", 
                                            "Pred-/Nutr+", "Pred+/Nutr+"))

turn_plot <- ggplot(pred_turn_df, aes(x = as.numeric(days_since_start), y = predicted, 
                                      color = treatment, group = treatment)) +
  geom_line(linewidth = 1, alpha =0.8) +
  geom_point(size = 2) +
  labs(x = "", y = "Model Estimate",
       title = "Turnover") +
  scale_color_manual(values = custom_colors, name = "Treatment") +
  scale_x_continuous(breaks = seq(1, 13), labels = month_labels[seq(1, 13, by = 1)]) +
  theme_jg_right +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1), 
        plot.margin = margin(5, 5, 5, 5))

turn_plot

# Nestedness Plot with legend on the right
pred_nest_df$treatment <- factor(pred_nest_df$treatment, 
                                 levels = c("Pred-/Nutr-", "Pred+/Nutr-", 
                                            "Pred-/Nutr+", "Pred+/Nutr+"))

nest_plot <- ggplot(pred_nest_df, aes(x = as.numeric(days_since_start), y = predicted, 
                                      color = treatment, group = treatment)) +
  geom_line(linewidth = 1, alpha = 0.8) +
  geom_point(size = 2) +
  labs(x = "", y = "Model Estimate",
       title = "Nestedness") +
  scale_color_manual(values = custom_colors, name = "Treatment") +
  scale_x_continuous(breaks = seq(1, 13), labels = month_labels[seq(1, 13, by = 1)]) +
  theme_jg_right +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1), 
        plot.margin = margin(5, 5, 5, 5))


# Combine plots with collected legend on the right
Fig2 <- tot_plot / turn_plot / nest_plot +
  plot_layout(heights = c(1, 1, 1), 
              guides = "collect") &  # Collect all legends
  theme(legend.position = "bottom",   # Position on right
        legend.justification = "center",  # Center the legend vertically
        legend.key.height = unit(0, "cm"),
        legend.key.width = unit(0, "cm"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))


Fig2
ggsave("figure/R2_Fig2_beta.pair_throughtime.png", Fig2, width = 6.5, height = 9, dpi = 600)
ggsave("figure/R2_Fig2_beta.pair_throughtime.pdf", Fig2, width = 6.5, height = 9, dpi = 600)

ggsave("figure/Total_beta.pair_throughtime.png", 
       tot_plot, width = 6, height = 5, dpi = 600)

ggsave("figure/Turn_beta.pair_throughtime.png", 
       turn_plot, width = 6, height = 5, dpi = 600)

ggsave("figure/Nest_beta.pair_throughtime.png", 
       nest_plot, width = 6, height = 5, dpi = 600)


## 12.2 FIGURE 3: Between phases unified model ####
### 12.2.1 Panel A: Unified model Three-way Interaction Effects Plot (From Unified Models) ####
# Extract predictions from the unified models

# Function to create interaction plot for a single component
create_interaction_plot <- function(model, component_name, tag, show_y_label = TRUE) {
  emm <- emmeans(model, ~ predator * nutrients * Phase)
  emm_df <- as.data.frame(emm)
  
  # Format for plotting
  emm_df <- emm_df %>%
    mutate(
      treatment = paste0(
        ifelse(predator == 1, "Pred+", "Pred-"),
        "/",
        ifelse(nutrients == 1, "Nutr+", "Nutr-")
      )
    )
  
  emm_df$treatment <- factor(emm_df$treatment, 
                             levels = c("Pred-/Nutr-", "Pred+/Nutr-", "Pred-/Nutr+", "Pred+/Nutr+"))
  
  levels(emm_df$Phase) <- c("Period1", "Period2", "Period3", "Period4")
  
  ggplot(emm_df, aes(x = Phase, y = emmean, color = treatment, group = treatment)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.1, 
                  size = 1) +
    labs(
      title = component_name,
      x = "",
      tag= tag,
      y = if(show_y_label) "Model Estimates (logit)" else "",
      color = "Treatment"
    ) +
    scale_color_manual(values = custom_colors)+
    theme_jg_bottom +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      title = element_text(size =12),
      axis.title.y = element_text(margin = margin(r = 10))
    )
}


# Create plots for each component
p_total <- create_interaction_plot(unified_model_tot, "Total Beta", "A", show_y_label = T)
p_turn <- create_interaction_plot(unified_model_turn, "Turnover", "", show_y_label = F)
p_nest <- create_interaction_plot(unified_model_nest, "Nestedness", "", show_y_label = F)

# Combine using patchwork
# Create Panel A with all three plots side by side
Fig3_panelA <- (p_total | p_turn | p_nest) +
  plot_layout(guides = 'collect') &
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.margin = margin(t = -20), 
    plot.margin = unit(c(1, 1, 1, 1), "pt")
  )

print(Fig3_panelA)

### 12.2.2 Panel B: Significance Overview Plot Unified model ####
# Extract p-values from ANOVA results
extract_anova_results <- function(model_anova) {
  data.frame(
    Effect = rownames(model_anova),
    p_value = model_anova[["Pr(>Chisq)"]]
  )
}

# Assuming you have the ANOVA results
tot_anova <- Anova(unified_model_tot)
turn_anova <- Anova(unified_model_turn)
nest_anova <- Anova(unified_model_nest)

# Convert p-values to significance levels
get_significance <- function(p) {
  case_when(
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    p < 0.1 ~ ".",
    TRUE ~ "ns"
  )
}

# Create significance data frame
significance_unified <- bind_rows(
  extract_anova_results(tot_anova) %>% mutate(Component = "Total Beta"),
  extract_anova_results(turn_anova) %>% mutate(Component = "Turnover"),
  extract_anova_results(nest_anova) %>% mutate(Component = "Nestedness")
) %>%
  mutate(
    Significance = get_significance(p_value),
    # Create a text label for displaying in the middle of the tiles
    SignificanceText = Significance,
    Effect = case_when(
      Effect == "predator" ~ "Pred",
      Effect == "nutrients" ~ "Nutr",
      Effect == "Phase" ~ "Period",
      Effect == "predator:nutrients" ~ "Pred × Nutr",
      Effect == "predator:Phase" ~ "Pred × Period",
      Effect == "nutrients:Phase" ~ "Nutr × Period",
      Effect == "predator:nutrients:Phase" ~ "Pred × Nutr × Period",
      TRUE ~ Effect
    )
  )

# Replace "ns" with empty string for display
significance_unified$SignificanceText[significance_unified$SignificanceText == "ns"] <- ""

significance_unified$Effect <- factor(significance_unified$Effect, 
                                      levels = c("Pred", "Nutr", 
                                                 "Period","Pred × Period", "Nutr × Period",
                                                 "Pred × Nutr", "Pred × Nutr × Period"
                                      ))

significance_unified$Component <- factor(significance_unified$Component, 
                                         levels = c("Nestedness", "Turnover", "Total Beta"))

significance_unified$Significance <- factor(significance_unified$Significance, 
                                            levels = c("***", "**", "*", ".", "ns"))

# Create the plot with stars in the middle and cleaned legend
Fig3_panelB <- ggplot(significance_unified, aes(x = Effect, y = Component, fill = Significance)) +
  geom_tile(color = "white", size = 1) +
  # Add text labels in the middle
  geom_text(aes(label = SignificanceText), 
            fontface = "bold", size = 4) +  # Adjust size as needed
  scale_fill_manual(
    values = c(
      "***" = "#bd0026",
      "**" = "#f03b20",
      "*" = "#fd8d3c",
      "." = "#fecc5c",
      "ns" = "#d9d9d9"
    ),
    # Clean up legend labels - remove stars and brackets
    labels = c(
      "p<0.001",
      "p<0.01",
      "p<0.05",
      "p<0.1",
      "n.s."
    ),
    drop = FALSE
  ) +
  labs(
    title = "Type II Wald chi-square tests",
    tag = "B",
    x = "",
    y = "",
    fill = "Significance Level"
  ) +
  theme_jg_right +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(hjust = 1),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.margin = margin(t = -20), 
    legend.box.spacing = unit(0.5, "pt"),  # space between heatmap and legend
    plot.margin = unit(c(1, 1, 1, 1), "pt"),  # Add some margin
    plot.title = element_text(size = 14),
    aspect.ratio = -0.25
  )


print(Fig3_panelB)

### 12.2.3 FIGURE 3: Between phases plot ####

# Update Figure 3 with the fixed dot plot
Fig3 <- Fig3_panelA / Fig3_panelB +
  plot_layout(heights = c(3, 1))


Fig3

ggsave("figure/R3_Fig3_beta.pair_BETWEEN_PHASE.png", Fig3, width = 8.1, height = 7, dpi = 300)
ggsave("figure/R3_Fig3_beta.pair_BETWEEN_PHASE.pdf", Fig3, width = 8.1, height = 7, dpi = 300)
ggsave("figure/R3_Fig3_beta.pair_BETWEEN_PHASE.svg", Fig3, width = 8.1, height = 7, dpi = 300)


## 12.3 FIGURE S9: Within phase Interaction plots treatment effects ####
# Combine all emmeans data
all_emm <- rbind(
  cbind(as.data.frame(emmeans(p1_tot_m1, ~predator*nutrients)), Phase = "P1", Component = "Total Beta"),
  cbind(as.data.frame(emmeans(p2_tot_m1, ~predator*nutrients)), Phase = "P2", Component = "Total Beta"),
  cbind(as.data.frame(emmeans(p3_tot_m1, ~predator*nutrients)), Phase = "P3", Component = "Total Beta"),
  cbind(as.data.frame(emmeans(p4_tot_m1, ~predator*nutrients)), Phase = "P4", Component = "Total Beta"),
  
  cbind(as.data.frame(emmeans(p1_turn_m1, ~predator*nutrients)), Phase = "P1", Component = "Turnover"),
  cbind(as.data.frame(emmeans(p2_turn_m1, ~predator*nutrients)), Phase = "P2", Component = "Turnover"),
  cbind(as.data.frame(emmeans(p3_turn_m1, ~predator*nutrients)), Phase = "P3", Component = "Turnover"),
  cbind(as.data.frame(emmeans(p4_turn_m1, ~predator*nutrients)), Phase = "P4", Component = "Turnover"),
  
  cbind(as.data.frame(emmeans(p1_nest_m1, ~predator*nutrients)), Phase = "P1", Component = "Nestedness"),
  cbind(as.data.frame(emmeans(p2_nest_m1, ~predator*nutrients)), Phase = "P2", Component = "Nestedness"),
  cbind(as.data.frame(emmeans(p3_nest_m1, ~predator*nutrients)), Phase = "P3", Component = "Nestedness"),
  cbind(as.data.frame(emmeans(p4_nest_m1, ~predator*nutrients)), Phase = "P4", Component = "Nestedness")
)

# Add significance indicators based on model results
# Create functions to extract p-values for all effects
get_predator_pvalue <- function(model) {
  # Extract Anova results for the predator main effect
  anova_results <- car::Anova(model)
  # Get p-value for predator main effect
  p_value <- anova_results["predator", "Pr(>Chisq)"]
  return(p_value)
}

get_nutrient_pvalue <- function(model) {
  # Extract Anova results for the nutrients main effect
  anova_results <- car::Anova(model)
  # Get p-value for nutrients main effect
  p_value <- anova_results["nutrients", "Pr(>Chisq)"]
  return(p_value)
}

get_interaction_pvalue <- function(model) {
  # Extract Anova results for the interaction term
  anova_results <- car::Anova(model)
  # Get p-value for predator:nutrients interaction
  p_value <- anova_results["predator:nutrients", "Pr(>Chisq)"]
  return(p_value)
}

# Helper function to format p-values with significance symbols
format_p_value <- function(p) {
  if(is.na(p)) return("")
  if(p < 0.001) return("***")
  if(p < 0.01) return("**") 
  if(p < 0.05) return("*")
  if(p < 0.1) return(".")
  return("n.s.")
}

# Create significance data frame
sig_data <- data.frame(
  Phase = rep(c("P1", "P2", "P3", "P4"), 3),
  Component = c(
    rep("Total Beta", 4),
    rep("Turnover", 4),
    rep("Nestedness", 4)
  ),
  pred_p = c(
    get_predator_pvalue(p1_tot_m1),
    get_predator_pvalue(p2_tot_m1),
    get_predator_pvalue(p3_tot_m1),
    get_predator_pvalue(p4_tot_m1),
    
    get_predator_pvalue(p1_turn_m1),
    get_predator_pvalue(p2_turn_m1),
    get_predator_pvalue(p3_turn_m1),
    get_predator_pvalue(p4_turn_m1),
    
    get_predator_pvalue(p1_nest_m1),
    get_predator_pvalue(p2_nest_m1),
    get_predator_pvalue(p3_nest_m1),
    get_predator_pvalue(p4_nest_m1)
  ),
  nutr_p = c(
    get_nutrient_pvalue(p1_tot_m1),
    get_nutrient_pvalue(p2_tot_m1),
    get_nutrient_pvalue(p3_tot_m1),
    get_nutrient_pvalue(p4_tot_m1),
    
    get_nutrient_pvalue(p1_turn_m1),
    get_nutrient_pvalue(p2_turn_m1),
    get_nutrient_pvalue(p3_turn_m1),
    get_nutrient_pvalue(p4_turn_m1),
    
    get_nutrient_pvalue(p1_nest_m1),
    get_nutrient_pvalue(p2_nest_m1),
    get_nutrient_pvalue(p3_nest_m1),
    get_nutrient_pvalue(p4_nest_m1)
  ),
  int_p = c(
    get_interaction_pvalue(p1_tot_m1),
    get_interaction_pvalue(p2_tot_m1),
    get_interaction_pvalue(p3_tot_m1),
    get_interaction_pvalue(p4_tot_m1),
    
    get_interaction_pvalue(p1_turn_m1),
    get_interaction_pvalue(p2_turn_m1),
    get_interaction_pvalue(p3_turn_m1),
    get_interaction_pvalue(p4_turn_m1),
    
    get_interaction_pvalue(p1_nest_m1),
    get_interaction_pvalue(p2_nest_m1),
    get_interaction_pvalue(p3_nest_m1),
    get_interaction_pvalue(p4_nest_m1)
  )
)

# Format p-values with significance symbols
sig_data$pred_sig <- sapply(sig_data$pred_p, format_p_value)
sig_data$nutr_sig <- sapply(sig_data$nutr_p, format_p_value)
sig_data$int_sig <- sapply(sig_data$int_p, format_p_value)

# Create text for significance - plain text without label box
sig_data$sig_text <- paste(
  "P: ", sig_data$pred_sig, "\n",
  "N: ", sig_data$nutr_sig, "\n",
  "P×N: ", sig_data$int_sig,
  sep = ""
)

# Set component order
all_emm$Component <- factor(all_emm$Component, 
                            levels = c("Total Beta", "Turnover", "Nestedness"))
sig_data$Component <- factor(sig_data$Component, 
                             levels = c("Total Beta", "Turnover", "Nestedness"))

# Set phase order
all_emm$Phase <- factor(all_emm$Phase, levels = c("P1", "P2", "P3", "P4"))
sig_data$Phase <- factor(sig_data$Phase, levels = c("P1", "P2", "P3", "P4"))

# Fixed position for significance text in each panel
# This will place the significance text at fixed locations relative to plot coordinates
sig_data$x_pos <- 0.45  # Position text at a fixed x position (closer to 1 on x-axis)
sig_data$y_pos <- -Inf  # Position at bottom of panel
sig_data$vjust <- -0.2  # Slightly above the bottom (negative vjust moves text up)

# Ensure numeric x values are treated as discrete factors in plotting
all_emm$predator <- factor(all_emm$predator, levels = c("0", "1"), 
                           labels = c("Pred-", "Pred+"))

# Create a combined faceting variable for Phase and Component
all_emm$facet_var <- interaction(all_emm$Phase, all_emm$Component, sep = "_")
sig_data$facet_var <- interaction(sig_data$Phase, sig_data$Component, sep = "_")

# Define phase labels (more descriptive)
phase_labels <- c(
  "P1" = "Period 1",
  "P2" = "Period 2", 
  "P3" = "Period 3",
  "P4" = "Period 4"
)


# Create the improved plot
FigS9 <- ggplot(all_emm, aes(x = predator, y = emmean, group = nutrients, color = nutrients)) +
  # Base layers
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.2) +
  
  # Add significance text with fixed positioning
  geom_text(data = sig_data, 
            aes(x = x_pos, y = y_pos, label = sig_text, group = NULL), 
            size = 3, color = "black", vjust = sig_data$vjust, hjust = 0) +
  
  # Use facet_grid but with free y scales
  facet_grid(Phase ~ Component, 
             scales = "free_y",
             labeller = labeller(Phase = phase_labels)) +  # Add descriptive phase labels
  
  # Customize labels and title
  labs(
    x = "",
    y = "Estimated marginal mean (logit scale)",
    caption = "Significance codes: *** p<0.001, ** p<0.01, * p<0.05, · p<0.1\nP: Predator effect, N: Nutrient effect, P×N: Interaction"
  ) +
  
  # Apply color scheme
  scale_color_manual(values = nutrient_colors, 
                     labels = nutrient_labels,
                     name = "") +
  
  # Apply theme
  theme_jg_bottom +
  theme(
    # Add borders around each facet
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    
    # Improve facet titles
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 1),
    
    # Place strips outside the plot area
    strip.placement = "outside",
    
    # Make phase labels more readable
    strip.text.y = element_text(angle = 270, hjust = 0.5),
    
    # Add caption and title formatting
    plot.caption = element_text(size = 9, hjust = 0),
    
    # Position legend at bottom
    legend.box.margin = margin(t = 10))

# Print the plot
print(FigS9)

# Save with appropriate dimensions
ggsave("figure/S9_FigS9_interaction_plot_per_phase.png", 
       plot = FigS9,
       width = 12, 
       height = 10, 
       units = "in", 
       dpi = 300)

ggsave("figure/S9_FigS9_interaction_plot_per_phase.pdf", 
       plot = FigS9,
       width = 12, 
       height = 10, 
       units = "in", 
       dpi = 300)


## 12.4 FIGURE S10: Within phase time effect per treatment (Modified) ####
# Function to extract and format post-hoc comparison results
format_pairs_data <- function(pairs_object, phase_name, metric_name) {
  # Convert to data frame
  pairs_df <- as.data.frame(pairs_object)
  
  # Add phase and metric columns
  pairs_df$Phase <- phase_name
  pairs_df$Metric <- metric_name
  
  # Parse the days_since_start_pairwise column to get individual days
  pairs_df <- pairs_df %>%
    tidyr::separate(days_since_start_pairwise, into = c("day1", "day2"), sep = " - ") %>%
    mutate(
      day1 = as.numeric(as.character(day1)),
      day2 = as.numeric(as.character(day2)),
      # Flag significant comparisons
      significant = p.value < 0.05,
      # Create significance indicator
      sig_level = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01 ~ "**",
        p.value < 0.05 ~ "*",
        p.value < 0.1 ~ ".",
        TRUE ~ ""
      )
    )
  
  # Create the appropriate comparison label based on phase
  if(phase_name == "Phase 1") {
    pairs_df$comparison <- factor(paste(pairs_df$day1, "-", pairs_df$day2), 
                                  levels = c("0 - 28", "0 - 63", "28 - 63"))
  } else if(phase_name == "Phase 2") {
    pairs_df$comparison <- factor(paste(pairs_df$day1, "-", pairs_df$day2), 
                                  levels = c("91 - 119", "91 - 154", "119 - 154"))
  } else if(phase_name == "Phase 3") {
    pairs_df$comparison <- factor(paste(pairs_df$day1, "-", pairs_df$day2), 
                                  levels = c("183 - 211", "183 - 238", "211 - 238"))
  } else if(phase_name == "Phase 4") {
    pairs_df$comparison <- factor(paste(pairs_df$day1, "-", pairs_df$day2), 
                                  levels = c("279 - 306", "279 - 341", "306 - 341"))
  }
  
  return(pairs_df)
}

# Total Beta Diversity
total_beta_data <- bind_rows(
  format_pairs_data(pairs_p1_tot, "Phase 1", "Total Beta"),
  format_pairs_data(pairs_p2_tot, "Phase 2", "Total Beta"),
  format_pairs_data(pairs_p3_tot, "Phase 3", "Total Beta"),
  format_pairs_data(pairs_p4_tot, "Phase 4", "Total Beta")
)

# Turnover Component
turnover_data <- bind_rows(
  format_pairs_data(pairs_p1_turn, "Phase 1", "Turnover"),
  format_pairs_data(pairs_p2_turn, "Phase 2", "Turnover"),
  format_pairs_data(pairs_p3_turn, "Phase 3", "Turnover"),
  format_pairs_data(pairs_p4_turn, "Phase 4", "Turnover")
)

# Nestedness Component
nestedness_data <- bind_rows(
  format_pairs_data(pairs_p1_nest, "Phase 1", "Nestedness"),
  format_pairs_data(pairs_p2_nest, "Phase 2", "Nestedness"),
  format_pairs_data(pairs_p3_nest, "Phase 3", "Nestedness"),
  format_pairs_data(pairs_p4_nest, "Phase 4", "Nestedness")
)

# Combine all metrics
all_combined_data <- bind_rows(total_beta_data, turnover_data, nestedness_data)

# Add cleaner month labels based on days
all_combined_data <- all_combined_data %>%
  mutate(
    month1 = case_when(
      day1 == 0 ~ "Apr/May",
      day1 == 28 ~ "May/Jun",
      day1 == 63 ~ "Jun/Jul",
      day1 == 91 ~ "Jul/Aug",
      day1 == 119 ~ "Aug/Sep",
      day1 == 154 ~ "Sep/Oct",
      day1 == 183 ~ "Oct/Nov",
      day1 == 211 ~ "Nov/Dec",
      day1 == 238 ~ "Dec/Jan",
      day1 == 279 ~ "Jan/Feb",
      day1 == 306 ~ "Feb/Mar",
      day1 == 341 ~ "Mar/Apr",
      TRUE ~ as.character(day1)
    ),
    month2 = case_when(
      day2 == 0 ~ "Apr/May",
      day2 == 28 ~ "May/Jun",
      day2 == 63 ~ "Jun/Jul",
      day2 == 91 ~ "Jul/Aug",
      day2 == 119 ~ "Aug/Sep",
      day2 == 154 ~ "Sep/Oct",
      day2 == 183 ~ "Oct/Nov",
      day2 == 211 ~ "Nov/Dec",
      day2 == 238 ~ "Dec/Jan",
      day2 == 279 ~ "Jan/Feb",
      day2 == 306 ~ "Feb/Mar",
      day2 == 341 ~ "Mar/Apr",
      TRUE ~ as.character(day2)
    ),
    # Create a new comparison label with month names
    month_comparison = paste(month1, "-", month2)
  )


all_combined_data <- all_combined_data %>% 
  mutate(
    treatment = paste0(
      ifelse(predator == 1, "Pred+", "Pred-"),
      "/",
      ifelse(nutrients == 1, "Nutr+", "Nutr-")
    )
  )


# Ensure treatments are in the correct order
all_combined_data$treatment <- factor(all_combined_data$treatment, 
                                      levels = c("Pred+/Nutr+", "Pred-/Nutr+", 
                                                 "Pred+/Nutr-", "Pred-/Nutr-"))

# Ensure metrics are in correct order
all_combined_data$Metric <- factor(all_combined_data$Metric, 
                                   levels = c("Total Beta", "Turnover", "Nestedness"))

# Fix the month comparison ordering - make sure it's in chronological order
all_combined_data$month_comparison <- factor(
  all_combined_data$month_comparison,
  levels = c(
    "Apr/May - May/Jun", "Apr/May - Jun/Jul", "May/Jun - Jun/Jul",
    "Jul/Aug - Aug/Sep", "Jul/Aug - Sep/Oct", "Aug/Sep - Sep/Oct",
    "Oct/Nov - Nov/Dec", "Oct/Nov - Dec/Jan", "Nov/Dec - Dec/Jan",
    "Jan/Feb - Feb/Mar", "Jan/Feb - Mar/Apr", "Feb/Mar - Mar/Apr"
  )
)

all_combined_data$Phase <- factor(all_combined_data$Phase, 
                                  levels = c("Phase 1", "Phase 2", "Phase 3", "Phase 4"), 
                                  labels = c("Period 1", "Period 2", "Period 3", "Period 4"))

# Set consistent color palette
effect_colors <- c("Positive" = "#4daf4a", "Negative" = "#e41a1c")

# Create the dot plot - fixing the geom_rect issue
FigS10 <- ggplot(all_combined_data, aes(x = month_comparison, y = treatment)) +
  # Add circles for effect sizes
  geom_point(aes(size = abs(estimate), 
                 color = ifelse(estimate > 0, "Positive", "Negative"),
                 alpha = significant)) +
  
  # Add significance indicators
  geom_text(aes(label = sig_level), vjust = -0.7, size = 4) +
  
  # Set colors for positive/negative effects
  scale_color_manual(values = effect_colors, name = "Direction") +
  
  # Set transparency for significance
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.3), name = "Significant") +
  
  # Set sizes for dots based on effect size
  scale_size_continuous(range = c(1, 8), name = "Effect Size \n(|log odds|)") +
  
  # Facet by metric and phase
  facet_grid(Metric ~ Phase, scales = "free_x", space = "free_x") +
  
  # Add labels
  labs(
    title = "",
    x = "",
    y = ""
  ) +
  
  # Use your custom theme
  theme_jg_right +
  theme(
    strip.background = element_rect(fill = "gray90", color = NA),
    strip.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.box = "vertical",
    panel.spacing.y = unit(0.8, "cm"), 
    legend.margin = margin(t = 15, r = 10, b = 15, l = 10)
  )

# Display the plot
print(FigS10)


ggsave("figure/S10_FigS10_dotplot_per_phase_post-hoc.png", FigS10, width = 8.5, height = 7.6, units = "in", dpi = 600)
ggsave("figure/S10_FigS10_dotplot_per_phase_post-hoc.pdf", FigS10, width = 8.2, height = 9, units = "in", dpi = 600)



## 12.5 FIGURE 4: Period 2 detailes analysis ####
### 12.5.1 Interaction plot period 2 ####
# Filter emmeans data for Period 2 only
emm_p2 <- all_emm %>%
  filter(Phase == "P2") %>%
  mutate(
    predator = factor(predator, levels = c("Pred-", "Pred+"))
  )

# Filter significance data for Period 2
sig_p2 <- sig_data %>%
  filter(Phase == "P2")

# Create Panel A
Fig4_panelA <- ggplot(emm_p2, 
                      aes(x = predator, y = emmean, 
                          group = nutrients, color = nutrients)) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), 
                width = 0.2, linewidth = 0.8) +
  geom_text(data = sig_p2,
            aes(x = x_pos, y = y_pos, label = sig_text, group = NULL),
            size = 3, color = "black", vjust = sig_p2$vjust, hjust = 0) +
  facet_wrap(~ Component, scales = "free_y", nrow = 1) +
  scale_color_manual(values = nutrient_colors,
                     labels = nutrient_labels,
                     name = "") +
  labs(x = "",
       y = "Estimated marginal mean\n(logit scale)",
       tag = "A") +
  theme_jg_bottom +
  theme(
    panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5),
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    legend.position = "right"
  )

print(Fig4_panelA)

### 12.5.2 Dot plot for Period 2 ####

# Filter dot plot data for Period 2 only
dotplot_p2 <- all_combined_data %>%
  filter(Phase == "Period 2")

# Define colour palette for direction
effect_colors <- c("Positive" = "#4daf4a", "Negative" = "#e41a1c")

# Create Panel B
Fig4_panelB <- ggplot(dotplot_p2, 
                      aes(x = month_comparison, y = treatment)) +
  geom_point(aes(size = abs(estimate),
                 color = ifelse(estimate > 0, "Positive", "Negative"),
                 alpha = significant)) +
  geom_text(aes(label = sig_level), vjust = -0.8, size = 3.5) +
  scale_color_manual(values = effect_colors, name = "Direction") +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.3),
                     name = "Significant") +
  scale_size_continuous(range = c(1, 8),
                        name = "Effect size\n(|log odds|)") +
  facet_wrap(~ Metric, scales = "free_x", nrow = 1) +
  labs(x = "",
       y = "",
       tag = "B") +
  theme_jg_right +
  theme(
    panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5),
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    legend.box = "vertical"
  )

print(Fig4_panelB)

### 12.5.3 FIGURE 4: Assemble Panel A and B ####

Fig4 <- Fig4_panelA / Fig4_panelB +
  plot_layout(heights = c(1.2, 1.5))

print(Fig4)

ggsave("figure/R4_Fig4_Period2_detailed.png", Fig4, 
       width = 10, height = 9, dpi = 300)
ggsave("figure/R4_Fig4_Period2_detailed.pdf", Fig4, 
       width = 10, height = 9, dpi = 300)
