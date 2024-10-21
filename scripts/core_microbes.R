### Setup env
library(tidyverse)

source('../MFD_colors.R')

setwd("/mfd_core")


### Import data
## Load data and create "complexes" correponding to full ontology strings
metadata <- data.table::fread('../section2/output/2024-05-10_samples-grid-10km.csv', na.strings = "") %>%
  mutate(mfd_areatype = str_c(mfd_sampletype, mfd_areatype, sep = ", "),
         complex1 = str_c(mfd_areatype, mfd_hab1, sep = ", "),
         complex2 = str_c(complex1, mfd_hab2, sep = ", "),
         complex3 = str_c(complex2, mfd_hab3, sep = ", ")) %>%
  group_by(complex1) %>%
  mutate(complex1_size = n()) %>%
  group_by(complex2) %>%
  mutate(complex2_size = n()) %>%
  group_by(complex3) %>%
  mutate(complex3_size = n()) %>%
  mutate(label1 = str_c(complex1, ", n = ", complex1_size, sep = ""),
         label2 = str_c(complex2, ", n = ", complex2_size, sep = ""),
         label3 = str_c(complex3, ", n = ", complex3_size, sep = "")) %>%
  ungroup() %>%
  filter(!complex1 %in% c("Water, Subterranean, Freshwater",
                         "Water, Urban, Sandfilter",
                         "Soil, Urban, Other",
                         "Sediment, Urban, Other",
                         "Soil, Urban, Roadside",
                         "Soil, Subterranean, Urban")) %>%
  mutate(across(mfd_sampletype, ~factor(., levels = names(sampletype.palette))),
         across(mfd_areatype, ~factor(., levels = names(areatype.palette))),
         across(complex1, ~factor(., levels = names(mfdo1.palette))),
         across(complex2:complex3, ~factor(.))) %>%
  filter(!is.na(complex1)) %>%
  select(-project_id)

## Pull sampel IDs
samples.subset <- metadata %>%
  pull(fieldsample_barcode)

## Load genus-aggregated observational table of 16S derived from the metagenomes 
data <- data.table::fread('../format_data/output/2024-03-07_arcbac-rarefaction-rel.csv') %>%
  select(any_of(samples.subset), Kingdom:Genus) %>%
  filter(rowSums(across(where(is.numeric)))!=0)

## Create summaries across the ontology levels
summary.type <- metadata %>%
  group_by(mfd_sampletype) %>%
  summarise(group_size = n(), .groups = "drop") %>%
  mutate(`Core` = ceiling(group_size*0.5)) %>%
  filter(!is.na(mfd_sampletype))

summary.area <- metadata %>%
  group_by(mfd_areatype) %>%
  summarise(group_size = n(), .groups = "drop") %>%
  mutate(`Core` = ceiling(group_size*0.5)) %>%
  filter(!is.na(mfd_areatype))

summary.mfdo1 <- metadata %>%
  group_by(complex1) %>%
  summarise(group_size = n(), .groups = "drop") %>%
  mutate(`Core` = ceiling(group_size*0.5)) %>%
  filter(!is.na(complex1))

summary.mfdo2 <- metadata %>%
  group_by(complex2) %>%
  summarise(group_size = n(), .groups = "drop") %>%
  mutate(`Core` = ceiling(group_size*0.5)) %>%
  filter(!is.na(complex2))

summary.mfdo3 <- metadata %>%
  group_by(complex3) %>%
  summarise(group_size = n(), .groups = "drop") %>%
  mutate(`Core` = ceiling(group_size*0.5)) %>%
  filter(!is.na(complex3))


### Core community members 
## Write function to identify core community members
core.by.group <- function(data, metadata, summary) {
  data.out <- list()
  
  levels <- summary %>%
    filter(group_size >= 6) %>%
    pull(1) %>%
    droplevels() %>%
    levels()
  
  for (i in 1:length(levels)) {
    filter <- metadata %>%
      select(colnames(summary[1]), fieldsample_barcode) %>%
      filter(if_any(1, ~ . == levels[i])) %>%
      pull(fieldsample_barcode)
    
    data.filt <- data %>%
      select(Genus, any_of(filter)) %>%
      filter(!Genus == "Unclassified") %>%
      filter(rowSums(across(where(is.numeric)))!=0)
    
    data.long <- data.filt %>%
      pivot_longer(!Genus, names_to = "fieldsample_barcode", values_to = "abundance")
    
    df <- data.long %>%
      group_by(Genus, fieldsample_barcode) %>%
      select(Genus, fieldsample_barcode, abundance) %>%
      ungroup() %>%
      group_by(Genus) %>%
      summarise(n_observations = sum(abundance > 0),
                n_abundant = sum(abundance >= 0.1),
                median_abundance = median(abundance),
                mean_abundance = mean(abundance),
                sd_abundance = sd(abundance)) %>%
      ungroup() %>%
      arrange(desc(mean_abundance)) %>%
      as_tibble() # triggers dtplyr
    
    df.core <- df %>%
      filter(n_abundant >= ceiling(length(filter)*0.5)) %>%
      group_by(Genus) %>%
      arrange(desc(mean_abundance)) %>%
      mutate(Genus_type = "Core") %>%
      as_tibble()
    
    genus_core <- rbind(df.core) %>%
      mutate(group = levels[[i]])
    
    data.out[[levels[i]]] <- genus_core
  }
  return(data.out)
}

## Run core by ontology level
core.community.type <- core.by.group(data, metadata, summary.type)
core.community.area <- core.by.group(data, metadata, summary.area)
core.community.mfdo1 <- core.by.group(data, metadata, summary.mfdo1)
core.community.mfdo2 <- core.by.group(data, metadata, summary.mfdo2)
core.community.mfdo3 <- core.by.group(data, metadata, summary.mfdo3)

## Bind results into single dataframes for each ontology level
df.core.type <- bind_rows(core.community.type) %>%
  rename(mfd_sampletype = group) %>%
  left_join(summary.type) %>%
  mutate(prevalence = n_observations/group_size) %>%
  group_by(Genus) %>%
  mutate(fidelity = n()) %>%
  ungroup()

df.core.area <- bind_rows(core.community.area) %>%
  rename(mfd_areatype = group) %>%
  left_join(summary.area) %>%
  mutate(prevalence = n_observations/group_size) %>%
  group_by(Genus) %>%
  mutate(fidelity = n()) %>%
  ungroup()

df.core.mfdo1 <- bind_rows(core.community.mfdo1) %>%
  rename(complex1 = group) %>%
  left_join(summary.mfdo1) %>%
  mutate(prevalence = n_observations/group_size) %>%
  group_by(Genus) %>%
  mutate(fidelity = n()) %>%
  ungroup()

df.core.mfdo2 <- bind_rows(core.community.mfdo2) %>%
  rename(complex2 = group) %>%
  left_join(summary.mfdo2) %>%
  mutate(prevalence = n_observations/group_size) %>%
  group_by(Genus) %>%
  mutate(fidelity = n()) %>%
  ungroup()

df.core.mfdo3 <- bind_rows(core.community.mfdo3) %>%
  rename(complex3 = group) %>%
  left_join(summary.mfdo3) %>%
  mutate(prevalence = n_observations/group_size) %>%
  group_by(Genus) %>%
  mutate(fidelity = n()) %>%
  ungroup()

## Evaluate the number of core genera acroos ontology levels
df.core.type %>%
  pull(Genus) %>% 
  unique %>%
  length()

df.core.area %>%
  pull(Genus) %>% 
  unique %>%
  length()

df.core.mfdo1 %>%
  pull(Genus) %>% 
  unique %>%
  length()

df.core.mfdo2 %>%
  pull(Genus) %>% 
  unique %>%
  length()

df.core.mfdo3 %>%
  pull(Genus) %>% 
  unique %>%
  length()

## Write results to output
data.table::fwrite(df.core.type, "output/2024-03-07_core-genera-type.csv")
data.table::fwrite(df.core.area, "output/2024-03-07_core-genera-area.csv")
data.table::fwrite(df.core.mfdo1, "output/2024-03-07_core-genera-mfdo1.csv")
data.table::fwrite(df.core.mfdo2, "output/2024-03-07_core-genera-mfdo2.csv")
data.table::fwrite(df.core.mfdo3, "output/2024-03-07_core-genera-mfdo3.csv")

## Same image
save.image(file = 'output/core_microbes.RData')
