### Setup env
library(tidyverse)

setwd("/mfd_core")


## Improt data
load("output/2025-04-23_core_microbes.RData")

## Mutate ontology level specific dataframes
df.core.type <- df.core.type %>%
  rename(complex = mfd_sampletype) %>%
  mutate(complexity = "mfd_sampletype")
df.core.area <- df.core.area %>%
  rename(complex = mfd_areatype) %>%
  mutate(complexity = "mfd_areatype")
df.core.mfdo1 <- df.core.mfdo1 %>%
  rename(complex = complex1) %>%
  mutate(complexity = "mfd_hab1") %>%
  mutate(across(complex, ~str_replace(., "Bogs, mires and fens", "Bogs mires and fens")),
         across(complex, ~str_replace(., "Poales, ", "Poales ")))
df.core.mfdo2 <- df.core.mfdo2 %>%
  rename(complex = complex2) %>%
  mutate(complexity = "mfd_hab2") %>%
  mutate(across(complex, ~str_replace(., "Bogs, mires and fens", "Bogs mires and fens")),
         across(complex, ~str_replace(., "Poales, ", "Poales ")))
df.core.mfdo3 <- df.core.mfdo3 %>%
  rename(complex = complex3) %>%
  mutate(complexity = "mfd_hab3") %>%
  mutate(across(complex, ~str_replace(., "Bogs, mires and fens", "Bogs mires and fens")),
         across(complex, ~str_replace(., "Poales, ", "Poales ")))

## Get number of distinct genera per ontology level
df.core.type %>%
  pull(Genus) %>%
  n_distinct()

df.core.area %>%
  pull(Genus) %>%
  n_distinct()

df.core.mfdo1 %>%
  pull(Genus) %>%
  n_distinct()

df.core.mfdo2 %>%
  pull(Genus) %>%
  n_distinct()

df.core.mfdo3 %>%
  pull(Genus) %>%
  n_distinct()

## Set levels
levels <- c("mfd_sampletype", "mfd_areatype", "mfd_hab1", "mfd_hab2", "mfd_hab3")

## Combine results
df <- rbind(df.core.type, df.core.area, df.core.mfdo1, df.core.mfdo2, df.core.mfdo3) %>%
  mutate(across(complexity, ~factor(., levels = levels)))

## Evaluate proportion of de novo genera
df.taxonomy <- df %>%
  select(Genus) %>%
  distinct() %>%
  reframe(total = n(),
          de_novo = sum(str_count(Genus, "MFD_g_")),
          do_novo_pct = round(de_novo/total*100, 1),
          known = total-de_novo,
          known_pct = round(known/total*100, 1))

## Create long format
data.long <- data %>%
  filter(!Genus == "Unclassified") %>%
  select(Genus, starts_with("MFD")) %>%
  pivot_longer(!Genus, names_to = "fieldsample_barcode", values_to = "abundance")

## Evaluate total number of observed genera across ontology
genera.sum <- data.long %>%
  filter(!abundance == 0) %>%
  left_join(metadata %>% select(fieldsample_barcode, mfd_sampletype:mfd_hab3)) %>%
  mutate(mfd_hab1 = str_c(mfd_areatype, mfd_hab1, sep = ", "),
         mfd_hab2 = str_c(mfd_hab1, mfd_hab2, sep = ", "),
         mfd_hab3 = str_c(mfd_hab2, mfd_hab3, sep = ", ")) %>%
select(-c(fieldsample_barcode, abundance)) %>%
  pivot_longer(cols = mfd_sampletype:mfd_hab3, names_to = "complexity", values_to = "complex") %>%
  mutate(across(complex, ~str_replace(., "Bogs, mires and fens", "Bogs mires and fens")),
         across(complex, ~str_replace(., "Poales, ", "Poales ")),
         across(complexity, ~factor(., levels = levels))) %>%
  filter(!is.na(complex)) %>%
  group_by(complexity, complex) %>%
  reframe(genera_total = n_distinct(Genus))

## Create group summaries
df.sum <- df %>%
  mutate(denovo_abundance = case_when(str_detect(Genus, "MFD_g_") ~ mean_abundance,
                                      TRUE~0),
         uniq_abundance = case_when(fidelity == 1 ~mean_abundance,
                                    TRUE~0)) %>%
  group_by(complex, complexity, group_size) %>%
  reframe(core_size = n(),
          uniq_core = sum(fidelity == 1),
          core_abundance = round(sum(mean_abundance), 1),
          uniq_abundance = round(sum(uniq_abundance), 1),
          denovo_abundance = round(sum(denovo_abundance), 1)) %>%
  rename(samples = group_size) %>%
  left_join(genera.sum) %>%
  relocate(genera_total, .after = "samples")

## Agri MFDO2 crop type
core.agri.crop <- df.sum %>%
  filter(complexity %in% c("mfd_hab1", "mfd_hab2"),
         str_detect(complex, "Agriculture,")) %>%
  group_by(complexity) %>%
  reframe(median = median(core_size),
          mean = mean(core_size),
          sd = sd(core_size))

## MFDO1 soil specific 
df.sum.soil.mfdo1 <- df.sum %>%
  filter(complexity == "mfd_hab1",
         str_detect(complex, "Soil")) %>%
  group_by(complexity) %>%
  reframe(sum_core = sum(core_size),
          sum_uniq = sum(uniq_core),
          median = round(median(uniq_core), 0), 
          mean = round(mean(uniq_core), 0),
          sd = round(sd(uniq_core), 0))

## Summary of smmaries
df.sum.sum <- df.sum %>%
  group_by(complexity) %>%
  reframe(size_core = sum(core_size),
          min_core = min(core_size),
          max_core = max(core_size),
          median_core = round(median(core_size), 0),
          mean_core = round(mean(core_size), 0),
          sd_core = round(sd(core_size), 0),
          min_abuncance = round(min(core_abundance), 1),
          max_abuncance = round(max(core_abundance), 1),
          median_abuncance = round(median(core_abundance), 1),
          mean_abuncance = round(mean(core_abundance), 1),
          sd_abuncance = round(sd(core_abundance), 1))

## Top genera
top.genera <- df %>%
  filter(complexity == "mfd_hab1") %>%
  group_by(complexity, Genus) %>%
  reframe(prevalence = n(), 
          # sum = round(sum(mean_abundance), 1),
          median = round(median(mean_abundance), 1),
          mean = round(mean(mean_abundance), 1),
          sd = round(sd(mean_abundance), 1))

## Top 10 prevalent
top.prevalence <- top.genera %>% 
  arrange(desc(prevalence)) %>%
  slice_head(n = 10)

## Top 10 abundant
top.abundance <- top.genera %>% 
  arrange(desc(mean)) %>%
  slice_head(n = 10)

## Combine
tmp.comb <- rbind(top.prevalence, top.abundance) %>%
  distinct() %>%
  arrange(desc(prevalence))

## Ontology summaries
df.sum.type <- df.core.type %>%
  group_by(complex) %>%
  rename(sample_type = complex) %>%
  reframe(core_type = n())

df.sum.area <- df.core.area %>%
  group_by(complex) %>%
  reframe(core_area = n()) %>%
  rename(area_type = complex) %>%
  mutate(sample_type = str_extract(area_type, "^([^,]*,){0}[^,]*"), .before = area_type)

df.sum.mfdo1 <- df.core.mfdo1 %>%
  group_by(complex) %>%
  reframe(core_mfdo1 = n()) %>%
  rename(mfd_hab1 = complex) %>%
  mutate(area_type = str_extract(mfd_hab1, "^([^,]*,){1}[^,]*"), .before = mfd_hab1)

df.sum.mfdo2 <- df.core.mfdo2 %>%
  group_by(complex) %>%
  reframe(core_mfdo2 = n()) %>%
  rename(mfd_hab2 = complex) %>%
  mutate(mfd_hab1 = str_extract(mfd_hab2, "^([^,]*,){2}[^,]*"), .before = mfd_hab2)

df.sum.mfdo3 <- df.core.mfdo3 %>%
  group_by(complex) %>%
  reframe(core_mfdo3 = n()) %>%
  rename(mfd_hab3 = complex) %>%
  mutate(mfd_hab2 = str_extract(mfd_hab3, "^([^,]*,){3}[^,]*"), .before = mfd_hab3)

## Combine
df.sum.combined <- df.sum.type %>%
  left_join(df.sum.area) %>%
  left_join(df.sum.mfdo1) %>%
  left_join(df.sum.mfdo2) %>%
  left_join(df.sum.mfdo3)

## Write results to output
data.table::fwrite(df, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_core_members.csv"))
data.table::fwrite(df.sum, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_core_summary_long.csv"))
data.table::fwrite(df.sum.combined, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_core_summary_wide.csv"))
















