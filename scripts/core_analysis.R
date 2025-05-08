### Setup env
library(tidyverse)
library(UpSetR)
library(ComplexUpset)
library(patchwork)

source('scripts/MFD_colors.R')

setwd("/mfd_core")

metadata <- data.table::fread("output/2025-04-23_MFD_core_metadata.tsv")

## Pull sample IDs
samples.subset <- metadata %>%
  pull(fieldsample_barcode)

## Load genus-aggregated observational table of 16S derived from the metagenomes 
data <- data.table::fread('data/2025-02-13_MFD_arcbac_genus_rarefaction_rel.csv') %>%
  select(any_of(samples.subset), Kingdom:Genus) %>%
  filter(rowSums(across(where(is.numeric)))!=0)

tax <- data %>%
  select(Kingdom:Genus)

df.core.mfdo1 <- data.table::fread("output/2025-04-23_MFD_core_genera_mfdo1.tsv") %>%
  rename(complex = complex1) %>%
  mutate(across(complex, ~factor(., levels = mfdo1.levels)))

gc()


### Core abundance
# Wide format
core.genera <- df.core.mfdo1 %>%
  select(Genus) %>%
  distinct() %>%
  pull(Genus)

df.core.mfdo1 %>%
  group_by(complex) %>%
  reframe(count = n()) %>%
  arrange(count)

df.core.long <- df.core.mfdo1 %>%
  select(Genus, complex, Genus_type) %>%
  complete(Genus, complex) %>%
  mutate(res = if_else(is.na(Genus_type), FALSE, TRUE)) %>%
  select(-Genus_type)

data.long <- data %>%
  select(Genus, starts_with("MFD")) %>%
  pivot_longer(!Genus, names_to = "fieldsample_barcode", values_to = "rel.abund") %>%
  filter(!rel.abund == 0) %>%
  left_join(tax) %>%
  select(Kingdom:Family, Genus, rel.abund, fieldsample_barcode)

data.long.complex <- data %>%
  select(Genus, starts_with("MFD")) %>%
  pivot_longer(!Genus, names_to = "fieldsample_barcode", values_to = "rel.abund") %>%
  left_join(metadata) %>%
  group_by(Genus, complex) %>%
  summarise(across(rel.abund, ~mean(.))) %>%
  ungroup() %>%
  filter(!rel.abund == 0) %>%
  left_join(tax) %>%
  select(Kingdom:Family, Genus, rel.abund, complex)

gc()

# data.mfdo1 <- data.long.complex %>%
#   group_by(Genus) %>%
#   pivot_wider(names_from = "complex", values_from = "rel.abund")

reads.summary <- data.long.complex %>%
  left_join(df.core.mfdo1 %>% select(complex, Genus, Genus_type)) %>%
  mutate(across(Genus_type, ~if_else(Genus == "Unclassified", Genus, Genus_type)),
         across(Genus_type, ~replace_na(., "Non-core"))) %>%
  group_by(complex, Genus_type) %>%
  reframe(RA = sum(rel.abund))

reads.summary %>%
  filter(!Genus_type == "Unclassified") %>%
  group_by(complex, Genus_type) %>%
  summarise(Sum_RA = sum(RA)) %>%
  ungroup() %>%
  rstatix::shapiro_test(Sum_RA)

reads.summary %>%
  filter(!Genus_type == "Unclassified") %>%
  summarise(Sum_RA = sum(RA)) %>%
  ungroup() %>%
  summarise(median_RA = median(Sum_RA),
            mean_RA = mean(Sum_RA),
            sd_RA = sd(Sum_RA))

summary.mfdo1 <- metadata %>%
  group_by(complex) %>%
  summarise(group_size = n(), .groups = "drop") %>%
  mutate(`Core` = ceiling(group_size*0.5)) %>%
  filter(!is.na(complex))

hab.levels <- summary.mfdo1 %>%
  select(complex) %>%
  droplevels() %>%
  pull(complex)

mfdo1.palette.sub <- mfdo1.palette[hab.levels]

plot.sum <- reads.summary %>%
  mutate(across(complex, ~factor(., levels = mfdo1.levels))) %>%
  mutate(across(Genus_type, ~factor(., levels = rev(c("Core", "Non-core", "Unclassified"))))) %>%
  ggplot(aes(x = complex, 
             y = RA, 
             fill = complex,
             alpha = Genus_type)) +
  geom_bar(color = "black",
           stat = "identity", 
           width = 0.6) +
  guides(alpha = guide_legend(title = "Genus type", override.aes = list(shape = 22, size = 8)),
         fill = "none") + 
  scale_fill_manual(values = mfdo1.palette.sub) +
  coord_flip() +
  scale_alpha_discrete(range=c(0,1))+
  scale_y_continuous(position = "left") +
  scale_x_discrete(lim = rev) +
  #ggtitle("Classified 16S rRNA gene fragments",
  #        subtitle = "") +
  labs(x = "",
       y = "Relative abundance (%)") +
  theme_minimal() +
  theme(axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_blank())

plot.sum

### Core upset plot
# Upset
upset.df <- df.core.long %>%
  group_by(Genus) %>%
  pivot_wider(names_from = "complex", values_from = "res") %>%
  ungroup() %>%
  select(Genus, any_of(hab.levels))

sets = colnames(upset.df[hab.levels])

plot.upset <- upset(
  upset.df,
  sets,
  set_sizes = upset_set_size(geom = geom_bar(width = 0.6, color = "black"), position = "right") +
    theme(axis.title.x = element_text(size = 12)),
  name = "",
  guides = 'over',
  sort_sets = FALSE,
  sort_intersections = 'ascending',
  # intersections = 'all',
  min_size = 2,
  sort_intersections_by=c('degree', 'cardinality'),
  queries = list(
    upset_query(set = (names(mfdo1.palette.sub[1])), fill = mfdo1.palette.sub[1]),
    upset_query(set = (names(mfdo1.palette.sub[2])), fill = mfdo1.palette.sub[2]),
    upset_query(set = (names(mfdo1.palette.sub[3])), fill = mfdo1.palette.sub[3]),
    upset_query(set = (names(mfdo1.palette.sub[4])), fill = mfdo1.palette.sub[4]),
    upset_query(set = (names(mfdo1.palette.sub[5])), fill = mfdo1.palette.sub[5]),
    upset_query(set = (names(mfdo1.palette.sub[6])), fill = mfdo1.palette.sub[6]),
    upset_query(set = (names(mfdo1.palette.sub[7])), fill = mfdo1.palette.sub[7]),
    upset_query(set = (names(mfdo1.palette.sub[8])), fill = mfdo1.palette.sub[8]),
    upset_query(set = (names(mfdo1.palette.sub[9])), fill = mfdo1.palette.sub[9]),
    upset_query(set = (names(mfdo1.palette.sub[10])), fill = mfdo1.palette.sub[10]),
    upset_query(set = (names(mfdo1.palette.sub[11])), fill = mfdo1.palette.sub[11]),
    upset_query(set = (names(mfdo1.palette.sub[12])), fill = mfdo1.palette.sub[12]),
    upset_query(set = (names(mfdo1.palette.sub[13])), fill = mfdo1.palette.sub[13]),
    upset_query(set = (names(mfdo1.palette.sub[14])), fill = mfdo1.palette.sub[14]),
    upset_query(set = (names(mfdo1.palette.sub[15])), fill = mfdo1.palette.sub[15]),
    upset_query(set = (names(mfdo1.palette.sub[16])), fill = mfdo1.palette.sub[16]),
    upset_query(set = (names(mfdo1.palette.sub[17])), fill = mfdo1.palette.sub[17]),
    upset_query(set = (names(mfdo1.palette.sub[18])), fill = mfdo1.palette.sub[18]),
    upset_query(set = (names(mfdo1.palette.sub[19])), fill = mfdo1.palette.sub[19])),
  matrix = (intersection_matrix(geom = geom_point(shape = 'circle filled', size = 3)) +
              scale_y_discrete(position = "left")),
  base_annotations=list(
    'Intersection size'=intersection_size(
      counts=FALSE,
      color = "black",
      text=list(
        hjust = -1,
        vjust = 0.5,
        angle = 90),
      text_colors = c(on_background = 'black', on_bar = 'black')) +
      ylab('Number of genera') + 
      scale_y_continuous(breaks = seq(0, 40, 10), limits = c(0, 40)) +
      theme(axis.title.y = element_text(size = 12))))

plot.upset 


p.tmp <- plot.upset[[3]] + theme(axis.text.x = element_text(angle = 90, hjust = 1))

plot.upset[[1]] + p.tmp + plot_layout(nrow = 2, heights = c(1,3))


# Custom intersections
intersections <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19",
                   "2-4", "3-5", "4-5", "4-18", "7-17", "13-15", "16-17", "3-5-8", "6-8-9", "7-16-17", "13-14-15", "7-9-16-17",
                   "10-11-14-15", "10-11-13-14-15", "6-7-9-12-16-17", "10-11-12-13-14-15", "3-5-6-7-8-9-12-16-17",
                   "6-7-9-10-11-12-13-14-16-17", "6-7-8-10-11-12-13-14-15-16", "6-7-8-10-11-12-13-14-15-16-17")

# Customize upset plot
data.upset1 <- plot.upset[[1]][["data"]] %>%
  filter(intersection %in% intersections) %>%
  # filter(intersection %in% seq(1, 18, 1)) %>%
  select(intersection, in_exclusive_intersection) %>%
  mutate(across(intersection, ~factor(., levels = intersections))) %>%
  complete(., intersection, fill=list(in_exclusive_intersection = 0)) %>%
  group_by(intersection) %>%
  summarise(sum = sum(in_exclusive_intersection))
  # mutate(complex = sets,
  #        across(complex, ~factor(., levels = mfdo1.levels)))

# Bar plot
p1 <- data.upset1 %>%
  ggplot(aes(x = intersection, y = sum)) +
  geom_col(color = "black", position = "stack") +
  xlab("") +
  ylab("Number of genera") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title = element_text(size = 14))

p1

# Intersection matrix
data.upset2 <- plot.upset[[3]][["data"]] %>%
  filter(intersection %in% intersections) %>%
  filter(value == "TRUE") %>%
  mutate(across(intersection, ~factor(., levels = intersections))) %>%
  mutate(across(group, ~factor(., levels = seq(1, 19, 1)))) %>%
  complete(., group, intersection) %>%
  mutate(across(group:intersection, ~as.character(.))) %>%
  mutate(value = replace(value, is.na(value) & intersection == group, TRUE)) %>%
  filter(!is.na(value)) %>%
  mutate(across(intersection, ~factor(., levels = intersections))) %>%
  mutate(across(group, ~factor(., levels = seq(1, 19, 1)))) %>%
  mutate(start = str_extract(intersection, "^."),
         end = str_extract(intersection, "[^-]+$"))

data.upset2.1 <- plot.upset[[4]][["data"]] %>%
  filter(intersection %in% intersections) %>%
  mutate(across(intersection, ~factor(., levels = intersections))) %>%
  mutate(across(group, ~factor(., levels = seq(1, 19, 1)))) %>%
  complete(., group, intersection)

p2 <- data.upset2 %>%
  ggplot(aes(x = intersection, y = group)) +
  ggstats::geom_stripped_rows(odd = "grey95", even = "white") +
  geom_point(data = data.upset2.1, fill = "white", color = "grey70", pch = 21, size = 4) +
  geom_segment(aes(x = intersection, xend = intersection, y = group, yend = end)) +
  geom_point(aes(fill = group), pch = 21, size = 4) +
  scale_fill_manual(values = unname(mfdo1.palette.sub)) +
  scale_y_discrete(limits=rev, labels = rev(names(mfdo1.palette.sub)), position = "left") +
  xlab("") +
  ylab("") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12),
        legend.position = "none")

p2

# Set sizes
p3 <- plot.upset[[4]] + 
  scale_x_discrete(limits=rev) +
  xlab("") +
  ylab("Number of core genera") +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14))

p3


### Managed vs. Natural
metadata.sub <- metadata %>%
  filter(str_detect(complex, "Soil"),
         !str_detect(complex, "Subterranean"))

samples.sub <- metadata.sub %>%
  pull(fieldsample_barcode)

df.core.mfdo1.sub <- df.core.mfdo1 %>%
  filter(str_detect(complex, "Soil")) %>%
  mutate(manage = case_when(str_detect(complex, "Natural") ~ "Natural",
                            TRUE ~ "Disturbed"))

# data.table::fwrite(tmp, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE,
#                    paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_core_managed_natural.csv"))

df.disturbed <- df.core.mfdo1.sub %>%
  filter(manage == "Disturbed")

df.natural <- df.core.mfdo1.sub %>%
  filter(manage == "Natural")

df.core.mfdo1.sub %>%
  pull(Genus) %>%
  unique() %>%
  length()

df.disturbed %>%
  pull(Genus) %>%
  unique() %>%
  length()

df.natural %>%
  pull(Genus) %>%
  unique() %>%
  length()

setdiff(str_c(df.disturbed$Family, df.disturbed$Genus, sep = " - "), str_c(df.natural$Family, df.natural$Genus, sep = " - "))

setdiff(str_c(df.natural$Family, df.natural$Genus, sep = " - "), str_c(df.disturbed$Family, df.disturbed$Genus, sep = " - "))

intersect(str_c(df.disturbed$Family, df.disturbed$Genus, sep = " - "), str_c(df.natural$Family, df.natural$Genus, sep = " - "))

genera.int <- c(setdiff(df.disturbed$Genus, df.natural$Genus), setdiff(df.natural$Genus, df.disturbed$Genus)) %>%
  unique()

df.core.mfdo1.sum <- df.core.mfdo1.sub %>%
  select(Genus, manage, median_abundance, mean_abundance, sd_abundance) %>%
  filter(Genus %in% genera.int) %>%
  group_by(Genus) %>%
  reframe(mean = mean(mean_abundance),
          sd = mean(sd_abundance)) %>%
  ungroup()

## Evaluate total number of observed genera across ontology
# data.long.sub <- data.long %>%
#   filter(Genus %in% genera.int)

data.long.sub <- data.long %>%
  filter(fieldsample_barcode %in% samples.sub) %>%
  filter(Genus %in% c("MFD_g_198", "MFD_g_4907", "Nitrospira", "mle1-7", "Ellin6067"))

genera.sub.sum <- data.long.sub %>%
  left_join(metadata.sub %>% select(fieldsample_barcode, complex)) %>%
  mutate(Tax = str_c(Family, ",\n", Genus)) %>%
  group_by(Genus, Tax, complex) %>%
  reframe(n_observations = sum(rel.abund > 0),
          n_abundant = sum(rel.abund >= 0.1),
          median_abundance = median(rel.abund),
          mean_abundance = mean(rel.abund),
          sd_abundance = sd(rel.abund)) %>%
  left_join(df.core.mfdo1.sub %>% select(complex, Genus, Genus_type)) %>%
  mutate(across(Genus_type, ~replace_na(., "Non-core"))) %>%
  left_join(summary.mfdo1) %>%
  group_by(Genus, complex) %>%
  mutate(prevalence = n_observations/group_size) %>%
  mutate(manage = case_when(str_detect(complex, "Natural") ~ "Natural",
                            TRUE ~ "Disturbed"))

p4 <- genera.sub.sum %>%
  left_join(distinct(metadata.sub %>% select(complex, mfd_areatype, mfd_hab1))) %>%
  mutate(complex = str_c(mfd_areatype, ",\n", mfd_hab1)) %>%
  ggplot(aes(x = Tax, y = complex)) +
  geom_point(aes(fill = mean_abundance, size = prevalence, color = Genus_type), shape = 21) +
  scale_fill_viridis_c() +
  # scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(3, "RdYlBu"))) +
  scale_color_manual(values = c("black", "grey")) +
  facet_grid(cols = vars(manage), scales = "free", switch = "y", space = "free") +
  scale_x_discrete(position = "top") +
  scale_y_discrete(position = "right") +
  scale_size(range = c(2,10)) +
  guides(fill = guide_colorbar(title = "Relative\nabundance %", theme(legend.title.position = "top")),
         size = guide_legend(title = "Prevalence", theme(legend.title.position = "top")),
         color = guide_legend(title = "Genus type", theme(legend.title.position = "top"),
                              override.aes = list(size = 4))) +
  coord_flip() +
  xlab("") +
  ylab("") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 0),
        axis.text.y = element_text(size = 12))


p4

### Link 16S reference seqeunces and 16S sequences binned in MAGs
## Load USEARCH results
usearch <- data.table::fread("data/shallow_mags_16S_vs_MFD_ssu_database_v1.3_NR987.txt.b6",
                             col.names = c("Query", "Target", "Percent_ID",
                                           "Alignment_lengt", "#mismatch",
                                           "#gap", "start_query", "end_query",
                                           "start_target", "end_target",
                                           "E-value", "Bit")) %>%
  select(-"E-value", -"Bit") %>%
  mutate(bin = str_remove(Query, "_16S+.*"), .before = "Query")


## Import list MAGs
mags <- data.table::fread("data/mags_shallow_all.tsv")

## Import laboratory metadata
lab.meta <- data.table::fread("data/2025-05-07_MFD_laboratory_metadata_collapsed.csv") %>%
  select(fieldsample_barcode, flat_name) %>%
  mutate(across(flat_name, ~str_remove(., ".fastq.gz")))

## Import MFD metadata
full.meta <- readxl::read_excel("data/2025-04-14_mfd_db.xlsx") %>%
  select(fieldsample_barcode, starts_with("mfd"))

## Combine files
comb <- usearch %>%
  left_join(mags) %>%
  mutate(flat_name = str_remove(bin, "\\.+.*"), .before = "bin") %>%
  left_join(lab.meta) %>%
  left_join(full.meta) %>%
  select(-c(bin, flat_name))

## Filter for genera of interest
comb.filt <- comb %>%
  filter(str_detect(Target, "Nitrospira|MFD_g_198|MFD_g_4907|Ellin6067|mle1-7"))

## Write output
data.table::fwrite(comb.filt, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_MFD_core_to_mags.tsv"))


## Customise plot layout
design <- "
11122
33345"

## Arrange plots
ggarranged <- p1 + p4 + p2 + p3 + plot.sum +
  plot_layout(design = design, guides = "keep")


ggarranged


## Write arranged plots to output
png(file = 'output/core.png',
    width = 1900,
    height = 1200) 
ggarranged
dev.off()

pdf(file = 'output/core.pdf',
    width = 19,
    height = 12) 
ggarranged
dev.off()

tiff(file = 'output/core.tiff',
     width = 1900,
     height = 1200) 
ggarranged
dev.off()

jpeg(file = 'output/core.jpeg',
     width = 1900,
     height = 1200) 
ggarranged
dev.off()

ggsave("output/core.svg", 
       plot = ggarranged, width = 19, height = 12, 
       units = "in", dpi = "retina")

