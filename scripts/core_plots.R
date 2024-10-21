### Setup env
library(tidyverse)
library(UpSetR)
library(ComplexUpset)
library(patchwork)
library(ampvis2)
library(ggtree)

color_vector_heat <- viridisLite::viridis(n = 4)

setwd("/mfd_core")

load("output/core_microbes.RData")

load('data/2024-03-07_MFD-ampvis-arcbac-data.RData')

mfd.ampvis.arcbac.ra[["metadata"]] <- mfd.ampvis.arcbac.ra[["metadata"]] %>%
  mutate(across(mfd_hab1, ~str_replace(., "Sclerophyllous scrub", "Temperate heath and scrub")),
         across(mfd_hab2, ~str_replace(., "Scrub", "Sclerophyllous scrub"))) %>%
  mutate(complex = str_c(mfd_sampletype, mfd_areatype, mfd_hab1, sep = ", "))

seq.meta <- data.table::fread('data/2024-03-07_combined_metadata.csv') %>%
  filter(fieldsample_barcode %in% samples.subset)

lab.meta <- data.table::fread('data/2023-09-22_corrected_combined_metadata.csv') %>%
  filter(fieldsample_barcode %in% samples.subset)

full.meta <- seq.meta %>%
  left_join(lab.meta)

tax <- data %>%
  select(Kingdom:Genus)

# De novo taxonomy
de.novo <- data %>%
  select(Genus) %>%
  distinct() %>%
  reframe(De_novo = sum(str_detect(Genus, "MFD_g_")),
          Known = nrow(.)-De_novo)

de.novo.core <- df.core.mfdo1 %>%
  select(Genus) %>%
  distinct() %>%
  reframe(De_novo = sum(str_detect(Genus, "MFD_g_")),
          Known = nrow(.)-De_novo)


### Core abundance
# Wide format
core.genera <- df.core.mfdo1 %>%
  select(Genus) %>%
  distinct() %>%
  pull(Genus)

df.core.long <- df.core.mfdo1 %>%
  select(Genus, complex1, Genus_type) %>%
  complete(Genus, complex1) %>%
  mutate(res = if_else(is.na(Genus_type), FALSE, TRUE)) %>%
  select(-Genus_type)

data.long <- data %>%
  select(Genus, starts_with("MFD")) %>%
  pivot_longer(!Genus, names_to = "fieldsample_barcode", values_to = "abundance")

data.long.complex <- data %>%
  select(Genus, starts_with("MFD")) %>%
  pivot_longer(!Genus, names_to = "fieldsample_barcode", values_to = "abundance") %>%
  left_join(metadata) %>%
  group_by(Genus, complex1) %>%
  summarise(across(abundance, ~mean(.))) %>%
  ungroup()

data.mfdo1 <- data.long.complex %>%
  group_by(Genus) %>%
  pivot_wider(names_from = "complex1", values_from = "abundance")

reads.summary <- data.long.complex %>%
  left_join(df.core.mfdo1 %>% select(complex1, Genus, Genus_type)) %>%
  mutate(across(Genus_type, ~if_else(Genus == "Unclassified", Genus, Genus_type)),
         across(Genus_type, ~replace_na(., "Other"))) %>%
  group_by(complex1, Genus_type) %>%
  summarise(RA = sum(abundance))

reads.summary %>%
  filter(!Genus_type == "Unclassified") %>%
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

plot.sum <- reads.summary %>%
  mutate(across(Genus_type, ~factor(., levels = rev(c("Core", "Other", "Unclassified"))))) %>%
  ggplot(aes(x = complex1, 
             y = RA, 
             fill = Genus_type)) +
  geom_bar(color = "black",
           stat = "identity", 
           width = 0.8) +
  scale_fill_viridis_d() +
  guides(fill = guide_legend(override.aes = list(shape = 22, size = 8))) +  
  scale_y_continuous(position = "left") +
  ggtitle("Classified 16S fragments",
          subtitle = "") +
  labs(x = "",
       y = "Relative abundance (%)") +
  guides(fill = guide_legend(title = "Genus type")) +
  theme_minimal(base_size = 18) +
  theme(title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

plot.sum

png(file = 'output/core-sum-relabundance.png',
    width = 1900,
    height = 1200)
plot.sum
dev.off()

pdf(file = 'output/core-sum-relabundance.pdf',
    width = 1900,
    height = 1200)
plot.sum
dev.off()

tiff(file = 'output/core-sum-relabundance.tiff',
     width = 1900,
     height = 1200)
plot.sum
dev.off()

ggsave("output/core-sum-relabundance.svg", plot = plot.sum, width = 19, height = 12, units = "in", dpi = "retina")

df.core.mfdo1 %>%
  group_by(complex1) %>% 
  summarise(n = n())

reads.summary %>%
  filter(Genus_type == "Core") %>%
  summarise(Sum_RA = sum(RA)) %>%
  ungroup() %>%
  rstatix::shapiro_test(Sum_RA)

reads.summary %>%
  filter(Genus_type == "Core") %>%
  summarise(Sum_RA = sum(RA)) %>%
  ungroup() %>%
  summarise(median_RA = median(Sum_RA),
            mean_RA = mean(Sum_RA),
            sd_RA = sd(Sum_RA))

# Percentage of genera
(df.core.mfdo1$Genus %>% unique() %>% length())/(nrow(data)-1)*100

pct.classified <- reads.summary %>%
  filter(!Genus_type == "Unclassified") %>%
  summarise(Sum_RA = sum(RA)) %>%
  ungroup() %>%
  summarise(mean_RA = mean(Sum_RA)) %>%
  pull(mean_RA)

pct.core <- reads.summary %>%
  filter(Genus_type == "Core") %>%
  summarise(Sum_RA = sum(RA)) %>%
  ungroup() %>%
  summarise(mean_RA = mean(Sum_RA)) %>%
  pull(mean_RA)

pct.core/pct.classified*100


### Prevalence
data.pa <- data %>%
  filter(!Genus == "Unclassified") %>%
  mutate(across(where(is.numeric), ~+as.logical(.x)))

prevalence.samples <- data.pa %>%
  select(Genus, starts_with("MFD")) %>%
  mutate(prevalence = rowSums(pick(where(is.numeric), -Genus)), .keep = "unused") %>%
  mutate(de_novo = case_when(str_detect(Genus, "MFD") ~ TRUE,
                             TRUE ~ FALSE))

p.prevalence.samples <- prevalence.samples %>%
  ggplot() +
  geom_bar(aes(x = prevalence), fill = "grey40", color = "black") +
  scale_x_continuous(trans = "reverse", breaks = rev(seq(1, 25, 1)), limits = c(26,0), name = "Number of samples") +
  # scale_y_log10() +
  scale_y_continuous(breaks = seq(0,7000,1000), limits = c(0,7000), name = "Number of genera", position = "right") +
  coord_flip() +
  theme_minimal(base_size = 18) +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1),
        aspect.ratio = 1)

p.prevalence.samples

p.prevalence.samples.bin <- prevalence.samples %>%
  mutate(bin = cut_interval(prevalence, n = 25)) %>%
  ggplot() +
  geom_bar(aes(x = bin, fill = de_novo)) +
  xlab('Number of samples') +
  ylab('Number of genera') +
  theme_minimal(base_size = 18)

p.prevalence.samples.bin

prevalence.mfdo1 <- data.mfdo1 %>%
  filter(!Genus == "Unclassified") %>%
  ungroup() %>%
  mutate(across(where(is.numeric), ~+as.logical(.x))) %>%
  mutate(prevalence = rowSums(pick(where(is.numeric), -Genus)), .keep = "unused") %>%
  mutate(de_novo = case_when(str_detect(Genus, "MFD") ~ TRUE,
                             TRUE ~ FALSE))

# Random groups
data.long.complex.random <- data %>%
  select(Genus, starts_with("MFD")) %>%
  pivot_longer(!Genus, names_to = "fieldsample_barcode", values_to = "abundance") %>%
  left_join(metadata) %>%
  group_by(complex1) %>%
  mutate(group2 = as.character(cur_group_id())) %>%
  ungroup()

random.summary <- function(data) {
  
  tmp <- data %>%
    mutate(group3 = sample(group2)) %>%
    group_by(Genus, group3) %>%
    summarise(across(abundance, ~mean(.))) %>%
    ungroup() %>%
    group_by(Genus) %>%
    pivot_wider(names_from = "group3", values_from = "abundance") %>%
    filter(!Genus == "Unclassified") %>%
    ungroup() %>%
    mutate(across(where(is.numeric), ~+as.logical(.x))) %>%
    mutate(prevalence = rowSums(pick(where(is.numeric), -Genus)), .keep = "unused") %>%
    mutate(de_novo = case_when(str_detect(Genus, "MFD") ~ TRUE,
                               TRUE ~ FALSE)) %>%
    group_by(prevalence) %>%
    summarise(n = n())
  
  return(tmp)
}

set.seed(123)
list <- lst()

for (i in 1:1000) {
  res <- random.summary(data.long.complex.random)
  
  list[[i]] <- res
}

tmp <- bind_rows(list) %>%
  group_by(prevalence) %>%
  summarise(mean = round(mean(n), 1),
            sd = round(sd(n), 1))

p.prevalence.mfdo1.random <- ggplot() +
  geom_bar(data = prevalence.mfdo1, aes(x = prevalence), fill = "grey40", color = "black") +
  geom_errorbar(data = tmp, aes(x = prevalence, y = mean, ymin = mean-sd, ymax = mean+sd), width = 0.2) +
  # geom_point(data = tmp, aes(x = prevalence, y = mean), fill = "grey40", color = "black") +
  # geom_smooth(data = tmp, aes(x = prevalence, y = mean), color = "black", linetype = 2,
  #             method = "lm", formula = y ~ exp(x), se = T) +
  geom_line(data = tmp, aes(x = prevalence, y = mean)) +
  scale_x_continuous(trans = "reverse", breaks = rev(seq(1, 18, 1)), limits = c(19,0), name = "Prevalence") +
  # scale_y_log10() +
  scale_y_continuous(breaks = seq(0,10000,1000), limits = c(0,10000), name = "Number of genera", position = "right") +
  coord_flip() +
  theme_minimal(base_size = 18) +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1),
        aspect.ratio = 1)

p.prevalence.mfdo1.random


p.prevalence.mfdo1 <- prevalence.mfdo1 %>%
  ggplot() +
  geom_bar(aes(x = prevalence), fill = "grey40", color = "black") +
  scale_x_continuous(trans = "reverse", breaks = rev(seq(1, 18, 1)), limits = c(19,0), name = "Prevalence") +
  # scale_y_log10() +
  scale_y_continuous(breaks = seq(0,10000,1000), limits = c(0,10000), name = "Number of genera", position = "right") +
  coord_flip() +
  theme_minimal(base_size = 18) +
  theme(axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1),
        plot.background = element_rect(colour = "black", fill=NA, linewidth=1),
        aspect.ratio = 1)

p.prevalence.mfdo1

ggsave("output/prevalent-genera.svg", plot = p.prevalence.mfdo1, width = 8, height = 8, units = "in", dpi = "retina", bg = "white")

prevalent.genera.mfdo1 <- prevalence.mfdo1 %>%
  filter(prevalence == 18) %>% 
  left_join(tax) %>%
  relocate(Genus, .after = "Family")


# Single inhabitants
prevalence.mfdo1 %>%
  filter(prevalence == 1) %>% 
  nrow()

prevalence.samples %>%
  filter(prevalence == 1) %>% 
  nrow()

# Check prevalence across samples
prevalent.genera.samples <- prevalence.samples %>%
  filter(prevalence >= length(samples.subset)*0.75) 

intersect(prevalent.genera.mfdo1$Genus, prevalent.genera.samples$Genus)

intersect(prevalent.genera.mfdo1$Genus, df.core.mfdo1$Genus)

setdiff(prevalent.genera.mfdo1$Genus, df.core.mfdo1$Genus)

data.table::fwrite(prevalent.genera.mfdo1, 'output/prevalent-genera.csv')

prevalence.mfdo1 %>%
  reframe(de_novo = sum(str_count(de_novo, "TRUE")),
          known = n()-de_novo,
          de_novo_pct = de_novo/(de_novo+known)*100)

# Heatmap of prevalent genera
heatmap.data.prevalent <- mfd.ampvis.arcbac.ra %>%
  amp_subset_samples(fieldsample_barcode %in% samples.subset) %>%
  amp_subset_taxa(tax_vector = prevalent.genera.mfdo1$Genus, normalise = TRUE, remove = FALSE) %>%
  amp_heatmap(group_by = "complex", 
              tax_show = nrow(prevalent.genera.mfdo1), 
              min_abundance = 0.0001,
              tax_aggregate = "Genus", 
              tax_add = c("Phylum"),
              normalise = FALSE,
              textmap = TRUE) %>%
  rownames_to_column(var = "Genus") %>%
  pivot_longer(!Genus, names_to = "complex", values_to = "abund")

p.heatmap.prevalent <- heatmap.data.prevalent %>%
  ggplot(aes(x = reorder(Genus, abund, decreasing = T), y = complex, fill = log10(abund))) +
  geom_tile(color = "black") + 
  scale_x_discrete(position = "bottom") +
  scale_y_discrete(limits=rev) +
  scale_fill_gradientn(colors = color_vector_heat, na.value = color_vector_heat[1],
                       # limits = c(0, 2.5),
                       # labels = rev(c("6.25", "4.00", "2.25", "1.00", "0.25", "0.00")),
                       name = "Relative \nabundance %") +
  guides(fill = guide_colourbar(title.position = "top", ticks.colour = "black")) +
  theme_minimal() +
  theme(axis.text.x.bottom = element_text(angle = 60, hjust = 1),
        axis.title.x.bottom = element_blank(),
        # axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        legend.direction = "vertical",
        legend.position = "left",
        axis.title = element_text(size = 12, face = "bold"))

p.heatmap.prevalent

# Save plot 
ggsave("output/prevalent-heatmap.svg", plot = p.heatmap.prevalent, width = 16, height = 10, units = "in", dpi = "retina")


### Core upset plot
# Upset
hab.levels <- summary.mfdo1 %>%
  select(complex1) %>%
  droplevels() %>%
  pull(complex1)

mfdo1.palette.sub <- mfdo1.palette[as.character(hab.levels)]

upset.df <- df.core.long %>%
  group_by(Genus) %>%
  pivot_wider(names_from = "complex1", values_from = "res") %>%
  ungroup() %>%
  select(Genus, any_of(hab.levels))

sets = colnames(upset.df[hab.levels])

plot.upset <- upset(
  upset.df,
  sets,
  set_sizes = upset_set_size(geom = geom_bar(width = 0.6, color = "black"), position = "left") +
    theme(axis.title.x = element_text(size = 18, face = "bold")),
  name = "",
  guides = 'over',
  sort_sets = FALSE,
  sort_intersections = 'ascending',
  # intersections = 'all',
  sort_intersections_by=c('degree'),
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
    upset_query(set = (names(mfdo1.palette.sub[18])), fill = mfdo1.palette.sub[18])),
  matrix = (intersection_matrix(geom = geom_point(shape = 'circle filled', size = 3)) +
              scale_y_discrete(position = "right")),
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
      theme(axis.title.y = element_text(size = 18, face = "bold"))))

# plot.upset

# Custom intersections
intersections <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18",
                      "1-4", "2-3", "2-5", "2-15", "6-11", "6-14", "7-14", "1-4-7", "6-8-14", "6-7-14", "6-11-14",
                      "9-10-13", "9-10-12-13", "6-7-8-11-14", "9-10-11-12-13", "7-9-10-11-12-13", "6-7-9-10-11-12-13-14")

# Customize upset plot
data.upset1 <- plot.upset[[2]][["data"]] %>%
  filter(intersection %in% intersections) %>%
  # filter(intersection %in% seq(1, 18, 1)) %>%
  select(intersection, in_exclusive_intersection) %>%
  mutate(across(intersection, ~factor(., levels = intersections))) %>%
  complete(., intersection, fill=list(in_exclusive_intersection = 0)) %>%
  group_by(intersection) %>%
  summarise(sum = sum(in_exclusive_intersection))

# Bar plot
p1 <- data.upset1 %>%
  ggplot(aes(x = intersection, y = sum)) +
  geom_col(color = "black", position = "stack") +
  xlab("") +
  ylab("Number of genera") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title = element_text(size = 12, face = "bold"))

p1

# Intersection matrix
data.upset2 <- plot.upset[[4]][["data"]] %>%
  filter(intersection %in% intersections) %>%
  filter(value == "TRUE") %>%
  mutate(across(intersection, ~factor(., levels = intersections))) %>%
  mutate(across(group, ~factor(., levels = seq(1, 18, 1)))) %>%
  complete(., group, intersection) %>%
  mutate(across(group:intersection, ~as.character(.))) %>%
  mutate(value = replace(value, is.na(value) & intersection == group, TRUE)) %>%
  filter(!is.na(value)) %>%
  mutate(across(intersection, ~factor(., levels = intersections))) %>%
  mutate(across(group, ~factor(., levels = seq(1, 18, 1)))) %>%
  mutate(start = str_extract(intersection, "^."),
         end = str_extract(intersection, "[^-]+$"))

data.upset2.1 <- plot.upset[[4]][["data"]] %>%
  filter(intersection %in% intersections) %>%
  mutate(across(intersection, ~factor(., levels = intersections))) %>%
  mutate(across(group, ~factor(., levels = seq(1, 18, 1)))) %>%
  complete(., group, intersection)

p2 <- data.upset2 %>%
  ggplot(aes(x = intersection, y = group)) +
  ggstats::geom_stripped_rows(odd = "grey95", even = "white") +
  geom_point(data = data.upset2.1, fill = "white", color = "grey70", pch = 21, size = 3.5) +
  geom_segment(aes(x = intersection, xend = intersection, y = group, yend = end)) +
  geom_point(aes(fill = group), pch = 21, size = 3.5) +
  scale_fill_manual(values = unname(mfdo1.palette.sub)) +
  scale_y_discrete(limits=rev, labels = rev(names(mfdo1.palette.sub)), position = "left") +
  xlab("") +
  ylab("") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        legend.position = "none")

p2

# Set sizes
plot.upset.tmp <- plot.upset

plot.upset.tmp[[3]][["data"]] <- plot.upset.tmp[[3]][["data"]] %>%
  filter(intersection %in% seq(1, 18, 1)) %>%
  mutate(across(intersection, ~factor(., levels = seq(1, 18, 1))))

p3 <- plot.upset.tmp[[3]] + 
  scale_x_discrete(limits=rev) +
  xlab("") +
  ylab("Number of core genera") +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.title = element_text(size = 12, face = "bold"))

p3

# Combiend plot
plot.upset.new <- plot_spacer() + p1 + p3 + p2 + plot_layout(ncol = 2, widths = c(1,2))

plot.upset.new

# Save plot 
ggsave("output/core-bar-plot.svg", plot = p1, width = 16, height = 10, units = "in", dpi = "retina")
ggsave("output/core-intersection-matrix.svg", plot = p2, width = 16, height = 10, units = "in", dpi = "retina")
ggsave("output/core-set-sizes.svg", plot = p3, width = 16, height = 10, units = "in", dpi = "retina")
ggsave("output/core-upset.svg", plot = plot.upset.new, width = 16, height = 10, units = "in", dpi = "retina")



### Core genera intersects
intersects <- upset.df %>%
  mutate(count = across(2:19, ~str_count(., "TRUE"))) %>%
  rowwise() %>%
  mutate(count_yes = across(.cols = contains("count"), .fns = sum)) %>%
  left_join(tax) %>%
  relocate(Kingdom:Family, .before = "Genus")

# Non-overlapping genera
genera.intersect.unique <- intersects %>%
  select(1:24, count_yes) %>%
  filter(count_yes$count == 1) %>%
  select(-count_yes) %>%
  pivot_longer(!Kingdom:Genus, names_to = "complex", values_to = "res") %>%
  filter(res == TRUE) %>%
  select(-res)

data.table::fwrite(genera.intersect.unique, "output/2024-03-07_unique-intersects.csv")

# Disturbed
genera.intersect.disturbed <- intersects %>%
  select(Kingdom:Genus, 12, 20, count_yes) %>%
  filter(if_all(where(is.logical), ~ . %in% TRUE)) %>%
  filter(count_yes$count == 2) %>%
  select(-count_yes) %>%
  mutate(Intersect = paste(colnames(select(., where(is.logical))), collapse = ' ; ')) %>%
  select(Intersect, Kingdom:Genus)

# Natural soil
genera.intersect.natural1 <- intersects %>%
  select(Kingdom:Genus, 15, 16, 18, 19, count_yes) %>%
  filter(if_all(where(is.logical), ~ . %in% TRUE)) %>%
  filter(count_yes$count == 4) %>%
  select(-count_yes) %>%
  mutate(Intersect = paste(colnames(select(., where(is.logical))), collapse = ' ; ')) %>%
  select(Intersect, Kingdom:Genus)

# Natural soil 2
genera.intersect.natural2 <- intersects %>%
  select(Kingdom:Genus, 15:19, count_yes) %>%
  filter(if_all(where(is.logical), ~ . %in% TRUE)) %>%
  filter(count_yes$count == 5) %>%
  select(-count_yes) %>%
  mutate(Intersect = paste(colnames(select(., where(is.logical))), collapse = ' ; ')) %>%
  select(Intersect, Kingdom:Genus)

# Natural soil 3
genera.intersect.natural3 <- intersects %>%
  select(Kingdom:Genus, c(12,13,15:20), count_yes) %>%
  filter(if_all(where(is.logical), ~ . %in% TRUE)) %>%
  filter(count_yes$count == 8) %>%
  select(-count_yes) %>%
  mutate(Intersect = paste(colnames(select(., where(is.logical))), collapse = ' ; ')) %>%
  select(Intersect, Kingdom:Genus)

# Freshwater
genera.intersect.freshwater <- intersects %>%
  select(Kingdom:Genus, 7, 10, count_yes) %>%
  filter(if_all(where(is.logical), ~ . %in% TRUE)) %>%
  filter(count_yes$count == 2) %>%
  select(-count_yes) %>%
  mutate(Intersect = paste(colnames(select(., where(is.logical))), collapse = ' ; ')) %>%
  select(Intersect, Kingdom:Genus)

# Freshwater and Bogs, mires and fens
genera.intersect.freswater.bogs <- intersects %>%
  select(Kingdom:Genus, 7, 10, 13, count_yes) %>%
  filter(if_all(where(is.logical), ~ . %in% TRUE)) %>%
  filter(count_yes$count == 3) %>%
  select(-count_yes) %>%
  mutate(Intersect = paste(colnames(select(., where(is.logical))), collapse = ' ; ')) %>%
  select(Intersect, Kingdom:Genus)




## Core heatmap
genus.list <- mfd.ampvis.arcbac.ra %>%
  amp_subset_samples(fieldsample_barcode %in% samples.subset) %>%
  amp_subset_taxa(tax_vector = core.genera, normalise = TRUE) %>%
  amp_merge_replicates(merge_var = "complex") %>%
  amp_export_long(metadata_vars = c("mfd_sampletype", "mfd_areatype", "mfd_hab1", "complex"), tax_levels = "Genus") %>%
  select(-OTU) %>%
  group_by(Genus, complex) %>%
  mutate(cumabund = sum(count)) %>%
  filter(cumabund > 0) %>%
  ungroup() %>%
  arrange(desc(cumabund)) %>%
  select(complex, Genus, cumabund) %>%
  group_by(complex) %>%
  slice_head(n = 4) %>%
  pull(Genus) %>%
  unique()

## Soil
genus.soil <- rbind(genera.intersect.disturbed, genera.intersect.natural1, 
                    genera.intersect.natural2, genera.intersect.natural3) %>%
  pull(Genus)

ampvis.core.soil <- mfd.ampvis.arcbac.ra %>%
  amp_subset_samples(fieldsample_barcode %in% samples.subset) %>%
  amp_subset_taxa(tax_vector = genus.soil, normalise = TRUE) %>%
  amp_merge_replicates(merge_var = "complex")

heatmap.data.core.soil <- ampvis.core.soil %>% 
  amp_heatmap(group_by = "complex", 
              tax_show = length(genus.soil), 
              min_abundance = 0.0001,
              tax_aggregate = "Genus", 
              # tax_add = c("Phylum"),
              normalise = FALSE,
              textmap = TRUE) %>%
  rownames_to_column(var = "Genus") %>%
  mutate(across(Genus, ~factor(., levels = genus.soil))) %>%
  pivot_longer(!Genus, names_to = "complex", values_to = "abund") %>%
  mutate(facet = case_when(Genus %in% genera.intersect.disturbed$Genus ~ "Disturbed",
                           Genus %in% genera.intersect.natural1$Genus ~ "Natural",
                           Genus %in% genera.intersect.natural2$Genus ~ "Natural", 
                           Genus %in% genera.intersect.natural3$Genus ~ "Disturbed + Natural"),
         across(facet, ~factor(., levels = c("Disturbed", "Natural", "Disturbed + Natural"))))

p.heatmap.core.soil <- heatmap.data.core.soil %>%
  ggplot(aes(x = complex, y = Genus, fill = log10(abund))) +
  geom_tile(color = "black") + 
  # scale_y_discrete(position = "bottom") +
  scale_x_discrete(limits=rev) +
  scale_fill_gradientn(colors = color_vector_heat, na.value = color_vector_heat[1],
                       # limits = c(0, 4),
                       # labels = rev(c("16", "9", "4", "1", "0.0")),
                       name = "Relative \nabundance %") +
  guides(fill = guide_colourbar(title.position = "top", ticks.colour = "black")) +
  facet_grid(rows = vars(facet), scales = "free", space = "free") +
  theme_minimal() +
  theme(axis.text.x.bottom = element_text(angle = 60, hjust = 1),
        axis.title.x.bottom = element_blank(),
        # axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        legend.direction = "vertical",
        legend.position = "left",
        axis.title = element_text(size = 12, face = "bold"))

p.heatmap.core.soil

## Sediment
genus.snf <- genera.intersect.unique %>% 
  filter(complex == "Sediment, Natural, Freshwater") %>%
  pull(Genus)

genus.suf <- genera.intersect.unique %>% 
  filter(complex == "Sediment, Urban, Freshwater") %>%
  pull(Genus)

genus.snb <- genera.intersect.unique %>% 
  filter(complex == "Soil, Natural, Bogs, mires and fens") %>%
  pull(Genus)

genus.sediment <- rbind(genera.intersect.freshwater, genera.intersect.freswater.bogs) %>%
  rbind(genera.intersect.unique %>% rename(Intersect = complex) %>%
          filter(Intersect %in% c("Sediment, Natural, Freshwater", "Sediment, Urban, Freshwater", "Soil, Natural, Bogs, mires and fens"))) %>%
  pull(Genus)

ampvis.core.sediment <- mfd.ampvis.arcbac.ra %>%
  amp_subset_samples(fieldsample_barcode %in% samples.subset) %>%
  amp_subset_taxa(tax_vector = genus.sediment, normalise = TRUE) %>%
  amp_merge_replicates(merge_var = "complex")

heatmap.data.core.sediment <- ampvis.core.sediment %>% 
  amp_heatmap(group_by = "complex", 
              tax_show = length(genus.sediment), 
              min_abundance = 0.0001,
              tax_aggregate = "Genus", 
              # tax_add = c("Phylum"),
              normalise = FALSE,
              textmap = TRUE) %>%
  rownames_to_column(var = "Genus") %>%
  mutate(across(Genus, ~factor(., levels = genus.sediment))) %>%
  pivot_longer(!Genus, names_to = "complex", values_to = "abund") %>%
  mutate(facet = case_when(Genus %in% genera.intersect.freshwater$Genus ~ "Freshwater",
                           Genus %in% genera.intersect.freswater.bogs$Genus ~ "Freshwater + BMF",
                           Genus %in% genus.snf ~ "Natural",
                           Genus %in% genus.suf ~ "Urban",
                           Genus %in% genus.snb ~ "BMF"),
         across(facet, ~factor(., levels = c("Natural", "Urban", "Freshwater", "BMF", "Freshwater + BMF"))))

p.heatmap.core.sediment <- heatmap.data.core.sediment %>%
  ggplot(aes(x = Genus, y = complex, fill = log10(abund))) +
  geom_tile(color = "black") + 
  scale_x_discrete(position = "bottom") +
  scale_y_discrete(limits=rev) +
  scale_fill_gradientn(colors = color_vector_heat, na.value = color_vector_heat[1],
                       # limits = c(0, 4),
                       # labels = rev(c("16", "9", "4", "1", "0.0")),
                       name = "Relative \nabundance %") +
  guides(fill = guide_colourbar(title.position = "top", ticks.colour = "black")) +
  facet_grid(cols = vars(facet), scales = "free", space = "free") +
  theme_minimal() +
  theme(axis.text.x.bottom = element_text(angle = 60, hjust = 1),
        axis.title.x.bottom = element_blank(),
        # axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        legend.direction = "vertical",
        legend.position = "left",
        axis.title = element_text(size = 12, face = "bold"))

p.heatmap.core.sediment




# Selection
sample.tmp <- metadata %>%
  # filter(str_detect(complex, "Soil")) %>%
  filter(str_detect(complex, "Soil|Sediment, Natural, Freshwater|Sediment, Urban, Freshwater")) %>%
  pull(fieldsample_barcode)

sample.hab <- metadata %>%
  # filter(str_detect(complex, "Soil")) %>%
  filter(str_detect(complex, "Soil|Sediment, Natural, Freshwater|Sediment, Urban, Freshwater")) %>%
  pull(complex) %>%
  unique() 
  # droplevels()

sample.genus <- df.core.mfdo1 %>%
  # filter(str_detect(complex, "Soil")) %>%
  filter(str_detect(complex1, "Soil|Sediment, Natural, Freshwater|Sediment, Urban, Freshwater")) %>%
  select(Genus) %>% 
  pull(Genus) %>% 
  unique()

sample.data <- data.mfdo1 %>%
  ungroup() %>%
  select(Genus, any_of(sample.hab)) %>%
  filter(Genus %in% sample.genus)

sample.genus.select <- sample.data %>% 
  rowwise(Genus) %>%
  mutate(range = max(c_across(where(is.numeric)))-min(c_across(where(is.numeric))),
         mean = mean(c_across(where(is.numeric))),
         sd = sd(c_across(where(is.numeric))), .keep = "unused") %>%
  ungroup() %>%
  arrange(desc(sd))

genus.select.top <- sample.genus.select %>%
  slice_head(n = 30) %>%
  pull(Genus)

reduced.genus.hell <- sample.data %>%
  filter(Genus %in% sample.genus.select$Genus) %>%
  column_to_rownames(var = "Genus") %>%
  t() %>%
  vegan::decostand(., method = "hellinger") %>%
  as.data.frame()

BC.wrap <- function(counts, lab){
  mat <- parallelDist::parallelDist(counts,
                       method = "bray",
                       binary = F,
                       diag = T,
                       upper = T,
                       threads = 10) %>%
    as.matrix()
  mat.to_return <- get_av_dist(mat, lab) %>%
    as.dist()
  return(mat.to_return)
}

## hclust example
createHclustObject <- function(x)hclust(vegan::vegdist(x, "bray"), "ave")

## bootstrap
set.seed(123)
b.small <- bootstrap::bootstrap(reduced.genus.hell, fun=createHclustObject, n = 100L, mc.cores = 10)

## plot
hc.small <- createHclustObject(reduced.genus.hell)
plot(hc.small)

## draw bootstrap values to corresponding node
bootstrap::bootlabels.hclust(hc.small, b.small, col="blue")

tmpa <- hc.small$order

tmpb <- hc.small$labels

tmpc <- tmpb[tmpa]

tmpc



# tmp <- hclust(dist(t(reduced.genus.hell), method = "euclidean"))
tmp <- hclust(vegan::vegdist(t(reduced.genus.hell), method = "bray"))

tmp2 <- tmp$order

tmp3 <- tmp$labels %>% str_replace(., "Burkholderia-Caballeronia-Paraburkholderia", "BCP")

tmp4 <- tmp3[tmp2]

tmp4

heatmap.data.sample <- sample.data %>%
  mutate(across(Genus, ~str_replace(., "Burkholderia-Caballeronia-Paraburkholderia", "BCP"))) %>%
  filter(Genus %in% genus.select.top) %>%
  mutate(across(Genus, ~factor(., levels = tmp4))) %>%
  pivot_longer(!Genus, names_to = "complex", values_to = "abund") %>%
  mutate(across(complex, ~factor(., levels = rev(tmpc))))

p.heatmap.core.sample <- heatmap.data.sample %>%
  ggplot(aes(x = complex, y = Genus, fill = sqrt(abund))) +
  geom_tile(color = "black") + 
  # scale_y_discrete(position = "bottom") +
  scale_x_discrete(limits=rev, position = "bottom") +
  scale_y_discrete(position = "right") +
  scale_fill_gradientn(colors = color_vector_heat, na.value = color_vector_heat[1],
                       limits = c(0, 4.28),
                       labels = rev(c("16", "9", "4", "1", "0")),
                       name = "Relative \nabundance %") +
  guides(fill = guide_colourbar(title.position = "top", ticks.colour = "black")) +
  theme_minimal() +
  # geom_rect(mapping = aes(xmin = 0.5, xmax = 4.5, ymin = 0.5, ymax = 5.5),
  #           fill = NA, col = "white") + 
  # geom_rect(mapping = aes(xmin = 5.5, xmax = 9.5, ymin = 5.5, ymax = 10.5),
  #           fill = NA, col = "white") + 
  # geom_rect(mapping = aes(xmin = 0.5, xmax = 4.5, ymin = 12.5, ymax = 18.5),
  #           fill = NA, col = "white") + 
  theme(axis.text.x = element_blank(),
        # axis.text.x = element_text(angle = 60, hjust = 1),
        axis.title.x = element_blank(),
        # axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        legend.direction = "vertical",
        legend.position = "right",
        axis.title = element_text(size = 12, face = "bold"))

p.heatmap.core.sample

## TREE
p.rep <- hc.small %>%
  ggtree(ladderize = F) +
  layout_dendrogram()

p.rep

b.df <- data.frame(x = -(round(hc.small$height, 6)),
                   bootstrap = b.small) %>%
  left_join((p.rep$data %>% select(x, node) %>% as.data.frame() %>% mutate(across(x, ~round(., 6)))),
            by = "x")

p.colors.rep <- p.rep %<+% data.frame(label = colnames(reduced.genus.hell)) +
  geom_tiplab(align = T, linesize = .5, as_ylab = T, angle = 90) +
  xlim(0, -0.45) +
  # xlim(0, -0.35) +
  geom_tippoint(aes(fill = label), shape = 21, size = 5) +
  scale_fill_manual(values = mfdo1.palette) +
  guides(fill = "none") +
  ggnewscale::new_scale_fill() +
  theme(axis.text.x = element_text(hjust = 0))

p.colors.rep

p.bootstrap.rep <- p.colors.rep %<+% b.df +
  geom_point(data = (. %>% filter(!isTip)), aes(fill=bootstrap), shape = 23, size = 3, color = "black") +
  scale_fill_distiller(type = "seq", palette = "Greys", direction = 1,  na.value = NA) +
  scale_fill_stepsn(colors = rev(grey.colors(8))[c(1, 8)], na.value = NA,
                    limits = c(0.4, 1),
                    breaks = c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0),
                    labels = c("0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0"),
                    name = "Bootstrap") +
  #scale_fill_binned(type = rev(grey.colors(8))[c(2, 4, 6, 8)], na.value = NA) +
  scale_color_continuous(na.value = NA) +
  ggnewscale::new_scale_fill() +
  ggnewscale::new_scale_color()

p.bootstrap.rep

plot.tmp <- p.bootstrap.rep + p.heatmap.core.sample + plot_layout(nrow = 2, heights = c(1,2)) 

plot.tmp

## Arrange panel
design1 <- "
122
344"

ggarranged1 <- plot_spacer() + p1 + p3 + p2 + plot.tmp + plot_layout(ncol = 2, widths = c(1,2)) + plot_layout(design = design1)

ggarranged1


png(file = 'output/core-panel-v1.png',
    width = 1000,
    height = 600) 
ggarranged1
dev.off()

ggsave("output/core-panel-v1.svg", plot = ggarranged1, width = 16, height = 10, units = "in", dpi = "retina")


design2 <- "
1225
3445"

ggarranged2 <- plot_spacer() + p1 + p3 + p2 + plot.tmp + plot_layout(ncol = 2, widths = c(1,2)) + plot_layout(design = design2)

ggarranged2


png(file = 'output/core-panel-v2.png',
    width = 1000,
    height = 600) 
ggarranged2
dev.off()

ggsave("output/core-panel-v2.svg", plot = ggarranged2, width = 16, height = 10, units = "in", dpi = "retina")
