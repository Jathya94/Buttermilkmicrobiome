library(tidyverse)
library(qiime2R)
library(ggpubr)

setwd("/Users/jathya/Desktop/Masters/Spring_2026/Molecular_microbiome_analysis/Project/PA7/")
list.files()

# Create output folder if it does not exist
if (!dir.exists("output")) {
  dir.create("output")
}

meta <- read_q2metadata("metadata.tsv") %>%
  select(SampleID, Group)
str(meta)
head(meta)

meta$SampleID <- as.character(meta$SampleID)

observed_features <- read_qza("core-metrics-results/observed_features_vector.qza")$data %>%
  rownames_to_column("SampleID")

shannon <- read_qza("core-metrics-results/shannon_vector.qza")$data %>%
  rownames_to_column("SampleID")

evenness <- read_qza("core-metrics-results/evenness_vector.qza")$data %>%
  rownames_to_column("SampleID")

faith_pd <- read_qza("core-metrics-results/faith_pd_vector.qza")$data %>%
  rownames_to_column("SampleID")

observed_features$SampleID <- as.character(observed_features$SampleID)
shannon$SampleID <- as.character(shannon$SampleID)
evenness$SampleID <- as.character(evenness$SampleID)
faith_pd$SampleID <- as.character(faith_pd$SampleID)

alpha_diversity <- observed_features %>%
  merge(shannon, by = "SampleID") %>%
  merge(evenness, by = "SampleID") %>%
  merge(faith_pd, by = "SampleID")

meta_alpha <- merge(meta, alpha_diversity, by = "SampleID")

nrow(meta_alpha)
head(meta_alpha)

str(meta_alpha)
head(meta_alpha)

meta_alpha$observed_features <- as.numeric(meta_alpha$observed_features)
meta_alpha$shannon <- as.numeric(meta_alpha$shannon)
meta_alpha$pielou_e <- as.numeric(meta_alpha$pielou_e)
meta_alpha$faith_pd <- as.numeric(meta_alpha$faith_pd)

meta_alpha$Group <- factor(meta_alpha$Group,
                           levels = c("Control diet", "Buttermilk diet"))

p1 <- ggplot(meta_alpha, aes(x = Group, y = shannon, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2) +
  theme_q2r() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  ) +
  labs(
    fill = "Diet Group",
    title = "observed_features",
    x = "Group",
    y = "observed_features"
  )

p1

p2 <- ggplot(meta_alpha, aes(x = Group, y = shannon, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(title = "Shannon Diversity",
       x = "Group",
       y = "Shannon")

p2

p3 <- ggplot(meta_alpha, aes(x = Group, y = shannon, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2) +
  theme_q2r() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  ) +
  labs(
    fill = "Diet Group",
    title = "pielou Evenness",
    x = "Group",
    y = "pielou Evenness"
  )

p3

p4 <- ggplot(meta_alpha, aes(x = Group, y = shannon, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2) +
  theme_q2r() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  ) +
  labs(
    fill = "Diet Group",
    title = "faith_pd",
    x = "Group",
    y = "faith_pd"
  )

p4

ggsave("observed_features_plot2.png", p1, width = 5, height = 4)
ggsave("shannon_plot2.png", p2, width = 5, height = 4)
ggsave("pielou_plot2.png", p3, width = 5, height = 4)
ggsave("faith_pd_plot2.png", p4, width = 5, height = 4)

meta <- read_q2metadata("metadata.tsv") %>%
  select(SampleID, Group)

meta_alpha$Group <- factor(meta_alpha$Group)








