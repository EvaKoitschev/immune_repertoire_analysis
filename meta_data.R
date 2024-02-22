library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(gridExtra)

# Load the data
all_repertoires <- fread("/mnt/susanne/datasets/all_repertoires.tsv")

# print all subject_id's for which study_id is not available
missing_studies <- all_repertoires %>% group_by(subject_id) %>% summarise(study_id = unique(study_id),sra_id = unique(sra_id))

# insert missing study_id's manually
all_repertoires[nchar(as.character(all_repertoires$subject_id)) == 5, "study_id"] <- "PRJNA248411"
all_repertoires[startsWith(as.character(all_repertoires$subject_id), "sample_"), "study_id"] <- "PRJNA549712"
all_repertoires[startsWith(as.character(all_repertoires$subject_id), "M"), "study_id"] <- "PRJNA248475"
all_repertoires[startsWith(as.character(all_repertoires$subject_id), "HD-"), "study_id"] <- "PRJNA738368"
all_repertoires[startsWith(as.character(all_repertoires$subject_id), "Patient"), "study_id"] <- "PRJNA675463"

fwrite(all_repertoires,"/mnt/susanne/datasets/all_repertoires.tsv", sep = "\t")

df_plot <- all_repertoires %>% select(study_id, sample_id, repertoire_id, status, subject_id)

# group by status and study_id and count the number of subjects
df_plot <- df_plot %>% group_by(study_id, status) %>% summarise(n = n_distinct(subject_id))
# plot one pie chart with number of subjects for each study_id and study id for each status

# make two pie charts, one for MS and one for healthy
df_plot_MS <- df_plot %>% filter(status == "MS")
df_plot_healthy <- df_plot %>% filter(status == "healthy")

# plot the pie charts onto one figure
par(mfrow=c(1,2))
p1 <- ggplot(df_plot_MS, aes(x = "", y = n, fill = study_id)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  labs(title = "MS", x = NULL, y = NULL) +
  theme_void() +
  theme(legend.position = "left") +
  scale_fill_brewer(palette="Greens") +
  # add labels with subject numbers
  geom_text(aes(label = n), position = position_stack(vjust = 0.5))

p2 <- ggplot(df_plot_healthy, aes(x = "", y = n, fill = study_id)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  labs(title = "Healthy", x = NULL, y = NULL) +
  theme_void() +
  theme(legend.position = "right") +
  scale_fill_brewer(palette = "Blues") +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5))

p_combined <- grid.arrange(p1, p2, ncol = 2) # add title to the combined plot
p_combined_with_title <- ggdraw(p_combined) + draw_label("Number of subjects per study", fontface = "bold", x = 0.5, hjust = 0.5, vjust = -13, size = 14)
p_combined_with_title
# save p_combined as pdf
ggsave("/mnt/susanne/datasets/number_of_subjects_per_study.pdf", p_combined_with_title, width = 8.63, height = 4.51)

