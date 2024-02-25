library("dplyr")
library("alakazam")
library("airr")
library("data.table")

input_MS <- "/mnt/susanne/datasets/MS_total"
input_healthy <- "/mnt/susanne/datasets/healthy_total"

if (dir.exists(input_MS)) {
  input_files <- list.files(input_MS, pattern = "*.tsv", full.names = T)
  input_rearr_MS <- lapply(input_files, read_rearrangement)
}

if (dir.exists(input_healthy)) {
  input_files <- list.files(input_healthy, pattern = "*.tsv", full.names = T)
  input_rearr_healthy <- lapply(input_files, read_rearrangement)
}

all_repertoires_MS <- bind_rows(input_rearr_MS)
all_repertoires_healthy <- bind_rows(input_rearr_healthy)
all_repertoires_MS$status <- "MS"
all_repertoires_healthy$status <- "healthy"
all_repertoires <- bind_rows(all_repertoires_MS,all_repertoires_healthy)

fwrite(all_repertoires,"/mnt/susanne/datasets/all_repertoires.tsv", sep = "\t")

### if starting from here, load the all_repertoires.tsv file ####

all_repertoires <- fread("/mnt/susanne/datasets/all_repertoires.tsv", sep = "\t")
all_repertoires_heavy <- all_repertoires %>% filter(locus=="IGH")

families_all_subjects_by_sequence <- countGenes(all_repertoires_heavy,gene="v_call",groups=c("subject_id","status"), mode="family")
fwrite(families_all_subjects_by_sequence,"/mnt/susanne/datasets/gene_usage/families_all_subjects_by_sequence.tsv", sep = "\t")
# get the gene usage for each subject and status by sequence
genes_all_subjects_by_sequence <- countGenes(all_repertoires_heavy,gene="v_call",groups=c("subject_id","status"),mode="gene")
fwrite(genes_all_subjects_by_sequence,"/mnt/susanne/datasets/gene_usage/genes_all_subjects_by_sequence.tsv", sep = "\t")

# get the gene usage for each subject and status by clone
genes_all_subjects_by_clone <- countGenes(all_repertoires_heavy,gene="v_call",groups=c("subject_id","status"),clone="clone_id",mode="gene")
fwrite(genes_all_subjects_by_clone,"/mnt/susanne/datasets/gene_usage/genes_all_subjects_by_clone.tsv", sep = "\t")

################### Clonally Expanded Sequences ################################

# extract clonally expanded sequences
all_repertoires_clonally_expanded <- all_repertoires_heavy[all_repertoires_heavy$clone_size_count >= 2,]

# get the gene usage for each subject and status by sequence for clonally expanded sequences
genes_clonally_expanded_by_sequence <- countGenes(all_repertoires_clonally_expanded,gene="v_call",groups=c("subject_id","status"),mode="gene")
fwrite(genes_clonally_expanded_by_sequence,"/mnt/susanne/datasets/gene_usage/genes_clonally_expanded_by_sequence.tsv", sep = "\t")

# get the gene usage for each subject and status by clone for clonally expanded sequences
genes_clonally_expanded_by_clone <- countGenes(all_repertoires_clonally_expanded,gene="v_call",groups=c("subject_id","status"),clone="clone_id",mode="gene")
fwrite(genes_clonally_expanded_by_clone,"/mnt/susanne/datasets/gene_usage/genes_clonally_expanded_by_clone.tsv", sep = "\t")


# function to get genes from v_call column (some sequences have multiple genes in the v_call column)
get_genes <- function(dataframe){
  return(unique(dataframe$gene))
}

get_stars <- function(p_value) {
  if (p_value < 0.001) {
    return("****")
  } else if (p_value < 0.01) {
    return("***")
  } else if (p_value < 0.05) {
    return("**")
  } else {
    return("")
  }
}
################## Plot family usage for MS and healthy subjects #####################
families_all_subjects_by_sequence <- fread("/mnt/susanne/datasets/gene_usage/families_all_subjects_by_sequence.tsv", sep = "\t")
family_colors <- c("IGHV1" = "coral", "IGHV2" = "red", "IGHV3" = "hotpink", "IGHV4" = "purple", "IGHV5" = "dodgerblue", "IGHV6" = "green", "IGHV7" = "darkgreen", "Other" = "black")

family_usage = data.frame()
df = families_all_subjects_by_sequence
families = get_genes(families_all_subjects_by_sequence)

for (current_family in families){
  family_MS <- df %>% filter(gene == current_family & status == "MS") %>% select(gene,status,seq_freq)
  if(nrow(family_MS) < 7){
    family_MS <- rbind(family_MS,data.frame(gene=current_family, status="MS",seq_freq=rep(0,7-nrow(family_MS))))
  }
  family_healthy <- df %>% filter(gene == current_family & status == "healthy") %>% select(gene,status,seq_freq)
  if(nrow(family_healthy) < 7){
    family_healthy <- rbind(family_healthy,data.frame(gene=current_family, status="healthy",seq_freq=rep(0,7-nrow(family_healthy))))
  }
  p.value <- round(wilcox.test(family_MS$seq_freq,family_healthy$seq_freq,alternative="greater")$p.value,4)
  family_both <- rbind(family_MS,family_healthy)
  family_both$p.value <- p.value
  family_usage <- rbind(family_usage,family_both)
}

# apply the Benjamini-Hochberg correction for multiple testing to the family_usage dataframe (one p-value per family)
family_usage$adjusted_p.value <- p.adjust(family_usage$p.value, method = "BH")

colours <- sapply(sort(unique(family_usage$gene)), function(x) family_colors[x])

# plot the family usage for MS and healthy subjects for each family as two boxplots using ggplot(gene~status*seq_freq)
p0 <- ggplot(family_usage, aes(x = gene,y=seq_freq, fill = status)) +
  geom_boxplot(width = 0.8, outlier.shape = 1, outlier.size = 1) +
  scale_fill_manual(values = c("MS" = "lightcoral", "healthy" = "lavender")) +
  coord_cartesian(ylim = c(-0.01, 0.6)) +
  labs(title = "Family usage for MS and healthy subjects", x = "IGHV family", y = "Frequency") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,colour = colours),plot.title = element_text(hjust = 0.5)) + 
  # add ** to the plot if the adjusted p-value is less than 0.05 and *** if it is less than 0.01 and **** if it is less than 0.001
  geom_text(aes(label = "*", x = gene, y = 0.61), size = 5) +
  theme(legend.position = c(0.91,0.78))
p0

# save with high resolution as svg and in A5 format
ggsave("/mnt/susanne/datasets/gene_usage/family_usage_MS_vs_healthy_significant.pdf", p0,width = 8.27, height = 5.83, units = "in", dpi = 500)

############## Plot gene usage for MS and healthy subjects #####################

## summarise and plot the gene usage for MS and healthy subjects
genes_all_subjects_by_sequence <- fread("/mnt/susanne/datasets/gene_usage/genes_all_subjects_by_sequence.tsv", sep = "\t")
genes_all_subjects_by_clone <- fread("/mnt/susanne/datasets/gene_usage/genes_all_subjects_by_clone.tsv",sep="\t")



# function: for each gene in the list, get the list of percentages of MS and healthy subjects (into new dataframe) and compute the mann-whitney u test for each gene between MS and healthy subjects
genes <- get_genes(genes_all_subjects_by_sequence)
length(genes)

gene_usage <- data.frame()
df <- genes_all_subjects_by_sequence

for (current_gene in genes){
  gene_MS <- df %>% filter(gene == current_gene & status == "MS") %>% select(gene,status,seq_freq)
  if(nrow(gene_MS) < 27){
    gene_MS <- rbind(gene_MS,data.frame(gene=current_gene, status="MS",seq_freq=rep(0,27-nrow(gene_MS))))
  }
  gene_healthy <- df %>% filter(gene == current_gene & status == "healthy") %>% select(gene,status,seq_freq)
  if(nrow(gene_healthy) < 27){
    gene_healthy <- rbind(gene_healthy,data.frame(gene=current_gene, status="healthy",seq_freq=rep(0,27-nrow(gene_healthy))))
  }
  # compute the Mann-Whitney U test for the current gene between MS and healthy subjects and store the gene and the p-value in gene_usage
  p.value <- round(wilcox.test(gene_MS$seq_freq,gene_healthy$seq_freq,alternative="greater")$p.value,4)
  gene_both <- rbind(gene_MS,gene_healthy)
  gene_both$p.value <- p.value
  gene_usage <- rbind(gene_usage,gene_both)

}


# apply the Benjamini-Hochberg correction for multiple testing to the gene_usage dataframe (one p-value per gene)
gene_usage$adjusted_p.value <- p.adjust(gene_usage$p.value, method = "BH")

colours <- sapply(sort(unique(gene_usage$gene)), function(x) {
  gene_family <- substr(x, 1, 5)
  return(family_colors[gene_family])})

# plot the gene usage for MS and healthy subjects for each gene as two boxplots using ggplot(gene~status*seq_freq)
p1 <- ggplot(gene_usage, aes(x = gene,y=seq_freq, fill = status)) +
  geom_boxplot(width = 0.9, outlier.shape = 1, outlier.size = 1) +
  scale_fill_manual(values = c("MS" = "lightcoral", "healthy" = "lavender")) +
  coord_cartesian(ylim = c(0, 0.18)) +
  labs(title = "Gene usage for MS and healthy subjects", x = "IGHV gene", y = "Frequency") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,colour=colours),plot.title = element_text(hjust = 0.5)) + 
  # add ** to the plot if the adjusted p-value is less than 0.05 and *** if it is less than 0.01 and **** if it is less than 0.001
  geom_text(aes(label = sapply(adjusted_p.value,get_stars), x = gene, y = 0.18),size=2.7) +
  theme(legend.position = "top")
p1


ggsave("/mnt/susanne/datasets/gene_usage/gene_usage_MS_vs_healthy_full.pdf", p1, width = 10.41, height = 3.96)


# plot significant genes only
gene_usage_significant <- gene_usage %>% filter(adjusted_p.value < 0.05)
colours <- sapply(sort(unique(gene_usage_significant$gene)), function(x) {
  gene_family <- substr(x, 1, 5)
  return(family_colors[gene_family])})

p2 <- ggplot(gene_usage_significant, aes(x = gene,y=seq_freq, fill = status)) +
  geom_boxplot(width = 0.8, outlier.shape = 1, outlier.size = 1) +
  scale_fill_manual(values = c("MS" = "lightcoral", "healthy" = "lavender")) +
  coord_cartesian(ylim = c(0, 0.11)) +
  labs(title = "Gene usage for MS and healthy subjects", x = "IGHV gene", y = "Frequency") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,colour = colours),plot.title = element_text(hjust = 0.5)) + 
  # add ** to the plot if the adjusted p-value is less than 0.05 and *** if it is less than 0.01 and **** if it is less than 0.001
  geom_text(aes(label = sapply(adjusted_p.value,get_stars), x = gene, y = 0.11), size = 5) +
  # legend inside the plot
  theme(legend.position = c(0.93,0.75))

p2
# save the plot as pdf
ggsave("/mnt/susanne/datasets/gene_usage/gene_usage_MS_vs_healthy_significant.pdf", p2, width = 9, height = 4.51)

##################### PERIPHERY VS CNS ANALYSIS ############################

# create a column compartment with tissue type (periphery or CNS) based on tissue column (blood -> periphery, CSF -> CNS, cervical lymph nodes -> periphery, brain white matter -> CNS)
all_repertoires_heavy$compartment <- NA  # Initialize compartment column

# assign compartment based on tissue type
all_repertoires_heavy$compartment[all_repertoires_heavy$tissue %in% c("blood", "cervical lymph node","PBMC","peripheral blood")] <- "periphery"
all_repertoires_heavy$compartment[all_repertoires_heavy$tissue %in% c("CSF","brain white matter")] <- "CNS"
# filter out the sequences that are not in the periphery or CNS compartment
all_repertoires_filtered_CNS <- all_repertoires_heavy[!is.na(all_repertoires_heavy$compartment),]

# filter subjects with less than 100 sequences in the CNS compartment
subjects_CNS <- all_repertoires_filtered_CNS %>% filter(compartment =="CNS") %>% group_by(subject_id) %>% filter(n() >= 100) %>% ungroup()
subjects_CNS <- unique(subjects_CNS$subject_id)
all_repertoires_filtered_CNS <- all_repertoires_filtered_CNS %>% filter(subject_id %in% subjects_CNS)

fwrite(all_repertoires_filtered_CNS,"/mnt/susanne/datasets/all_repertoires_filtered_CNS.tsv", sep = "\t")
# get the gene usage for each subject and status by sequence for blood vs. CNS compartment
families_tissue_by_sequence <- countGenes(all_repertoires_filtered_CNS, gene="v_call", groups=c("subject_id","compartment"), mode="family")
fwrite(families_tissue_by_sequence,"/mnt/susanne/datasets/gene_usage/families_periphery_CNS_by_sequence.tsv", sep = "\t")

families_tissue_by_sequence <- fread("/mnt/susanne/datasets/gene_usage/families_periphery_CNS_by_sequence.tsv", sep = "\t")

families <- get_genes(families_tissue_by_sequence)
length(families)

gene_usage = data.frame()
df = families_tissue_by_sequence

for (current_family in families){
  family_periphery <- df %>% filter(gene == current_family & compartment == "periphery") %>% select(gene,compartment,seq_freq)
  if(nrow(family_periphery) < 15){
    family_periphery <- rbind(family_periphery,data.frame(gene=current_family, compartment="periphery",seq_freq=rep(0,15-nrow(family_periphery))))
  }
  family_CNS <- df %>% filter(gene == current_family & compartment == "CNS") %>% select(gene,compartment,seq_freq)
  if(nrow(family_CNS) < 15){
    family_CNS <- rbind(family_CNS,data.frame(gene=current_family, compartment="CNS",seq_freq=rep(0,15-nrow(family_CNS))))
  }
  p.value <- round(wilcox.test(family_periphery$seq_freq, family_CNS$seq_freq, paired = TRUE, alternative = "two.sided")$p.value, 4)
  family_both <- rbind(family_periphery,family_CNS)
  family_both$p.value <- p.value
  gene_usage <- rbind(gene_usage,family_both)
}

# apply the Benjamini-Hochberg correction for multiple testing
gene_usage$adjusted_p.value <- p.adjust(gene_usage$p.value, method = "BH")

# Generate colors based on sorted gene names
colours <- sapply(sort(unique(gene_usage$gene)), function(x) {
  gene_family <- substr(x, 1, 5)
  return(family_colors[gene_family])
})

# plot the gene usage for periphery and CNS compartment for all gene families (significant or not)
p3 <- ggplot(gene_usage, aes(x = gene,y=seq_freq, fill = compartment)) +
  scale_fill_manual(values = c("periphery" = "seagreen", "CNS" = "skyblue")) +
  geom_boxplot(width = 0.8, outlier.shape = 1, outlier.size = 1) +
  coord_cartesian(ylim = c(-0.01,1)) +
  labs(title = "Family usage for periphery and CNS compartment", x = "IGHV family", y = "Frequency") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,colour=colours),plot.title = element_text(hjust = 0.5)) + 
  # add ** to the plot if the adjusted p-value is less than 0.05 and *** if it is less than 0.01 and **** if it is less than 0.001
  geom_text(aes(label = sapply(adjusted_p.value,get_stars), x = gene, y = 1), size = 5) +
  theme(legend.position = c(0.9,0.8))

p3
# save the plot as pdf
ggsave("/mnt/susanne/datasets/gene_usage/gene_usage_periphery_CNS.pdf", myplot_combined, width = 10.41, height = 3.96)





#### Gene usage for individual genes

# get the gene usage for each subject and status by sequence for blood vs. CNS compartment
genes_tissue_by_sequence <- countGenes(all_repertoires_filtered_CNS, gene="v_call", groups=c("subject_id","compartment"), mode="gene")
fwrite(genes_tissue_by_sequence,"/mnt/susanne/datasets/gene_usage/genes_periphery_CNS_by_sequence.tsv", sep = "\t")


genes_tissue_by_sequence <- fread("/mnt/susanne/datasets/gene_usage/genes_periphery_CNS_by_sequence.tsv", sep = "\t")

genes <- get_genes(genes_tissue_by_sequence)
length(genes)

gene_usage = data.frame()
df = genes_tissue_by_sequence

for (current_gene in genes){
  gene_periphery <- df %>% filter(gene == current_gene & compartment == "periphery") %>% select(gene,compartment,seq_freq)
  if(nrow(gene_periphery) < 15){
    gene_periphery <- rbind(gene_periphery,data.frame(gene=current_gene, compartment="periphery",seq_freq=rep(0,15-nrow(gene_periphery))))
  }
  gene_CNS <- df %>% filter(gene == current_gene & compartment == "CNS") %>% select(gene,compartment,seq_freq)
  if(nrow(gene_CNS) < 15){
    gene_CNS <- rbind(gene_CNS,data.frame(gene=current_gene, compartment="CNS",seq_freq=rep(0,15-nrow(gene_CNS))))
  }
  p.value <- round(wilcox.test(gene_periphery$seq_freq, gene_CNS$seq_freq, paired = TRUE, alternative = "two.sided")$p.value, 4)
  gene_both <- rbind(gene_periphery,gene_CNS)
  gene_both$p.value <- p.value
  gene_usage <- rbind(gene_usage,gene_both)
}

# apply the Benjamini-Hochberg correction for multiple testing
gene_usage$adjusted_p.value <- p.adjust(gene_usage$p.value, method = "BH")


# Generate colors based on sorted gene names
colours <- sapply(sort(unique(gene_usage$gene)), function(x) {
  gene_family <- substr(x, 1, 5)
  return(family_colors[gene_family])
})

# Plot the gene usage for periphery and CNS compartment for all genes
p4 <- ggplot(gene_usage, aes(x = gene,y=seq_freq, fill = compartment)) +
  scale_fill_manual(values = c("periphery" = "seagreen", "CNS" = "skyblue")) +
  geom_boxplot(width = 0.8, outlier.shape = 1, outlier.size = 1) +
  coord_cartesian(ylim = c(0,0.27)) +
  labs(title = "Gene usage for periphery and CNS compartment", x = "IGHV gene", y = "Frequency") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,colour=colours),plot.title = element_text(hjust = 0.5)) + 
  # add ** to the plot if the adjusted p-value is less than 0.05 and *** if it is less than 0.01 and **** if it is less than 0.001
  geom_text(aes(label = sapply(adjusted_p.value,get_stars), x = gene, y = 0.27), size = 3.25) +
  theme(legend.position = "top")

p4
# save the plot as pdf
ggsave("/mnt/susanne/datasets/gene_usage/gene_usage_periphery_CNS_genes.pdf", p4, width = 8.27, height = 11.69)

# plot significant genes only
gene_usage_significant <- gene_usage %>% filter(adjusted_p.value < 0.05)

colours <- sapply(sort(unique(gene_usage_significant$gene)), function(x) {
  gene_family <- substr(x, 1, 5)
  return(family_colors[gene_family])
})

p5 <- ggplot(gene_usage_significant, aes(x = gene,y=seq_freq, fill = compartment)) +
  geom_boxplot(width = 0.8, outlier.shape = 1, outlier.size = 1) +
  scale_fill_manual(values = c("periphery" = "seagreen", "CNS" = "skyblue")) +
  coord_cartesian(ylim = c(0,0.065)) +
  labs(title = "Gene usage for periphery and CNS compartment", x = "IGHV gene", y = "Frequency") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,colour=colours),plot.title = element_text(hjust = 0.5)) + 
  # add ** to the plot if the adjusted p-value is less than 0.05 and *** if it is less than 0.01 and **** if it is less than 0.001
  geom_text(aes(label = sapply(adjusted_p.value,get_stars), x = gene, y = 0.065), size = 5) +
  theme(legend.position = c(0.9,0.8))

p5

# save the plot as pdf
ggsave("/mnt/susanne/datasets/gene_usage/gene_usage_periphery_CNS_genes_significant.pdf", p5, width = 10.41, height = 3.96)



