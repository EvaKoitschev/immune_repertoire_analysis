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









############## Plot gene usage for MS and healthy subjects #####################

## summarise and plot the gene usage for MS and healthy subjects
genes_all_subjects_by_sequence <- fread("/mnt/susanne/datasets/gene_usage/genes_all_subjects_by_sequence.tsv", sep = "\t")
genes_all_subjects_by_clone <- fread("/mnt/susanne/datasets/gene_usage/genes_all_subjects_by_clone.tsv",sep="\t")

# function to get genes from v_call column (some sequences have multiple genes in the v_call column)
get_genes <- function(dataframe){
  return(unique(dataframe$gene))
  }

# function: for each gene in the list, get the list of percentages of MS and healthy subjects (into new dataframe) and compute the mann-whitney u test for each gene between MS and healthy subjects
genes <- get_genes(genes_all_subjects_by_sequence)
length(genes)

gene_usage_list <- list()
df <- genes_all_subjects_by_sequence

for (current_gene in genes){
  gene_MS <- df %>% filter(gene == current_gene & status == "MS") %>% select(seq_freq)
  if(nrow(gene_MS) < 27){
    gene_MS <- rbind(gene_MS,data.frame(seq_freq=rep(0,27-nrow(gene_MS))))
  }
  gene_healthy <- df %>% filter(gene == current_gene & status == "healthy") %>% select(seq_freq)
  if(nrow(gene_healthy) < 27){
    gene_healthy <- rbind(gene_healthy,data.frame(seq_freq=rep(0,27-nrow(gene_healthy))))
  }
  # compute the Mann-Whitney U test for the current gene between MS and healthy subjects and store the gene and the p-value in gene_usage
  p.value <- round(wilcox.test(gene_MS$seq_freq,gene_healthy$seq_freq,alternative="greater")$p.value,4)
  MS_freqs = as.list(gene_MS$seq_freq)
  healthy_freqs = as.list(gene_healthy$seq_freq)
  # store current_gene and p_value in gene_usage as well as the dataframes gene_MS and gene_healthy
  gene_usage_list[[current_gene]] <- list(
    gene = current_gene,
    p.value = p.value,
    MS_freqs = MS_freqs,
    healthy_freqs = healthy_freqs
  )

}
gene_usage <- do.call(rbind, gene_usage_list)
gene_usage <- as.data.frame(gene_usage)

# apply the Benjamini-Hochberg correction for multiple testing
gene_usage$adjusted_p.value <- p.adjust(gene_usage$p.value, method = "BH")

# filter out the genes that are significantly different between MS and healthy subjects
gene_usage_MS_vs_healthy_filtered <- gene_usage[gene_usage$adjusted_p.value < 0.05,]
myplot <- list()

# Iterate through each row of gene_usage_MS_vs_healthy_filtered
for (i in 1:nrow(gene_usage_MS_vs_healthy_filtered)) {
  gene <- gene_usage_MS_vs_healthy_filtered$gene[i]
  MS_freqs <- gene_usage_MS_vs_healthy_filtered$MS_freqs[[i]]
  healthy_freqs <- gene_usage_MS_vs_healthy_filtered$healthy_freqs[[i]]
  
  # Create a data frame with both MS and healthy frequencies in long format
  plot_data <- data.frame(
    Status = rep(c("MS", "healthy"), each = length(MS_freqs)),
    Frequency = c(unlist(MS_freqs), unlist(healthy_freqs))
  )
  
  # Create the boxplot for the current gene
  myplot[[i]] <- ggplot(plot_data, aes(x = Status, y = Frequency, fill = Status)) +
    geom_boxplot() +
    labs(title=gene,x=NULL,y=NULL)
}

# combine plots into one figure with shared y axis
myplot_combined <- ggarrange(plotlist = myplot, ncol = 7, nrow = 1, common.legend = TRUE, legend = "bottom",align = "hv") + theme(axis.text.y = element_blank(), axis.title.y = element_blank())
# save the combined plot as pdf
ggsave("/mnt/susanne/datasets/gene_usage/gene_usage_MS_vs_healthy.pdf", myplot_combined, width = 10.41, height = 3.96)




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



families <- get_genes(families_tissue_by_sequence)
length(families)

gene_usage_list <- list()
for (current_gene in families){
  gene_periphery <- families_tissue_by_sequence %>% filter(gene == current_gene & compartment == "periphery") %>% select(seq_freq)
  if(nrow(gene_periphery) < 15){
    gene_periphery <- rbind(gene_periphery,data.frame(compartment="periphery",seq_freq=rep(0,15-nrow(gene_periphery))))
  }
  gene_CNS <- families_tissue_by_sequence %>% filter(gene == current_gene & compartment == "CNS") %>% select(seq_freq)
  if(nrow(gene_CNS) < 15){
    gene_CNS <- rbind(gene_CNS,data.frame(compartment="CNS",seq_freq=rep(0,15-nrow(gene_CNS))))
  }
  
  # compute the Mann-Whitney U test for the current gene between periphery and CNS compartment and store the gene and the p-value in gene_usage
  p.value <- round(wilcox.test(gene_periphery$seq_freq,gene_CNS$seq_freq,alternative="two.sided")$p.value,4)
  periphery_freqs = as.list(gene_periphery$seq_freq)
  CNS_freqs = as.list(gene_CNS$seq_freq)
  # store current_gene and p_value in gene_usage as well as the dataframes gene_periphery and gene_CNS
  gene_usage_list[[current_gene]] <- list(
    family = current_gene,
    p.value = p.value,
    periphery_freqs = periphery_freqs,
    CNS_freqs = CNS_freqs
  )
}
gene_usage <- do.call(rbind, gene_usage_list)
gene_usage <- as.data.frame(gene_usage)

# apply the Benjamini-Hochberg correction for multiple testing
gene_usage$adjusted_p.value <- p.adjust(gene_usage$p.value, method = "BH")

# plot the gene usage for periphery and CNS compartment for all gene families (significant or not)
myplot <- list()
for (i in 1:nrow(gene_usage)) {
  gene <- gene_usage$family[i]
  periphery_freqs <- gene_usage$periphery_freqs[[i]]
  CNS_freqs <- gene_usage$CNS_freqs[[i]]
  
  # Create a data frame with both periphery and CNS frequencies in long format
  plot_data <- data.frame(
    Compartment = rep(c("periphery", "CNS"), each = length(periphery_freqs)),
    Frequency = c(unlist(periphery_freqs), unlist(CNS_freqs))
  )
  
  # Create the boxplot for the current gene and add the significance level
  
  myplot[[i]] <- ggplot(plot_data, aes(x = Compartment, y = Frequency, fill = Compartment)) +
    geom_boxplot() +
    labs(title=gene,x=NULL,y=NULL) + 
    #y lim 0,1
    coord_cartesian(ylim = c(0, 1)) +
    # add p-value to the plot at x = 0.5 and y = 1 (in the plot)
    annotate("text", x = 0.5, y = 1, label = paste("p =", gene_usage$adjusted_p.value[i]))
    }
# sort plots based on title
myplot <- myplot[order(sapply(myplot, function(x) x$layers[[1]]$mapping$title))]

# combine plots into one figure with shared y axis
myplot_combined <- ggarrange(plotlist = myplot, ncol = 7, nrow = 1, common.legend = TRUE, legend = "bottom",align = "hv") + theme(axis.text.y = element_blank(), axis.title.y = element_blank())
# save the combined plot as pdf
ggsave("/mnt/susanne/datasets/gene_usage/gene_usage_periphery_CNS.pdf", myplot_combined, width = 10.41, height = 3.96)
