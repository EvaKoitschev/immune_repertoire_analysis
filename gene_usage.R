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

# get the gene usage for each subject and status by sequence
genes_all_subjects_by_sequence <- countGenes(all_repertoires,gene="v_call",groups=c("subject_id","status"),mode="gene")
fwrite(genes_all_subjects_by_sequence,"/mnt/susanne/datasets/gene_usage/genes_all_subjects_by_sequence.tsv", sep = "\t")
# get the gene usage for each subject and status by clone
genes_all_subjects_by_clone <- countGenes(all_repertoires,gene="v_call",groups=c("subject_id","status"),clone="clone_id",mode="gene")
fwrite(genes_all_subjects_by_clone,"/mnt/susanne/datasets/gene_usage/genes_all_subjects_by_clone.tsv", sep = "\t")

# extract clonally expanded sequences (clones with more than 5 sequences, clone_size_count >= 5)
all_repertoires_clonally_expanded <- all_repertoires[all_repertoires$clone_size_count >= 5,]
# get the gene usage for each subject and status by sequence for clonally expanded sequences
genes_clonally_expanded_by_sequence <- countGenes(all_repertoires_clonally_expanded,gene="v_call",groups=c("subject_id","status"),mode="gene")
fwrite(genes_clonally_expanded_by_sequence,"/mnt/susanne/datasets/gene_usage/genes_clonally_expanded_by_sequence.tsv", sep = "\t")
# get the gene usage for each subject and status by clone for clonally expanded sequences
genes_clonally_expanded_by_clone <- countGenes(all_repertoires_clonally_expanded,gene="v_call",groups=c("subject_id","status"),clone="clone_id",mode="gene")
fwrite(genes_clonally_expanded_by_clone,"/mnt/susanne/datasets/gene_usage/genes_clonally_expanded_by_clone.tsv", sep = "\t")

## summarise and plot the gene usage for MS and healthy subjects
genes_all_subjects_by_sequence <- fread("/mnt/susanne/datasets/gene_usage/genes_all_subjects_by_sequence.tsv", sep = "\t")

# function to get genes from v_call column (some sequences have multiple genes in the v_call column)
get_genes <- function(df){
  return(unique(df$gene) %>% filter(starts_with("IGHV")))
  }

# function: for each gene in the list, get the list of percentages of MS and healthy subjects (into new dataframe) and compute the mann-whitney u test for each gene between MS and healthy subjects


genes <- get_genes(genes_all_subjects_by_sequence)

gene_usage <- data.frame()
df <- genes_all_subjects_by_sequence

for (current_gene in genes){
  gene_MS <- df %>% filter(gene == current_gene & status == "MS") %>% select(seq_freq)
  if(nrow(gene_MS) < 27){
    gene_MS <- rbind(gene_MS,data.frame(seq_freq= rep(0,27-nrow(gene_MS))))
  }
  gene_healthy <- df %>% filter(gene == current_gene & status == "healthy") %>% select(seq_freq)
  if(nrow(gene_healthy) < 27){
    gene_healthy <- rbind(gene_healthy,data.frame(seq_freq= rep(0,27-nrow(gene_healthy))))
  }
  gene_usage <- rbind(gene_usage, data.frame(gene=current_gene,p.value=round(wilcox.test(gene_MS$seq_freq,gene_healthy$seq_freq)$p.value,3)))
}

# apply the Benjamini-Hochberg correction for multiple testing
gene_usage$p.value <- p.adjust(gene_usage$p.value, method = "BH")
# filter out genes with p.value > 0.05
gene_usage_MS_vs_healthy_filtered <- gene_usage[gene_usage$p.value < 0.05,]
# plot the gene usage for MS and healthy subjects in two boxplots per gene (one for MS and one for healthy subjects)
gene_usage_MS_vs_healthy_filtered <- melt(gene_usage_MS_vs_healthy_filtered, id.vars = "gene")
ggplot(gene_usage_MS_vs_healthy_filtered, aes(x=gene, y=value, fill=variable)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(title="Gene usage for MS and healthy subjects", x="Gene", y="Sequence frequency") + scale_fill_manual(values=c("red","blue"))








##################### PERIPHERY VS CNS ANALYSIS ############################

# create a column compartment with tissue type (periphery or CNS) based on tissue column (blood -> periphery, CSF -> CNS, cervical lymph nodes -> periphery, brain white matter -> CNS)
all_repertoires$compartment <- NA  # Initialize compartment column

if (any(all_repertoires$tissue %in% c("blood", "cervical lymph node","PBMC","peripheral blood"))) {
  all_repertoires$compartment <- "periphery"
}
if (any(all_repertoires$tissue %in% c("CSF","brain white matter"))) {
  all_repertoires$compartment <- "CNS"
}

# filter out the sequences that are not in the periphery or CNS compartment
all_repertoires_filtered_CNS <- all_repertoires[!all_repertoires$compartment == "NA",]
# filter subjects with less than 100 sequences in the CNS compartment
all_repertoires_filtered_CNS <- all_repertoires_filtered_CNS[!all_repertoires_filtered_CNS$subject_id %in% all_repertoires_filtered_CNS[all_repertoires_filtered_CNS$compartment == "CNS",.N < 100,by=subject_id]$subject_id,]

# get the gene usage for each subject and status by sequence for blood vs. CNS compartment
families_tissue_by_sequence <- countGenes(all_repertoires_filtered_CNS, gene="v_call", groups=c("subject_id","status","compartment"), mode="family")
fwrite(families_tissue_by_sequence,"/mnt/susanne/datasets/gene_usage/families_periphery_CNS_by_sequence.tsv", sep = "\t")



