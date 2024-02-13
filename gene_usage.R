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

##################### PERIPHERY VS CNS ANALYSIS ############################

# create a column compartment with tissue type (periphery or CNS) based on tissue column (blood -> periphery, CSF -> CNS, cervical lymph nodes -> periphery, brain white matter -> CNS)
all_repertoires$compartment <- if(all_repertoires$tissue %in% c("blood","cervical lymph nodes"), "periphery")
all_repertoires$compartment <- if(all_repertoires$tissue %in% c("CSF","brain white matter"), "CNS")

# filter out the sequences that are not in the periphery or CNS compartment
all_repertoires_filtered_CNS <- all_repertoires[!is.na(all_repertoires$compartment),]
# filter subjects with less than 100 sequences in the CNS compartment
all_repertoires_filtered_CNS <- all_repertoires_filtered_CNS[!all_repertoires_filtered_CNS$subject_id %in% all_repertoires_filtered_CNS[all_repertoires_filtered_CNS$compartment == "CNS",.N < 100,by=subject_id]$subject_id,]

# get the gene usage for each subject and status by sequence for CNS compartment
families_tissue_by_sequence <- countGenes(all_repertoires_filtered_CNS, gene="v_call", groups=c("subject_id","status","compartment"), mode="family")

