suppressPackageStartupMessages(library("DT"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("airr"))
suppressPackageStartupMessages(library("alakazam"))
suppressPackageStartupMessages(library("shazam"))
suppressPackageStartupMessages(library("scoper"))
suppressPackageStartupMessages(library("stringr"))
library(ggplot2)
library(data.table)
library(dplyr)
library(alakazam)
library(grid)
library(mltools)
library(pheatmap)
library(forcats)

# read in .tsv files for MS patients and HC (healthy controls)
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

# create dataframe of all repertoires with status = MS for MS patients or status = healthy for HC
all_repertoires_MS <- bind_rows(input_rearr_MS)
all_repertoires_healthy <- bind_rows(input_rearr_healthy)
all_repertoires_MS$status <- "MS"
all_repertoires_healthy$status <- "healthy"
all_repertoires <- bind_rows(all_repertoires_MS,all_repertoires_healthy)

# filter dataframe to only keep heavy chain sequences
all_heavy <- all_repertoires %>% filter(locus=="IGH")

# print number of subjects and samples
n_distinct(all_heavy$sample_id)
n_distinct(all_heavy$subject_id)

#save dataframe to file for consecutive analyses
write.csv(all_heavy, "/mnt/susanne/datasets/db_heavy.csv", row.names = FALSE)
# clean up
rm(all_repertoires, all_repertoires_MS, all_repertoires_healthy, input_rearr_MS, input_rearr_healthy)
gc()

all_heavy <- read.csv("/mnt/susanne/datasets/db_heavy.csv")
# store number of sequences per subject for further analyses
sequences_per_subject <- all_heavy %>%
  select(subject_id, status, sequence_id) %>%
  group_by(subject_id, status) %>%
  summarise(n_sequences = n())
write.table(sequences_per_subject, "/mnt/susanne/datasets/number_of_sequences_per_subject.tsv", sep="\t", quote = F, row.names = F)


####################### CONVERGENCE ANALYSIS #############################################

Vgene_freq <- alakazam::countGenes(all_heavy, gene = "v_call", groups = c("subject_id","status"), mode="gene")

unique_genes <- sortGenes(unique(Vgene_freq$gene))

all_heavy$convergent_clone_id <- NULL
# create list to store individual gene dataframes
dfs <- list()

# loop over genes and identify clones
for (v_gene_subset in unique_genes) {
  all_repertoires_vsubset <- all_heavy %>% filter(grepl(v_gene_subset, v_call))
  
  df_convergence <- hierarchicalClones(all_repertoires_vsubset,
                                       threshold=0.2,
                                       method="aa",
                                       linkage="single",
                                       normalize="len",
                                       junction="junction",
                                       v_call="v_call",
                                       j_call="j_call",
                                       clone="convergent_clone_id",
                                       fields=NULL,
                                       cell_id=NULL,
                                       locus="locus",
                                       only_heavy=TRUE,
                                       split_light=FALSE,
                                       first=FALSE,
                                       cdr3=FALSE, mod3=FALSE,
                                       max_n=0,
                                       nproc=16,
                                       verbose=T, log=NULL,
                                       summarize_clones=FALSE)
  df_convergence$v_gene_subset_convergence <- v_gene_subset
  # write individual clones to file
  write.table(df_convergence, paste0("/mnt/susanne/datasets/convergence/convergent_clones_",str_replace_all(v_gene_subset, "/","-"),".tsv"),
              sep = "\t", quote = F, row.names = F)
  # create new column with overall unique clone id
  df_convergence$unique_convergent_clone_id <- paste(df_convergence$v_gene_subset_convergence, df_convergence$convergent_clone_id, sep = "_")
  # append df to list
  dfs <- c(dfs, list(df))
}
# concatenate all dataframes into one
df_convergence_all_genes <- rbindlist(dfs)
write.table(df_convergence_all_genes, "/mnt/susanne/datasets/convergence/convergent_clones_all_genes.tsv", sep = "\t", quote = F, row.names = F)

# count all convergent clones 
count_convergent_clones_all <- countClones(df_convergence_all_genes, clone = "unique_convergent_clone_id")
count_convergent_clones_subject <- countClones(df_convergence_all_genes, groups="subject_id",clone="unique_convergent_clone_id")
count_convergent_clones_status <- countClones(df_convergence_all_genes, groups="status", clone="unique_convergent_clone_id")
write.table(count_convergent_clones_all, "/mnt/susanne/datasets/convergent_clones_all.tsv", sep = "\t", quote = F, row.names = F)
write.table(count_convergent_clones_subject, "/mnt/susanne/datasets/convergent_clones_subject.tsv", sep = "\t", quote = F, row.names = F)
write.table(count_convergent_clones_status, "/mnt/susanne/datasets/convergent_clones_status.tsv", sep = "\t", quote = F, row.names = F)


