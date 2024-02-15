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


####################### IDENTIFYING CONVERGENT CLONES #############################################

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
fwrite(df_convergence_all_genes, "/mnt/susanne/datasets/convergence/convergent_clones_all_genes.tsv", sep = "\t", quote = F, row.names = F)


################################# COLLAPSING CONVERGENT CLUSTERS #############################################
df_convergence_all_genes <- fread("/mnt/susanne/datasets/convergence/convergent_clones_all_genes.tsv", sep = "\t")

# get dataframe columns sequence_id, unique_convergent_clone_id
clusters <- df_convergence_all_genes %>% select(sequence_id, unique_convergent_clone_id)
fwrite(clusters,"/mnt/susanne/datasets/convergence/clusters.tsv",sep="\t")


################################ python used to identify connected components in the graph of convergent clones

# read in the output of the python script
cc <- fread("/mnt/susanne/datasets/convergence/connected_components.tsv", sep = "\t", header = F)
df_convergence_all_genes$collapsed_convergent_clone_id <- NA

#create a mapping from each unique convergent clone id to the representative collapsed convergent clone id (first clone_id)
unique_convergent_clone_ids <- lapply(strsplit(as.character(cc$V1), ","), unique)
representatives <- sapply(unique_convergent_clone_ids, `[`, 1)
flattened_ids <- unlist(unique_convergent_clone_ids)
mapping <- setNames(rep(representatives, lengths(unique_convergent_clone_ids)), flattened_ids)
#map each unique convergent clone id to the representative collapsed convergent clone id
df_convergence_all_genes$collapsed_convergent_clone_id <- mapping[df_convergence_all_genes$unique_convergent_clone_id]


################################### ANALYSIS OF CLONES ###################################

# count all convergent clones 
count_convergent_clones_all <- countClones(df_convergence_all_genes, clone = "collapsed_convergent_clone_id")
count_convergent_clones_subject_status <- countClones(df_convergence_all_genes, groups=c("subject_id","status"),clone="collapsed_convergent_clone_id")

fwrite(count_convergent_clones_all, "/mnt/susanne/datasets/convergent_clones_all.tsv", sep = "\t", quote = F, row.names = F)
fwrite(count_convergent_clones_subject_status, "/mnt/susanne/datasets/convergent_clones_subject_status.tsv", sep = "\t", quote = F, row.names = F)



########################################## ANALYSIS OF SEQUENCES ###############################################################

conv_clone_freq_all <- as.data.frame(table(count_convergent_clones_subject$unique_convergent_clone_id)) %>% arrange(desc(Freq))
conv_clone_freq_all_top <- as.character(conv_clone_freq_all$Var1[1:10])
conv_clone_freq_all_shared <- conv_clone_freq_all %>% filter(Freq > 1)# %>% pull(Var1)
n_distinct(conv_clone_freq_all_shared)


count_convergent_clones_all_shared <- count_convergent_clones_all %>% 
  filter(unique_convergent_clone_id %in% conv_clone_freq_all_shared) %>%
  left_join(conv_clone_freq_all, by = join_by(unique_convergent_clone_id == Var1)) %>%
  rename(subject_frequency = Freq)
count_convergent_clones_all_shared$rank <- c(1:nrow(count_convergent_clones_all_shared))

p1 <- ggplot(count_convergent_clones_all_shared, aes(x=rank, y=seq_count)) +
  geom_line() +
  geom_point() +
  scale_x_log10() +
  theme_bw() +
  labs(x="Convergent clusters by rank",
       y="Number of sequences",
       title=paste0("Convergent cluster size (N=",nrow(count_convergent_clones_all_shared),")"))
p1

count_convergent_clones_all_shared <- count_convergent_clones_all_shared %>% 
  arrange(desc(subject_frequency)) %>%
  mutate(rank_subj=c(1:nrow(count_convergent_clones_all_shared)))

p2 <- ggplot(count_convergent_clones_all_shared, aes(x=rank_subj, y=subject_frequency)) +
  geom_line() +
  geom_point() +
  scale_x_log10() +
  theme_bw() +
  labs(x="Convergent clusters by rank",
       y="Number of subjects per cluster",
       title=paste0("Convergent clusters (N=",nrow(count_convergent_clones_all_shared),")"))
p2

conv_clone_freq_all_top_100 <- as.character(conv_clone_freq_all$Var1[1:100])

count_convergent_clones_sampleid_foronehot <- count_convergent_clones_subject %>% 
  select(unique_convergent_clone_id, subject_id) %>%
  filter(unique_convergent_clone_id %in% conv_clone_freq_all_top_100) %>%
  mutate_at(c("subject_id"), as.factor)

count_convergent_clones_onehot <- one_hot(as.data.table(count_convergent_clones_sampleid_foronehot), cols = "subject_id") %>%
  aggregate(.~unique_convergent_clone_id, function(x)sum(x)) %>%
  tibble::column_to_rownames("unique_convergent_clone_id")
colnames(count_convergent_clones_onehot) <- sub('^subject_id_', '', colnames(count_convergent_clones_onehot))

col_annot_project <- all_heavy %>% select(subject_id, status) %>% distinct() %>%
  filter(subject_id %in% colnames(count_convergent_clones_onehot)) %>%
  arrange(desc(status)) %>%
  tibble::column_to_rownames("subject_id")

# prepare colors for annotations
mycolors <- c("#FC4E07", "#00AFBB")
names(mycolors) <- c("MS", "healthy")
mycolors <- list(status = mycolors)

count_convergent_clones_onehot_ordered <- count_convergent_clones_onehot[conv_clone_freq_all_top_100,]
count_convergent_clones_onehot_ordered <- count_convergent_clones_onehot_ordered[,rownames(col_annot_project)]

pheat<-pheatmap(as.matrix(count_convergent_clones_onehot_ordered),
                cellwidth = 3,
                cellheight = 10,
                cluster_rows = F,
                cluster_cols = F,
                color = colorRampPalette(c("white", "black"))(2),
                border_color = "gray80",
                scale="none",
                show_colnames = F,
                annotation_col = col_annot_project,
                annotation_colors = mycolors
)
pheat
