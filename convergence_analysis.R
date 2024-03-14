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
library(forcats)


# file paths
input_MS <- "/mnt/susanne/datasets/MS_total"
input_healthy <- "/mnt/susanne/datasets/healthy_total"

# read in data and make dataframes in AIRR format
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
# add diagnosis in status column
all_repertoires_MS$status <- "MS"
all_repertoires_healthy$status <- "healthy"
all_repertoires <- bind_rows(all_repertoires_MS,all_repertoires_healthy)

# filter dataframe to only keep heavy chain sequences
all_heavy <- all_repertoires %>% filter(locus=="IGH")

#save dataframe to file for consecutive analyses
write.csv(all_heavy, "/mnt/susanne/datasets/db_heavy.csv", row.names = FALSE)

# clean up
rm(all_repertoires, all_repertoires_MS, all_repertoires_healthy, input_rearr_MS, input_rearr_healthy)
gc()

####################### IDENTIFYING CONVERGENT CLONES #############################################

# get all unique V genes
Vgene_freq <- alakazam::countGenes(all_heavy, gene = "v_call", groups = c("subject_id","status"), mode="gene")
unique_genes <- sortGenes(unique(Vgene_freq$gene))

# initialise convergent_clone_id column
all_heavy$convergent_clone_id <- NULL

# create list to store individual gene dataframes
dfs <- list()

# loop over genes and identify convergent clusters using Scoper
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

# get dataframe columns sequence_id, unique_convergent_clone_id
clusters <- df_convergence_all_genes %>% select(sequence_id, unique_convergent_clone_id)
# save clusters to file
fwrite(clusters,"/mnt/susanne/datasets/convergence/clusters.tsv",sep="\t")


################################ !python used to identify connected components in the graph of convergent clones! ################################

# read in the output of the python script (connected components)
cc <- fread("/mnt/susanne/datasets/convergence/connected_components.tsv", sep = "\t", header = F)

# initialise new column collapsed_convergent_clone_id
df_convergence_all_genes$collapsed_convergent_clone_id <- NA

#create a mapping from each unique convergent clone id to the representative collapsed convergent clone id (first clone_id)
unique_convergent_clone_ids <- lapply(strsplit(as.character(cc$V1), ","), unique)
representatives <- sapply(unique_convergent_clone_ids, `[`, 1)
flattened_ids <- unlist(unique_convergent_clone_ids)
mapping <- setNames(rep(representatives, lengths(unique_convergent_clone_ids)), flattened_ids)
#map each unique convergent clone id to the representative collapsed convergent clone id
df_convergence_all_genes$collapsed_convergent_clone_id <- mapping[df_convergence_all_genes$unique_convergent_clone_id]

# count all convergent clusters and count convergent clusters per subject and status 
count_convergent_clones_all <- countClones(df_convergence_all_genes, clone = "collapsed_convergent_clone_id")
count_convergent_clones_subject_status <- countClones(df_convergence_all_genes, groups=c("subject_id","status"),clone="collapsed_convergent_clone_id")

fwrite(count_convergent_clones_all, "/mnt/susanne/datasets/convergence/convergent_clones_all.tsv", sep = "\t", quote = F, row.names = F)
fwrite(count_convergent_clones_subject_status, "/mnt/susanne/datasets/convergence/convergent_clones_subject_status.tsv", sep = "\t", quote = F, row.names = F)



########################################## ANALYSIS OF SEQUENCES ###############################################################

# filter count_convergent_clones_subject_status to only keep clones present in at least two MS patients
conv_clone_freq_MS <- count_convergent_clones_subject_status %>% filter(status == "MS") %>% select(subject_id,status,collapsed_convergent_clone_id)
conv_clone_freq_MS <- as.data.frame(table(conv_clone_freq_MS$collapsed_convergent_clone_id)) %>% arrange(desc(Freq))
conv_clone_freq_MS_shared <- conv_clone_freq_MS %>% filter(Freq > 1) %>% pull(Var1)
n_distinct(conv_clone_freq_MS_shared)
count_convergent_clones_MS_shared <- count_convergent_clones_subject_status %>% filter(collapsed_convergent_clone_id %in% conv_clone_freq_MS_shared)

# group by collapsed_convergent_clone_id, count number of MS patients and rank the clones by number of MS patients
conv_clone_freq_MS_shared_sorted <- count_convergent_clones_MS_shared %>% group_by(collapsed_convergent_clone_id) %>% summarise(n_MS_subjects = sum(status == "MS", na.rm = TRUE)) %>% arrange(desc(n_MS_subjects))
conv_clone_freq_MS_shared_sorted$rank <- c(1:nrow(conv_clone_freq_MS_shared_sorted))

# join with the original dataframe to get all columns
count_convergent_clones_MS_shared <- count_convergent_clones_MS_shared %>% left_join(conv_clone_freq_MS_shared_sorted, by = "collapsed_convergent_clone_id")

# plot the number of sequences per convergent cluster sorted by rank
p1 <- ggplot(conv_clone_freq_MS_shared_sorted, aes(x=rank, y=n_MS_subjects)) +
  geom_line() +
  geom_point(size=0.8) +
  scale_x_log10() +
  scale_y_continuous(breaks = seq(0, 27, by = 3)) +
  theme_bw() +
  labs(x="Convergent clusters by rank",
        y="Number of MS subjects per cluster",
        title=paste0("Convergent clusters (N=",nrow(conv_clone_freq_MS_shared_sorted),")")) +
  theme(plot.title = element_text(hjust = 0.5))
p1
ggsave(p1, filename = "/mnt/susanne/datasets/convergence/conv_clusters_ranked.png", width = 10, height = 6, units = "in", dpi = 300)


# group by collapsed_convergent_clone_id, for each group get number of MS and number of healthy subjects and calulate a Fisher's exact test
conv_clone_freq_MS_healthy <- count_convergent_clones_MS_shared %>% group_by(collapsed_convergent_clone_id) %>% summarise(n_MS_subjects=sum(status == "MS"),n_healthy_subjects = sum(status == "healthy", na.rm = TRUE)) %>% arrange(desc(n_MS_subjects))
conv_clone_freq_MS_healthy$p.value <- apply(conv_clone_freq_MS_healthy, 1, function(x) fisher.test(matrix(c(as.numeric(x[2]),as.numeric(x[3]),27-as.numeric(x[2]),27-as.numeric(x[3])), nrow=2),alternative = "greater" )$p.value)
conv_clone_freq_MS_healthy$adjusted_p.value <- p.adjust(conv_clone_freq_MS_healthy$p.value, method = "BH")

# number of significant clones
n_significant_clones <- sum(conv_clone_freq_MS_healthy$adjusted_p.value < 0.05)
# number of clusters with unadjusted p-value below 0.05
n_promising_clones <- sum(conv_clone_freq_MS_healthy$p.value < 0.05)

# filter for promising clones
conv_clone_freq_MS_healthy_promising <- conv_clone_freq_MS_healthy %>% filter(p.value < 0.05)

# get the full dataframe (df_convergence) of all significant clones (collapsed_convergent_clone_id)
promising_clones <- df_convergence_all_genes %>% filter(collapsed_convergent_clone_id %in% conv_clone_freq_MS_healthy_promising$collapsed_convergent_clone_id)

# get annotations for promising clusters
table_annotations <- promising_clones %>% group_by(collapsed_convergent_clone_id) %>% summarise(IGHV_genes = list(unique(substr(v_call, 1,8))), junction_length = list(unique(junction_length)), MS_subjects = n_distinct(subject_id[status == "MS"]),healthy_subjects = n_distinct(subject_id[status=="healthy"]) , CSF_present = ifelse(sum(tissue %in% c("CSF", "brain white matter", "choroid plexus", "pia mater")) >= 1, "YES", "NO"))

############ sequence distribution in convergent clusters ######## 

#group by collapsed_convergent_clone_id and subject_id and get the number of rows for each group
promising_clones_sequences <- promising_clones %>% group_by(collapsed_convergent_clone_id, subject_id,status) %>% summarise(n = n())

# plot the number of sequences per promising clone for each subject with MS or healthy status
p2 <- ggplot(promising_clones_sequences, aes(x=collapsed_convergent_clone_id, y=n, fill=status)) +
  geom_dotplot(binaxis='y', stackdir="center", dotsize=0.7, position=position_dodge(width = 0.7),binwidth=0.1) +
  scale_fill_manual(values=c("lavender", "lightcoral"), name="Status", labels=c("healthy donor", "MS patient")) +
  theme_bw() +
  scale_y_log10() +
  labs(x="Convergent clusters", y="Number of sequences", title=paste0("Sequences per subject in convergent clusters")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  #make legend position top
  theme(legend.position = "top") +
  annotation_logticks(base=10, sides="l") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

p2

ggsave(p2, filename = "/mnt/susanne/datasets/convergence/sequences_per_subject.png", width = 10, height = 6, units = "in", dpi = 300)

#get most prevalent junction amino acid sequences of promising clones
junctions <- promising_clones %>% group_by(collapsed_convergent_clone_id,junction_aa) %>% summarise(n = n())
most_frequent_junctions <- junctions %>% filter(n==max(n)) %>% select(collapsed_convergent_clone_id,junction_aa)


