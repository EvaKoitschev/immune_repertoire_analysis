library(data.table)
library(dplyr)
library(ggplot2)

# file paths
input_MS <- "/mnt/susanne/datasets/MS_total"
input_healthy <- "/mnt/susanne/datasets/healthy_total"

# read data and make dataframes in AIRR format
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

# add diagnosis as status variable
all_repertoires_MS$status <- "MS"
all_repertoires_healthy$status <- "healthy"

# combine cohorts
all_repertoires <- bind_rows(all_repertoires_MS,all_repertoires_healthy)

# get subjects for each cohort
MS_subjects <- unique(all_repertoires %>% filter(status=="MS") %>% select(subject_id))
healthy_subjects <- unique(all_repertoires %>% filter(status=="healthy") %>% select(subject_id))


# plot the histogram of junction_length for each MS subject in one figure
p0 <- all_repertoires %>% filter(status=="MS") %>% group_by(subject_id) %>% ggplot(aes(x = junction_length, y=after_stat(density))) +
  geom_histogram(breaks=c(seq(0,150,by=1),Inf),color="black") +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(breaks = seq(0, 150, by = 20)) +
  facet_wrap(~ subject_id, scales = "free") +
  labs(title = "Junction length distribution per MS subject", x = "Junction length (nucleotides)", y = "Frequency") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 4),
        axis.text.y = element_text(size = 4),
        axis.title = element_text(size = 8),
        strip.text = element_text(size = 6),
        plot.title = element_text(hjust = 0.5, size = 10))
p0
# save as portrait A4 pdf
ggsave("/mnt/susanne/datasets/junction_region_length/junction_length_per_MS_subject.pdf",p0, width = 8.27, height = 11.69)


# plot the histogram of junction_length for each healthy subject in one figure
p1 <- all_repertoires %>% filter(status=="healthy") %>% group_by(subject_id) %>% ggplot(aes(x = junction_length,y=after_stat(density))) +
  geom_histogram(breaks=c(seq(0,150,by=1),Inf), closed="right",colour="black") +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(breaks = seq(0, 150, by = 20)) +
  facet_wrap(~subject_id, scales = "free") +
  labs(title = "Junction length distribution per healthy subject", x = "Junction length (nucleotides)", y = "Frequency") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 4),
        axis.text.y = element_text(size = 4),
        axis.title = element_text(size = 8),
        strip.text = element_text(size = 6),
        plot.title = element_text(hjust = 0.5, size = 10))
p1
# save as portrait A4 pdf
ggsave("/mnt/susanne/datasets/junction_region_length/junction_length_per_healthy_subject.pdf", width =11.69, height = 8.27)



###### plot the mean histogram for MS and healthy subjects ######

# specify breaks
breaks <- seq(0, 150, by = 1)


# function to return histogram counts for a given subject
get_hist <- function(current_subject,data){
  subject_data <- data.frame()
  subject_data <- data %>% filter(subject_id == current_subject)
  hist <- hist(subject_data$junction_length, breaks = c(breaks,Inf),plot = FALSE,freq=TRUE,right=TRUE)
  print(hist$counts)
  print(hist$density)
  return(hist$density)
}

# create dataframe for MS subjects
hist_data_MS <- data.frame(breaks=c(seq(1,150,by=1),Inf))

# apply the get_hist function to each subject and store the results in a dataframe
for (current_subject in MS_subjects$subject_id){
  print(current_subject)
  hist_data_MS[[current_subject]] <- get_hist(current_subject,all_repertoires)
}

# get mean of each bin across all subjects
hist_data_MS$mean <- rowMeans(hist_data_MS[,2:ncol(hist_data_MS)],na.rm=TRUE)

# plot the mean histogram for breaks 1-150
hist_data_MS_filtered <- hist_data_MS %>% filter(breaks <= 150)

p1 <- hist_data_MS_filtered %>% ggplot(aes(x = breaks, y = mean)) +
  geom_bar(stat = "identity", width = 1, color = "black", fill = "lightcoral") +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(breaks = seq(0, 150, by = 5)) +
  labs(title = "Junction length distribution for MS subjects", x = "Junction length", y = "Frequency") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 12))
p1

# save as landscape A5 pdf
ggsave("/mnt/susanne/datasets/junction_region_length/mean_junction_length_MS.pdf", height= 5.83, width = 8.27)

# create dataframe for healthy subjects
hist_data_healthy <- data.frame(breaks=c(seq(1,150,by=1),Inf))

# apply the get_hist function to each healthy subject and store the results in a dataframe
for (current_subject in healthy_subjects$subject_id){
  print(current_subject)
  hist_data_healthy[[current_subject]] <- get_hist(current_subject,all_repertoires)
}

# compute the mean of each bin across all subjects
hist_data_healthy$mean <- rowMeans(hist_data_healthy[,2:ncol(hist_data_healthy)],na.rm=TRUE)

# plot the mean histogram for breaks 1-150
hist_data_healthy_filtered <- hist_data_healthy %>% filter(breaks <= 150)

p2 <- hist_data_healthy_filtered %>% ggplot(aes(x = breaks, y = mean)) +
  geom_bar(stat = "identity", width = 1, color = "black", fill = "lavender") +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(breaks = seq(0, 150, by = 5)) +
  labs(title = "Mean junction length distribution for healthy subjects", x = "Junction length", y = "Frequency") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 12))
p2


# save as landscape A5 pdf
ggsave("/mnt/susanne/datasets/junction_region_length/mean_junction_length_healthy.pdf", height= 5.83, width = 8.27)


# combine mean columns for MS and healthy subjects
hist_data_MS_filtered <- hist_data_MS_filtered %>% select(breaks,mean) %>% mutate(status = "MS")
hist_data_healthy_filtered <- hist_data_healthy_filtered %>% select(breaks,mean) %>% mutate(status = "healthy")
hist_combined <- rbind(hist_data_MS_filtered,hist_data_healthy_filtered)

# function to format y labels as percentages
format_y_labels <- function(x) {
  ifelse(x < 0, scales::percent(abs(x)), scales::percent(x))
}

# plot the combined mean histogram
p3 <- hist_combined %>% filter(breaks <= 110) %>% ggplot(aes(x = breaks, y = mean, fill = status)) +
  geom_bar(stat = "identity", position = "dodge", width = 2.3, color = "black") +
  scale_fill_manual(values = c("MS" = "lightcoral", "healthy" = "#D1D1F6")) +
  scale_y_continuous(labels = format_y_labels, breaks = seq(0, 1, by = 0.05)) +
  scale_x_continuous(breaks = seq(0, 110, by = 3)) +
  labs(title = "Junction length distribution for MS and healthy subjects", x = "Junction length", y = "Frequency", fill = "Status") +
  theme_bw() +
  # add legend
  theme(legend.position = c(0.9, 0.9), legend.background = element_rect(fill = "white", colour = "black")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 12))
p3


# make a plot of sample_2 and sample 28 next to each other
sample_2 <- all_repertoires %>% filter(subject_id == "sample_2")
sample_28 <- all_repertoires %>% filter(subject_id == "sample_28")

# plot the histogram of junction_length for sample_2
p4 <- sample_2 %>% ggplot(aes(x = junction_length, y = after_stat(density))) +
  geom_histogram(breaks = c(seq(9, 100, by = 1), Inf),color = "black",size=0.05,fill="#D1D1F6") +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(breaks = seq(9, 100, by = 3)) +
  labs(title = "sample_2", x = "Junction length", y = "Frequency") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 20),
        plot.title = element_text(hjust = 0.5, size = 20))

p4

# plot the histogram of junction_length for sample_28
p5 <- sample_28 %>% ggplot(aes(x = junction_length, y = after_stat(density))) +
  geom_histogram(breaks = c(seq(9, 100, by = 1), Inf), color = "black",size=0.2,fill="lightcoral") +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(breaks = seq(9, 100, by = 3)) +
  labs(title = "sample_28", x = "Junction length", y = "Frequency") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 20),
        plot.title = element_text(hjust = 0.5, size = 20))

p5






############################# CNS vs periphery junction length comparison #############################

# Initialize compartment column
all_repertoires$compartment <- NA  

# assign compartment based on tissue type
all_repertoires$compartment[all_repertoires$tissue %in% c("blood", "cervical lymph node","PBMC","peripheral blood")] <- "periphery"
all_repertoires$compartment[all_repertoires$tissue %in% c("CSF","brain white matter","pia mater","choroid plexus")] <- "CNS"

# filter out subjects with no CNS samples
subjects_CNS <- all_repertoires %>% filter(compartment =="CNS") %>% group_by(subject_id) %>% filter(n() >0) %>% ungroup()
subjects_CNS <- unique(subjects_CNS$subject_id)
all_repertoires_CNS <- all_repertoires %>% filter(subject_id %in% subjects_CNS) %>% filter(!is.na(compartment))


# calculate the mean junction_length per compartment and plot a boxplot for each compartment
mean_junction_lengths <- all_repertoires_CNS %>% group_by(subject_id, compartment) %>% summarise(mean_junction_length = mean(junction_length))

p1 <- mean_junction_lengths %>% ggplot(aes(x = compartment, y = mean_junction_length, fill = compartment)) +
  geom_boxplot(outlier.shape = 11, width = 0.5, position = position_dodge(width = 0.5),alpha=0.8) +
  scale_fill_manual(values = c("CNS" = "#2976A6", "periphery" =  "#BF3F57")) +
  labs(title = "Mean junction lengths\n across MS subjects", x = "Compartment", y = "Average Junction length") +
  theme_bw() +
  theme(legend.position = c(0.5,0.1), legend.background = element_rect(fill = "white", colour = "black")) +
  coord_cartesian(ylim = c(35, 70))+
  theme(axis.text.x = element_text(size=10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 12),
        legend.text = element_text(size = 12))

# compute Wilcoxon signed rank test to test increased junction length in CNS
CNS_values <- mean_junction_lengths %>% filter(compartment == "CNS") %>% select(subject_id,mean_junction_length) %>% arrange(subject_id)
periphery_values <- mean_junction_lengths %>% filter(compartment == "periphery") %>% select(subject_id,mean_junction_length) %>% arrange(subject_id)
test <- wilcox.test(CNS_values$mean_junction_length,periphery_values$mean_junction_length,alternative="greater",paired=TRUE)
p_value <- test$p.value

# annotate plot
p1 <- p1 + annotate("text", x = 1.5, y = 65, label = "*", size = 8)

p1


# save as landscape A5 pdf
ggsave("/mnt/susanne/datasets/junction_region_length/junction_length_per_compartment.pdf", height= 5.83, width = 8.27)


######## plot junction length distribution for CNS and periphery ########

# get average histogram for CNS and periphery
hist_data_CNS <- data.frame(breaks=c(seq(1,150,by=1),Inf))
hist_data_periphery <- data.frame(breaks=c(seq(1,150,by=1),Inf))

# get histogram for each subject and tissue and store in dataframe
for (current_subject in subjects_CNS){
  print(current_subject)
  CNS_data <- all_repertoires_CNS %>% filter(subject_id == current_subject) %>% filter(compartment == "CNS")
  hist_data_CNS[[current_subject]] <- get_hist(current_subject,CNS_data)
  periphery_data <- all_repertoires_CNS %>% filter(subject_id == current_subject) %>% filter(compartment == "periphery")
  hist_data_periphery[[current_subject]] <- get_hist(current_subject,periphery_data)
}

# compute the mean of each bin across all subjects
hist_data_CNS$mean <- rowMeans(hist_data_CNS[,2:ncol(hist_data_CNS)],na.rm=TRUE)
hist_data_periphery$mean <- rowMeans(hist_data_periphery[,2:ncol(hist_data_periphery)],na.rm=TRUE)

# combine mean columns for CNS and periphery
hist_data_CNS_filtered <- hist_data_CNS %>% filter(breaks <= 150) %>% mutate(compartment = "CNS")
hist_data_periphery_filtered <- hist_data_periphery %>% filter(breaks <= 150) %>% mutate(compartment = "periphery")
hist_combined <- rbind(hist_data_CNS_filtered,hist_data_periphery_filtered)

# plot the combined mean histogram
p2 <- hist_combined %>% filter(breaks <= 102) %>% filter(breaks >=15) %>%
  ggplot(aes(x = breaks, y = mean, fill = compartment)) +
  geom_bar(stat="identity",position="dodge",width = 2.3, color = "black",alpha=0.8) +
  scale_fill_manual(values = c("CNS" = "#2976A6", "periphery" =  "#BF3F57")) +
  scale_y_continuous(labels=format_y_labels,breaks = seq(-1, 1, by = 0.05)) +
  scale_x_continuous(breaks = seq(15, 102, by = 3)) +
  labs(title = "Junction length distribution for CNS and periphery", 
       x = "Junction length (nucleotides)", y = "Frequency", fill="Compartment") +
  theme_bw() +
  theme(legend.position = c(0.9,0.9),legend.background = element_rect(fill = "white", colour = "black")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 12))
p2

# save as landscape A5 pdf
ggsave("/mnt/susanne/datasets/junction_region_length/mean_junction_length_CNS_periphery.pdf", height= 5.83, width = 8.27)
