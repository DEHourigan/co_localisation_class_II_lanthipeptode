# Load packages
.libPaths()
.libPaths(c("/data/san/data0/users/david/rstudio/packages", .libPaths()))
newlib <- "/data/san/data0/users/david/rstudio/packages"
.libPaths()

writeFasta <- function(data, filename) {
  fastaLines <- c()
  for (rowNum in 1:nrow(data)) {
    fastaLines <- c(fastaLines, as.character(paste(">", data[rowNum, "name"], sep = "")))
    fastaLines <- c(fastaLines, as.character(data[rowNum, "seq"]))
  }
  fileConn <- file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

## C:\Users\dhourigan\OneDrive - University College Cork (2)\Desktop\main_results.csv
if (!require(data.table)) {
  install.packages("data.table", repos = "http://cran.us.r-project.org", lib = newlib)
  library(data.table)
}
if (!require(multidplyr)) {
  install.packages("multidplyr", repos = "http://cran.us.r-project.org", lib = newlib)
  library(multidplyr)
}
library(readr)
if (!require(formattable)) {
  install.packages("formattable", repos = "http://cran.us.r-project.org", lib = newlib)
  library(formattable)
}
if (!require(fs)) {
  install.packages("fs", repos = "http://cran.us.r-project.org", lib = newlib)
  library(fs)
}
if (!require(dplyr)) {
  install.packages("dplyr", repos = "http://cran.us.r-project.org", lib = newlib)
  library(dplyr)
}
if (!require(ggplot2)) {
  install.packages("ggplot2", repos = "http://cran.us.r-project.org", lib = newlib)
  library(ggplot2)
}
if (!require(purrr)) {
  install.packages("purrr", repos = "http://cran.us.r-project.org", lib = newlib)
  library(purrr)
}
if (!require(ggthemes)) {
  install.packages("ggthemes", repos = "http://cran.us.r-project.org", lib = newlib)
  library(ggthemes)
}
if (!require(BiocManager)) {
  install.packages("BiocManager", repos = "http://cran.us.r-project.org", lib = newlib)
  library(BiocManager)
}
if (!require(gplots)) {
  install.packages("gplots", repos = "http://cran.us.r-project.org", lib = newlib)
  library(gplots)
}
if (!require(gridExtra)) {
  install.packages("gridExtra", repos = "http://cran.us.r-project.org", lib = newlib)
  library(gridExtra)
}
if (!require(grid)) {
  install.packages("grid", repos = "http://cran.us.r-project.org", lib = newlib)
  library(grid)
}
if (!require(ggplot2)) {
  install.packages("ggplot2", repos = "http://cran.us.r-project.org", lib = newlib)
  library(ggplot2)
}

if (!require(forcats)) {
  install.packages("forcats", repos = "http://cran.us.r-project.org", lib = newlib)
  library(forcats)
}

if (!require(tidyr)) {
  install.packages("tidyr", repos = "http://cran.us.r-project.org", lib = newlib)
  library(tidyr)
}
if (!require(dtplyr)) {
  install.packages("dtplyr", repos = "http://cran.us.r-project.org", lib = newlib)
  library(dtplyr)
}
if (!require(topGO)) {
  BiocManager::install("topGO")
  library(topGO)
}
if (!require(SparseM)) {
  install.packages("SparseM")
  library(SparseM)
}

### for topGO ensure 4.3
#Sys.getenv("PATH")
#BiocManager::install(version = "3.17")
# Function to parse the GenBank file and find the proteins
library(Biostrings)
library(GenomicRanges)
library(seqinr)
library(stringr)
library(dplyr)
library(readxl)
library(thacklr)
library(gggenomes)







# ~~~~~~~~~~~~~~~~~~~~~~~~~
# concatenate all the BGCs into one file
# ~~~~~~~~~~~~~~~~~~~~~~~~~

LanM_40k_window = fread("/data/san/data0/users/david/intelligence/lanthipeptide_rodeo_2020paper/LanM_paper_40k_window.tsv")
LanM_40k_window = LanM_40k_window %>% mutate(bgc_type = "lanthipeptide_2") %>% as_tibble() %>%
    mutate(unique_window = paste(Nucleotide_acc, start_window, sep = "_"))

DUF95_paper_40k_window = fread("/data/san/data0/users/david/intelligence/lanthipeptide_rodeo_2020paper/LanM_paper_40k_window.tsv")
DUF95_paper_40k_window$unique_window <- paste(DUF95_paper_40k_window$Nucleotide_acc, DUF95_paper_40k_window$start_window, sep = "_")
DUF95_paper_40k_window = DUF95_paper_40k_window %>% mutate(bgc_type = "circular") %>% as_tibble()


# ~~~~~~~~~~~~~~~~~~~~~~~~~
# CHECKPOINT
# ~~~~~~~~~~~~~~~~~~~~~~~~~
colnames(DUF95_paper_40k_window)

LanM_40k_window = fread("/data/san/data0/users/david/intelligence/lanthipeptide_rodeo_2020paper/LanM_paper_40k_window.tsv")
LanM_40k_window$unique_window <- paste(LanM_40k_window$Nucleotide_acc, LanM_40k_window$start_window, sep = "_")
DUF95_phyla <- LanM_40k_window %>% 
  dplyr::select(phylum) 
  # %>% distinct() %>%
  # filter(phylum %in% c("Actinomycetota","Bacillota" ))

#  
# Merge PF00365, PF00905 and PF01225 with pfam annotation
#
refseq_tsv_wPF00365_40kb = fread("/data/san/data0/users/david/intelligence/negative_dataset/databases/refseq_tsv_wPF00365_40kb_window.tsv")
refseq_tsv_wPF00905_40kb = fread("/data/san/data0/users/david/intelligence/negative_dataset/databases/refseq_tsv_wPF00905_40kb_window.tsv")
refseq_tsv_wPF01225_40kb = fread("/data/san/data0/users/david/intelligence/negative_dataset/databases/refseq_tsv_wPF01225_40kb_window.tsv")


## databases from mmseqs (clustering of all proteins within the refseq whole genome). 
refseq_pfam_annot <- fread("/data/san/data1/users/david/intelligence_database/rep_v_pfam_searchout.m8")
colnames(refseq_pfam_annot) <- c("locus_tag", "subject", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")


#left join this with the refseq_tsv_wPF00365_40kb_window.tsv
# and we can get pfam per protein in the window
PF00365_df <- left_join(refseq_tsv_wPF00365_40kb, refseq_pfam_annot, by = "locus_tag", relationship = "many-to-many") 
PF00905_df <- left_join(refseq_tsv_wPF00905_40kb, refseq_pfam_annot, by = "locus_tag", relationship = "many-to-many") 
PF01225_df <- left_join(refseq_tsv_wPF01225_40kb, refseq_pfam_annot, by = "locus_tag", relationship = "many-to-many") 


library(taxonomizr)
# wget was used in terminal ## getAccession2taxid(baseUrl='https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/')
setwd("/data/san/data0/users/david/taxonomizer") #SET UP DATABASE ONCE #
sqlFile <- "/data/san/data0/users/david/taxonomizer/accessionTaxa.sql"
  
  
# define the function
append_taxonomy <- function(df, sqlFile) {
  # append taxid
  df$taxid <- taxonomizr::accessionToTaxa(df$subject.x, sqlFile)
  
  # define the taxonomy levels
  taxa <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
  
  # append each taxonomic level
  for(taxon in taxa) {
    df[[taxon]] <- getTaxonomy(df$taxid, sqlFile, desiredTaxa = taxon)
  }
  
  return(df)
}

# control blast hits with annotation, 40kb region and taxonomy
PF00365_df <- append_taxonomy(PF00365_df, sqlFile) %>% ungroup
PF00905_df <- append_taxonomy(PF00905_df, sqlFile)
PF01225_df <- append_taxonomy(PF01225_df, sqlFile)



n_distinct(PF00365_df$Nucleotide_acc)
n_distinct(LanM_40k_window$Nucleotide_acc)



print_num_distinct <- function(df, filter_df, filter_column) {
  num_distinct_original <- df %>% dplyr::select(Nucleotide_acc) %>% n_distinct()
  num_distinct_filtered <- df %>% 
    filter(phylum %in% filter_df[[filter_column]]) %>%
    dplyr::select(Nucleotide_acc) %>%
    n_distinct()
  
  cat("\nNumber of distinct Nucleotide_acc in the original dataframe:", num_distinct_original)
  cat("\nNumber of distinct Nucleotide_acc in the filtered dataframe:", num_distinct_filtered, "\n")
}

# apply the function to each data frame
print_num_distinct(PF00365_df, DUF95_phyla, "phylum")
print_num_distinct(PF00905_df, DUF95_phyla, "phylum")
print_num_distinct(PF01225_df, DUF95_phyla, "phylum")

# 
# Just get phylum that have Lans
#
create_table <- function(df, filter_df, filter_column) { # Define the function
  result_table <- df %>% 
    filter(phylum %in% filter_df[[filter_column]]) %>%
    dplyr::select(Nucleotide_acc, subject.y) %>% 
    distinct() 

  colnames(result_table) <- c("Nucleotide_acc", "pfam")
  result_table$group <- "control"
  result_table$pfam = sub("\\..*", "", result_table$pfam)

  return(result_table)
}

# Abi from paper Diep PF02517
# 
# Tables that have including phyla 
# 
PF00365_table <- create_table(PF00365_df, DUF95_phyla, "phylum")
PF00905_table <- create_table(PF00905_df, DUF95_phyla, "phylum")
PF01225_table <- create_table(PF01225_df, DUF95_phyla, "phylum")

PF00365_table$id <- "PF00365"
PF00905_table$id <- "PF00905"
PF01225_table$id <- "PF01225"
DUF95_table <- LanM_40k_window %>% dplyr::select(Nucleotide_acc, pfam) %>% distinct()
DUF95_table$group <- "lanII bacteriocin"
DUF95_table$id <- "lan"


ctrl_v_circ = rbind(PF00365_table, DUF95_table) 
ctrl_df = rbind(PF00365_table)
write_tsv(ctrl_v_circ, "/data/san/data1/users/david/intelligence_database/ctrl_v_LanM_vs_Pks.tsv")
write_tsv(ctrl_df, "/data/san/data1/users/david/intelligence_database/LanM_ctrl_df_pfams_nucs.tsv")

pfam_descs = fread("/data/san/data0/users/david/intelligence/tables/pfam_desc.tsv",
  col.names = c("pfam","CL","desc_CL","desc_pfam","name"))


### Plot what the control df looks like
library(ggsci)
set.seed(1)
(PF00365_df %>%
filter(Nucleotide_acc %in% PF00365_table$Nucleotide_acc) %>%
  group_by(Nucleotide_acc) %>%
  dplyr::summarise(unique = n_distinct(locus_tag)) %>%
  ggplot(aes(x = "All Genomes", y = unique)) +
  geom_jitter(width = 0.3, size = 0.3, alpha = 0.3) +
  geom_violin(alpha = 0.7, fill="goldenrod") +
  theme_bw() +
  labs(x = "", y = "Protein Count", title = "PF00365 control") +
   theme(
    plot.title = element_text(size = 12),  # Adjust size as needed
    axis.title.x = element_text(size = 10),  # Adjust size as needed
    axis.title.y = element_text(size = 10),  # Adjust size as needed
    axis.text.x = element_text(size = 10),  # Adjust size as needed for x-axis labels
    axis.text.y = element_text(size = 10)   # Adjust size as needed for y-axis labels
  ) ) %>%
ggsave(filename = "/data/san/data0/users/david/intelligence/figures/PF00365_control_window_protein_count_LanM_september23.png", 
    width = 6, height = 10, units = "cm")

(PF00905_df %>%
filter(Nucleotide_acc %in% PF00365_table$Nucleotide_acc) %>%
  group_by(Nucleotide_acc) %>%
  dplyr::summarise(unique = n_distinct(locus_tag)) %>%
  ggplot(aes(x = "All Genomes", y = unique)) +
  geom_jitter(width = 0.3, size = 0.3, alpha = 0.3) +
  geom_violin(alpha = 0.7, fill="goldenrod") +
  theme_bw() +
  labs(x = "", y = "Protein Count", title = "PF00905 (PBP) control") +
   theme(
    plot.title = element_text(size = 12),  # Adjust size as needed
    axis.title.x = element_text(size = 10),  # Adjust size as needed
    axis.title.y = element_text(size = 10),  # Adjust size as needed
    axis.text.x = element_text(size = 10),  # Adjust size as needed for x-axis labels
    axis.text.y = element_text(size = 10)   # Adjust size as needed for y-axis labels
  ) ) %>%
ggsave(filename = "/data/san/data0/users/david/intelligence/figures/PF00905_df_control_window_protein_count_LanM_september23.png", 
    width = 6, height = 10, units = "cm")

(PF01225_df %>%
filter(Nucleotide_acc %in% PF00365_table$Nucleotide_acc) %>%
  group_by(Nucleotide_acc) %>%
  dplyr::summarise(unique = n_distinct(locus_tag)) %>%
  ggplot(aes(x = "All Genomes", y = unique)) +
  geom_jitter(width = 0.3, size = 0.3, alpha = 0.3) +
  geom_violin(alpha = 0.7, fill="goldenrod") +
  theme_bw() +
  labs(x = "", y = "Protein Count", title = "PF01225 (MurD) control") +
   theme(
    plot.title = element_text(size = 12),  # Adjust size as needed
    axis.title.x = element_text(size = 10),  # Adjust size as needed
    axis.title.y = element_text(size = 10),  # Adjust size as needed
    axis.text.x = element_text(size = 10),  # Adjust size as needed for x-axis labels
    axis.text.y = element_text(size = 10)   # Adjust size as needed for y-axis labels
  ) ) %>%
ggsave(filename = "/data/san/data0/users/david/intelligence/figures/PF01225_df_control_window_protein_count_LanM_september23.png", 
    width = 6, height = 10, units = "cm")



##
# CHECKPOINT 2
##
ctrl_v_circ = fread("/data/san/data1/users/david/intelligence_database/ctrl_v_LanM_vs_Pks.tsv")


# SET UP DUF95/Circular Bacteriocin TABLE
create_contingency_table <- function(df, DUF95_df, control_name) {
  DUF95_table <- DUF95_df %>% 
    dplyr::select(Nucleotide_acc, pfam) %>% 
    distinct()
  DUF95_table$group <- "lanthipeptide II"

  # Select the relevant columns from df and add a group column
  df <- df %>%
    dplyr::select(Nucleotide_acc, pfam) %>%
    distinct()
  df$group <- control_name

  # make contingency table for pfams
  circ_v_control_df <- rbind(DUF95_table, df) %>% 
    distinct() 
  # Filter out rows with complete cases
  circ_v_control_df <- circ_v_control_df[complete.cases(circ_v_control_df), ]  
  # Replace NA values
  circ_v_control_df[is.na(circ_v_control_df)] <- "Unknown"

  contingency_table <- table(circ_v_control_df$group, circ_v_control_df$pfam)
  total_circular_bacteriocin <- sum(contingency_table[1, ])
  total_control <- sum(contingency_table[2, ])
  
  cat("\nTotal for lanII bacteriocin:", total_circular_bacteriocin)
  cat("\nTotal for", control_name, ":", total_control, "\n")
  
  return(contingency_table)
}

# Apply the function to each data frame
cont_control_PF00365 <- create_contingency_table(PF00365_table, DUF95_paper_40k_window, "PF00365")
cont_control_PF00905 <- create_contingency_table(PF00905_table, DUF95_paper_40k_window, "PF00905")
cont_control_PF01225 <- create_contingency_table(PF01225_table, DUF95_paper_40k_window, "PF01225")

head(cont_control_PF01225)

# ~~
# Counts vs a Random distribution of counts i.e. is there more 
# ~~

fwrite(tablex, "/data/san/data0/users/david/intelligence/tables/lanthipeptide_pfam_table.tsv")


tablex = fread("/data/san/data0/users/david/intelligence/tables/lanthipeptide_pfam_table.tsv")

LanMin = DUF95_paper_40k_window %>%
	dplyr::select(Nucleotide_acc, pfam) %>%
	mutate(group = "lanthipeptide",
		id = "lanthipeptide II") %>%
	filter(grepl("^PF", pfam))
tablex = rbind(LanMin,PF00365_table) %>% dplyr::select(-id)

plot_size_of_pfam_pools_df = tablex %>% 
  group_by(group) %>%
  summarize(total = n())



# Calculate total counts and counts per group
counts <- tablex %>%
  group_by(pfam) %>%
  summarize(
    total_genes = n(),
    genes_in_lanthipeptide = sum(group == "lanthipeptide"),
    genes_in_control = sum(group == "control")
  ) %>%
  filter(pfam != "") %>%
  arrange(desc(total_genes))

# Calculate fractions for each group
N <- sum(counts$total_genes)
n_lanthipeptide <- sum(counts$genes_in_lanthipeptide)
n_control <- sum(counts$genes_in_control)
f_lanthipeptide <- n_lanthipeptide / N
f_control <- n_control / N

# Calculate expected numbers and standard deviations for each group
counts <- counts %>%
  mutate(
    m_prime_lanthipeptide = total_genes * f_lanthipeptide,
    s_prime_lanthipeptide = sqrt(total_genes * f_lanthipeptide * (1 - f_lanthipeptide)),
    m_prime_control = total_genes * f_control,
    s_prime_control = sqrt(total_genes * f_control * (1 - f_control))
  )

# Calculate Z-scores for both groups
significant_threshold <- 0.05

counts <- counts %>%
  mutate(
    Z_score_lanthipeptide = (genes_in_lanthipeptide - m_prime_lanthipeptide) / s_prime_lanthipeptide,
    Z_score_control = (genes_in_control - m_prime_control) / s_prime_control,
    p_value = pnorm(-abs(Z_score_lanthipeptide)),  # two-tailed test
    bonferroni_p_value = pmin(1, p_value * n()),
    p_value_control = pnorm(-abs(Z_score_control)),  # two-tailed test
    bonferroni_p_value_control = pmin(1, p_value_control * n()),
    significant = bonferroni_p_value < significant_threshold
  ) %>% 
  filter(!is.na(pfam))

# View results
print(counts)

counts <- counts %>%
  left_join(., pfam_descs) 

# Assuming Z-score and significance are appropriately calculated for both groups
histogram_z_scores <- ggplot(counts, 
    aes(x = Z_score_lanthipeptide, fill = significant)) +
  geom_histogram(bins = 100, alpha = 0.6, color = "black", size=0.1) +
  #geom_histogram(aes(x = Z_score_control, fill = significant), bins = 100, alpha = 0.6, color = "darkgreen") +
  labs(x = "Z-score", y = "Frequency", title = "Histogram of Z-scores for Pfam Families") +
  scale_fill_manual(values = c("grey", "red", "darkgreen")) +
  theme_classic(base_size = 8)

ggsave(filename = "/data/san/data0/users/david/intelligence/figures/histogram_z_scores.png", 
  plot = histogram_z_scores, width = 10, height = 6, units = "cm")



expected_vs_observed_plot = ggplot(counts) +
   geom_point(aes(x = m_prime_control, y = m_prime_control,
   fill = "Expected",
    ), 
    alpha = 0.3,
    shape=21) +
  geom_point(aes(x = m_prime_control, y = genes_in_control,
    fill = "Control"), shape=21,
    alpha = 0.3) +
  geom_point(aes(x = m_prime_lanthipeptide, y = genes_in_lanthipeptide, 
    fill = "Lanthipeptide"), shape=21,) +
  scale_fill_manual("", values = c("Expected" = "darkgrey", 
                                "Control" = "goldenrod", 
                                "Lanthipeptide" = "red")) +               
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
    coord_cartesian(clip = "off") +
  labs(x = "Expected Number of Pfams", y = "Observed Number of Pfams",
       title = "Expected vs. Observed Pfams in 'lanthipeptide'") +
  theme_classic(base_size = 8) +
  geom_text(data = subset(counts, (name == "ABC transporter" & significant == TRUE) | 
    name == "Lanthionine synthetase C-like protein"),
             aes(x = m_prime_lanthipeptide, y = genes_in_lanthipeptide, label = name),
             color = "black",
             fontface = "bold",
             size=2,
             vjust = -0.8,
             position = position_nudge(y = 1.5))  + # Labels for "ABC transporter"
    geom_text(data = subset(counts, (name == "ABC transporter" & significant == TRUE) | 
    name == "Lanthionine synthetase C-like protein"), 
             aes(x = m_prime_control, y = genes_in_control, label = name),
             color = "black",
             fontface = "bold",
             size=2,
             vjust = -0.8,
             position = position_nudge(y = 1.5)) 

ggsave(filename = "/data/san/data0/users/david/intelligence/figures/expected_vs_observed_plot.png", 
  plot = expected_vs_observed_plot, width = 10, height = 7, units = "cm")


subset_count = counts %>%
  filter(m_prime_lanthipeptide < 500 &  m_prime_control < 500) %>%
  filter(genes_in_lanthipeptide > m_prime_lanthipeptide) %>%
  filter(bonferroni_p_value < 0.05)

expected_vs_observed_plot_subset = ggplot(subset_count) +
  geom_point(aes(x = m_prime_lanthipeptide, y = genes_in_lanthipeptide),
    shape=21,
    color = "black",
    # alpha = (1- p_value),
    fill="red") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  labs(x = "Expected Number of Pfams", y = "Observed Number of Pfams", 
       title = "Expected vs. Observed Pfams in 'lanthipeptide'")  +
  theme_classic(base_size = 8) +
  geom_label(data = subset(subset_count, (name == "Competence-damaged protein") | 
    desc_pfam == "Methylase_S" | desc_pfam == "Methyltransf_25" |
    desc_pfam == "Mersacidin"), 
             aes(x = m_prime_lanthipeptide, y = genes_in_lanthipeptide, label = desc_pfam),
             color = "black", 
            #  fontface = "italic",
             size=2,
            #  hjust=-.011,
            #  vjust = -0.8,
             position = position_nudge(y = 1.5)) +
  geom_label(data = subset(subset_count, desc_pfam == "AbiJ_NTD3"), 
             aes(x = m_prime_lanthipeptide, y = genes_in_lanthipeptide, label = desc_pfam),
             color = "black", 
             size=2,
             position = position_nudge(y = 120)) +
  geom_label(data = subset(subset_count, desc_pfam == "MqsA_antitoxin"), 
             aes(x = m_prime_lanthipeptide, y = genes_in_lanthipeptide, label = desc_pfam),
             color = "black", 
             size=2,
             position = position_nudge(y = 180)) +
  geom_label(data = subset(subset_count, desc_pfam == "Eco57I"), 
             aes(x = m_prime_lanthipeptide, y = genes_in_lanthipeptide, label = desc_pfam),
             color = "black", 
             size=2,
             position = position_nudge(y = 260))

ggsave(filename = "/data/san/data0/users/david/intelligence/figures/expRed circcted_vs_observed_plot_subset.png", 
  plot = expected_vs_observed_plot_subset, width = 10, height = 7, units = "cm")





# ~~
# Q-Q plot
# ~~
# Filter data for lanthipeptide group and count Pfam families
lanthipeptide_counts <- table(tablex %>% filter(group == "lanthipeptide") %>% .$pfam)
control_counts <- table(tablex %>% filter(group == "control") %>% .$pfam)
lanthipeptide_counts <- as.numeric(lanthipeptide_counts)
control_counts <- as.numeric(control_counts)


qqplot(lanthipeptide_counts, control_counts, main = "Q-Q Plot", xlab = "Lanthipeptide", ylab = "Control")
abline(0, 1, col = "red")
length(lanthipeptide_counts)
length(control_counts)

# If the lengths are not the same, pad the shorter one with zeros
max_length <- max(length(lanthipeptide_counts), length(control_counts))

lanthipeptide_counts <- c(lanthipeptide_counts, rep(0, max_length - length(lanthipeptide_counts)))
control_counts <- c(control_counts, rep(0, max_length - length(control_counts)))

# Perform the Q-Q plot
qqplot(lanthipeptide_counts, control_counts, main = "Q-Q Plot", xlab = "Lanthipeptide", ylab = "Control")
abline(0, 1, col = "red")


# ~~
# Pearson correlation
# ~~
# null hypothesis is that there is no correlation between the two groups
contingency_table <- table(tablex$group, tablex$pfam)
chi_results_mc <- chisq.test(contingency_table, simulate.p.value = TRUE, B = 2000)
print(chi_results_mc)
colnames(tablex) 
# reject the null hypothesis meaning there is a correlation between pfams


cramers_v <- sqrt(chi_results_mc$statistic / (sum(contingency_table) * (min(dim(contingency_table)) - 1)))

# Print Cramér's V
print(cramers_v)
print("this value of 0.6892283 means there is a  association between group and Pfam")







# ~~
# stats on data
# ~~
def = fread("/data/san/data0/users/david/intelligence/tables/phage_defense.txt", header=FALSE, 
  col.names = c("pfam"))
pfams_proteins = LanM_40k_window %>% 
  dplyr::select(-Query) %>%
  distinct()

# what are the highest pfams
pfams_highest = pfams_proteins %>%
  filter(pfam %in% def$pfam) %>%
  dplyr::select(pfam,  desc,Protein_acc) %>%
  group_by(pfam, desc) %>%
  summarise(n = n()) %>%
  arrange(desc(n))
fwrite(pfams_highest, "/data/san/data0/users/david/intelligence/tables/pfams_highest.tsv")

def_red = tail(pfams_highest, 104)

pfams_proteins %>%
  filter(pfam %in% def_red$pfam) %>%
  distinct() %>%
  group_by(Nucleotide_acc, genus,species,"Genus/Species", unique_window) %>%
  summarize(n = n()) %>%
  arrange(desc(n)) %>% n_distinct()


pfams_def_genus_counts = pfams_proteins %>%
  filter(pfam %in% def_red$pfam) %>%
  distinct() %>%
  group_by(pfam, desc, genus) %>%
  summarize(n = n()) %>%
  arrange(desc(n))
fwrite(pfams_def_genus_counts, "/data/san/data0/users/david/intelligence/tables/pfams_def_genus_counts.tsv")


pfams_proteins %>% 
  distinct() %>%
  filter(Nucleotide_acc == "NZ_CP007513.1") %>%
  filter(pfam %in% def_red$pfam)

#### what have more 2 or more pfams from defence
pfams_proteins %>%
  filter(pfam %in% def_red$pfam) %>%
  distinct() %>%
  group_by(Nucleotide_acc, genus,species,"Genus/Species") %>%
  summarize(n = n()) %>%
  filter(n > 4) %>%
  arrange(desc(n))  %>%
  distinct(Nucleotide_acc) 

pfams_proteins %>%
  filter(pfam %in% def_red$pfam) %>%
  distinct() %>%
  group_by(genus) %>%
  summarize(n = n()) %>%
  filter(n > 4) %>%
  distinct(genus) 

coloc_df = pfams_proteins %>% 
  dplyr::select(Nucleotide_acc, genus, species, Protein_acc, pfam) %>% 
  distinct()

coloc_df %>%
  filter(pfam %in% def_red$pfam) %>%
  # merge pfams if they on the same protein
  group_by(Nucleotide_acc, Protein_acc) %>%
  summarise(pfams = paste(pfam, collapse = ";")) %>%
  group_by(Nucleotide_acc) %>%
  summarise(n = n()) %>%
  arrange(desc(n))


#### Bacteriocin
bacteriocins = pfams_proteins %>% 
  filter(pfam == "PF10439") %>%
  dplyr::select(Protein_acc) %>%
  distinct()

fwrite(bacteriocins, "/data/san/data0/users/david/intelligence/tables/bacteriocins.tsv")

#### COMPETENCEW
competence = c("PF07508","PF07508","PF05952","PF00154")
pfams_proteins %>%
  filter(pfam %in% competence) %>%
  distinct() %>%
  group_by(Nucleotide_acc, genus,species,"Genus/Species", unique_window) %>%
  summarize(n = n()) %>%
  arrange(desc(n)) %>% n_distinct()

pfams_competence_genus_counts = pfams_proteins %>%
  filter(pfam %in% competence) %>%
  distinct() %>%
  group_by(pfam, desc, genus) %>%
  summarize(n = n()) %>%
  arrange(desc(n))
fwrite(pfams_competence_genus_counts, "/data/san/data0/users/david/intelligence/tables/pfams_competence_genus_counts.tsv")

sum(is.na(tablex))
# Check unique values
length(unique(tablex$Nucleotide_acc))
length(unique(tablex$pfam))
length(unique(tablex$group))

# View structure of the dataframe
str(tablex)

# Frequency counts
table(tablex$pfam)
table(tablex$group)

# Cross-tabulation
table(tablex$group, tablex$pfam)

# Cross-tabulation
tbl <- table(tablex$group, tablex$pfam)

# Chi-Square Test of Independence with Monte Carlo simulation
chi_results_mc <- chisq.test(tbl, simulate.p.value = TRUE, B = 2000)
print(chi_results_mc)
print("These counts are not independent and reject the null hypotheis")

# Using expected frequencies from the Monte Carlo Chi-Square test
expected_mc <- chi_results_mc$expected

# Printing expected frequencies
print(expected_mc)

# Over-representation check (observed > expected)
over_represented <- tbl > expected_mc
print(over_represented)

# Reshape the data for visualization and further analysis
melted_data <- melt(tbl)
colnames(melted_data) <- c("group", "pfam", "count")
melted_expected <- melt(expected_mc)
colnames(melted_expected) <- c("group", "pfam", "expected")

# Merging datasets
comparison <- merge(melted_data, melted_expected, by=c("group", "pfam"))
comparison <- comparison %>%
  mutate(ratio = count / expected)

# remove pfams with NA
comparison <- comparison %>%
	left_join(., pfam_descs) %>%
	filter(!is.na(desc)) 

# Plotting the ratio of observed to expected counts, with colors by group
ratio_plot_observed_vs_expected <- ggplot(comparison, aes(x=reorder(pfam, -ratio), y=ratio, fill=group)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_manual(values = c("goldenrod","red")) +
  theme(axis.text.x = element_blank(),
    axis.text.x.top = element_blank(),
  	axis.text.x.bottom = element_blank(),
    axis.title.x = element_text()) +  # Remove Pfam names on x-axis but keep the label
  theme_classic() + 
  xlab("Pfam") +
  labs(title="Ratio of Observed to Expected Counts")

ggsave(filename = "/data/san/data0/users/david/intelligence/figures/ratio_plot_observed_vs_expected.png", 
	plot = ratio_plot_observed_vs_expected, 
  width = 15, height = 10, units = "cm")

comparison_filter = comparison %>% filter(pfam %in% PF00365_fishers_LanM_rodeo$pfam)

comparison %>% 
	filter(ratio > 1) %>%
	filter(group == "lanthipeptide") %>%
	arrange(desc(ratio)) %>% 
	filter(pfam %in% c("PF05016", "PF01420","PF05147"))



# Plotting the dot plot with four regressions
# dot_plot <- ggplot(comparison, aes(x = ratio, y = count, group = group)) +
#   geom_smooth(method = "lm", formula = y ~ x, color = "#ff0000") +
#   geom_smooth(aes(y = expected), method = "lm", formula = y ~ x, color = "black") +
#   scale_x_continuous(name = "Ratio") +
#   scale_y_continuous(name = "Count") +
#   labs(title = "Dot Plot with Four Regressions") +
#   theme_classic() +
#   facet_wrap(~group, scales = "free")

# # Save the plot as an image
# ggsave(filename = "/data/san/data0/users/david/intelligence/figures/dot_plot_test.png",
#        plot = dot_plot,
#        width = 20, height = 10, units = "cm")
comparison[comparison$group == "lanthipeptide" & comparison$count > 50 & comparison$expected > 50, ]
dotplotdf = comparison[comparison$count > 50 & comparison$expected > 50, ]

dot_plot2 <- ggplot(dotplotdf ,aes(x = expected, y = count, group = group)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_smooth(method = "lm", formula = y ~ x, color = "#ff0000") +
  geom_smooth(aes(y = expected), method = "lm", formula = y ~ x, color = "black") +
  scale_x_continuous(name = "expected") +
  scale_y_continuous(name = "observed") +
  theme_classic() 

ggsave(filename = "/data/san/data0/users/david/intelligence/figures/dot_plot_test2.png",
       plot = dot_plot2,
       width = 10, height = 10, units = "cm")



strep1= fread("/data/san/data1/users/david/db/streptococcus_pangenome/refseq/streptococcus_pirate_out/PIRATE.gene_families.tsv")
strepdf = strep1[,0:20]

aggregated_counts <- strepdf %>%
  group_by(gene_family) %>%
  summarize(count = sum(number_genomes), .groups = 'drop')

# Print the aggregated counts
print(aggregated_counts)

# Pivot the data to wide format
pivoted_data <- strepdf %>%
  pivot_wider(names_from = gene_family, values_from = number_genomes, values_fill = list(number_genomes = 0))

# Print the column names of the pivoted data
print(colnames(pivoted_data))

pa_st = fread("/data/san/data1/users/david/db/streptococcus_pangenome/refseq/streptococcus_pirate_out/presence-absence.tsv")
# Load necessary libraries
# Load necessary libraries
library(data.table)
library(gplots)

# Assuming you have your data in a CSV file
# df <- fread("your_file.csv")
# For the sake of example, I'm assuming your data is already loaded into df

# Remove the "Gene" column and convert the rest of the data to a matrix
gene_presence_matrix <- as.matrix(pa_st[, -1, with=FALSE])

# Transpose the matrix so that genes are in the columns and genomes are in the rows
transposed_gene_matrix <- t(gene_presence_matrix)

# Calculate the co-occurrence matrix for genes
gene_co_occurrence_matrix <- transposed_gene_matrix %*% gene_presence_matrix

# Set the diagonal to NA to ignore self-co-occurrences
diag(gene_co_occurrence_matrix) <- NA

# Visualize the co-occurrence matrix using a heatmap
heatmap.2(gene_co_occurrence_matrix, 
          trace = "none", 
          col = colorRampPalette(c("white", "blue"))(100), 
          margins = c(10, 10),
          cexRow = 0.5, cexCol = 0.5,
          main = "Gene Co-occurrence Matrix",
          na.color = "grey")

# If you want to extract specific pairs of co-occurring genes, you can do so by filtering the matrix
threshold <- 50 # Set a threshold for co-occurrence count
co_occurring_genes <- which(gene_co_occurrence_matrix > threshold, arr.ind=TRUE)

# Convert to a data frame for easier reading
co_occurring_genes_df <- data.frame(
  Gene1 = rownames(gene_co_occurrence_matrix)[co_occurring_genes[, 1]],
  Gene2 = colnames(gene_co_occurrence_matrix)[co_occurring_genes[, 2]],
  CoOccurrenceCount = gene_co_occurrence_matrix[co_occurring_genes]
)

# Display the co-occurring genes
print(head(co_occurring_genes_df))





# ~~
# Fisher's exact test
# ~~
perform_test <- function(df, control_group) {
  
  # Calculating total counts for each group
  total_circular_bacteriocin <- sum(df["lanthipeptide II",])
  total_control <- sum(df[control_group,])
  
  # Defining a new function to use with `apply`
  test_func <- function(column) {
    present_counts <- as.numeric(column)
    absent_counts <- c(total_circular_bacteriocin, total_control) - present_counts
    variable_table <- rbind(present_counts, absent_counts)
  
    # Perform Fisher's exact test
    test_result <- fisher.test(variable_table)
  
    # Calculate difference in proportions
    proportion_diff <- present_counts[1]/total_circular_bacteriocin - present_counts[2]/total_control
  
    return(data.frame(variable = names(column), p_value = test_result$p.value, proportion_diff = proportion_diff))
  }
  
  # Apply the function to all variables and save the p-values and proportion differences
  results <- apply(df, 2, test_func)
  
  return(results)
}

# Running the function on each dataframe
results_PF00365_fishers <- perform_test(cont_control_PF00365, "PF00365")


# Transform each list
results_PF00365_fishers_df <- do.call(rbind, lapply(results_PF00365_fishers, function(x) {
  x$control <- "PF00365"
  x
}))



all_fishers_df <- results_PF00365_fishers_df %>%
  rownames_to_column(var = "pfam")
 
all_fishers_df_clean = all_fishers_df %>% 
  # filter(grepl("^PF", pfam)) %>%
  dplyr::select(pfam, p_value, proportion_diff, control) %>% 
  distinct()

all_fishers_df_clean$pfam <- gsub("\\.\\d+$", "", all_fishers_df_clean$pfam)
all_fishers_df_clean = all_fishers_df_clean %>% distinct()
# Apply Bonferroni correction
all_fishers_df_clean$adjusted_p_value <- p.adjust(all_fishers_df_clean$p_value, method = "bonferroni")

pfam_descs = DUF95_paper_40k_window %>% dplyr::select(pfam, desc) %>% distinct() 

PF00365_fishers_LanM_rodeo = left_join(all_fishers_df_clean, pfam_descs) %>%
  distinct() %>% 
  arrange(adjusted_p_value) %>%
  filter(grepl("^PF", pfam)) %>%
  filter(adjusted_p_value < 0.05) %>%
  filter(proportion_diff > 0) %>%
  filter(control == "PF00365")
fwrite(PF00365_fishers_LanM_rodeo, "/data/san/data0/users/david/intelligence/tables/PF00365_fishers_LanM_rodeo.tsv", sep = "\t")



# ~~~~~~~~~~~~~~~~~~~~~~~~
# Volcano Plot
# ~~~~~~~~~~~~~~~~~~~~~~~~
library(ggrepel)

vplot1 = PF00365_fishers_LanM_rodeo %>% 
  filter(adjusted_p_value < 0.05) %>% 
  filter(proportion_diff > 0) %>%
  left_join(., (dplyr::select(DUF95_paper_40k_window, desc, pfam) %>% 
  distinct())) %>%
  arrange(desc(adjusted_p_value)) %>% 
  filter(control == "PF00365" | control == "PF01225") %>%
  ggplot(aes(x = log10(proportion_diff), y = -log10(adjusted_p_value))) +
  geom_jitter(alpha = 0.9, size = 0.6, aes(color=control)) +
  geom_label_repel(data = . %>% filter(desc %in% c(
		"Lanthionine synthetase C-like protein",
		"Type-A lantibiotic",
     "Competence-damaged protein",
     "Bacillus competence pheromone ComX",
		"ABC transporter")), 
            aes(label = desc), vjust = 1, hjust = 1, size = 2) +
  xlab("log10(Proportion difference)") +
  ylab("-log10(Adjusted P-Value)") +
  ggtitle("Volcano Plot: LanM-associated Pfams") +
  scale_color_manual(values = c("black", "#a22d2d")) +
  theme_classic(base_size = 8, base_family = "arial") +
  theme(legend.position = "none") +
           theme(
    plot.title = element_text(size = 10),  # Adjust size as needed
    axis.title.x = element_text(size = 10),  # Adjust size as needed
    axis.title.y = element_text(size = 10),  # Adjust size as needed
    axis.text.x = element_text(size = 10),  # Adjust size as needed for x-axis labels
    axis.text.y = element_text(size = 10)   # Adjust size as needed for y-axis labels
  ) +
  coord_cartesian(clip="off") 

ggsave("/data/san/data0/users/david/intelligence/figures/LanM_controlPk_volcano_fishers_prop_diff_2controls_pfam.png", 
	vplot1, width = 10, height = 7, units = "cm")

vplot2 = PF00365_fishers_LanM_rodeo %>% 
  filter(proportion_diff < 0.0015) %>%
  filter(adjusted_p_value < 0.05) %>% 
  filter(proportion_diff > 0) %>%
  ggplot(aes(x = -log10(proportion_diff), y = -log10(adjusted_p_value))) +
  geom_jitter(alpha = 0.9, size = 0.6, aes(color=control)) +
    geom_label_repel(data = . %>% filter(desc %in% c(
      "Zeta toxin",
      "Competence-damaged protein",
      "Type I restriction modification DNA specificity domain",
      "Type I restriction and modification enzyme - subunit R C terminal",
      "Competence-damaged protein")), 
              aes(label = desc), vjust = 1, hjust = 1, size = 2, 
              position = position_nudge(x = 0.2, y = 0.2)) +
  xlab("-log10(Proportion difference)") +
  ylab("-log10(Adjusted P-Value)") +
  ggtitle("Volcano Plot: LanM-associated Pfams") +
  scale_color_manual(values = c("black", "#a22d2d")) +
  theme_classic(base_size = 8, base_family = "arial") +
  theme(legend.position = "none") +
           theme(
    plot.title = element_text(size = 10),  # Adjust size as needed
    axis.title.x = element_text(size = 10),  # Adjust size as needed
    axis.title.y = element_text(size = 10),  # Adjust size as needed
    axis.text.x = element_text(size = 10),  # Adjust size as needed for x-axis labels
    axis.text.y = element_text(size = 10)   # Adjust size as needed for y-axis labels
  ) +
  coord_cartesian(clip="off") 

ggsave("/data/san/data0/users/david/intelligence/figures/LanM_controlPk_volcano_fishers_prop_diff_pfam_zoomed.png", 
	vplot2, width = 10, height = 7, units = "cm")

# ~~~~~~~~~~~~~~~~~~~~~~~~
# Random Forrest for PFAMs vs PF00365
# ~~~~~~~~~~~~~~~~~~~~~~~~
install.packages("ranger")
install.packages("caret")
install.packages("ROCR")
install.packages("pROC")
install.packages("randomForest")
library(randomForest)
library(ranger)
library(caret)
library(ROCR)
library(pROC)


# ~~~~~~~~~~~~~~~~~~~~~~~~

PF00905_table
PF01225_table
DUF95_table


pfam_descs = DUF95_paper_40k_window %>% 
	dplyr::select(pfam, desc) %>% 
	group_by(pfam, desc) %>%
	summarise(n = n()) %>% arrange(desc(n)) %>%
	as.data.frame()

rm_these = pfam_descs %>%
	dplyr::slice(1:50) %>%
	dplyr::select(pfam) %>%
	as.vector()


top_50_pfams_per_group <- tablex %>%
  group_by(group, pfam) %>%
  summarise(count = n(), .groups = 'drop') %>%
  arrange(group, desc(count)) %>%
  group_by(group) %>%
  slice_head(n = 50) %>%
  ungroup()

# Save the top 100 Pfams to a dataframe called "remove_these"
remove_these <- top_50_pfams_per_group %>%
  arrange(desc(count)) %>%
  slice_head(n = 100)

# Print the "remove_these" dataframe
print(remove_these)

# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Run RF Model for 5 seeds 
# ~~~~~~~~~~~~~~~~~~~~~~~~~
library(ggrepel)

perform_iterations <- function(circ_v_control_wo_query, num_iterations, seed_values) {
  all_roc_data <- list()
  feature_importance_list <- list()
  for (seed in seed_values) {
    set.seed(seed)
    
    training_rows <- caret::createDataPartition(circ_v_control_wo_query$group, p = 0.7, list = FALSE)
    training_data <- circ_v_control_wo_query[training_rows, ]
    testing_data <- circ_v_control_wo_query[-training_rows, ]
    
    rf_model <- ranger::ranger(group ~ pfam, data = training_data,
                               importance = "impurity", num.threads = 32,
                               probability = TRUE)
    
    predictions <- predict(rf_model, data = testing_data, num.threads = 32, type = "response")
    prob_control <- predictions$predictions[, "control"]
    
    roc_obj <- pROC::roc(testing_data$group, prob_control)
    roc_df <- data.frame(
      FPR = roc_obj$specificities,
      TPR = roc_obj$sensitivities,
      Seed = seed
    )
    auc_val <- pROC::auc(roc_obj)
    
    all_roc_data[[as.character(seed)]] <- list(roc_df = roc_df, auc_val = auc_val)
    feature_importance_list[[as.character(seed)]] <- rf_model$variable.importance
  }
  
  return(all_roc_data)
}

circ_v_control_df <- tablex
circ_v_control_df <- na.omit(tablex)
circ_v_control_df <- circ_v_control_df %>% filter(!pfam %in% remove_these$pfam)
circ_v_control_df$pfam <- as.factor(circ_v_control_df$pfam)
circ_v_control_df$group <- as.factor(circ_v_control_df$group)

values_to_filter <- c("hypothetical", "PF00365", "PF00005",
                      "PF01225", "PF08245", "PF02875", "PF01225",
                      "PF13241", "PF04101", "PF03033", "PF00905",
                      "PF03717", "PF06271", "PF08245", "PF02875",
                      "PF03033", "PF00905", "PF03717", "PF01795",
                      "PF01098", "PF00091", "PF12327", "PF00953",
                      "PF08478", "PF07478", "PF13535", "PF14450",
                      "PF06723", "PF01820", "PF09221", rm_these)

circ_v_control_wo_query <- circ_v_control_df %>%
  filter(!(pfam %in% values_to_filter))

seed_values <- c(123, 456, 789, 987, 654)  # Add your desired seed values
roc_data <- perform_iterations(circ_v_control_wo_query, num_iterations = 5, seed_values = seed_values)
# Create a combined ROC curve plot
all_roc_df <- do.call(rbind, lapply(roc_data, function(x) x$roc_df))
average_auc <- data.frame(Seed = names(roc_data), AUC = sapply(roc_data, function(x) x$auc_val))
average_auc_subset <- average_auc[average_auc$Seed == 123, ]
plot_ROC = ggplot(all_roc_df, aes(x = FPR, y = TPR, color = as.factor(Seed))) +
  geom_line() +
  labs(x = "False Positive Rate", y = "True Positive Rate", color = "Seed") +
  ggtitle("ROC Curves - Pks") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 10),  # Adjust size as needed
    axis.title.x = element_text(size = 10),  # Adjust size as needed
    axis.title.y = element_text(size = 10),  # Adjust size as needed
    axis.text.x = element_text(size = 10),  # Adjust size as needed for x-axis labels
    axis.text.y = element_text(size = 10),   # Adjust size as needed for y-axis labels
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10)
  ) +
  scale_x_reverse() +
  geom_text_repel(data = average_auc_subset, aes(x = 0.8, y = 0.1 * as.numeric(as.factor(Seed)), 
                                              label = paste("Avg. AUC = ", round(AUC, 3))), 
                size = 3, hjust = -0.4,
                vjust=-20,
                segment.color = NA, 
                color = "black") +
  scale_color_brewer(palette = "Reds")

ggsave("/data/san/data0/users/david/intelligence/figures/LanM_bacteriocin_RF_predict_PKS_newDF.png", 
       plot_ROC, width = 10, height = 10, units = "cm")


# Extract and combine feature importance from all iterations
feature_importance <- do.call(rbind, lapply(results$feature_importance, function(x) as.data.frame(t(x))))
feature_importance <- feature_importance %>%
  rownames_to_column(var = "pfam") %>%
  gather(key = "Seed", value = "Importance", -pfam)

# Calculate average importance for each pfam
average_importance <- feature_importance %>%
  group_by(pfam) %>%
  summarise(Avg_Importance = mean(Importance, na.rm = TRUE)) %>%
  arrange(desc(Avg_Importance))

# Print or visualize the top contributing Pfams
top_contributing_pfams <- average_importance %>%
  slice_head(n = 10)  # Adjust the number as needed to see more or fewer Pfams

print(top_contributing_pfams)




# ~~~~~~~~~~~~~~~~~~~~~~~~
# VOLCANO PLOT
# ~~~~~~~~~~~~~~~~~~~~~~~~
cont_control_PF00365


# assuming `contingency_table` is your table
diff_pfam_PF00365 <- cont_control_PF00365[1,] - cont_control_PF00365[2,]
abs_diff_PF00365 <- abs(diff_pfam_PF00365)

# sorting to find the most significant differences
top_diff_PF00365 <- sort(abs_diff_PF00365, decreasing = TRUE)[1:100]


# ~~~~~~~~~~~~~~~~~~~~~~~~
# PFAM TO GO
# ~~~~~~~~~~~~~~~~~~~~~~~~
ctrl_v_circ
library(ragp)
# Specify the Bioconductor version when installing
install.packages("https://bioconductor.org/packages/release/bioc/src/contrib/topGO_2.52.0.tar.gz", repos = NULL, type = "source")
library(topGO)

ctrl_v_circ 
circular_DF_GO <- pfam2go(data_pfam = ctrl_v_circ, pfam = "pfam")
contingency_table_GO <- table(circular_DF_GO$id, circular_DF_GO$GO_acc)


# Number of top GO terms to select
top_n <- 30

# Find top N GO terms for each group
top_GO_terms <- circular_DF_GO %>%
  filter(!is.na(GO_name)) %>%
  filter(GO_name != "hypothetical") %>%
  filter(GO_name != "") %>%
  dplyr::select(GO_name, group, Nucleotide_acc, id) %>%
  group_by(id, GO_name) %>%
  summarise(total = n(), .groups = "drop") %>%
  group_by(id) %>%
  arrange(desc(total)) %>%
  top_n(top_n, total) 


# ~~~~~~~~~~~~~~~~~~~~~~~~
# GO (1) Biological Process
# ~~~~~~~~~~~~~~~~~~~~~~~~

library(GO.db)

go_ids <- as.character(circular_DF_GO$GO_acc)
# Get the GO categories (ontologies) corresponding to these IDs
go_categories <- select(GO.db, keys = go_ids, columns = "ONTOLOGY", keytype = "GOID")
colnames(go_categories) <- c("GO_acc", "category")
go_categories = distinct(go_categories)

# filter our proteins with pfams
values_to_filter

circular_DF_GO_filtered <- circular_DF_GO %>% filter(!pfam %in% values_to_filter)
circular_DF_GO_filtered = circular_DF_GO_filtered %>% left_join(., go_categories, by = "GO_acc") 

# Find top N GO terms for each group
top_GO_terms_BP <- circular_DF_GO_filtered %>%
  filter(category == "BP") %>%
  filter(!is.na(GO_name)) %>%
  filter(GO_name != "hypothetical") %>%
  filter(GO_name != "") %>%
  dplyr::select(GO_name, id, Nucleotide_acc) %>%
  group_by(id, GO_name) %>%
  summarise(total = n(), .groups = "drop") %>%
  group_by(id) %>%
  arrange(desc(total)) %>%
  top_n(top_n, total) 

# Plot
GO_BP_duf95vsPks = ggplot(top_GO_terms_BP, aes(x = reorder(GO_name, total), y = total, fill = id)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  labs(x = "", y = "Count", fill = "Group") +
  theme_bw() +
  scale_fill_manual(values=c("lan" = "#d02727", "PF00365" = "goldenrod", "PF0090" = "#686a68", "PF01225" = "#4c9f7b")) +
  ggtitle(paste("Top", top_n, "GO: Biological Process"))
ggsave("/data/san/data0/users/david/intelligence/figures/LanM_GO_BP_v_pks.png", GO_BP_duf95vsPks, width = 20, height = 26, units = "cm")




# ~~~~~~~~~~~~~~~~~~~~~~~~
# correlational networking
# ~~~~~~~~~~~~~~~~~~~~~~~~
install.packages("cooccur")
library(igraph)
library(ggraph)
library(cooccur)
library(reshape2)


matrix <- dcast(tablex, Nucleotide_acc~pfam)
matrix <- column_to_rownames(matrix, var = "Nucleotide_acc")

co <- cooccur(matrix, spp_names = TRUE)
matrix[1:5, 1:5]



# ~~~~~~~~~~~~~~~~~~~~~~~~
# correlational networking
# ~~~~~~~~~~~~~~~~~~~~~~~~
install.packages("GGally")
library(GGally)

tablex_transformed <- tablex %>%
  dplyr::count(Nucleotide_acc, pfam) %>%
  pivot_wider(names_from = pfam, values_from = n, values_fill = list(n = 0))





# ~~~~~~~~~~~~~~~~~~~~~~~~
# Look for patterns of systems    1) PHAGE DEFENSE SYSTEMS
# ~~~~~~~~~~~~~~~~~~~~~~~~
ctrl_v_circ = fread("/data/san/data1/users/david/intelligence_database/ctrl_v_circ.tsv")


library("readxl")
rsphagepaper1 <- read_excel("/data/san/data0/users/david/intelligence/rotem_sorek/aar4120_tabless1-s5.xlsx", sheet = 1)
rsphagepaper2 <- read_excel("/data/san/data0/users/david/intelligence/rotem_sorek/aar4120_tabless1-s5.xlsx", sheet = 2)
rsphagepaper3 <- read_excel("/data/san/data0/users/david/intelligence/rotem_sorek/aar4120_tabless1-s5.xlsx", sheet = 3)
rsphagepaper4 <- read_excel("/data/san/data0/users/david/intelligence/rotem_sorek/aar4120_tabless1-s5.xlsx", sheet = 4)
rsphagepaper5 <- read_excel("/data/san/data0/users/david/intelligence/rotem_sorek/aar4120_tabless1-s5.xlsx", sheet = 5)

# Create a vector with the 5 objects
list_of_tibbles  <- c(rsphagepaper1, rsphagepaper2, rsphagepaper3, rsphagepaper4, rsphagepaper5)

phage_defense_pfams1 =  rsphagepaper1 %>% unique() %>% filter(str_starts(Family, "pfam")) %>% 
    dplyr::select(Family) %>% distinct() %>% arrange(Family) %>%
    distinct() %>% 
    mutate(Family = str_replace_all(Family, "pfam", "PF"))

phage_defense_pfams2 =  rsphagepaper2 %>% unique() %>% filter(str_starts(Family, "pfam")) %>% 
    dplyr::select(Family) %>% distinct() %>% arrange(Family) %>%
    distinct() %>% 
    mutate(Family = str_replace_all(Family, "pfam", "PF"))

phage_defense_pfams3 =  rsphagepaper3 %>% unique() %>% filter(str_starts(`Anchor protein families`, "pfam")) %>% 
    dplyr::select(`Anchor protein families`) %>% distinct() %>% arrange(`Anchor protein families`) %>%
    distinct() %>% 
    mutate(Family = str_replace_all(`Anchor protein families`, "pfam", "PF")) %>% dplyr::select(Family)

phage_defense_pfams4 =  rsphagepaper4 %>% unique() %>% filter(str_starts(`Anchor protein families`, "pfam")) %>% 
    dplyr::select(`Anchor protein families`) %>% distinct() %>% arrange(`Anchor protein families`) %>%
    distinct() %>% 
    mutate(Family = str_replace_all(`Anchor protein families`, "pfam", "PF")) %>% dplyr::select(Family) %>%
    separate_rows(Family, sep = ";")


phage_defense_pfams5 =  rsphagepaper5 %>% unique() %>% filter(str_starts(`Associated domains`, "pfam")) %>% 
    dplyr::select(`Associated domains`) %>% distinct() %>% arrange(`Associated domains`) %>%
    distinct() %>% 
    mutate(`Associated domains` = str_replace_all(`Associated domains`, "pfam", "PF")) %>% 
    separate_rows(`Associated domains`, sep = ";")


phage_defense_pfams_df = rbind(phage_defense_pfams1, phage_defense_pfams2, phage_defense_pfams3, phage_defense_pfams4) 
# Create a data frame from the list
phage_defense_pfams_df = phage_defense_pfams_df %>% distinct() %>% arrange(Family) %>% 
  filter(str_detect(Family, "^PF")) %>% distinct()






########################################
# Align competence proteins
########################################
library(Biostrings)
library(msa)

input_file = "/data/san/data0/users/david/intelligence/lanthipeptide_rodeo_2020paper/PF05952.faa"
output_dir = "/data/san/data0/users/david/intelligence/figures"
proteins = readAAStringSet(input_file)

aligned_proteins = msa(proteins, type = "protein", method = "ClustalW", 
  verbose= TRUE, order= "aligned" )

# Generate output file path
output_base_name <- tools::file_path_sans_ext(basename(input_file))
output_tex = file.path(output_dir, paste0(output_base_name, "_alignment_output.tex"))

# Generate alignment output and convert to PDF
msaPrettyPrint(aligned_proteins, 
  file=output_tex, output="tex", verbose=FALSE, askForOverwrite=FALSE,
  showNames="left", showNumbering="none", showLogo="top",
  paperWidth = 12, paperHeight = 4,
               showConsensus="bottom", logoColors="hydropathy",
               shadingMode = "identical",shadingModeArg = 90)

tinytex::pdflatex(output_tex)


########################################
# Figure 4 
########################################
LanM_40k_window
defense_focus = c(#"NZ_CP019655.1",
"NZ_EQ973330.1","NZ_CAAIOF010000003.1","NZ_PIEU01000022.1","NZ_LN890331.1","NZ_PJTG01000002.1",
"NZ_LAWY01000037.1"
  )

lanthipeptide_pfam = c("PF05147","PF13575","PF14867","PF03412","PF00082","PF16934","PF14867")
dna_def_pfam = c("PF20473","PF01420","PF04851","PF02384","PF08463","PF08463","PF00270","PF01938","PF13649","PF07669","PF07669","PF20473",
  "PF02384","PF13337", "PF08665")
dna_def_prot = c("WP_110028332.1","WP_110028334.1")
library(dplyr)
library(tidyr)

fig_4_df  <- LanM_40k_window %>%
  filter(Nucleotide_acc %in% defense_focus) %>%
  mutate(color = "grey88") %>%
  dplyr::rename(seq_id = unique_window, 
                stand = dir,
                start = start,
                end = end,
                bin_id = species) %>%
  group_by(Protein_acc) %>%
  mutate(color = case_when(
    pfam %in% lanthipeptide_pfam ~ "goldenrod",
    pfam %in% dna_def_pfam ~ "red",
    Protein_acc %in% dna_def_prot ~ "red",
  )) %>%
  fill(color, .direction = "downup") %>%  # Fill the color for all rows within each group
  ungroup() %>%
  mutate(color = ifelse(is.na(color), "grey88", color)) 

fig_4_plot <- gggenomes(fig_4_df) +
  geom_seq() +
  geom_gene(aes(fill = color),
    size = 4) +
  scale_fill_identity() +
  coord_cartesian(clip = 'off') +
  scale_x_continuous(expand = c(0.2,0)) +
  geom_bin_label(size = 3) +
  geom_seq_label(size = 2, vjust = 1.5) +
  scale_x_continuous(expand=c(0.4,0.7,0.01,0.7)) +
  theme(legend.position="right") 

ggsave("/data/san/data0/users/david/intelligence/figures/fig_4_plot.png", 
  plot = fig_4_plot, 
  width = 18, height = 10, units = "cm")

########################################
# Figure 5
########################################
LanM_40k_window
competence_focus = c("NC_018081.1",
"NZ_CNVF02000013.1",
"NZ_NJFO02000003.1",
"NZ_JH792105.1",
"NZ_FUXA01000008.1",
# "NC_009328.1",
"NZ_PIJH01000018.1",
"NZ_LMBZ01000008.1"
)

lanthipeptide_pfam = c("PF05147","PF13575","PF14867","PF03412","PF00082","PF16934","PF14867","PF07730")
competence_pfam = c("PF18146","PF00154","PF05952","PF06133","PF12072","PF02464","PF12072","PF05389")


fig_5_df  <- LanM_40k_window %>%
  filter(Nucleotide_acc %in% competence_focus) %>%
  mutate(color = "grey88") %>%
  dplyr::rename(seq_id = unique_window, 
                stand = dir,
                start = start,
                end = end,
                bin_id = species) %>%
  group_by(Protein_acc) %>%
  mutate(color = case_when(
    pfam %in% lanthipeptide_pfam ~ "goldenrod",
    pfam %in% competence_pfam ~ "red",
    Protein_acc %in% competence_focus ~ "red",
  )) %>%
  fill(color, .direction = "downup") %>%  # Fill the color for all rows within each group
  ungroup() %>%
  mutate(color = ifelse(is.na(color), "grey88", color)) %>%
  mutate(bin_id = case_when(Nucleotide_acc== "NZ_LMBZ01000008.1" ~ "Solibacillus cecembensis",
    TRUE ~ bin_id))

fig_5_plot <- gggenomes(fig_5_df) %>%
  focus(pfam %in% lanthipeptide_pfam, .expand = c(10000,12000)) +
  geom_seq() +
  geom_gene(aes(fill = color),
    size = 4) +
  scale_fill_identity() +
  coord_cartesian(clip = 'off') +
  scale_x_continuous(expand = c(0.2,0)) +
  geom_bin_label(size = 3) +
  geom_seq_label(size = 2, vjust = 1.5) +
  scale_x_continuous(expand=c(0.4,0.7,0.01,0.7)) +
  theme(legend.position="right") 

ggsave("/data/san/data0/users/david/intelligence/figures/fig_5_plot.png", 
  plot = fig_5_plot, 
  width = 18, height = 10, units = "cm")


########################################
# Figure 6
########################################
LanM_40k_window
sugarfocus = c(
  "NZ_CP028837.1",
  "NC_019042.1",
  "NZ_CP008926.1",
  "NZ_AODG01000003.1"

)

lanthipeptide_pfam = c("PF05147","PF13575","PF14867","PF03412","PF00082","PF16934","PF14867","PF07730")

sugar_pfam2 = LanM_40k_window %>%
  filter(grepl("PTS", desc))  %>% 
  dplyr::select(pfam, desc) %>% distinct
sugar_pfam = c("PF00834","PF02502","PF00294","PF00294", sugar_pfam2$pfam)
sugar_prot = c("WP_086950606.1",
"WP_086950605.1")

fig_6_df  <- LanM_40k_window %>%
  filter(Nucleotide_acc %in% sugarfocus) %>%
  mutate(color = "grey88") %>%
  dplyr::rename(seq_id = unique_window, 
                stand = dir,
                start = start,
                end = end,
                bin_id = species) %>%
  group_by(Protein_acc) %>%
  mutate(color = case_when(
    pfam %in% lanthipeptide_pfam ~ "goldenrod",
    pfam %in% sugar_pfam ~ "red",
    Protein_acc %in% sugar_prot ~ "red",
  )) %>%
  fill(color, .direction = "downup") %>%  # Fill the color for all rows within each group
  ungroup() %>%
  mutate(color = ifelse(is.na(color), "grey88", color)) %>%
  mutate(bin_id = case_when(Nucleotide_acc== "NZ_LMBZ01000008.1" ~ "Solibacillus cecembensis",
    TRUE ~ bin_id))

fig_6_plot <- gggenomes(fig_6_df) %>%
  focus(pfam %in% lanthipeptide_pfam, .expand = c(10000,12000)) +
  geom_seq() +
  geom_gene(aes(fill = color),
    size = 4) +
  scale_fill_identity() +
  coord_cartesian(clip = 'off') +
  scale_x_continuous(expand = c(0.2,0)) +
  geom_bin_label(size = 3) +
  geom_seq_label(size = 2, vjust = 1.5) +
  scale_x_continuous(expand=c(0.4,0.7,0.01,0.7)) +
  theme(legend.position="right") 

ggsave("/data/san/data0/users/david/intelligence/figures/fig_6_plot.png", 
  plot = fig_6_plot, 
  width = 18, height = 10, units = "cm")
