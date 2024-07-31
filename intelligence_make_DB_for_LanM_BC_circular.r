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
  dplyr::select(phylum) %>% distinct()

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




### Plot what the control df looks like
library(ggsci)
set.seed(1)
(PF00365_df %>%
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



# Tables that have including phyla that have CIRC bacs
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

# Create a list of your data frames for easier iteration
list_of_tables <- list("PF00365" = cont_control_PF00365, "PF00905" = cont_control_PF00905, "PF01225" = cont_control_PF01225)
#list_of_tables <- list("PF00365" = cont_control_PF00365) %>% as.list()

# Initialize an empty data frame to hold results
results_df <- data.frame()

# # For each table, perform chi-square test and Spearman correlation
# for (control in names(list_of_tables)) {
#   # Calculate total Pfams for the group
#   total_pfams <- sum(list_of_tables[[control]]["lanthipeptide II",])
  
#   # For each column in the data frame, excluding the first (row names)
#   for (j in 2:ncol(list_of_tables[[control]])) {
#     # Check if the count for "circular bacteriocin" is more than 10
#     if (list_of_tables[[control]]["lanthipeptide II",j] > 10) {
#       # Extract the column and row names for current test
#       pfam <- colnames(list_of_tables[[control]])[j]
      
#       # Calculate pfam count for the group
#       pfam_count <- list_of_tables[[control]]["lanthipeptide II",j]
      
#       # Perform the chi-square test
#       test_result <- chisq.test(list_of_tables[[control]][,j])
      
#       # Perform Spearman correlation
#       spearm_result <- cor.test(list_of_tables[[control]]["lanthipeptide II",],
#                                 list_of_tables[[control]]["control",],
#                                 method = "spearman")
      
#       # Determine the group with higher frequency
#       group <- ifelse(list_of_tables[[control]][1,j] > list_of_tables[[control]][2,j], 
#                   rownames(list_of_tables[[control]])[1], 
#                   rownames(list_of_tables[[control]])[2])
      
#       # Append the result to the data frame
#       results_df <- rbind(results_df, 
#                           data.frame("Control" = control, 
#                                      "Pfam" = pfam, 
#                                      "Group" = group, 
#                                      "PfamCount" = pfam_count,
#                                      "TotalPfams" = total_pfams,
#                                      "P.Value" = test_result$p.value,
#                                      "Spearman.Rho" = spearm_result$estimate,
#                                      "Spearman.P.Value" = spearm_result$p.value))
#     }
#   }
# }

# # Add a new column "Significant" based on the p-value
# results_df <- results_df %>%
#   mutate(Significant = ifelse(P.Value < 0.05, "Yes", "No"))

# results_df %>% filter("Spearman.Rho" > 0.3) %>% View()
#   filter(Group == "circular bacteriocin") %>%
#   arrange(desc("Spearman.Rho")) %>%
#   View()

#   head(results_df %>% filter(Spearman.Rho < 0), 10)

# a low p-value and a high positive Spearman's rho could 
# suggest that there's a statistically significant positive 
# correlation between the counts of individual Pfams in the "circular bacteriocin" 
# group and their overall ranking.



######### 

# ctrl_v_circ$group 

# count_pfams <- function(df, id1, id2) {
#   count_table <- df %>%
#     filter(id %in% c(id1, id2)) %>%
#     group_by(pfam, id) %>%
#     summarise(count = n_distinct(Nucleotide_acc), .groups = "drop")
  
#   return(count_table)
# }

# chisq_test_pfams <- function(count_table, id1, id2) {
#   test_results <- data.frame()

#   for(pfam in unique(count_table$pfam)) {
    
#     count_id1 <- ifelse(pfam %in% count_table$pfam[count_table$id == id1], count_table$count[count_table$pfam == pfam & count_table$id == id1], 0)
#     count_id2 <- ifelse(pfam %in% count_table$pfam[count_table$id == id2], count_table$count[count_table$pfam == pfam & count_table$id == id2], 0)
    
#     contingency_table <- matrix(c(count_id1, count_id2), nrow = 2)
    
#     test_result <- tryCatch({
#       chisq.test(contingency_table)
#     }, warning = function(w) {
#       return(NULL)
#     }, error = function(e) {
#       return(NULL)
#     }, finally = {
#       NULL
#     })
    
#     if (!is.null(test_result)) {
#       test_results <- rbind(test_results, data.frame("Pfam" = pfam, "ChiSquare" = test_result$statistic, "P.Value" = test_result$p.value))
#     }
#   }
  
#   return(test_results)
# }


# # Usage:

# PF00365_circ_counts <- count_pfams(ctrl_v_circ, "lanII", "PF00365")
# PF00365_circ_chisq <- chisq_test_pfams(PF00365_circ_counts, "lanII", "PF00365")

# # Repeat for the other control groups:
# PF00905_circ_counts <- count_pfams(ctrl_v_circ, "lanII", "PF00905")
# PF00905_circ_chisq <- chisq_test_pfams(PF00905_circ_counts, "lanII", "PF00905")

# PF01225_circ_counts <- count_pfams(ctrl_v_circ, "lanII", "PF01225")
# PF01225_circ_chisq <- chisq_test_pfams(PF01225_circ_counts, "lanII", "PF01225")


# # Add an 'id' column to each data frame
# PF00365_circ_chisq$id <- "PF00365"
# PF00905_circ_chisq$id <- "PF00905"
# PF01225_circ_chisq$id <- "PF01225"

# # Combine all the results into a single data frame
# all_results_chi <- rbind(PF00365_circ_chisq, PF00905_circ_chisq, PF01225_circ_chisq)

# # Apply Bonferroni correction
# all_results_chi$P.Value.Bonferroni <- p.adjust(all_results_chi$P.Value, method = "bonferroni")

# all_results_chi %>% filter(P.Value.Bonferroni < 0.05) %>% 
#   arrange(P.Value.Bonferroni) 

# bonferoni_table_pfams = all_results_chi



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
results_PF00905_fishers <- perform_test(cont_control_PF00905, "PF00905")
results_PF01225_fishers <- perform_test(cont_control_PF01225, "PF01225")


# Transform each list
results_PF00365_fishers_df <- do.call(rbind, lapply(results_PF00365_fishers, function(x) {
  x$control <- "PF00365"
  x
}))
results_PF00905_fishers_df <- do.call(rbind, lapply(results_PF00905_fishers, function(x) {
  x$control <- "PF00905"
  x
}))
results_PF01225_fishers_df <- do.call(rbind, lapply(results_PF01225_fishers, function(x) {
  x$control <- "PF01225"
  x
}))


# Combine into a single dataframe
library(tibble)
all_fishers_df <- rbind(results_PF00365_fishers_df, results_PF00905_fishers_df, results_PF01225_fishers_df) %>%
  rownames_to_column(var = "pfam")
 
all_fishers_df_clean = all_fishers_df %>% 
  filter(grepl("^PF", pfam)) %>%
  dplyr::select(pfam, p_value, proportion_diff, control) %>% 
  distinct()

all_fishers_df_clean$pfam <- gsub("\\.\\d+$", "", all_fishers_df_clean$pfam)
all_fishers_df_clean = all_fishers_df_clean %>% distinct()
# Apply Bonferroni correction
all_fishers_df_clean$adjusted_p_value <- p.adjust(all_fishers_df_clean$p_value, method = "bonferroni")

pfam_descs = DUF95_paper_40k_window %>% dplyr::select(pfam, desc) %>% distinct() 

PF00365_fishers_LanM_rodeo = left_join(all_fishers_df_clean, pfam_descs) %>%
  distinct() %>% arrange(adjusted_p_value) %>%
  filter(grepl("^PF", pfam)) %>%
  filter(adjusted_p_value < 0.05) %>%
  filter(proportion_diff > 0) %>%
  filter(control == "PF00365") 
  fwrite(PF00365_fishers_LanM_rodeo, "/data/san/data0/users/david/intelligence/tables/PF00365_fishers_LanM_rodeo.tsv", sep = "\t")


library(formattable)
library(htmlwidgets)
library(htmltools)
library(pandoc)


options(scipen=4)
fmt_DUF95 <- formattable(PF00365_fishers_LanM_rodeo)

#htmlwidgets::saveWidget(fmt_DUF95, "/data/san/data0/users/david/intelligence/tables/LanM_rodeo_vs_control_fishers.html", selfcontained = TRUE)
widget_obj <- formattable::as.htmlwidget(fmt_DUF95)
htmlwidgets::saveWidget(widget_obj, 
           file = "/data/san/data0/users/david/intelligence/tables/LanM_rodeo_vs_control_fishers.html", 
           selfcontained = TRUE)

library(ragp)
pfam2go(PF00365_fishers_LanM_rodeo, "pfam") %>% 
  group_by(GO_name) %>% 
  #distinct(pfam) %>%
  dplyr::summarise(n = n()) %>%
  mutate(percentage = n / length(unique(PF00365_fishers_LanM_rodeo$pfam)) * 100) %>% 
  arrange(desc(percentage)) %>% View()

left_join(DUF95_PF01944_percentage_present, 
          (pfam2go(putative_circular_with_DUF95, "pfam") 
           %>% select(pfam, GO_name, GO_acc))) %>%
  unique() %>% 
  View()




# Get min and max values to normalize proportion_diff
min_value <- min(all_fishers_df_clean$proportion_diff)
max_value <- max(all_fishers_df_clean$proportion_diff)
# Normalize to 0-1
all_fishers_df_clean$normalized_proportion_diff <- (all_fishers_df_clean$proportion_diff - min_value) / (max_value - min_value)
# Scale to -1 to +1
all_fishers_df_clean$normalized_proportion_diff <- all_fishers_df_clean$normalized_proportion_diff * 2 - 1

all_fishers_df_clean %>% left_join(., pfam_descs) %>% 
  distinct() %>% 
  arrange(adjusted_p_value) %>% 
  filter(grepl("^PF", pfam)) %>%
  filter(adjusted_p_value < 0.05) %>%
  filter(proportion_diff > 0) %>%
  filter(control == "PF00365") %>%
  arrange(adjusted_p_value) %>% 
  fwrite(., "/data/san/data0/users/david/intelligence/tables/PF00365_fishers_LanM_rodeo_nrmalised.tsv", sep = "\t")

nrml_fish_lan =fread("/data/san/data0/users/david/intelligence/tables/PF00365_fishers_LanM_rodeo_nrmalised.tsv") 
nrml_fish_lan  

install.packages("htmlTable")
library(htmlTable)
html_output <- htmlTable(nrml_fish_lan)
writeLines(html_output, "/data/san/data0/users/david/intelligence/tables/PF00365_fishers_LanM_rodeo_nrmal.html")



(all_fishers_df_clean %>% 
  filter(adjusted_p_value < 0.05) %>% 
  filter(proportion_diff > 0) %>%
  left_join(., (dplyr::select(DUF95_paper_40k_window, desc, pfam) %>% 
  distinct())) %>%
  arrange(desc(adjusted_p_value)) %>% 
  filter(control == "PF00365" | control == "PF01225") %>%
  ggplot(aes(x = normalized_proportion_diff, y = -log10(adjusted_p_value))) +
  geom_jitter(alpha = 0.9, size = 0.6, aes(color=control)) +
  xlab("Normalized Proportion Difference") +
  ylab("-log10(Adjusted P-Value)") +
  ggtitle("Volcano Plot: Pfam (DUF95)") +
  scale_color_manual(values = c("black", "#a22d2d")) +
  theme_bw() +
  theme(legend.position = "none") +
           theme(
    plot.title = element_text(size = 10),  # Adjust size as needed
    axis.title.x = element_text(size = 10),  # Adjust size as needed
    axis.title.y = element_text(size = 10),  # Adjust size as needed
    axis.text.x = element_text(size = 10),  # Adjust size as needed for x-axis labels
    axis.text.y = element_text(size = 10)   # Adjust size as needed for y-axis labels
  )) %>%
ggsave("/data/san/data0/users/david/intelligence/figures/LanM_controlPk_volcano_fishers_prop_diff_2controls_pfam.png", ., width = 7, height = 7, units = "cm")




# ~~~~~~~~~~~~~~
# Find whats within red regions        5 K  region
# ~~~~~~~~~~~~~~


refseq_pfam_annot <- fread("/data/san/data1/users/david/intelligence_database/rep_v_pfam_searchout.m8")
colnames(refseq_pfam_annot) <- c("locus_tag", "subject", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
k5_circ = fread("/data/san/data0/users/david/intelligence/circular_bacteriocin_paper/DUF95_paper_5k_window.tsv")
k5_PF00365 = fread("/data/san/data0/users/david/intelligence/negative_dataset/databases/refseq_tsv_wPF00365_5kb_window.tsv")
k5_PF00365 = k5_PF00365 %>% dplyr::rename(pfam = subject.y )
# avove, renamae as need pfam in both

k5_window_PF00365 <- create_contingency_table(k5_PF00365, k5_circ, "PF00365")
k5_circ_PF00365_fishers <- perform_test(k5_window_PF00365, "PF00365")

      k5_circ_PF00365_fishers_df <- do.call(rbind, lapply(k5_circ_PF00365_fishers, function(x) {
        x$control <- "PF00365"
        x
      }))
      library(tibble)
      k5_circ_PF00365_fishers_df = k5_circ_PF00365_fishers_df %>% rownames_to_column(var = "pfam")

      k5_circ_PF00365_fishers_df_clean = k5_circ_PF00365_fishers_df %>% 
        dplyr::select(pfam, p_value, proportion_diff, control) %>% 
        distinct()

      k5_circ_PF00365_fishers_df_clean$pfam <- gsub("\\.\\d+$", "", k5_circ_PF00365_fishers_df_clean$pfam)
      k5_circ_PF00365_fishers_df_clean = k5_circ_PF00365_fishers_df_clean %>% distinct()
      # Apply Bonferroni correction
      k5_circ_PF00365_fishers_df_clean$adjusted_p_value <- p.adjust(k5_circ_PF00365_fishers_df_clean$p_value, method = "bonferroni")
      

pfam_descs = k5_circ %>% dplyr::select(pfam, desc ) %>% distinct()

left_join(k5_circ_PF00365_fishers_df_clean, pfam_descs) %>%
  distinct() %>% arrange(adjusted_p_value) %>%
  filter(adjusted_p_value < 0.05) %>%
  View()

phage_integrase_pfams = c("PF02899", "PF13495", "PF13102")

  pi5kcircnucs = k5_circ %>% 
      group_by(Nucleotide_acc) %>% 
      filter(pfam %in% phage_integrase_pfams) %>%
      dplyr::select(Nucleotide_acc) %>% distinct()

    
    
    gggenomes() +
    geom_seq(color = "black") +
    geom_seq_label(color = "black", label="") +
    geom_gene(aes(fill = color), size=4) +     
    theme(panel.background = element_rect(fill = "#EAEAEA"),
        legend.background = element_rect(fill = "EAEAEA")) +                            
    scale_fill_identity()

library(gggenomes)
library(thacklr)

k5plot = k5_circ %>% group_by(Nucleotide_acc) %>%
filter(Nucleotide_acc %in% pi5kcircnucs$Nucleotide_acc) %>%
mutate(seq_id = paste(Nucleotide_acc, start_window, sep="_")) %>%
sample_n_groups(30) %>%
ungroup() %>%
group_by(pfam) %>%
dplyr::mutate(
    type = "CDS",
    color = case_when(
      any(pfam %in% phage_integrase_pfams) ~ 'darkgreen', # phage integrase
      any(pfam %in% circular_bacteriocin_pfams) ~ '#a22a15',
      is.na(pfam) ~ "grey"
    )
  ) %>%
  ungroup() 


gggenomes(k5plot) +
    geom_seq(color = "black") +
    geom_seq_label(color = "black", label="") +
    geom_gene(size=5, aes(fill=color)) +
    theme(panel.background = element_rect(fill = "#EAEAEA"),
      legend.background = element_rect(fill = "#EAEAEA")) +
    scale_fill_identity()




# Ferric uptake regulator


# ~~~~~~~~~~~~~~~~~~~~~~~~
# Random Forrest for PFAMs vs PF00365
# ~~~~~~~~~~~~~~~~~~~~~~~~

install.packages("randomForest")
library(randomForest)
install.packages("ranger")
library(ranger)
install.packages("caret")
library(caret)
install.packages("ROCR")
library(ROCR)
install.packages("pROC")
library(pROC)


# ~~~~~~~~~~~~~~~~~~~~~~~~

PF00905_table
PF01225_table
DUF95_table


pfam_descs = DUF95_paper_40k_window %>% dplyr::select(pfam, desc) %>% distinct() 



# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Run RF Model for 5 seeds 
# ~~~~~~~~~~~~~~~~~~~~~~~~~
library(ggrepel)

perform_iterations <- function(circ_v_control_wo_query, num_iterations, seed_values) {
  all_roc_data <- list()
  
  for (seed in seed_values) {
    set.seed(seed)
    
    training_rows <- caret::createDataPartition(circ_v_control_wo_query$group, p = 0.8, list = FALSE)
    training_data <- circ_v_control_wo_query[training_rows, ]
    testing_data <- circ_v_control_wo_query[-training_rows, ]
    
    rf_model <- ranger::ranger(group ~ pfam, data = training_data,
                       importance = 'impurity', num.threads = 32, 
                       probability = TRUE)
    
    predictions <- predict(rf_model, data = testing_data, num.threads = 32, type = "response")
    prob_control <- predictions$predictions[, "control"]
    prob_circ <- predictions$predictions[, "circular bacteriocin"]
    
    roc_obj <- pROC::roc(testing_data$group, prob_control)
    roc_df <- data.frame(
      FPR = roc_obj$specificities,
      TPR = roc_obj$sensitivities,
      Seed = seed
    )
    auc_val <- pROC::auc(roc_obj)
    
    all_roc_data[[as.character(seed)]] <- list(roc_df = roc_df, auc_val = auc_val)
  }
  
  return(all_roc_data)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~
# PF00365_table
# ~~~~~~~~~~~~~~~~~~~~~~~~
rbind(DUF95_table, PF00365_table, PF00905_table, PF01225_table) %>%
  group_by(id, pfam) %>%
  summarise(count = n()) %>%
  arrange(id, desc(count)) %>%
  mutate(rank = row_number()) %>%
  filter(rank <= 40) %>% 
  left_join(., pfam_descs, by = c("pfam" = "pfam")) %>%
  View()


circ_v_control_df <- rbind(DUF95_table, PF00365_table) %>% distinct()
circ_v_control_df <- na.omit(circ_v_control_df)
circ_v_control_df$pfam <- as.factor(circ_v_control_df$pfam )
circ_v_control_df$group <- as.factor(circ_v_control_df$group )

DUF95_paper_40k_window %>% filter(genus == "Enterococcus") %>% 
View()

## Remove DUf95 and ABC transporter from the training and test dataset
  values_to_filter <- c("hypothetical","PF00365", "PF00005", 
  "PF01225", "PF08245", "PF02875", "PF01225", 
  "PF13241", "PF04101", "PF03033",  "PF00905", 
  "PF03717", "PF06271", "PF08245", "PF02875", 
  "PF03033", "PF00905", "PF03717", "PF01795",
  "PF01098",  "PF00091", "PF12327", "PF00953",
  "PF08478", "PF07478",  "PF13535", "PF14450",
  "PF06723", "PF01820",  "PF09221")

  circ_v_control_wo_query = circ_v_control_df %>%
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

ggsave("/data/san/data0/users/david/intelligence/figures/Circular_bacteriocin_RF_predict_PKS_newDF.png", 
       plot_ROC, width = 10, height = 10, units = "cm")


PF00905_table
PF01225_table
DUF95_table
# ~~~~~~~~~~~~~~~~~~~~~~~~
# PF01225_table               Number 2
# ~~~~~~~~~~~~~~~~~~~~~~~~
circ_v_control_df <- rbind(DUF95_table, PF01225_table) %>% distinct()
circ_v_control_df <- na.omit(circ_v_control_df)
circ_v_control_df$pfam <- as.factor(circ_v_control_df$pfam )
circ_v_control_df$group <- as.factor(circ_v_control_df$group )

## Remove DUf95 and ABC transporter from the training and test dataset
## Remove DUf95 and ABC transporter from the training and test dataset
  values_to_filter <- c("hypothetical","PF00365", "PF00005", 
  "PF01225", "PF08245", "PF02875", "PF01225", 
  "PF13241", "PF04101", "PF03033",  "PF00905", 
  "PF03717", "PF06271", "PF08245", "PF02875", 
  "PF03033", "PF00905", "PF03717", "PF01795",
  "PF01098",  "PF00091", "PF12327", "PF00953",
  "PF08478", "PF07478",  "PF13535", "PF14450",
  "PF06723", "PF01820",  "PF09221")

  circ_v_control_wo_query = circ_v_control_df %>%
          filter(!(pfam %in% values_to_filter))

roc_data <- perform_iterations(circ_v_control_wo_query, num_iterations = 1, seed_values = seed_values)

# Create a combined ROC curve plot
  all_roc_df <- do.call(rbind, lapply(roc_data, function(x) x$roc_df))
  average_auc <- data.frame(Seed = names(roc_data), AUC = sapply(roc_data, function(x) x$auc_val))
  average_auc_subset <- average_auc[average_auc$Seed == 123, ]

plot_ROC = ggplot(all_roc_df, aes(x = FPR, y = TPR, color = as.factor(Seed))) +
  geom_line() +
  labs(x = "False Positive Rate", y = "True Positive Rate", color = "Seed") +
  ggtitle("ROC Curves") +
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

ggsave("/data/san/data0/users/david/intelligence/figures/Circular_bacteriocin_RF_predict_PF01225_newDF.png", 
       plot_ROC, width = 10, height = 10, units = "cm")



# ~~~~~~~~~~~~~~~~~~~~~~~~
# PF01225_table               Number 3
# ~~~~~~~~~~~~~~~~~~~~~~~~
circ_v_control_df <- rbind(DUF95_table, PF00905_table) %>% distinct()
circ_v_control_df <- na.omit(circ_v_control_df)
circ_v_control_df$pfam <- as.factor(circ_v_control_df$pfam )
circ_v_control_df$group <- as.factor(circ_v_control_df$group )

## Remove DUf95 and ABC transporter from the training and test dataset
  values_to_filter <- c("hypothetical","PF00365", "PF00005", 
  "PF01225", "PF08245", "PF02875", "PF01225", 
  "PF13241", "PF04101", "PF03033",  "PF00905", 
  "PF03717", "PF06271", "PF08245", "PF02875", 
  "PF03033", "PF00905", "PF03717", "PF01795",
  "PF01098",  "PF00091", "PF12327", "PF00953",
  "PF08478", "PF07478",  "PF13535", "PF14450",
  "PF06723", "PF01820",  "PF09221")

  circ_v_control_wo_query = circ_v_control_df %>%
          filter(!(pfam %in% values_to_filter))

roc_data <- perform_iterations(circ_v_control_wo_query, num_iterations = 5, seed_values = seed_values)

# Create a combined ROC curve plot
all_roc_df <- do.call(rbind, lapply(roc_data, function(x) x$roc_df))
average_auc <- data.frame(Seed = names(roc_data), AUC = sapply(roc_data, function(x) x$auc_val))
average_auc_subset <- average_auc[average_auc$Seed == 123, ]

plot_ROC = ggplot(all_roc_df, aes(x = FPR, y = TPR, color = as.factor(Seed))) +
  geom_line() +
  labs(x = "False Positive Rate", y = "True Positive Rate", color = "Seed") +
  ggtitle("ROC Curves: PF00905") +
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
                                              label = paste("Avg. AUC = ", round(AUC, 2))), 
                size = 3, hjust = -0.4,
                vjust=-20,
                segment.color = NA, 
                color = "black") +
  scale_color_brewer(palette = "Reds")

ggsave("/data/san/data0/users/david/intelligence/figures/Circular_bacteriocin_RF_predict_PF00905_newDF.png", 
       plot_ROC, width = 10, height = 10, units = "cm")











# ~~~~~~~~~~~~~~~~~~~~~~~~
# VOLCANO PLOT
# ~~~~~~~~~~~~~~~~~~~~~~~~
cont_control_PF00365
cont_control_PF00905 
cont_control_PF01225


# assuming `contingency_table` is your table
diff_pfam_PF00365 <- cont_control_PF00365[1,] - cont_control_PF00365[2,]
diff_pfam_PF00905 <- cont_control_PF00905[1,] - cont_control_PF00905[2,]
diff_pfam_PF01225 <- cont_control_PF01225[1,] - cont_control_PF01225[2,]

# absolute differences for comparison
abs_diff_PF00365 <- abs(diff_pfam_PF00365)
abs_diff_PF00905 <- abs(diff_pfam_PF00905)
abs_diff_PF01225 <- abs(diff_pfam_PF01225)

# sorting to find the most significant differences
top_diff_PF00365 <- sort(abs_diff_PF00365, decreasing = TRUE)[1:100]
top_diff_PF00905 <- sort(abs_diff_PF00905, decreasing = TRUE)[1:100]
top_diff_PF01225 <- sort(abs_diff_PF01225, decreasing = TRUE)[1:100]



# Assuming `contingency_table` is your data
## FISHERS RULE OF 5
# Assuming your 'circular bacteriocin' data is in the first row of the contingency table
contingency_table_filtered <- cont_control_PF00365[, cont_control_PF00365[1,] > 100]
contingency_table_filtered <- cont_control_PF00905 [, cont_control_PF00905 [1,] > 100]

# Calculate fold change
fold_change <- log2(contingency_table_filtered[1,]/contingency_table_filtered[2,])
# Calculate p-values
apply(contingency_table_filtered, 2, function(x) length(x) >= 2)
p_values <- apply(contingency_table_filtered, 2, function(x) chisq.test(x)$p.value)
p_adjusted <- p.adjust(p_values, method = "BH")
# Create a vector of pfams that you want to highlight
highlight_pfams <- c("PF01944", "PF00005")

# Adjust for multiple testing
volcano_data <- data.frame(pfam = names(fold_change),
                           log2FoldChange = fold_change,
                           negLog10PValue = -log10(p_adjusted))
# Add a column to your dataframe indicating whether each pfam should be highlighted
volcano_data$highlight <- ifelse(volcano_data$pfam %in% highlight_pfams, "Highlight", "Don't highlight")
volcano_data$negLog10PValue <- pmin(volcano_data$negLog10PValue, 300)
volcano_data$log2FoldChange <- pmin(volcano_data$log2FoldChange, 20)
(ggplot(volcano_data, aes(x = log2FoldChange, y = negLog10PValue, color = highlight)) +
  geom_point(alpha = 0.4, size = 0.2) +
  theme_classic() +
  xlab("Log2 Fold Change") +
  ylab("-Log10 Adjusted P-value") +
  ggtitle("Volcano plot of Pfam Counts") +
  scale_color_manual(values = c("black", "red")) +
  theme_bw() +
  theme(legend.position = "none") +
           theme(
    plot.title = element_text(size = 10),  # Adjust size as needed
    axis.title.x = element_text(size = 10),  # Adjust size as needed
    axis.title.y = element_text(size = 10),  # Adjust size as needed
    axis.text.x = element_text(size = 10),  # Adjust size as needed for x-axis labels
    axis.text.y = element_text(size = 10)   # Adjust size as needed for y-axis labels
  )  )%>%
ggsave("/data/san/data0/users/david/intelligence/figures/Cicular_bacteriocin_volcano_new_P905_DB.png", ., width = 7, height = 7, units = "cm")

# ~~~~~~~~~~~~~~~~~~~~~~~~
# Plot difference
# ~~~~~~~~~~~~~~~~~~~~~~~~
abs_diff_PF00365
abs_diff_PF00905
abs_diff_PF01225


top_diff_pfam_df = top_diff_pfam %>% as.data.frame()
top_diff_pfam_df$pfam = rownames(top_diff_pfam_df)

# Order the data by the absolute difference
hypothetical <- top_diff_pfam_df[order(abs(top_diff_pfam_df$.)), ]
# Convert the vector to a dataframe
hypothetical_df <- data.frame(
  Protein = names(hypothetical),
  Difference = hypothetical
)

# Take the top 100
hypothetical <- head(hypothetical, 100)

# Make the difference negative for the control group
hypothetical$.[hypothetical$group == "Control"] <- - hypothetical$diff[hypothetical$group == "Control"]

# Create the plot
ggplot(hypothetical, aes(x = reorder(protein, diff), y = diff, fill = group)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  labs(x = "Protein", y = "Difference in Abundance")


# ~~~~~~~~~~~~~~~~~~~~~~~~
# Bar chart of top pfams per group
# ~~~~~~~~~~~~~~~~~~~~~~~~

## layout together
ggplot(top_GO_terms, aes(x = reorder(GO_name, total), y = total, fill = group)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  labs(x = "GO Term", y = "Count", fill = "Group") +
  theme_minimal() +
  ggtitle(paste("Top", top_n, "GO terms per group"))

# Number of top GO terms to select
top_n <- 30

# Find top N GO terms for each group
top_pfams_terms <- ctrl_v_circ %>%
  filter(!is.na(pfam)) %>%
  filter(pfam != "hypothetical") %>%
  filter(pfam != "") %>%
  dplyr::select(pfam, id, Nucleotide_acc) %>%
  group_by(id, pfam) %>%
  summarise(total = n(), .groups = "drop") %>%
  group_by(id) %>%
  arrange(desc(total)) %>%
  top_n(top_n, total) 

# Plot
pfam_duf95vsPks = ggplot(top_pfams_terms, aes(x = reorder(pfam, total), y = total, fill = id)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  labs(x = "Pfam", y = "Count", fill = "Group") +
  theme_bw() +
  scale_fill_manual(values=c("circ" = "#d02727", "PF00365" = "#4c4cbc", "PF0090" = "#686a68", "PF01225" = "#4c9f7b")) +
  ggtitle(paste("Top", top_n, "pfam per group"))
ggsave("/data/san/data0/users/david/intelligence/figures/Cicular_bacteriocin_pfams_v_top_newDB.png", pfam_duf95vsPks, width = 20, height = 26, units = "cm")




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
#  RF on GO iterations  (1) Biological Process
# ~~~~~~~~~~~~~~~~~~~~~~~~
# https://bioconductor.org/packages/release/data/annotation/src/contrib/GO.db_3.17.0.tar.gz
install.packages("https://bioconductor.org/packages/release/data/annotation/src/contrib/GO.db_3.17.0.tar.gz", repos = NULL, type = "source")
library(GO.db)
#library(topGO)
library(ranger)
library(caret)
library(ROCR)
library(pROC)
# DF = circular_DF_GO
seed_values <- c(123, 456, 789, 159, 753) # Use your desired seed values

circular_DF_GO$GO_acc
#categories <- unique(circ_df_GOont$category) # Get all unique categories

# Assuming "df" is your dataframe and "GO_acc" is the column of GO IDs
go_ids <- as.character(circular_DF_GO$GO_acc)
# Get the GO categories (ontologies) corresponding to these IDs
go_categories <- select(GO.db, keys = go_ids, columns = "ONTOLOGY", keytype = "GOID")
colnames(go_categories) <- c("GO_acc", "category")
go_categories = distinct(go_categories)

# filter our proteins with pfams
values_to_filter
colnames(circular_DF_GO_filtered)

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
  scale_fill_manual(values=c("circ" = "#d02727", "PF00365" = "#4c4cbc", "PF0090" = "#686a68", "PF01225" = "#4c9f7b")) +
  ggtitle(paste("Top", top_n, "GO: Biological Process"))
ggsave("/data/san/data0/users/david/intelligence/figures/Cicular_bacteriocin_GO_BP_v_pks_newDB.png", GO_BP_duf95vsPks, width = 20, height = 26, units = "cm")




# ~~~~~~~~~~~~~~~~~~~~~~
# RANDOM FORREST FOR GO TERMS
# ~~~~~~~~~~~~~~~~~~~~~~

# Initialize lists to store the results
roc_dfs <- list() 
auc_values <- list() 
categories = go_categories %>% filter(!is.na(category)) %>% 
  dplyr::select(category) %>%
  distinct()


circ_df_GOont_narm = circular_DF_GO_filtered %>% filter(!is.na(GO_acc))

colnames(circ_df_GOont_narm)
# Loop over categories
for(cat in categories$category){
  # Create empty data frames to store the results for each category
  roc_dfs[[cat]] <- data.frame()
  auc_values[[cat]] <- data.frame()

  # Subset the data for the current category

  circular_DF_GO2 <- circ_df_GOont_narm[circ_df_GOont_narm$category == cat,]
  circular_DF_GO2 = circular_DF_GO2 %>% filter(!is.na(GO_acc))

  circ_df_GOont_narm$GO_name = as.factor(circ_df_GOont_narm$GO_name)
  circ_df_GOont_narm$group = as.factor(circ_df_GOont_narm$group)
  # Loop over seed values
  for(i in 1:length(seed_values)){
    set.seed(seed_values[i]) 
    training_rows_go <- createDataPartition(circular_DF_GO2$group, p = 0.8, list = FALSE)
    training_data_go <- circular_DF_GO2[training_rows_go, ]
    testing_data_go <- circular_DF_GO2[-training_rows_go, ]
  
    rf_model_go <- ranger(group ~ GO_acc, data = circ_df_GOont_narm,
                          importance = 'impurity', num.threads = 12, 
                          probability = TRUE)
  
    predictions_go <- predict(rf_model_go, data = testing_data_go, num.threads = 12, type = "response")
    prob_control_go <- predictions_go$predictions[, "control"]
    prob_circ_go <- predictions_go$predictions[, "circular bacteriocin"]
  
    # Create a ROC curve
    roc_obj_go <- roc(testing_data_go$group, prob_control_go) 
    # Create a data frame for ggplot
    roc_df_go <- data.frame(
      Seed = seed_values[i],
      FPR = roc_obj_go$specificities,
      TPR = roc_obj_go$sensitivities
    )
  
    auc_val <- auc(roc_obj_go)
    auc_df_temp <- data.frame(Seed = seed_values[i], AUC = auc_val)
  
    # Save the data frame to the list
    roc_dfs[[cat]] <- rbind(roc_dfs[[cat]], roc_df_go)
    auc_values[[cat]] <- rbind(auc_values[[cat]], auc_df_temp)
  }
}

# Combine all the ROC data frames into one and add category as a column
all_roc_df <- do.call(rbind, lapply(names(roc_dfs), function(cat) {
  df <- roc_dfs[[cat]]
  df$Category <- cat
  df
}))

# Similarly for AUC values
all_auc_df <- do.call(rbind, lapply(names(auc_values), function(cat) {
  df <- auc_values[[cat]]
  df$Category <- cat
  df
}))

# Calculate average AUC per category
average_auc <- aggregate(AUC ~ Category, data = all_auc_df, FUN = mean)

# Plot ROC curves
# Define your color palette
color_palette <- c("Seed1" = "#FF0000", 
                   "Seed2" = "#FF3333", 
                   "Seed3" = "#FF6666", 
                   "Seed4" = "#FF9999",
                   "Seed5" = "#FFCCCC") # adjust these to your seeds names and preferred colors


library(ggrepel)
plot_ROC_GO = ggplot(all_roc_df, aes(x = FPR, y = TPR, color = as.factor(Seed), linetype = Category, group = interaction(Seed, Category))) +
  geom_line() +
  labs(x = "False Positive Rate", y = "True Positive Rate", color = "Seed", linetype = "Category") +
  ggtitle("ROC Curves") +
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
    geom_text_repel(data = average_auc, 
                aes(x = 0.8, y = 0.2, label = paste(Category, ": Avg. AUC = ", round(AUC, 2))), 
                size = 3, hjust = -0.2, inherit.aes = FALSE, segment.color = NA, nudge_y = 0.05) +
  scale_color_brewer(palette = "Reds") 
ggsave("/data/san/data0/users/david/intelligence/figures/Cicular_bacteriocin_RF_GO_category_newDB.png", plot_ROC_GO, width = 10, height = 10, units = "cm")


# ~~~~~~~~~~~~~~~~~~~~~~~~
# TREEMAP
# ~~~~~~~~~~~~~~~~~~~~~~~~
install.packages('treemap')

# Then, load the package into your workspace.
library(treemap)
colnames(circ_df_GOont)
treedata_control =  circ_df_GOont %>% #filter(group == "control") %>%
  group_by(GO_name) %>%
  summarise(n=sum(n())) %>%
  left_join(., circ_df_GOont) %>% 
  dplyr::select(GO_name, GO_acc, pfam, n, group)

treedata_control$n = as.double(treedata_control$n)

treemap(
    treedata_control,
    index = c("GO_name", "group"),
    vSize = "n",
    title = "Global Population by Continent and Country"
  )
dev.off()



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



ctrl_v_circ %>% 
  filter(pfam %in% phage_defense_pfams_df$Family) %>% 
  group_by(group) %>%
  summarise(count = n()) %>%
  arrange(desc(count))


pfam_descs = DUF95_paper_40k_window %>% dplyr::select(pfam, desc) %>% distinct() 
ctrl_v_circ %>% 
  filter(pfam %in% phage_defense_pfams_df$Family) %>% 
  left_join(., pfam_descs) %>%
  dplyr::select(pfam, desc) %>%
  group_by(desc) %>%
  summarise(count = n()) %>%
  arrange(desc(count)) %>%
  View()

## Find out he genus and unique regions that have multiple phage defense systems putatively
joined_df <- inner_join(DUF95_paper_40k_window, phage_defense_pfams_df, by = c("pfam" = "Family"))
grouped_df <- joined_df %>%
  group_by(unique_window) %>%
  summarise(n = n(), .groups = 'drop') %>%
  left_join(., (DUF95_paper_40k_window %>% dplyr::select(unique_window, phylum, genus, Nucleotide_acc)), by = "unique_window") %>%
  arrange(desc(n)) %>%
  distinct()

ordered_df <- grouped_df %>%
  arrange(desc(n)) %>% 
  filter(n > 15) %>%
  left_join(., (DUF95_paper_40k_window), by = "unique_window") %>%
  filter(pfam %in% phage_defense_pfams1$Family) %>%
  dplyr::select(unique_window, n) %>%
  distinct() %>%
  arrange(desc(n)) 

 
library()

# ~~~~~~~~~~~~~~~~~~~~
# Plotting the genomes related to phage
# ~~~~~~~~~~~~~~~~~~~~

  library(thacklr)
  library(gggenomes)
  library(crayon)
  library(ggtree)
  library(RcolorBrewer)
      sample_n_groups = function(grouped_df, size, replace = FALSE, weight=NULL) {
                      grp_var <- grouped_df %>% 
                      groups %>%
                      unlist %>% 
                      as.character
                    random_grp <- grouped_df %>% 
                      summarise() %>% 
                      sample_n(size, replace, weight) %>% 
                      mutate(unique_id = 1:NROW(.))
                    grouped_df %>% 
                      right_join(random_grp, by=grp_var) %>% 
                      group_by_(grp_var) 
      }
    window_size_threshold <- 50000
    # Function to split seq_id based on the condition
    split_seq_id <- function(seq_id, start, end) {
      if (length(start) == 0) {
        return(character(0))
      }

      seq_id_split <- character(length(start))
      current_group <- 1

      for (i in 1:length(start)) {
        if (i > 1 && (start[i] - end[i - 1]) > 30000) {
          current_group <- current_group + 1
        }
        seq_id_split[i] <- paste0(seq_id, "_", current_group)
      }

      return(seq_id_split)
    }
 
    sample_circ <- DUF95_paper_40k_window %>%
        filter(window_size > 20000) %>%
        dplyr::rename(seq_id = unique_window) %>%
        group_by(seq_id) %>%
        filter(n() > 10) %>%
        arrange(seq_id, start) %>%
        ungroup() %>% 
        group_by(seq_id) %>%
        sample_n_groups(20)


    top_pdefense <- DUF95_paper_40k_window %>%
        filter(unique_window %in% ordered_df$unique_window) %>%
        filter(window_size > 20000) %>%
        dplyr::rename(seq_id = unique_window) %>%
        group_by(seq_id) %>%
        filter(n() > 10) %>%
        arrange(seq_id, start) %>%
        ungroup() %>% 
        group_by(seq_id) 
    
DUF95_paper_40k_window %>% 
  dplyr::select(unique_window, pfam, desc) %>%
  filter(grepl("restriction", desc)) %>%
  filter(pfam %in% bonferoni_table_pfams$Pfam) %>%
  View()

statistical_rm_pfams = DUF95_paper_40k_window %>% 
  filter(grepl("restriction", desc)) %>%
  filter(pfam %in% bonferoni_table_pfams$Pfam) %>%
  dplyr::select(pfam, Query)  %>%
  distinct() 


## RM systems seem to statistically localise near circular bacteriocins
bonferoni_table_pfams %>% filter(Pfam %in% statistical_rm_pfams$pfam )
top_rms <- DUF95_paper_40k_window %>%
    dplyr::rename(seq_id = unique_window) %>%
    dplyr::group_by(seq_id) %>%
    dplyr::filter(any(pfam %in% unique(statistical_rm_pfams$pfam)) & n() > 10 & window_size > 20000)


top_rms <- DUF95_paper_40k_window %>%
    dplyr::rename(seq_id = unique_window) %>%
    dplyr::group_by(seq_id) %>%
    dplyr::filter(any(pfam %in% unique(bonferoni_table_pfams$Pfam)) & n() > 10 & window_size > 20000)


set.seed(123)  # for reproducibility
# Split the dataframe into a list of dataframes, one for each genus
grouped_top_rms <- top_rms %>% 
  group_by(genus) %>%
  group_split()

# Sample 15 dataframes from the list and bind them together
sampled_df_list <- sample(grouped_top_rms, size = )
top_rms_sample <- bind_rows(sampled_df_list)

circular_bacteriocin_pfams = c("PF01944", "PF00005", "PF16576", "PF16942", "PF09221")
top_rms_sample <- top_rms_sample %>%
  group_by(Protein_acc) %>%
  mutate(
    color = case_when(
      any(pfam %in% statistical_rm_pfams$pfam) ~ 'darkgreen', # phage defense
      any(pfam %in% circular_bacteriocin_pfams) ~ '#a22a15', # DUF95
      TRUE ~ 'grey' # default color
    )
  ) %>%
  ungroup()

  ###
    # GGenomes plot
  ###
( top_rms_sample %>%
  gggenomes() +
    geom_seq(color = "black") +
    geom_seq_label(color = "black", label="") +
    geom_gene(aes(fill = color), size=4) +     
    theme(panel.background = element_rect(fill = "#EAEAEA"),
        legend.background = element_rect(fill = "EAEAEA")) +                            
    scale_fill_identity() ) %>%
  ggsave("/data/san/data0/users/david/intelligence/figures/phage_defense_systems_near_circular_RMs.png", .,width = 8, height = 6, units = "in", dpi = 300, )


# Write all protein accessions to a file and download the sequences
top_rms %>%
  ungroup() %>%
  dplyr::select(Protein_acc) %>%
  distinct() %>% 
  write_csv("/data/san/data0/users/david/intelligence/phage_defense/phage_defense_systems_near_circular_RMs_protein_accs.csv", col_names = FALSE)
library(rentrez)

# download and cluster the sequences
# This was done in folder
#

# Load necessary library
library(igraph)
library(ggraph)
library(tidygraph)
rm_net <- read_tsv("/data/san/data0/users/david/intelligence/phage_defense/ava_blast.net", col_names = FALSE)
colnames(rm_net) <- c("Protein_acc", "Protein_acc2") 

# Create an igraph object
g_rm <- graph_from_data_frame(rm_net, directed=T)  # Set directed=TRUE if you want a directed graph

# Load your additional data for node attributes
node_data <- fread("/data/san/data0/users/david/intelligence/phage_defense/phage_defense_systems_near_circular_RMs_protein_accs_bakta.tsv", skip = 3 )
colnames(node_data)[1] <- c("Protein_acc")

# list the accessions in the largest nodes
node_degrees <- degree(g_rm)
df_num_nodes <- data.frame(
  node = names(node_degrees),
  degree = node_degrees
)

# Run the walktrap community detection algorithm
wc <- cluster_walktrap(g_rm)
# Get the membership vector (this tells you which cluster each node belongs to)
membership <- membership(wc)
# Count the number of nodes per cluster
cluster_sizes <- table(membership)
# Get the Protein_acc for each cluster
proteins <- sapply(1:max(membership), function(i) {
  # Get the names of the nodes in this cluster
  node_names <- names(membership)[membership == i]
  
  # Return the node names
  paste(node_names, collapse = ", ")
})
# Convert to a data frame
cluster_info <- data.frame(
  cluster = names(cluster_sizes),
  size = as.numeric(cluster_sizes),
  Protein_acc = proteins
)

# Assign the cluster IDs to the nodes in your graph
V(g_rm)$cluster <- membership(wc)
V(g_rm)$cluster  %>% sort()
  # Convert your data frame to long format
cluster_info_long <- cluster_info %>%
  separate_rows(Protein_acc, sep = ", ") %>%
  mutate(Protein_acc = trimws(Protein_acc))
cluster_info_joined <- left_join(cluster_info_long, node_data, by = "Protein_acc") 


# Now you can use the 'color' attribute for color in your plot
g_rm_protein_plot = ggraph(g_rm) +
  geom_edge_link(width = 0.1) +
  geom_node_point(size = 2, 
                      color = case_when(
                        V(g_rm)$cluster == 344 ~ "#8f1d1d",
                        V(g_rm)$cluster == 279 ~ "#8f1d1d",
                        V(g_rm)$cluster == 254 ~ "#8f1d1d",
                        V(g_rm)$cluster == 2 ~ "darkgreen",
                        V(g_rm)$cluster == 243 ~ "darkgreen",
                        V(g_rm)$cluster == 253 ~ "orange",
                        V(g_rm)$cluster == 342 ~ "black",
                        V(g_rm)$cluster == 4 ~ "black",
                        V(g_rm)$cluster == 328 ~ "black",
                        V(g_rm)$cluster == 338 ~ "black",
                        V(g_rm)$cluster == 339 ~ "black",
                        V(g_rm)$cluster == 340 ~ "black",
                        V(g_rm)$cluster == 336 ~ "black",
                        V(g_rm)$cluster == 337 ~ "black",
                        V(g_rm)$cluster == 343 ~ "black",
                        V(g_rm)$cluster == 341 ~ "darkblue",
                        V(g_rm)$cluster == 20 ~ "#8f1d1d",
                        V(g_rm)$cluster == 11 ~ "#1d8f8f",
  TRUE ~ "grey"  # default color for any other clusters
)) + 
  theme_void() +
  theme(panel.background = element_rect(fill = "#EAEAEA"),
        legend.background = element_rect(fill = "EAEAEA")) +
  scale_color_identity() # use 'color' for color
ggsave("/data/san/data0/users/david/intelligence/figures/phage_defense_systems_near_circular_RMs_network.png", g_rm_protein_plot,width = 14, height = 8, units = "in", dpi = 300 )

cluster_info_joined %>% left_join(., DUF95_paper_40k_window) %>%
  filter(cluster == c(342, 4, 328, 338, 339, 340, 336, 337, 343)) %>%
  select(Protein_acc, cluster, size, genus) %>%
  View()








## 

##


( gggenomes(top_pdefense) +
    geom_seq(color = "black") +
    geom_seq_label(color = "black") +
    geom_gene(aes(fill = dplyr::case_when(
              pfam %in% phage_defense_pfams_df$Family ~ 'darkgreen', # phage defense
              pfam == "PF01944" ~ '#a22a15', # DUF95
              pfam == "PF00005" ~ 'tomato', # ABC transport # nolint
              pfam == "PF16576" ~ 'tomato',
              pfam == "PF16942" ~ 'tomato', # bac
              pfam == "PF09221" ~ 'tomato'))) +                                             
    ggplot2::scale_fill_identity() ) %>%
  ggsave("/data/san/data0/users/david/intelligence/figures/phage_defense_systems_near_circular_test.png", .,width = 8, height = 8, units = "in", dpi = 300, )


top_pdefense$seq_id %>% unique() 

top_pdefense$pfam == "PF07751" 

top_pdefense %>% filter(seq_id == "AIXV01000001.1_155560") %>% 
distinct() %>% dplyr::select(pfam, desc, genus, seq_id, Protein_acc) %>%
filter(!is.na(pfam)) %>%
filter(pfam != "") %>%
filter(pfam != "hypothetical") %>%
distinct() %>%
View()

# ~~~~~~~~~~~~~~~~~~~~~~~~
# From PFAMs and GOs to systems  :  GRAPH network for PFAMs
# ~~~~~~~~~~~~~~~~~~~~~~~~

  library(ggraph)
  library(igraph)
  library(tidygraph)
  library(ggdark)

distinct_taxa_w_CIRC =  DUF95_paper_40k_window %>% 
  dplyr::select(taxid) %>%
  distinct()


gene_lists <- aggregate(pfam ~ Nucleotide_acc, data = DUF95_paper_40k_window, FUN = paste, collapse = ",")
  # Split the combined pfam values into individual pfams
  gene_lists$pfam <- strsplit(gene_lists$pfam, ",")
  # Create a table of gene co-occurrences
  gene_pairs <- lapply(gene_lists$pfam, combn, m = 2, simplify = FALSE)
  gene_pairs <- unlist(gene_pairs, recursive = FALSE)
  gene_pairs <- t(simplify2array(gene_pairs))
  gene_pairs <- data.frame(gene1 = gene_pairs[,1], gene2 = gene_pairs[,2])
  gene_pairs_3 <- gene_pairs %>% 
    filter(!grepl('hypothetical', gene1)) %>% filter(!grepl('hypothetical', gene2)) %>% 
    filter(!grepl('PF01944', gene2)) %>% filter(!grepl('PF00005', gene2)) %>%
    filter(!grepl('PF01944', gene1)) %>% filter(!grepl('PF00005', gene1)) %>%
    filter(!is.na(gene1)) %>% filter(!is.na(gene2)) %>% filter(gene1 != gene2)

  genes_count_noDUF <- as.data.frame(table(gene_pairs_3))

  colnames(genes_count_noDUF) <- c("gene1", "gene2", "Freq")
  desc_table1 <- (DUF95_paper_40k_window %>% select(pfam, desc)) %>% as.data.frame()
  colnames(desc_table1) <- c("gene1", "desc1")
  desc_table2 <- (DUF95_paper_40k_window %>% select(pfam, desc)) %>% as_tibble()
  colnames(desc_table2) <- c("gene2", "desc2")
  left_join(as.data.frame(genes_count_noDUF), as.data.frame(desc_table1), by="gene1") %>% distinct() %>%
    left_join(., as.data.frame(desc_table2), by="gene2") %>% distinct() %>% View()


  genes_count_500 <- genes_count_noDUF %>% 
  filter(gene1 != "") %>%
  filter(Freq > 500 & Freq < 2000)




## filter to remove
  circular_pfams_to_filter = c("PF09221","PF14089","PF00005", "", "PF13555", "PF01381",
  "PF02687", "PF00528", "PF13304", "PF02463", "PF01408", "PF01380")

  bacteriocin_pfams = c("PF14867", "PF13537", "PF13575", "PF12730", "PF12730", "PF02441", "PF01944",
                 "PF18218", "PF01320", "PF00072", "PF00486", "PF00005", "PF02518", "PF00512",
                 "PF03857", "PF08951", "PF05147", "PF00881", "PF06182", "PF00296", "PF05402",
                 "PF15565", "PF15586", "PF01719", "PF02794", "PF04055", "PF13165", "PF15007",
                 "PF11083", "PF13471", "PF00733", "PF02540", "PF05147", "PF02052", "PF04604",
                 "PF08130")

  
genomad_df = fread("/data/san/data0/databases/DRAM/amg_database.20221222.tsv")
genomad_df$ANNOTATION_ACCESSIONS %>% unique()

dram_amg_df = fread("/data/san/data1/users/david/amg_database.20221222.tsv")
dram_amg_df$PFAM

# just take phae defense systems
      genes_count_plot = genes_count_plot %>% 
        filter(gene1 %in% phage_defense_pfams_df$Family & gene2 %in% phage_defense_pfams_df$Family) %>%
        filter(!gene1 %in% circular_pfams_to_filter & !gene2 %in% circular_pfams_to_filter) 

# just take statistically relevant pfams
      genes_count_plot = genes_count_plot %>% 
        filter(gene1 %in% bonferoni_table_pfams$Pfam & gene2 %in% bonferoni_table_pfams$Pfam) %>%
        filter(!gene1 %in% circular_pfams_to_filter & !gene2 %in% circular_pfams_to_filter) 


genes_count_plot <- genes_count_500 %>% dplyr::filter(!gene1 %in% circular_pfams_to_filter) %>%
    dplyr::filter(!gene2 %in% circular_pfams_to_filter)

genes_count_plot %>% arrange(desc(Freq)) %>% 
  left_join(., pfam_descs, by=c("gene1" = "pfam")) %>% 
  left_join(., pfam_descs, by=c("gene2" = "pfam")) %>%
  arrange(desc(Freq)) %>% View()

  graph <- as_tbl_graph(genes_count_plot ,directed = FALSE)
  V(graph)$clu <- as.character(membership(cluster_louvain(graph)))
  V(graph)$size <- igraph::degree(graph)
  E(graph)$arrow.mode <- 0
  E(graph)$weight <- genes_count_plot$Freq
  E(graph)$width <- E(graph)$weight/100


  # plot using ggraph #### 
  ggraph(graph, layout = "kk") + 
  geom_edge_link(aes(width = Freq), 
                edge_colour = "black", 
                edge_alpha = 0.2) + 
  geom_node_point(aes(size = size),
                  colour = case_when(names(V(graph)) %in% phage_defense_pfams_df$Family ~ "#a6230c",
                                      names(V(graph)) %in% bacteriocin_pfams ~ "#0c5aa6",
                                      names(V(graph)) %in% dram_amg_df$PFAM ~ "green",
                                     TRUE ~ "grey")) + 
  geom_node_text(aes(label = names(V(graph))), 
                color = 'black', 
                size = 6, 
                repel = TRUE,
                family="ariel") + 
    ggforce::theme_no_axes()  +
    theme_void() +
    theme( 
          text = element_text(size = 20, color = "black")) +
    ggtitle("Co-occurence network of PFAMs within 40kb of putative circular bacteriocin gene clusters (over 500 occurences)")





















# ~~~~~~~~~~~~~~~~~~~~~~~~
# PCOA? 
# ~~~~~~~~~~~~~~~~~~~~~~~~
install.packages("ggbeeswarm")
install.packages("ade4")
install.packages("ggpubr")
library(ade4)
library(ggbeeswarm)
library(ggpubr)
library(vegan)
# Create a contingency table
table <- circ_v_control_df %>%
  group_by(group, pfam) %>%
  summarise(count = n(), .groups = 'drop')
# remove na and hypothetical pfams
filtered_table <- table %>% filter(!(pfam %in% c("", "hypothetical")))
filtered_table$group <- as.factor(filtered_table$group)
filtered_table$count <- as.numeric(filtered_table$count)




# Reshape the data
reshaped_data <- spread(filtered_table, key = group, value = count, fill = 0)
count_matrix <- data.matrix(reshaped_data[-1])
dist_matrix <- vegdist(count_matrix, method = "bray")

# Check the size of the distance matrix
print(dim(dist_matrix))


# Transpose the data frame to have Pfams as rows and groups as columns
wide_df_t <- t(wide_df[,-1])
colnames(wide_df_t) <- wide_df$group
colnames(wide_df_t) <- c("control", "circ")
wide_test = wide_df_t %>% as.data.frame() %>% 
    filter(control > 20 & circ > 20)


# Perform chi-square test
chisq_result <- chisq.test(wide_test)





# Run PCoA
pcoa <- cmdscale(dist_matrix, eig = TRUE, k = 2)
# Compute Bray-Curtis dissimilarity
dist_matrix <- vegdist(count_matrix, method = "bray")
# Check the size of the distance matrix
print(dim(dist_matrix))
# Compute Bray-Curtis dissimilarity
dist_matrix <- vegdist(matrix, method = "bray")
# Perform the PCoA
pcoa <- cmdscale(dist_matrix, eig = TRUE, k = 2) # k is the number of dimensions
# Convert PCoA results to a data frame
pcoa_df <- as.data.frame(pcoa$points)
# Add column names for easier plotting
names(pcoa_df) <- c("PCOA1", "PCOA2")
# Add 'pfam' column to data frame
pcoa_df$pfam <- rownames(pcoa_df)

ggplot(pcoa_df, aes(x = PCOA1, y = PCOA2)) +
  geom_point() +
  theme_minimal() +
  labs(x = "PCoA 1", y = "PCoA 2", title = "PCoA of Pfams") +
  theme(text = element_text(size = 18))


install.packages("FactoMineR")
library("FactoMineR")
install.packages("corrr")
library("corrr")
devtools::install_github("kassambara/ggcorrplot")
library(ggcorrplot)




rownames(reshaped_data) <- reshaped_data[1,]
numerical_data <- reshaped_data[-1,] 
numerical_data <- as.data.frame(numerical_data)
numerical_data <- numerical_data %>% drop_na(pfam)
rownames(numerical_data) <- numerical_data[,1]
numerical_data = numerical_data[,-1]


df_transposed <- as.data.frame(t(numerical_data))
data_normalized  <- scale(df_transposed)
head(data_normalized)

corr_matrix <- cor(data_normalized)
# Remove rows with non-finite values
corr_matrix <- corr_matrix[complete.cases(corr_matrix), ]
ggcorrplot(corr_matrix) 

all(is.finite(corr_matrix))
data.pca <- princomp(corr_matrix)
summary(data.pca)



set.seed(1)
DUF95_table <- DUF95_paper_40k_window %>% dplyr::select(Nucleotide_acc, pfam) %>% distinct()
DUF95_table$group <- "circular bacteriocin"

PF00365_table <- PF00365_df %>% dplyr::select(Nucleotide_acc, subject.y) %>% distinct() 
colnames(PF00365_table) <- c("Nucleotide_acc", "pfam")
PF00365_table$group <- "control"
PF00365_table$pfam = sub("\\..*", "", PF00365_table$pfam)














# Assuming df is your data frame, and it has columns 'group' and 'pfam'
# And group has levels 'control' and 'circular'
unique_pfams <- unique(df$pfam)

# Initialize a data frame to store the p-values
p_values <- data.frame(pfam=character(), p_value=numeric())

# Loop over the unique pfams
for (pfam in unique_pfams) {
  # Create a subset of the data for this pfam
  df_subset <- df[df$pfam == pfam, ]
  
  # Generate the contingency table
  contingency_table <- table(df_subset$group)
  
  # Check that the pfam appears at least once in each group
  if (all(contingency_table > 0)) {
    # Perform the chi-square test
    chi2_test <- chisq.test(contingency_table, correct=FALSE)
  
    # Store the p-value
    p_values <- rbind(p_values, data.frame(pfam=pfam, p_value=chi2_test$p.value))
  } else {
    # If the pfam doesn't appear in one or both of the groups, set the p-value to NA
    p_values <- rbind(p_values, data.frame(pfam=pfam, p_value=NA))
  }
}

# Adjust the p-values for multiple comparisons
p_values$p_adjusted <- p.adjust(p_values$p_value, method = "bonferroni")

# Print the p-values
print(p_values)



df = circ_v_control_df
unique_pfams <- unique(df$pfam)

 top_pfams  =   circ_v_control_df %>% 
        filter(group == "circular bacteriocin") %>% 
        group_by(pfam) %>%
        summarise(n = n()) %>%
        arrange(desc(n)) %>%
        dplyr::select(pfam)


ten_pfams = top_pfams[1:50,]


    contingency_tables <- list()
    # Loop over the unique pfams
    for (pfam in unique_pfams) {
      # Create a subset of the data for this pfam
      df_subset <- df[df$pfam == pfam, ]
      
      # Generate the contingency table
      contingency_table <- table(df_subset$group)
      
      # Add the contingency table to the list
      contingency_tables[[pfam]] <- contingency_table
      
      # Check that the pfam appears at least once in each group
      if (length(contingency_table) > 1 && all(contingency_table > 0)) {
        # Perform the chi-square test
        chi2_test <- chisq.test(contingency_table, correct=FALSE)
      
        # Store the p-value
        p_values <- rbind(p_values, data.frame(pfam=pfam, p_value=chi2_test$p.value))
      } else {
        # If the pfam doesn't appear in one or both of the groups, set the p-value to NA
        p_values <- rbind(p_values, data.frame(pfam=pfam, p_value=NA))
      }
    }




    # Initialize a data frame to store the p-values
    p_values_specific <- data.frame(pfam=character(), p_value=numeric())

    # Loop over the specific pfams
    for (pfam in top_pfams$pfam) {
      # Retrieve the contingency table from the list
      contingency_table <- contingency_tables[[pfam]]
      # Check that the pfam appears at least once in each group
      if (length(contingency_table) > 1 && all(contingency_table > 0)) {
        # Perform the chi-square test
        chi2_test <- chisq.test(contingency_table, correct=FALSE)
        # Store the p-value
        p_values_specific <- rbind(p_values_specific, data.frame(pfam=pfam, p_value=chi2_test$p.value))
      } else {
        # If the pfam doesn't appear in one or both of the groups, set the p-value to NA
        p_values_specific <- rbind(p_values_specific, data.frame(pfam=pfam, p_value=NA))
      }
    }
    # Adjust the p-values for multiple comparisons
    p_values_specific$p_adjusted <- p.adjust(p_values_specific$p_value, method = "bonferroni")
    # Print the p-values
    print(p_values_specific)


pfam_descs = DUF95_paper_40k_window %>% dplyr::select(pfam, desc ) %>% distinct()

p_values_specific %>% filter(p_adjusted < 0.0005) %>%
            arrange(p_adjusted) %>%
            left_join(., pfam_descs) %>% 
            filter("PF05147" %in% pfam) %>% 
            













# region proportion
  prop.table(table((circ_v_control_df %>% distinct(Nucleotide_acc, group))$group))

# pfam proportion
  prop.table(table(circ_v_control_df$group))





# apply binomial test to each row
p.values <- apply(contingency_table, 1, function(x) {
    binom.test(x[1], sum(x))$p.value
})

# create a data frame of the results
results <- data.frame(row.names(contingency_table), p.values)

# adjust for multiple testing (optional)
results$adj.p.values <- p.adjust(results$p.values, method = "bonferroni")

# view results
print(results)
View(results)

# find significant rows (optional)
significant_rows <- row.names(results)[results$adj.p.values < 0.05]
print(significant_rows)
















gene_lists <- aggregate(pfam ~ Nucleotide_acc, data = total_refseq_pfam_anot, FUN = paste, collapse = ",")
  # Split the combined pfam values into individual pfams
  gene_lists$pfam <- strsplit(gene_lists$pfam, ",")
  # Create a table of gene co-occurrences
  gene_pairs <- lapply(gene_lists$pfam, combn, m = 2, simplify = FALSE)
  gene_pairs <- unlist(gene_pairs, recursive = FALSE)
  gene_pairs <- t(simplify2array(gene_pairs))
  gene_pairs <- data.frame(gene1 = gene_pairs[,1], gene2 = gene_pairs[,2])
  gene_pairs_3 <- gene_pairs %>% 
    filter(!grepl('hypothetical', gene1)) %>% filter(!grepl('hypothetical', gene2)) %>% 
    filter(!grepl('PF01944', gene2)) %>% filter(!grepl('PF00005', gene2)) %>%
    filter(!grepl('PF01944', gene1)) %>% filter(!grepl('PF00005', gene1))

  genes_count_noDUF <- as.data.frame(table(gene_pairs_3))

  colnames(genes_count_noDUF) <- c("gene1", "gene2", "Freq")
  desc_table1 <- (putative_circular_with_DUF95 %>% select(pfam, desc)) %>% as.data.frame()
  colnames(desc_table1) <- c("gene1", "desc1")
  desc_table2 <- (putative_circular_with_DUF95 %>% select(pfam, desc)) %>% as_tibble()
  colnames(desc_table2) <- c("gene2", "desc2")
  left_join(as.data.frame(genes_count_100), as.data.frame(desc_table1), by="gene1") %>% distinct() %>%
    left_join(., as.data.frame(desc_table2), by="gene2") %>% distinct() %>% View()
  genes_count_10 <- genes_count_noDUF %>% filter(Freq > 10)