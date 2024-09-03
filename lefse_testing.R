#!/usr/bin/Rscript
library(optparse)
library(maaslin3)

# Define the command line options
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Path to input txt file", metavar = "character")
)

# Parse the command line options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("Input file must be provided", call. = FALSE)
}

run_maaslin_analysis <- function(input_file) {
  # Maaslin3_path <- ""
  # for (R_file in dir(Maaslin3_path, pattern = "*.R$")) {
  #   source(file.path(Maaslin3_path, R_file))
  # }
  
  taxa_table <- read.csv(input_file, sep = '\t', header = F)
  metadata <- t(taxa_table[1:3,])
  colnames(metadata) <- metadata[1,]
  metadata <- metadata[-1,]
  rownames(metadata) <- paste0("Sample", 1:nrow(metadata))
  metadata <- data.frame(metadata)
  metadata$oxygen_availability <- factor(metadata$oxygen_availability, levels = c('Low_O2', 'Mid_O2', 'High_O2'))
  
  taxa_table <- taxa_table[-c(1:3),]
  rownames_taxa_table <- taxa_table[,1]
  taxa_table <- apply(taxa_table[,-1], 2, function(x){x <- as.numeric(x); x <- ifelse(x == min(x), 0, x)})
  rownames(taxa_table) <- rownames_taxa_table
  colnames(taxa_table) <- paste0("Sample", 1:nrow(metadata))
  taxa_table <- data.frame(taxa_table)
  
fit_out <- maaslin3(input_data = taxa_table, 
                     input_metadata = metadata, 
                     min_abundance = 0, 
                     min_prevalence = 0, 
                     output = 'output/', 
                     min_variance = 0, 
                     normalization = 'NONE', 
                     transform = 'LOG',
                     formula = '~ oxygen_availability + (1 | subject_id)', 
                     plot_associations = FALSE, 
                     save_models = FALSE, 
                     plot_summary_plot = F,
                     max_significance = 0.1, 
                     augment = TRUE, 
                     cores = 6)
  maaslin_write_results_lefse_format(fit_out)
}

# Run the analysis with the provided input file
run_maaslin_analysis(opt$input)
