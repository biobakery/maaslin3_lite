#!/usr/bin/Rscript
library(optparse)
library(maaslin3)

# Define the command line options
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Path to input txt file", metavar = "character"),
  make_option(c("-n", "--normalize"), type = "logical", default = FALSE,
              help = "Whether to normalize the data (TRUE for TSS, FALSE for NONE)", metavar = "logical"),
  make_option(c("-c", "--class"), type = "character", default = "oxygen_availability",
              help = "Class for the fixed effect in the formula", metavar = "character"),
  make_option(c("-s", "--subclass"), type = "character",
              help = "Subclass for the fixed effect in the formula", metavar = "character"),
  make_option(c("-r", "--random_component"), type = "character", default = "subject_id",
              help = "Random component for the formula (e.g., subject_id)", metavar = "character"),
  make_option(c("-a", "--alpha_threshold"), type = "character", default = 0.1,
              help = "Maximum FDR corrected significance level", metavar = "numeric")
)

# Parse the command line options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("Input file must be provided", call. = FALSE)
}

run_maaslin_analysis <- function(input_file, normalize, class, subclass, random_component, alpha_threshold) {
  # Read input data
  taxa_table <- read.csv(input_file, sep = '\t', header = F)
  
  # Metadata setup
  metadata <- t(taxa_table[1:3,])
  colnames(metadata) <- metadata[1,]
  metadata <- metadata[-1,]
  rownames(metadata) <- paste0("Sample", 1:nrow(metadata))
  metadata <- data.frame(metadata)
  metadata[[class]] <- factor(metadata[[class]])
  metadata[[subclass]] <- factor(metadata[[subclass]])
  
  # Process taxa table
  taxa_table <- taxa_table[-c(1:3),]
  rownames_taxa_table <- taxa_table[,1]
  taxa_table <- apply(taxa_table[,-1], 2, function(x){x <- as.numeric(x); x <- ifelse(x == min(x), 0, x)})
  rownames(taxa_table) <- rownames_taxa_table
  colnames(taxa_table) <- paste0("Sample", 1:nrow(metadata))
  taxa_table <- data.frame(taxa_table)
  
  # Set normalization method based on user input
  normalization_method <- ifelse(normalize, 'TSS', 'NONE')
  
  # Create formula dynamically based on user input
  if (is.null(subclass)) {
    formula_str <- paste0("~ ", class, " + (1 | ", random_component, ")")
  } else {
    formula_str <- paste0("~ ", class, " + ", subclass, " + (1 | ", random_component, ")")
  }
  
  
  # Run Maaslin3 analysis
  fit_out <- maaslin3(input_data = taxa_table, 
                      input_metadata = metadata, 
                      min_abundance = 0, 
                      min_prevalence = 0, 
                      output = 'output/', 
                      min_variance = 0, 
                      normalization = normalization_method, 
                      transform = 'LOG',
                      formula = formula_str, 
                      plot_associations = FALSE, 
                      save_models = FALSE, 
                      plot_summary_plot = F,
                      max_significance = alpha_threshold, 
                      augment = TRUE, 
                      cores = 6)
  
  # Save results in LEfSe format
  maaslin_write_results_lefse_format('output', fit_out$fit_data_abundance, fit_out$fit_data_prevalence)
}

# Run the analysis with the provided options
run_maaslin_analysis(opt$input, opt$normalize, opt$class, opt$subclass, opt$random_component, as.numeric(opt$alpha_threshold))
