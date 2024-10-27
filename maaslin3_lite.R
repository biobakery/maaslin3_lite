#!/usr/bin/Rscript
library(optparse)
library(maaslin3)
library(dplyr)
library(ggplot2)

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
                      subtract_median = TRUE,
                      augment = TRUE, 
                      cores = 1)
  
  # Save results in LEfSe format
  maaslin_write_results_lefse_format('output', fit_out$fit_data_abundance, fit_out$fit_data_prevalence)

  return(fit_out)
}

run_make_coef_plot <- function(merged_results_sig, 
                           max_significance,
                           class) {
    coef_plot_vars <- class
    
    coef_plot_data <- merged_results_sig[merged_results_sig$metadata %in% coef_plot_vars,]
    
    # Limit plotted coefficients to median +/- 10 times distance to quartiles
    quantile_df <- coef_plot_data %>%
        dplyr::group_by(.data$full_metadata_name) %>%
        dplyr::summarise(
            lower_q = median(.data$coef) - 10 * 
                (median(.data$coef) - quantile(.data$coef, 0.25)),
            upper_q = median(.data$coef) + 10 * 
                (quantile(.data$coef, 0.75) - median(.data$coef))
        ) %>%
        data.frame()
    rownames(quantile_df) <- quantile_df$full_metadata_name
    
    # Make sure insignificant coefficients don't distort the plot
    coef_plot_data <-
        coef_plot_data[coef_plot_data$qval_individual < 
                           max_significance |
                           (coef_plot_data$coef > quantile_df[
                               coef_plot_data$full_metadata_name, 
                               'lower_q'] &
                                coef_plot_data$coef < quantile_df[
                                    coef_plot_data$full_metadata_name, 
                                    'upper_q']),]
    
    # Choose breaks for plot
    custom_break_fun <- function(n) {
        return(function(x) {
            extended_breaks <- scales::breaks_extended(n)(x)
            if (max(x) > 0) {
                extended_breaks <- extended_breaks[
                    extended_breaks <= max(x) * 0.9]
            } else {
                extended_breaks <- extended_breaks[
                    extended_breaks <= max(x) * 1.1]
            }
            if (min(x) > 0) {
                extended_breaks <- extended_breaks[
                    extended_breaks >= min(x) * 1.1]
            } else {
                extended_breaks <- extended_breaks[
                    extended_breaks >= min(x) * 0.9]
            }
            extended_breaks
        })
    }
    
    # Plot
    p1 <- ggplot(coef_plot_data, aes(x = feature, y = coef, fill = value, alpha = model)) +
        geom_bar(stat = "identity", 
                 position = position_dodge(width = 0.8), 
                 width = 0.7) +
        geom_errorbar(aes(ymin = coef - stderr, ymax = coef + stderr),
                      position = position_dodge(width = 0.8), 
                      width = 0.25) +
        scale_fill_brewer(palette = "Dark2") +
        labs(x = "Feature", y = "Coefficient", fill = "Class", alpha = "Model") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        theme(panel.spacing.x = unit(0.5, "lines")) + 
        coord_flip() + 
        facet_wrap(
            ~ value,
            scales = 'free_x',
            ncol = length(coef_plot_vars)
        ) + 
        scale_alpha_manual(values = c('Abundance' = 1, 'Prevalence' = 0.5)) + 
        scale_y_continuous(
            breaks = custom_break_fun(n = 6),
            limits = c(
                min(coef_plot_data$coef) - 
                    quantile(coef_plot_data$stderr, 0.8),
                max(coef_plot_data$coef) + 
                    quantile(coef_plot_data$stderr, 0.8)
            )
        )
    
    return(p1)
}


run_maaslin_plotting <- function(fit_out, class, alpha_threshold, first_n = 20) {
    fit_data_abundance <- fit_out$fit_data_abundance
    fit_data_prevalence <- fit_out$fit_data_prevalence
    
    if (is.null(fit_data_abundance$results)) {
        merged_results <- fit_data_prevalence$results
    } else if (is.null(fit_data_prevalence$results)) {
        merged_results <- fit_data_abundance$results
    } else {
        merged_results <- rbind(fit_data_abundance$results,
                                fit_data_prevalence$results)
    }
    
    merged_results <- maaslin3:::preprocess_merged_results(merged_results)
    
    # Subset associations for plotting
    merged_results_joint_only <-
        unique(merged_results[, c('feature', 'qval_joint')])
    merged_results_joint_only <-
        merged_results_joint_only[
            order(merged_results_joint_only$qval_joint),]
    if (length(unique(merged_results_joint_only$feature)) < first_n) {
        first_n <- length(unique(merged_results_joint_only$feature))
    }
    signif_taxa <-
        unique(merged_results_joint_only$feature)[seq(first_n)]
    
    merged_results_sig <- merged_results %>%
        dplyr::filter(.data$feature %in% signif_taxa)
    
    p1 <- run_make_coef_plot(merged_results_sig = merged_results_sig, 
                              max_significance = alpha_threshold, 
                       class = class)
    
    height_out <-
        5.5 + max(first_n / 5 - 5, 0) + nchar(as.character(class)) / 10
    width_out <-
        2 + max(nchar(merged_results$feature)) / 15
    
    ggplot2::ggsave(
        filename = 'output/summary_plot.png',
        plot = p1, 
        dpi = 600,
        width = width_out,
        height = height_out)
}

# Run the analysis with the provided options
fit_out <- run_maaslin_analysis(opt$input, opt$normalize, opt$class, opt$subclass, opt$random_component, as.numeric(opt$alpha_threshold))
run_maaslin_plotting(fit_out, opt$class, as.numeric(opt$alpha_threshold), first_n = 20)







