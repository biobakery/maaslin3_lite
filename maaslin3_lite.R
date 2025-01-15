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

find_indices <- function(vec) {
    result <- list()
    
    start_idx <- 1
    
    for (i in 2:length(vec)) {
        # Check if the current value is different from the previous one
        if (vec[i] != vec[i - 1]) {
            result[[as.character(vec[start_idx])]] <- c(start_idx, i)
            start_idx <- i
        }
    }
    
    result[[as.character(vec[start_idx])]] <- c(start_idx, length(vec) + 1)
    
    result <- lapply(result, FUN = function(x){x - 1})
    
    return(result)
}

# Function to create the named list
make_hierarchy <- function(v1, v2) {
    result <- list()
    
    current_name <- v1[1]
    current_values <- v2[1]
    names_vec <- c()
    
    for (i in 2:length(v1)) {
        # If the value in v1 changes, save the current list entry and start a new one
        if (v1[i] != current_name) {
            result[[current_name]] <- current_values
            names_vec <- c(names_vec, as.character(current_name))
            current_name <- v1[i]
            current_values <- v2[i]
        } else {
            # Otherwise, append the value to the current vector
            current_values <- c(current_values, v2[i])
        }
    }
    
    result[[current_name]] <- current_values
    names_vec <- c(names_vec, as.character(current_name))
    
    names(result) <- names_vec
    return(result)
}

list_to_json <- function(x) {
    json_str <- "{"
    names_x <- names(x)
    
    for (i in 1:length(x)) {
        name <- names_x[i]
        value <- x[[i]]
        json_str <- paste0(json_str, '"', name, '": ')
        
        if (is.list(value)) {
            # Recursively call list_to_json if it's a list
            json_str <- paste0(json_str, list_to_json(value))
        } else {
            if (name != 'norm') {
                # Otherwise, directly add the value
                if (!is.numeric(value)) {
                    value <- paste0('"', value, '"')
                    json_str <- paste0(json_str, '[', paste(value, collapse = ', '), ']')
                } else {
                    json_str <- paste0(json_str, '[', paste(format(value, scientific = FALSE), collapse = ', '), ']')
                }
            } else {
                json_str <- paste0(json_str, value)
            }
        }
        
        # Add a comma separator for the next item, unless it's the last one
        if (i < length(x)) {
            json_str <- paste0(json_str, ", ")
        }
    }
    
    json_str <- paste0(json_str, "}")
    
    return(json_str)
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
  
  if (!is.null(subclass)) {
      taxa_table <- taxa_table[,order(metadata[[class]], metadata[[subclass]])]
      metadata <- metadata[order(metadata[[class]], metadata[[subclass]]),]
  } else {
      taxa_table <- taxa_table[,order(metadata[[class]])]
      metadata <- metadata[order(metadata[[class]]),]
  }
  
  # Set normalization method based on user input
  normalization_method <- ifelse(normalize, 'TSS', 'NONE')
  
  # Create formula dynamically based on user input
  if (is.null(subclass)) {
    formula_str <- paste0("~ ", class)
  } else {
    formula_str <- paste0("~ ", class, " + ", subclass)
  }
  
  if (!is.null(random_component)) {
      formula_str <- paste0(formula_str, " + (1 | ", random_component, ")")
  }
  
  # Determine delimiter for taxonomic naming
  bar_count <- sapply(strsplit(rownames(taxa_table), "\\|"), length) - 1
  period_count <- sapply(strsplit(rownames(taxa_table), "\\."), length) - 1
  delimiter <- ifelse(mean(bar_count) > 1 | mean(period_count) > 1, 
                      ifelse(mean(bar_count) > mean(period_count), '|', '.'),
                      NA)
  if (!is.na(delimiter)) {
      tax_levels <- max(sapply(strsplit(rownames(taxa_table), delimiter, fixed = T), length))
      
      fit_out_growing <- list()
      fit_out_growing$fit_data_abundance$results <- data.frame()
      fit_out_growing$fit_data_prevalence$results <- data.frame()
      for (tax_level in 1:tax_levels) {
          taxa_table_tmp <- taxa_table[sapply(strsplit(rownames(taxa_table), delimiter, fixed = T), length) == tax_level,, drop=F]
          
          fit_out <- maaslin3(input_data = taxa_table_tmp, 
                              input_metadata = metadata, 
                              min_abundance = 0, 
                              min_prevalence = 0, 
                              output = 'output/', 
                              min_variance = -1, 
                              normalization = normalization_method, 
                              transform = 'LOG',
                              formula = formula_str, 
                              plot_associations = FALSE, 
                              save_models = FALSE, 
                              plot_summary_plot = F,
                              max_significance = alpha_threshold, 
                              subtract_median = TRUE,
                              augment = TRUE,
                              warn_prevalence = F,
                              cores = 1)
          
          fit_out_growing$fit_data_abundance$results <- rbind(fit_out_growing$fit_data_abundance$results,
                                                              fit_out$fit_data_abundance$results)
          fit_out_growing$fit_data_prevalence$results <- rbind(fit_out_growing$fit_data_prevalence$results,
                                                              fit_out$fit_data_prevalence$results)
      }
      fit_out <- fit_out_growing
  } else {
      tax_levels <- 1
      
      fit_out <- maaslin3(input_data = taxa_table, 
                          input_metadata = metadata, 
                          min_abundance = 0, 
                          min_prevalence = 0, 
                          output = 'output/', 
                          min_variance = -1, 
                          normalization = normalization_method, 
                          transform = 'LOG',
                          formula = formula_str, 
                          plot_associations = FALSE, 
                          save_models = FALSE, 
                          plot_summary_plot = F,
                          max_significance = alpha_threshold, 
                          subtract_median = TRUE,
                          augment = TRUE, 
                          warn_prevalence = F,
                          cores = 1)
  }
  
  # Save results in LEfSe format
  maaslin_write_results_lefse_format('output', fit_out$fit_data_abundance, fit_out$fit_data_prevalence)
  
  # Write second file for plotting
  feats <- apply(taxa_table, 1, function(row) unname(unlist(row)), simplify = F)
  names(feats) <- gsub("\\|", ".", rownames(taxa_table))
  cls_list <- list()
  if (!is.null(random_component)) {
      cls_list$subject <- as.character(metadata[[random_component]])
  }
  cls_list$class = as.character(metadata[[class]])
  if (!is.null(subclass)) {
      cls_list$subclass <- as.character(metadata[[subclass]])
  }
  class_sl = find_indices(metadata[[class]])
  
  output <- list(feats = feats,
       norm = ifelse(normalize, 1, max(colSums(taxa_table)) / tax_levels),
       cls = cls_list,
       class_sl = class_sl)
  
  if (!is.null(subclass)) {
      output$subclass_sl <- find_indices(metadata[[subclass]])
      output$class_hierarchy = make_hierarchy(metadata[[class]], metadata[[subclass]])
  }
  
  json_out <- list_to_json(output)
  writeLines(json_out, 'output/format_data.lefse_internal_for')
  
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







