# Script with helper functions for WT manuscript code
# 15.09.2026
# Franziska Oberhammer

#' Plot betas distinguishing different Infinium chemistries. Original function from sesame was updated to use ggplot!
#'
#' @param sdf SigDF
#' @param title main title in ggplots
#' @param ... additional options to plot
#' @return create a density plot

QC_plotBetaByDesign <- function(
    sdf, title="", ...) {
  
  # get info for plotting for density of beta value, as well as beta values for only Infinium I Red probes, only infinium I green probes and only infinium II probes
  df_data_dA <- data.frame(d = na.omit(sesame:::getBetas(sdf)), type = "All")
  df_data_dR <- data.frame(d = na.omit(sesame:::getBetas(sesame:::InfIR(sdf))), type = "Infinium-I Red")
  df_data_dG <- data.frame(d = na.omit(sesame:::getBetas(sesame:::InfIG(sdf))), type = "Infinium-I Grn")
  df_data_d2 <- data.frame(d = na.omit(sesame:::getBetas(sesame:::InfII(sdf))), type = "Infinium-II")
  
  # make df with all these metrix combined
  df_data <- rbind(
    df_data_dA,
    df_data_dR,
    df_data_dG,
    df_data_d2
  )
  
  # plot density of beta values
  ggplot(data = df_data, aes(x = d, y = ..density.., fill = type, colour = type)) +
    geom_density(alpha = 0.05) +
    guides(fill=guide_legend(title="Probe Type"),
           colour = guide_legend(title="Probe Type"))+
    ylab("Density")+
    xlab("ß-value")+
    ggtitle(title)+
    theme_minimal()
}

#' Plot red-green QQ-Plot using Infinium-I Probes. Original function from sesame was updated to use ggplot!
#'
#' @param sdf a \code{SigDF}
#' @param main plot title
#' @param ... additional options to ggplot
#' @return create a ggplot

QC_plotRedGrnQQ <- function(sdf, main="R-G QQ Plot", ...) {
  dG <- sesame:::InfIG(noMasked(sdf))
  dR <- sesame:::InfIR(noMasked(sdf))
  m <- max(c(dR$MR,dR$UR,dG$MG,dG$UG), na.rm=TRUE)
  
  quantiles = seq(0, 1, 0.00001)
  
  ggplot(mapping = aes(x = quantile(c(dR$MR, dR$UR), quantiles), 
                       y = quantile(c(dG$MG, dG$UG), quantiles))) + 
    geom_point() +
    geom_abline(aes(slope = 1, intercept = 0), linetype = 2) +
    ylim(0,m)+
    xlim(0,m)+
    labs(x = "Infinium-I Red Signal Quantiles", 
         y = "Infinium-I Green Signal Quantiles",
         title = main)+
    theme_minimal()
}



# Enrichment function -----------------------------------------------------
run_chiSquare_or_fisher_enrichment <- function(df, cluster_var, outcome_var,
                                               x_label = NULL,
                                               y_label = NULL,
                                               alpha = 0.05, n_sim = 10000, 
                                               seed = 1234567) {
  
  # Build contingency table dynamically
  tbl <- table(df[[cluster_var]], df[[outcome_var]])
  
  # Determine table size (number of cells)
  n_cells <- nrow(tbl) * ncol(tbl)
  
  do_detailed_tests <- FALSE  # flag for later
  
  if (n_cells <= 10) {
    # Sparse table → Fisher's exact test
    set.seed(seed)
    overall_test <- fisher.test(tbl)
    test_type <- "Fisher exact test"
    
    # Check if overall test is significant
    if (overall_test$p.value <= alpha) {
      do_detailed_tests <- TRUE
    }
    
  } else {
    # Larger table → Chi-squared test with simulation
    set.seed(seed)
    overall_test <- chisq.test(tbl, simulate.p.value = TRUE, B = n_sim)
    test_type <- "Chi-squared test with simulated p-value"
    
    # Check if overall test is significant
    if (overall_test$p.value <= alpha) {
      do_detailed_tests <- TRUE
    }
  }
  
  # Initialize empty outputs for detailed results
  enrichment_df <- NULL
  sig_df <- NULL
  p <- NULL
  per_cell_test <- NULL
  
  # Run detailed tests only if overall test is significant
  if (do_detailed_tests) {
    if (n_cells <= 10) {
      # Fisher per-cluster vs all others
      clusters <- rownames(tbl)
      pvals <- numeric(length(clusters))
      names(pvals) <- clusters
      
      for (i in seq_along(clusters)) {
        cluster <- clusters[i]
        row_cluster <- tbl[cluster,]
        row_other <- colSums(tbl[setdiff(clusters, cluster), , drop = FALSE])

        mini_table <- rbind(row_cluster, row_other)
        rownames(mini_table) <- c(cluster, "Other")
        fisher_result <- fisher.test(mini_table)
        pvals[i] <- fisher_result$p.value
      }
      
      p_adj <- p.adjust(pvals, method = "bonferroni")
      enrichment_df <- data.frame(
        Cluster = clusters,
        Outcome = rep(outcome_var, length(clusters)),
        Count = NA,
        StdResidual = NA,
        PValue = pvals,
        PAdj = p_adj,
        Significant = p_adj <= alpha,
        Direction = NA,
        signif_label = case_when(
          p_adj <= 0.001 ~ "***",
          p_adj <= 0.01  ~ "**",
          p_adj <= 0.05  ~ "*",
          TRUE ~ NA_character_
        )
      )
      
      sig_df <- enrichment_df %>% filter(Significant)
      per_cell_test <- "Fisher per-cluster tests"
      # No plot for Fisher
      p <- NULL
      
    } else {
      # Chi-squared per-cell enrichment
      std_residuals <- overall_test$stdres
      cell_p_values <- 2 * (1 - pnorm(abs(std_residuals)))
      p_adj <- p.adjust(cell_p_values, method = "bonferroni")
      significant_cells <- p_adj <= alpha
      
      enrichment_df <- as.data.frame(as.table(tbl)) %>%
        rename(Cluster = Var1, Outcome = Var2, Count = Freq) %>%
        mutate(
          StdResidual = as.vector(std_residuals),
          PValue = as.vector(cell_p_values),
          PAdj = p_adj,
          Significant = significant_cells,
          Direction = ifelse(StdResidual > 0, "Enriched", "Depleted"),
          signif_label = case_when(
            p_adj <= 0.001 ~ "***",
            p_adj <= 0.01  ~ "**",
            p_adj <= 0.05  ~ "*",
            TRUE           ~ NA_character_
          )
        )
      
      sig_df <- enrichment_df %>% filter(Significant)
      
      if (is.null(y_label)) y_label <- outcome_var
      if (is.null(x_label)) x_label <- cluster_var
      p <- ggplot(enrichment_df, aes(x = Cluster, y = Outcome, fill = StdResidual)) +
        geom_tile() +
        geom_tile(data = subset(enrichment_df, Significant),
                  aes(x = Cluster, y = Outcome),
                  fill = NA, color = "black", size = 0.75) +
        geom_text(aes(label = round(StdResidual, 2)), color = "black") +
        scale_fill_gradient2(
          low = "blue", mid = "white", high = "red", midpoint = 0,
          name = "Std. residual"
        ) +
        labs(x = x_label, y = y_label) +
        scale_y_discrete(limits = rev) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        theme(legend.key.width = unit(0.5, "cm")) +
        theme_Publication()
      
      per_cell_test <- "Chi-squared per-cell enrichment"
    }
  }
  
  # Return list in desired order:
  return(list(
    overall_test_type = test_type,
    overall_test_object = overall_test,
    overall_p_value = overall_test$p.value,
    overall_significant = do_detailed_tests,
    enrichment_df = enrichment_df,
    sig_df = sig_df,
    plot = p,
    per_cell_test = per_cell_test, 
    seed = seed
  ))
}



# ----------confusion matrix for Classifier validation ----------

# Generates a confusion-matrix summary table without grouping variable
make_conf_mat_simple <- function(data, original_col = "original", predicted_col = "predict", probability_col = "probability") {
  
  original_sym  <- rlang::sym(original_col)
  predicted_sym <- rlang::sym(predicted_col)
  probability_sym <- rlang::sym(probability_col)
  
  data %>%
    group_by(!!original_sym, !!predicted_sym) %>%
    summarise(
      n = n(),
      avg_prob = round(mean(!!probability_sym), 1),
      .groups = "drop"
    ) %>%
    group_by(!!original_sym) %>%
    mutate(freq = n / sum(n)) %>%
    ungroup() %>%
    mutate(
      label  = paste0(n, " (", avg_prob, "%)"),
      label2 = paste0(n, " (", scales::percent(freq, accuracy = 1), ")")
    )
}


# corresponding confusion matrix plot for confusion matrix of classifier validation
plot_conf_mat <- function(
    conf_mat,
    title = "Confusion Matrix",
    label_type = c("count", "prob", "relfreq", "count_prob", "count_relfreq"),
    high_colour = "darkorange3",
    original_col = "originalCall",
    predicted_col = "predict"
) {
  
  # Match argument
  label_type <- match.arg(label_type)
  
  # Subtitle text based on label_type
  subtitle_text <- dplyr::case_when(
    label_type == "count"          ~ "Displayed: Count per cell (n)",
    label_type == "prob"           ~ "Displayed: Average probability (%)",
    label_type == "relfreq"        ~ "Displayed: Relative frequency (%)",
    label_type == "count_prob"     ~ "Displayed: Count (n) and average probability (%)",
    label_type == "count_relfreq"  ~ "Displayed: Count (n) and relative frequency (%)"
  )
  
  # Compose label content
  conf_mat <- conf_mat %>%
    mutate(
      display_label = dplyr::case_when(
        label_type == "count"          ~ as.character(n),
        label_type == "prob"           ~ paste0(avg_prob, "%"),
        label_type == "relfreq"        ~ scales::percent(freq, accuracy = 1),
        label_type == "count_prob"     ~ paste0(n, " (", avg_prob, "%)"),
        label_type == "count_relfreq"  ~ paste0(n, " (", scales::percent(freq, accuracy = 1), ")")
      )
    )
  
  original_sym <- rlang::sym(original_col)
  predicted_sym <- rlang::sym(predicted_col)
  
  ggplot(conf_mat, aes(
    x = !!predicted_sym,
    y = !!original_sym,
    fill = freq * 100
  )) +
    geom_tile(color = "white") +
    geom_text(aes(label = display_label), size = 4) +
    scale_fill_gradient(low = "white", high = high_colour, limits = c(0, 100)) +
    labs(
      title = title,
      subtitle = subtitle_text,
      x = "Predicted Cluster",
      y = "Original Cluster",
      fill = "Relative Frequency"
    ) +
    theme_Publication() +
    theme(
      panel.background = element_rect(fill = "white", colour = "white"),
      legend.key.width = unit(1, 'cm')
    )
}


# summarising classifier prediction for multiple samples per patient -----------
summarise_patient_predictions <- function(df, 
                                          id_col = "M_Label",
                                          predict_col = "predict",
                                          prob_col = "probability",
                                          prob_threshold = 73) {
  
  dplyr::mutate(
    df,
    patient = stringr::str_remove(!!rlang::sym(id_col), ".$"),
    .prob_num = as.numeric(!!rlang::sym(prob_col))
  ) %>%
    dplyr::group_by(patient) %>%
    dplyr::group_modify(~ {
      patient_df <- .x
      
      # if only 1 sample → keep it
      if (nrow(patient_df) == 1) {
        used <- patient_df[[id_col]]
        discarded <- character(0)
        final_class <- patient_df[[predict_col]]
        avg_prob <- patient_df$.prob_num
        n_votes <- 1
        n_total <- 1
      } else {
        # multiple samples → filter by threshold
        above <- dplyr::filter(patient_df, .prob_num >= prob_threshold)
        
        if (nrow(above) == 0) {
          # no sample above threshold → keep all
          used <- patient_df[[id_col]]
          discarded <- character(0)
          filtered <- patient_df
        } else {
          used <- above[[id_col]]
          discarded <- setdiff(patient_df[[id_col]], used)
          filtered <- above
        }
        
        # majority vote
        votes <- dplyr::count(filtered, !!rlang::sym(predict_col), sort = TRUE)
        
        if (length(unique(votes$n)) == 1) {
          # tie → pick sample with highest probability
          final_class <- dplyr::slice(
            dplyr::arrange(filtered, dplyr::desc(.prob_num)),
            1
          ) %>%
            dplyr::pull(!!rlang::sym(predict_col))
        } else {
          # clear majority wins
          final_class <- votes[[predict_col]][1]
        }
        
        # number of votes supporting final prediction
        n_votes <- sum(filtered[[predict_col]] == final_class)
        n_total <- nrow(filtered)
        
        # average probability of used samples
        avg_prob <- mean(filtered$.prob_num)
      }
      
      # return one row per patient
      tibble::tibble(
        used_samples = paste(used, collapse = ", "),
        discarded_samples = paste(discarded, collapse = ", "),
        final_prediction = final_class,
        avg_probability = avg_prob,
        votes = paste0(n_votes, "/", n_total)
      )
    }) %>%
    dplyr::ungroup()
}

