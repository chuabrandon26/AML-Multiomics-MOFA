# Notes:
#   - The script assumes the four assignment files are in project_dir:
#       rna_count.tsv, gene_mutation.tsv, drug_auc.tsv, metadata.tsv
#   - It creates a local Python virtual environment named .venv_mofa if needed,
#     so MOFA2 can find mofapy2 without depending on a broken global Python setup.

project_dir <- "C:/Users/chuab/Desktop/Multi_Omics_folder/multi_omics_report"
qmd_name <- "MoBi_MultiOmic_4745778_Chua.qmd"
html_name <- "MoBi_MultiOmic_4745778_Chua.html"
qmd_file <- file.path(project_dir, qmd_name)
html_file <- file.path(project_dir, html_name)

# Set these to FALSE if you only want to regenerate the .qmd without rendering.
setup_mofa_python <- TRUE
render_report <- TRUE
overwrite_qmd <- TRUE

if (!dir.exists(project_dir)) {
  stop("Project directory does not exist: ", project_dir)
}
setwd(project_dir)

required_input_files <- c("rna_count.tsv", "gene_mutation.tsv", "drug_auc.tsv", "metadata.tsv")
missing_input_files <- required_input_files[!file.exists(file.path(project_dir, required_input_files))]
if (length(missing_input_files) > 0) {
  stop("Missing assignment input file(s): ", paste(missing_input_files, collapse = ", "))
}

find_python <- function() {
  candidates <- c(
    Sys.which("python"),
    "C:/Users/chuab/anaconda3/python.exe",
    "C:/Users/chuab/AppData/Local/Programs/Python/Python313/python.exe"
  )
  candidates <- unique(candidates[nzchar(candidates)])
  candidates[file.exists(candidates)][1]
}

run_command <- function(command, args = character(), label = command) {
  message("Running: ", label)
  output <- system2(command, args = args, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0) {
    cat(paste(output, collapse = "\n"), "\n")
    stop("Command failed: ", label, call. = FALSE)
  }
  invisible(output)
}

setup_mofa_backend <- function(project_dir) {
  venv_dir <- file.path(project_dir, ".venv_mofa")
  venv_python <- file.path(venv_dir, "Scripts", "python.exe")

  if (!file.exists(venv_python)) {
    base_python <- find_python()
    if (is.na(base_python) || !file.exists(base_python)) {
      stop("Could not find Python to create the MOFA virtual environment.")
    }
    run_command(base_python, c("-m", "venv", shQuote(venv_dir)), "create .venv_mofa")
    run_command(venv_python, c("-m", "pip", "install", "--upgrade", "pip"), "upgrade pip")
  }

  backend_check <- suppressWarnings(system2(
    venv_python,
    args = c("-c", shQuote("import numpy, scipy, pandas, h5py, mofapy2")),
    stdout = TRUE,
    stderr = TRUE
  ))
  backend_status <- attr(backend_check, "status")
  if (!is.null(backend_status) && backend_status != 0) {
    # Version pins keep the MOFA Python stack stable on Windows/RStudio.
    run_command(
      venv_python,
      c(
        "-m", "pip", "install",
        "numpy<2", "scipy<1.13", "pandas<3", "scikit-learn<1.6", "h5py<4", "mofapy2"
      ),
      "install MOFA Python backend"
    )
  } else {
    message("MOFA Python backend already available in: ", venv_python)
  }

  Sys.setenv(RETICULATE_PYTHON = normalizePath(venv_python, winslash = "/"))
  invisible(venv_python)
}

find_quarto <- function() {
  candidates <- c(
    Sys.which("quarto"),
    "C:/Program Files/RStudio/resources/app/bin/quarto/bin/quarto.exe",
    "C:/Program Files/Quarto/bin/quarto.exe"
  )
  candidates <- unique(candidates[nzchar(candidates)])
  quarto_bin <- candidates[file.exists(candidates)][1]
  if (is.na(quarto_bin) || !file.exists(quarto_bin)) {
    stop("Could not find Quarto. Install Quarto or open RStudio with bundled Quarto available.")
  }
  normalizePath(quarto_bin, winslash = "/")
}


qmd_text <- r'----------------------------------------(
---
title: "MoBi MultiOmic Assignment"
author: "Chua"
format:
  html:
    toc: true
    number-sections: true
    code-fold: show
---
```{r setup}
# Reproducibility and package setup
work_dir <- "C:/Users/chuab/Desktop/Multi_Omics_folder/multi_omics_report"
if (!dir.exists(work_dir)) {
  stop("Working directory was not found: ", work_dir)
}
setwd(work_dir)

seed_value <- 4745778
set.seed(seed_value)
options(stringsAsFactors = FALSE, scipen = 999, timeout = 600)

candidate_pythons <- c(
  file.path(work_dir, ".venv_mofa", "Scripts", "python.exe"),
  "C:/Users/chuab/anaconda3/python.exe"
)
candidate_python <- candidate_pythons[file.exists(candidate_pythons)][1]
if (!is.na(candidate_python)) {
  Sys.setenv(RETICULATE_PYTHON = candidate_python)
}

cran_packages <- c(
  "tidyverse", "ggplot2", "pheatmap", "matrixStats", "patchwork",
  "scales", "janitor", "readr", "tibble", "tidyr", "dplyr", "forcats",
  "uwot", "FNN", "igraph", "ggrepel", "broom"
)
bioc_packages <- c("DESeq2", "MOFA2")

missing_cran <- cran_packages[!vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran) > 0) {
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
missing_bioc <- bioc_packages[!vapply(bioc_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc) > 0) {
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(MOFA2)
  library(ggplot2)
  library(pheatmap)
  library(matrixStats)
  library(patchwork)
  library(scales)
  library(janitor)
  library(readr)
  library(tibble)
  library(tidyr)
  library(dplyr)
  library(forcats)
  library(uwot)
  library(FNN)
  library(igraph)
  library(ggrepel)
  library(broom)
})

theme_report <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 2),
      plot.subtitle = element_text(color = "grey30"),
      axis.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      strip.text = element_text(face = "bold"),
      plot.caption = element_text(color = "grey35", hjust = 0)
    )
}
theme_set(theme_report())

knitr::opts_chunk$set(
  warning = FALSE,
  message = FALSE,
  error = FALSE,
  fig.align = "center",
  fig.width = 8,
  fig.height = 5,
  dpi = 120
)

format_int <- function(x) comma(as.integer(x), accuracy = 1)
format_num <- function(x, digits = 2) {
  ifelse(is.na(x), "NA", format(round(x, digits), big.mark = ",", nsmall = digits, trim = TRUE))
}
format_pct <- function(x, accuracy = 0.1) {
  percent(x, accuracy = accuracy)
}
format_p <- function(p) {
  case_when(
    is.na(p) ~ "NA",
    p < 0.001 ~ "<0.001",
    TRUE ~ formatC(p, format = "f", digits = 3)
  )
}

empty_plot <- function(title, message) {
  ggplot() +
    annotate("text", x = 0, y = 0, label = message, size = 4, color = "grey30") +
    xlim(-1, 1) +
    ylim(-1, 1) +
    labs(title = title) +
    theme_void(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
}

row_center_scale <- function(mat) {
  mat <- as.matrix(mat)
  row_means <- rowMeans(mat, na.rm = TRUE)
  row_sds <- matrixStats::rowSds(mat, na.rm = TRUE)
  row_sds[is.na(row_sds) | row_sds == 0] <- 1
  sweep(sweep(mat, 1, row_means, "-"), 1, row_sds, "/")
}

impute_row_mean <- function(mat) {
  mat <- as.matrix(mat)
  for (i in seq_len(nrow(mat))) {
    missing_i <- is.na(mat[i, ])
    if (any(missing_i)) {
      replacement <- mean(mat[i, ], na.rm = TRUE)
      if (is.nan(replacement)) replacement <- 0
      mat[i, missing_i] <- replacement
    }
  }
  mat
}

run_pca <- function(feature_by_sample, scale_features = FALSE, n_pcs = 10) {
  mat <- as.matrix(feature_by_sample)
  mat <- mat[rowSums(!is.na(mat)) >= 2, , drop = FALSE]
  if (anyNA(mat)) {
    stop("PCA input contains missing values. Impute or filter before PCA.")
  }
  vars <- matrixStats::rowVars(mat)
  mat <- mat[is.finite(vars) & vars > 0, , drop = FALSE]
  if (nrow(mat) < 2 || ncol(mat) < 3) {
    return(list(ok = FALSE, reason = "Too few non-constant features or samples for PCA."))
  }
  pca <- prcomp(t(mat), center = TRUE, scale. = scale_features)
  variance <- tibble(
    PC = seq_along(pca$sdev),
    variance = pca$sdev^2 / sum(pca$sdev^2),
    variance_percent = 100 * variance
  )
  keep_pcs <- seq_len(min(n_pcs, ncol(pca$x)))
  scores <- as_tibble(pca$x[, keep_pcs, drop = FALSE], rownames = "sampleID")
  list(ok = TRUE, pca = pca, variance = variance, scores = scores, matrix = mat)
}

metadata_summary <- function(metadata) {
  metadata %>%
    summarise(across(
      -sampleID,
      list(
        non_missing = ~sum(!is.na(.x)),
        n_distinct = ~n_distinct(.x, na.rm = TRUE)
      ),
      .names = "{.col}__{.fn}"
    )) %>%
    pivot_longer(everything(), names_to = "metric", values_to = "value") %>%
    separate(metric, into = c("variable", "summary"), sep = "__") %>%
    pivot_wider(names_from = summary, values_from = value) %>%
    mutate(
      type = map_chr(variable, ~{
        x <- metadata[[.x]]
        if (is.numeric(x) && n_distinct(x, na.rm = TRUE) > 8) "numeric" else "categorical"
      }),
      example_values = map_chr(variable, ~{
        vals <- metadata[[.x]]
        vals <- vals[!is.na(vals)]
        vals <- unique(as.character(vals))
        paste(head(vals, 4), collapse = ", ")
      })
    ) %>%
    arrange(type, desc(n_distinct))
}

select_plot_variables <- function(metadata) {
  priorities <- c(
    "disease_stage", "response_to_induction_therapy", "flt3_itd", "npm1",
    "specimen_type", "vital_status", "sex", "ethnicity", "center_id"
  )
  available <- priorities[priorities %in% names(metadata)]
  usable <- available[vapply(available, function(v) {
    n <- n_distinct(metadata[[v]], na.rm = TRUE)
    n >= 2 && n <= 8
  }, logical(1))]
  color_var <- if (length(usable) > 0) usable[[1]] else NULL
  shape_candidates <- usable[usable != color_var]
  shape_var <- if (length(shape_candidates) > 0) {
    shape_usable <- shape_candidates[vapply(shape_candidates, function(v) {
      n_distinct(metadata[[v]], na.rm = TRUE) <= 6
    }, logical(1))]
    if (length(shape_usable) > 0) shape_usable[[1]] else NULL
  } else {
    NULL
  }
  list(color = color_var, shape = shape_var)
}

make_pca_plot <- function(pca_result, metadata, title, color_var = NULL, shape_var = NULL) {
  if (!isTRUE(pca_result$ok)) {
    return(empty_plot(title, pca_result$reason))
  }
  plot_df <- pca_result$scores %>% left_join(metadata, by = "sampleID")
  x_lab <- paste0("PC1 (", format_num(pca_result$variance$variance_percent[1], 1), "%)")
  y_lab <- paste0("PC2 (", format_num(pca_result$variance$variance_percent[2], 1), "%)")
  p <- ggplot(plot_df, aes(x = PC1, y = PC2))
  if (!is.null(color_var) && color_var %in% names(plot_df) &&
      !is.null(shape_var) && shape_var %in% names(plot_df)) {
    p <- p + geom_point(aes(color = .data[[color_var]], shape = .data[[shape_var]]), size = 2.6, alpha = 0.9)
  } else if (!is.null(color_var) && color_var %in% names(plot_df)) {
    p <- p + geom_point(aes(color = .data[[color_var]]), size = 2.6, alpha = 0.9)
  } else {
    p <- p + geom_point(color = "#2B6CB0", size = 2.6, alpha = 0.9)
  }
  p +
    labs(
      title = title,
      x = x_lab,
      y = y_lab,
      color = if (!is.null(color_var)) str_replace_all(color_var, "_", " ") else NULL,
      shape = if (!is.null(shape_var)) str_replace_all(shape_var, "_", " ") else NULL
    ) +
    theme_report()
}

make_scree_plot <- function(pca_result, title, n = 10) {
  if (!isTRUE(pca_result$ok)) {
    return(empty_plot(title, pca_result$reason))
  }
  pca_result$variance %>%
    slice_head(n = n) %>%
    mutate(PC = factor(paste0("PC", PC), levels = paste0("PC", PC))) %>%
    ggplot(aes(x = PC, y = variance_percent)) +
    geom_col(fill = "#2B6CB0", width = 0.75) +
    geom_text(aes(label = paste0(format_num(variance_percent, 1), "%")), vjust = -0.25, size = 3) +
    labs(title = title, x = NULL, y = "Variance explained (%)") +
    theme_report() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

top_pc_outliers <- function(pca_result, n = 3) {
  if (!isTRUE(pca_result$ok)) return("not available")
  pca_result$scores %>%
    mutate(pc_distance = sqrt(PC1^2 + PC2^2)) %>%
    arrange(desc(pc_distance)) %>%
    slice_head(n = n) %>%
    pull(sampleID) %>%
    paste(collapse = ", ")
}

choose_pc_cols <- function(pca_result, max_pcs = 10) {
  if (!isTRUE(pca_result$ok)) return(character())
  pc_cols <- names(pca_result$scores)[str_detect(names(pca_result$scores), "^PC[0-9]+$")]
  head(pc_cols, max_pcs)
}

association_tests <- function(scores, metadata, value_cols, id_col = "sampleID") {
  empty_result <- tibble(
    component = character(),
    metadata = character(),
    method = character(),
    estimate = numeric(),
    p_value = numeric(),
    p_adj = numeric()
  )
  if (!id_col %in% names(scores) || length(value_cols) == 0) return(empty_result)
  df <- scores %>%
    select(all_of(id_col), all_of(value_cols)) %>%
    dplyr::rename(sampleID = all_of(id_col)) %>%
    left_join(metadata, by = "sampleID")
  meta_cols <- setdiff(names(metadata), "sampleID")
  tests <- map_dfr(value_cols, function(value_col) {
    map_dfr(meta_cols, function(meta_col) {
      y <- df[[value_col]]
      x <- df[[meta_col]]
      keep <- !is.na(y) & !is.na(x)
      y <- y[keep]
      x <- x[keep]
      if (length(y) < 8 || n_distinct(x) < 2) return(tibble())
      if (is.numeric(x) && n_distinct(x) > 8) {
        test_df <- tibble(y = y, x = as.numeric(x))
        fit <- tryCatch(lm(y ~ x, data = test_df), error = function(e) NULL)
        if (is.null(fit)) return(tibble())
        aov_tab <- tryCatch(anova(fit), error = function(e) NULL)
        if (is.null(aov_tab) || nrow(aov_tab) < 1) return(tibble())
        tibble(
          component = value_col,
          metadata = meta_col,
          method = "ANOVA linear model",
          estimate = unname(coef(fit)[2]),
          p_value = aov_tab[["Pr(>F)"]][1]
        )
      } else {
        x <- droplevels(as.factor(x))
        if (nlevels(x) < 2 || nlevels(x) > 10 || any(table(x) < 2)) return(tibble())
        test_df <- tibble(y = y, x = x)
        fit <- tryCatch(aov(y ~ x, data = test_df), error = function(e) NULL)
        if (is.null(fit)) return(tibble())
        aov_tab <- tryCatch(summary(fit)[[1]], error = function(e) NULL)
        if (is.null(aov_tab) || nrow(aov_tab) < 1) return(tibble())
        eta_sq <- aov_tab[["Sum Sq"]][1] / sum(aov_tab[["Sum Sq"]], na.rm = TRUE)
        tibble(
          component = value_col,
          metadata = meta_col,
          method = "One-way ANOVA",
          estimate = eta_sq,
          p_value = aov_tab[["Pr(>F)"]][1]
        )
      }
    })
  })
  if (nrow(tests) == 0) return(empty_result)
  tests %>%
    mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
    arrange(p_adj, p_value)
}

format_top_association <- function(test_table) {
  if (nrow(test_table) == 0) {
    return("No suitable metadata association test could be computed.")
  }
  top <- test_table %>% dplyr::slice(1)
  paste0(
    top$component, " was most strongly associated with ",
    str_replace_all(top$metadata, "_", " "),
    " (", top$method, ", adjusted p = ", format_p(top$p_adj), ")."
  )
}

kruskal_association_tests <- function(scores, metadata, value_cols, id_col = "sampleID") {
  empty_result <- tibble(
    component = character(),
    metadata = character(),
    method = character(),
    statistic = numeric(),
    p_value = numeric(),
    p_adj = numeric()
  )
  if (!id_col %in% names(scores) || length(value_cols) == 0) return(empty_result)
  df <- scores %>%
    select(all_of(id_col), all_of(value_cols)) %>%
    dplyr::rename(sampleID = all_of(id_col)) %>%
    left_join(metadata, by = "sampleID")
  meta_cols <- setdiff(names(metadata), "sampleID")
  tests <- map_dfr(value_cols, function(value_col) {
    map_dfr(meta_cols, function(meta_col) {
      y <- df[[value_col]]
      x <- df[[meta_col]]
      if (is.numeric(x) && n_distinct(x, na.rm = TRUE) > 8) return(tibble())
      keep <- !is.na(y) & !is.na(x)
      y <- y[keep]
      x <- droplevels(as.factor(x[keep]))
      if (length(y) < 8 || nlevels(x) < 2 || nlevels(x) > 10 || any(table(x) < 2)) {
        return(tibble())
      }
      test <- tryCatch(kruskal.test(y ~ x), error = function(e) NULL)
      if (is.null(test)) return(tibble())
      tibble(
        component = value_col,
        metadata = meta_col,
        method = "Kruskal-Wallis",
        statistic = unname(test$statistic),
        p_value = test$p.value
      )
    })
  })
  if (nrow(tests) == 0) return(empty_result)
  tests %>%
    mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
    arrange(p_adj, p_value)
}

association_method_comparison <- function(scores, metadata, value_cols, id_col = "sampleID") {
  anova_tbl <- association_tests(scores, metadata, value_cols, id_col = id_col) %>%
    filter(method == "One-way ANOVA") %>%
    transmute(
      component,
      metadata,
      anova_eta_sq = estimate,
      anova_p_value = p_value,
      anova_p_adj = p_adj
    )
  kruskal_tbl <- kruskal_association_tests(scores, metadata, value_cols, id_col = id_col) %>%
    transmute(
      component,
      metadata,
      kruskal_statistic = statistic,
      kruskal_p_value = p_value,
      kruskal_p_adj = p_adj
    )
  full_join(anova_tbl, kruskal_tbl, by = c("component", "metadata")) %>%
    mutate(
      min_p_adj = pmin(anova_p_adj, kruskal_p_adj, na.rm = TRUE),
      min_p_adj = if_else(is.infinite(min_p_adj), NA_real_, min_p_adj),
      method_agreement = case_when(
        is.na(anova_p_adj) | is.na(kruskal_p_adj) ~ "Only one test available",
        anova_p_adj < 0.05 & kruskal_p_adj < 0.05 ~ "Both adjusted p < 0.05",
        anova_p_adj >= 0.05 & kruskal_p_adj >= 0.05 ~ "Both not adjusted-significant",
        TRUE ~ "Methods differ"
      )
    ) %>%
    arrange(min_p_adj, anova_p_adj, kruskal_p_adj)
}

format_top_method_comparison <- function(compare_table) {
  if (nrow(compare_table) == 0) {
    return("No categorical metadata comparison could be computed.")
  }
  top <- compare_table %>% dplyr::slice(1)
  paste0(
    top$component, " with ", str_replace_all(top$metadata, "_", " "),
    " had ANOVA adjusted p = ", format_p(top$anova_p_adj),
    " and Kruskal-Wallis adjusted p = ", format_p(top$kruskal_p_adj),
    "."
  )
}

format_comparison_for_table <- function(compare_table, n = 8) {
  compare_table %>%
    slice_head(n = n) %>%
    transmute(
      component,
      metadata = str_replace_all(metadata, "_", " "),
      anova_eta_sq = format_num(anova_eta_sq, 3),
      anova_p_adj = format_p(anova_p_adj),
      kruskal_statistic = format_num(kruskal_statistic, 2),
      kruskal_p_adj = format_p(kruskal_p_adj),
      method_agreement
    )
}

plot_association_heatmap <- function(test_table, title, max_log10 = 12) {
  if (nrow(test_table) == 0) {
    return(empty_plot(title, "No suitable metadata associations were available."))
  }
  numeric_components <- suppressWarnings(as.integer(str_remove(test_table$component, "^[^0-9]*")))
  component_levels <- if (all(is.na(numeric_components))) {
    unique(test_table$component)
  } else {
    test_table$component[order(numeric_components, na.last = TRUE)] %>% unique()
  }
  plot_df <- test_table %>%
    mutate(
      component = factor(component, levels = component_levels),
      metadata_label = str_replace_all(metadata, "_", " "),
      metadata_label = fct_reorder(metadata_label, -p_adj, .fun = min, .desc = FALSE),
      neg_log10_padj = pmin(-log10(p_adj), max_log10),
      label = case_when(
        is.na(p_adj) ~ "",
        p_adj < 0.001 ~ "***",
        p_adj < 0.01 ~ "**",
        p_adj < 0.05 ~ "*",
        TRUE ~ ""
      )
    )
  ggplot(plot_df, aes(x = component, y = metadata_label, fill = neg_log10_padj)) +
    geom_tile(color = "white", linewidth = 0.25) +
    geom_text(aes(label = label), size = 3, color = "black") +
    scale_fill_gradient(low = "grey95", high = "#B91C1C", na.value = "grey90") +
    labs(
      title = title,
      x = "Component",
      y = NULL,
      fill = "-log10 adjusted p"
    ) +
    theme_report(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_method_comparison_heatmap <- function(compare_table, title, max_log10 = 12) {
  if (nrow(compare_table) == 0) {
    return(empty_plot(title, "No categorical metadata comparisons were available."))
  }
  numeric_components <- suppressWarnings(as.integer(str_remove(compare_table$component, "^[^0-9]*")))
  component_levels <- if (all(is.na(numeric_components))) {
    unique(compare_table$component)
  } else {
    compare_table$component[order(numeric_components, na.last = TRUE)] %>% unique()
  }
  plot_df <- compare_table %>%
    mutate(metadata_label = str_replace_all(metadata, "_", " ")) %>%
    select(any_of("view"), component, metadata_label, anova_p_adj, kruskal_p_adj) %>%
    pivot_longer(
      c(anova_p_adj, kruskal_p_adj),
      names_to = "method",
      values_to = "p_adj"
    ) %>%
    mutate(
      method = recode(
        method,
        anova_p_adj = "One-way ANOVA",
        kruskal_p_adj = "Kruskal-Wallis"
      ),
      component = factor(component, levels = component_levels),
      neg_log10_padj = pmin(-log10(p_adj), max_log10),
      label = case_when(
        is.na(p_adj) ~ "",
        p_adj < 0.001 ~ "***",
        p_adj < 0.01 ~ "**",
        p_adj < 0.05 ~ "*",
        TRUE ~ ""
      )
    )
  p <- ggplot(plot_df, aes(x = component, y = metadata_label, fill = neg_log10_padj)) +
    geom_tile(color = "white", linewidth = 0.25) +
    geom_text(aes(label = label), size = 2.6) +
    scale_fill_gradient(low = "grey96", high = "#B91C1C", na.value = "grey90") +
    labs(title = title, x = "Component", y = NULL, fill = "-log10 adjusted p") +
    theme_report(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  if ("view" %in% names(plot_df)) {
    p + facet_grid(method ~ view, scales = "free_y", space = "free_y")
  } else {
    p + facet_wrap(~method, ncol = 1)
  }
}

extract_association <- function(test_table, component_name, metadata_name) {
  test_table %>%
    filter(component == component_name, metadata == metadata_name) %>%
    slice_head(n = 1)
}

run_umap_from_pca <- function(pca_result, max_pcs = 10, n_neighbors = 15, min_dist = 0.2) {
  if (!isTRUE(pca_result$ok)) {
    return(list(ok = FALSE, reason = pca_result$reason, embedding = tibble()))
  }
  pc_cols <- choose_pc_cols(pca_result, max_pcs = max_pcs)
  if (length(pc_cols) < 2 || nrow(pca_result$scores) < 5) {
    return(list(ok = FALSE, reason = "Too few samples or PCs for UMAP.", embedding = tibble()))
  }
  umap_input <- pca_result$scores %>%
    select(all_of(pc_cols)) %>%
    as.matrix()
  set.seed(seed_value)
  emb <- uwot::umap(
    umap_input,
    n_neighbors = min(n_neighbors, nrow(umap_input) - 1),
    min_dist = min_dist,
    metric = "euclidean",
    n_components = 2,
    verbose = FALSE,
    ret_model = FALSE
  )
  colnames(emb) <- c("UMAP1", "UMAP2")
  list(
    ok = TRUE,
    pc_cols = pc_cols,
    n_neighbors = min(n_neighbors, nrow(umap_input) - 1),
    embedding = bind_cols(pca_result$scores["sampleID"], as_tibble(emb))
  )
}

make_umap_plot <- function(umap_result, metadata, title, color_var = NULL, shape_var = NULL) {
  if (!isTRUE(umap_result$ok)) {
    return(empty_plot(title, umap_result$reason))
  }
  plot_df <- umap_result$embedding %>% left_join(metadata, by = "sampleID")
  p <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2))
  if (!is.null(color_var) && color_var %in% names(plot_df) &&
      !is.null(shape_var) && shape_var %in% names(plot_df)) {
    p <- p + geom_point(aes(color = .data[[color_var]], shape = .data[[shape_var]]), size = 2.5, alpha = 0.9)
  } else if (!is.null(color_var) && color_var %in% names(plot_df)) {
    p <- p + geom_point(aes(color = .data[[color_var]]), size = 2.5, alpha = 0.9)
  } else {
    p <- p + geom_point(color = "#2B6CB0", size = 2.5, alpha = 0.9)
  }
  p +
    labs(
      title = title,
      x = "UMAP1",
      y = "UMAP2",
      color = if (!is.null(color_var)) str_replace_all(color_var, "_", " ") else NULL,
      shape = if (!is.null(shape_var)) str_replace_all(shape_var, "_", " ") else NULL
    ) +
    theme_report()
}

umap_pc1_alignment <- function(umap_result, pca_result) {
  if (!isTRUE(umap_result$ok) || !isTRUE(pca_result$ok)) return(NA_real_)
  plot_df <- umap_result$embedding %>%
    left_join(pca_result$scores %>% select(sampleID, PC1), by = "sampleID")
  suppressWarnings(cor(plot_df$UMAP1, plot_df$PC1, use = "pairwise.complete.obs"))
}

build_mnn_graph <- function(pca_result, k = 10, max_pcs = 10) {
  empty <- list(ok = FALSE, reason = "MNN graph unavailable.", edges = tibble(), clusters = tibble())
  if (!isTRUE(pca_result$ok)) return(empty)
  pc_cols <- choose_pc_cols(pca_result, max_pcs = max_pcs)
  if (length(pc_cols) < 2 || nrow(pca_result$scores) < 4) return(empty)
  score_df <- pca_result$scores %>% select(sampleID, all_of(pc_cols))
  mat <- score_df %>% select(all_of(pc_cols)) %>% as.matrix()
  sample_ids <- score_df$sampleID
  k_use <- min(k, nrow(mat) - 1)
  knn <- FNN::get.knn(mat, k = k_use)
  directed_edges <- tibble(
    from = rep(sample_ids, each = k_use),
    to = sample_ids[as.vector(t(knn$nn.index))],
    rank = rep(seq_len(k_use), times = length(sample_ids)),
    distance = as.vector(t(knn$nn.dist))
  )
  reverse_edges <- tibble(from = directed_edges[["to"]], to = directed_edges[["from"]])
  mnn_edges <- directed_edges %>%
    semi_join(reverse_edges, by = c("from", "to")) %>%
    mutate(sample_a = pmin(from, to), sample_b = pmax(from, to)) %>%
    group_by(sample_a, sample_b) %>%
    summarise(distance = mean(distance), .groups = "drop") %>%
    transmute(from = sample_a, to = sample_b, distance)
  graph <- igraph::graph_from_data_frame(
    mnn_edges,
    directed = FALSE,
    vertices = tibble(name = sample_ids)
  )
  membership <- if (igraph::gsize(graph) > 0) {
    igraph::membership(igraph::cluster_louvain(graph))
  } else {
    igraph::components(graph)$membership
  }
  clusters <- tibble(
    sampleID = names(membership),
    drug_mnn_cluster = paste0("MNN", as.integer(membership)),
    mnn_degree = as.integer(igraph::degree(graph))
  )
  list(ok = TRUE, k = k_use, pc_cols = pc_cols, edges = mnn_edges, clusters = clusters, graph = graph)
}

make_mnn_pca_plot <- function(pca_result, mnn_result, metadata, color_var = "drug_mnn_cluster") {
  if (!isTRUE(pca_result$ok) || !isTRUE(mnn_result$ok)) {
    return(empty_plot("Drug MNN graph on PCA", "MNN graph unavailable."))
  }
  plot_df <- pca_result$scores %>%
    left_join(mnn_result$clusters, by = "sampleID") %>%
    left_join(metadata, by = "sampleID")
  edge_df <- mnn_result$edges %>%
    left_join(plot_df %>% select(sampleID, PC1, PC2), by = c("from" = "sampleID")) %>%
    dplyr::rename(x = PC1, y = PC2) %>%
    left_join(plot_df %>% select(sampleID, PC1, PC2), by = c("to" = "sampleID")) %>%
    dplyr::rename(xend = PC1, yend = PC2)
  ggplot(plot_df, aes(x = PC1, y = PC2)) +
    geom_segment(
      data = edge_df,
      aes(x = x, y = y, xend = xend, yend = yend),
      inherit.aes = FALSE,
      color = "grey55",
      alpha = 0.18,
      linewidth = 0.25
    ) +
    geom_point(aes(color = .data[[color_var]]), size = 2.4, alpha = 0.92) +
    labs(
      title = "Drug-response mutual nearest-neighbour graph",
      x = paste0("PC1 (", format_num(pca_result$variance$variance_percent[1], 1), "%)"),
      y = paste0("PC2 (", format_num(pca_result$variance$variance_percent[2], 1), "%)"),
      color = str_replace_all(color_var, "_", " ")
    ) +
    theme_report()
}
```

# Overview

This report analyses a multi-omics acute myeloid leukaemia (AML) dataset consisting of RNA-seq counts, targeted gene mutation calls, ex vivo drug-response AUC values, and clinical/sample metadata. The analysis first inspects and preprocesses each omics layer independently, then uses PCA to evaluate single-omics structure, and finally integrates the available samples with MOFA2 while leaving missing views as missing values. The objective is not only to produce plots, but to interpret whether latent molecular axes appear to reflect expression variation, mutation burden or genotype, drug sensitivity/resistance, or a combination of these signals.

The report is written to be reproducible from the four tab-separated files in the working directory. RNA-seq counts are normalised with DESeq2 size factors and variance-stabilised before PCA and MOFA. Mutation calls are encoded as binary sample-by-gene events and modelled with a Bernoulli likelihood in MOFA. Drug AUC values are treated as continuous response phenotypes. Missing values are inspected, and imputed only for PCA, and retained as missing values for MOFA where possible.

# Task 1: Load and Inspect the Data

```{r data-loading}
required_files <- c("rna_count.tsv", "gene_mutation.tsv", "drug_auc.tsv", "metadata.tsv")
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required file(s): ", paste(missing_files, collapse = ", "))
}

rna_raw <- readr::read_tsv("rna_count.tsv", show_col_types = FALSE, progress = FALSE)
mutation_raw <- readr::read_tsv("gene_mutation.tsv", show_col_types = FALSE, progress = FALSE) %>%
  janitor::clean_names()
drug_raw <- readr::read_tsv("drug_auc.tsv", show_col_types = FALSE, progress = FALSE) %>%
  janitor::clean_names()
metadata_raw <- readr::read_tsv("metadata.tsv", show_col_types = FALSE, progress = FALSE)

metadata <- metadata_raw %>%
  janitor::clean_names() %>%
  dplyr::rename(sampleID = sample_id) %>%
  mutate(
    sampleID = as.character(sampleID),
    across(where(is.character), ~na_if(.x, "NA")),
    center_id = if ("center_id" %in% names(.)) as.factor(center_id) else center_id
  ) %>%
  distinct(sampleID, .keep_all = TRUE)

rna_samples <- colnames(rna_raw)[-1]
mutation_samples <- sort(unique(as.character(mutation_raw$sample_id)))
drug_samples <- sort(unique(as.character(drug_raw$sample_id)))
metadata_samples <- sort(unique(metadata$sampleID))

raw_dimension_table <- tibble(
  dataset = c("RNA-seq count table", "Gene mutation table", "Drug AUC table", "Metadata table"),
  raw_rows = c(nrow(rna_raw), nrow(mutation_raw), nrow(drug_raw), nrow(metadata)),
  raw_columns = c(ncol(rna_raw), ncol(mutation_raw), ncol(drug_raw), ncol(metadata)),
  sample_count = c(length(rna_samples), length(mutation_samples), length(drug_samples), length(metadata_samples)),
  feature_count = c(
    nrow(rna_raw),
    n_distinct(mutation_raw$gene),
    n_distinct(drug_raw$inhibitor),
    ncol(metadata) - 1
  )
)

sample_sets <- list(
  RNAseq = rna_samples,
  Mutations = mutation_samples,
  Drug_AUC = drug_samples,
  Metadata = metadata_samples
)

overlap_table <- tibble(
  comparison = c(
    "RNA-seq and mutations",
    "RNA-seq and drug AUC",
    "Mutations and drug AUC",
    "All three omics layers",
    "All three omics layers plus metadata"
  ),
  shared_samples = c(
    length(intersect(rna_samples, mutation_samples)),
    length(intersect(rna_samples, drug_samples)),
    length(intersect(mutation_samples, drug_samples)),
    length(Reduce(intersect, sample_sets[c("RNAseq", "Mutations", "Drug_AUC")])),
    length(Reduce(intersect, sample_sets))
  )
)

all_samples <- sort(unique(unlist(sample_sets)))
sample_presence <- tibble(sampleID = all_samples) %>%
  mutate(
    RNAseq = sampleID %in% rna_samples,
    Mutations = sampleID %in% mutation_samples,
    Drug_AUC = sampleID %in% drug_samples,
    Metadata = sampleID %in% metadata_samples,
    omics_layers_present = RNAseq + Mutations + Drug_AUC
  )

sample_presence_long <- sample_presence %>%
  arrange(desc(omics_layers_present), sampleID) %>%
  mutate(sampleID = factor(sampleID, levels = sampleID)) %>%
  pivot_longer(c(RNAseq, Mutations, Drug_AUC, Metadata), names_to = "layer", values_to = "present")

id_consistency_table <- tibble(
  check = c(
    "RNA-seq samples missing from metadata",
    "Mutation samples missing from metadata",
    "Drug AUC samples missing from metadata",
    "Metadata samples without RNA-seq",
    "Metadata samples without mutation calls",
    "Metadata samples without drug AUC"
  ),
  n_samples = c(
    length(setdiff(rna_samples, metadata_samples)),
    length(setdiff(mutation_samples, metadata_samples)),
    length(setdiff(drug_samples, metadata_samples)),
    length(setdiff(metadata_samples, rna_samples)),
    length(setdiff(metadata_samples, mutation_samples)),
    length(setdiff(metadata_samples, drug_samples))
  )
)
```

```{r raw-dimensions-table}
knitr::kable(
  raw_dimension_table,
  caption = "Raw dimensions and feature/sample counts for the four input files."
)
```

The dimension table confirms that the four files have different structures. RNA-seq is already a feature-by-sample matrix with `r format_int(length(rna_samples))` samples, while the mutation and drug files are long tables that need to be reshaped into matrices before PCA and MOFA. The feature counts also show the scale difference between views: RNA-seq starts with `r format_int(nrow(rna_raw))` gene rows, whereas the mutation panel contains only `r format_int(n_distinct(mutation_raw$gene))` genes and the drug screen contains `r format_int(n_distinct(drug_raw$inhibitor))` inhibitors. This matters because the RNA view can capture broad continuous transcriptional variation, while the mutation view is much sparser and more targeted.

```{r overlap-table}
knitr::kable(
  overlap_table,
  caption = "Sample overlap across omics layers and metadata."
)
```

The overlap table shows that the limiting step for integration is not RNA-seq or mutation availability, but the overlap with drug-response measurements. The three omics layers share `r format_int(length(Reduce(intersect, sample_sets[c("RNAseq", "Mutations", "Drug_AUC")])))` samples. This is the practical sample size for MOFA, whereas single-omics plots can use larger assay-specific sets such as the `r format_int(length(rna_samples))` RNA samples or `r format_int(length(mutation_samples))` mutation samples.

```{r id-consistency-table}
knitr::kable(
  id_consistency_table,
  caption = "Consistency of sample identifiers between omics tables and metadata."
)
```

This identifier check separates true assay missingness from sample-name problems. The RNA, mutation, drug, and metadata tables all use the same AML-style sample IDs, but not every sample appears in every assay. For example, `r format_int(length(setdiff(metadata_samples, drug_samples)))` metadata samples have no drug AUC measurements, which explains why the integrated analysis has fewer samples than the metadata table.

```{r sample-overlap-plot, cache=FALSE, fig.width=8, fig.height=4, fig.cap="Sample presence/absence across RNA-seq, mutation, drug-response, and metadata tables. Each column is a sample, ordered by the number of available omics layers."}
ggplot(sample_presence_long, aes(x = sampleID, y = layer, fill = present)) +
  geom_tile(color = "white", linewidth = 0.15) +
  scale_fill_manual(
    values = c("FALSE" = "#E2E8F0", "TRUE" = "#2F855A"),
    breaks = c(FALSE, TRUE),
    labels = c("Absent", "Present")
  ) +
  labs(title = "Sample availability across data layers", x = "Samples", y = NULL, fill = "Status") +
  theme_report() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
```

The heatmap is a visual check of the same overlap pattern. Most columns have RNA, mutation, and metadata entries, but the drug-response row contains the main gaps. This means that missingness is not a subtle modelling detail: it determines how much information each sample contributes to the integrated model. The plot also reassures me that there is no obvious sample-ID formatting error, because absent measurements appear as structured assay gaps rather than as a completely empty layer.

## Task 1 interpretation: data dimensions and sample matching

The RNA-seq count file contains `r format_int(nrow(rna_raw))` gene rows and `r format_int(ncol(rna_raw))` columns, where the first column is the gene symbol and the remaining `r format_int(length(rna_samples))` columns are samples. The mutation file contains `r format_int(nrow(mutation_raw))` long-format sample-gene rows for `r format_int(n_distinct(mutation_raw$gene))` genes across `r format_int(length(mutation_samples))` samples. The drug-response file contains `r format_int(nrow(drug_raw))` long-format inhibitor/sample AUC measurements for `r format_int(n_distinct(drug_raw$inhibitor))` inhibitors across `r format_int(length(drug_samples))` samples. The metadata file contains `r format_int(nrow(metadata))` samples and `r format_int(ncol(metadata) - 1)` metadata variables.

The key integration constraint is sample overlap. There are `r format_int(length(Reduce(intersect, sample_sets[c("RNAseq", "Mutations", "Drug_AUC")])))` samples shared by RNA-seq, mutation, and drug-response layers, and the same `r format_int(length(Reduce(intersect, sample_sets)))` samples are also represented in the metadata table. A strict complete-case analysis would therefore be limited by drug availability. In the updated MOFA analysis, however, samples are aligned on the union of available IDs and missing whole views are left as `NA`, so MOFA can still learn from samples that have only one or two measured omics layers.

Additionally, the identifier check does not show evidence that sample identifiers follow incompatible naming conventions, but rather, the mismatch is caused by genuine assay availability differences. This is usual for multi-omics studies, where not every biospecimen yields every molecular or pharmacological measurement. The consequence is that single-omics PCA can use assay-specific sample sets, while MOFA uses all aligned samples and treats unmeasured views as missing rather than as biological zeroes.

## Metadata handling for plots and associations

```{r metadata-summary}
metadata_info <- metadata_summary(metadata)
plot_vars <- select_plot_variables(metadata)
metadata_color <- plot_vars$color
metadata_shape <- plot_vars$shape

knitr::kable(
  metadata_info,
  caption = "Metadata variables available for exploratory colouring and factor association testing."
)
```

The metadata include clinical and sample-level variables such as specimen type, vital status, disease stage, induction response, age, white blood cell count, FLT3-ITD status, and NPM1 status. For PCA and factor visualisation, the report automatically selects metadata columns that are neither all missing nor all unique. In this run, the primary colour variable is `r ifelse(is.null(metadata_color), "not available", str_replace_all(metadata_color, "_", " "))`, the shape variable is `r ifelse(is.null(metadata_shape), "not available", str_replace_all(metadata_shape, "_", " "))`. These choices are intentionally conservative as metadata with too many unique values are better tested statistically than used as a crowded plot legend.

Any metadata association reported below should be interpreted as exploratory. The analysis is unsupervised and does not adjust for all possible confounders simultaneously. In AML, disease stage, specimen type, mutation status, prior treatment, and sample handling can all influence expression and drug response, so apparent associations are best treated as hypotheses for follow-up rather than as definitive causal explanations.

## Statistical tests used for PC and factor associations

For the PC and MOFA-factor association screens I keep ANOVA as the main test and add Kruskal-Wallis beside it for categorical metadata. The response variable in these tests is always a continuous score, such as a PC score or a MOFA factor value. The metadata variable is then tested to see whether it explains part of that score.

One-way ANOVA is used when the metadata variable is categorical, such as centre, disease stage, FLT3-ITD status, or response category. It compares the mean score between groups. The F statistic is based on the ratio of between-group variation to within-group variation. If the group means are far apart compared with the spread inside each group, the p-value becomes small. In the tables, the ANOVA estimate for categorical variables is eta-squared. This is the fraction of score variance explained by the grouping variable in that one test.

The linear-model version of ANOVA is used when the metadata variable is continuous, such as age or white blood cell count. In this case the model is score ~ metadata. The test asks whether adding the metadata variable improves the model compared with an intercept-only model. This is still an ANOVA F test, but it is written as a simple linear regression because the predictor is continuous. The estimate shown in the table is the fitted slope. A positive slope means higher metadata values tend to have higher PC or factor scores, while a negative slope means the opposite.

Kruskal-Wallis is added as a non-parametric companion for categorical metadata. It first ranks all PC or factor scores and then tests whether the rank distributions differ between groups. This makes it less dependent on normally distributed residuals and less sensitive to outliers than one-way ANOVA. However, it tests differences in ranked distributions rather than differences in group means. It also does not give the same direct eta-squared effect size. Therefore, I use it as a sensitivity check beside ANOVA rather than as the only result.

Spearman correlation is not used as the main method here. Spearman works by ranking two ordered variables and testing whether they have a monotonic relationship. This is useful when both variables are continuous or ordinal and the expected relationship is consistently increasing or decreasing. It is not well matched to questions such as whether RNA PC2 differs by centre, whether a mutation PC differs by FLT3-ITD group, or whether a MOFA factor differs by response category. Those questions are group-comparison questions. ANOVA and Kruskal-Wallis answer them more directly. For continuous metadata, the linear-model ANOVA is also easier to compare with the rest of the report because it keeps all p-values in a variance-explained framework.

# Task 2: RNA-seq Preprocessing with DESeq2

```{r rna-preprocessing}
gene_column <- colnames(rna_raw)[1]
rna_counts_tbl <- rna_raw %>%
  dplyr::rename(gene = tidyselect::all_of(gene_column)) %>%
  mutate(gene = as.character(gene)) %>%
  filter(!is.na(gene), gene != "") %>%
  mutate(across(-gene, ~suppressWarnings(as.numeric(.x))))

rna_numeric_na_count <- sum(is.na(as.matrix(select(rna_counts_tbl, -gene))))
rna_duplicate_genes <- sum(duplicated(rna_counts_tbl$gene))

rna_counts_tbl <- rna_counts_tbl %>%
  mutate(across(-gene, ~replace_na(.x, 0))) %>%
  group_by(gene) %>%
  summarise(across(everything(), ~sum(.x, na.rm = TRUE)), .groups = "drop")

rna_counts <- rna_counts_tbl %>%
  column_to_rownames("gene") %>%
  as.matrix()
storage.mode(rna_counts) <- "integer"
rna_counts <- round(rna_counts)

rna_library_sizes <- tibble(
  sampleID = colnames(rna_counts),
  library_size = colSums(rna_counts)
) %>%
  left_join(metadata, by = "sampleID")

min_samples_expressed <- ceiling(0.05 * ncol(rna_counts))
rna_keep <- rowSums(rna_counts >= 10) >= min_samples_expressed
rna_counts_filtered <- rna_counts[rna_keep, , drop = FALSE]

dds <- DESeqDataSetFromMatrix(
  countData = rna_counts_filtered,
  colData = data.frame(sampleID = colnames(rna_counts_filtered), row.names = colnames(rna_counts_filtered)),
  design = ~ 1
)
dds <- estimateSizeFactors(dds)
vsd <- vst(dds, blind = TRUE)
rna_vst <- assay(vsd)

rna_gene_variance <- matrixStats::rowVars(rna_vst)
rna_hvg_n <- min(2000, length(rna_gene_variance))
rna_hvg_genes <- names(sort(rna_gene_variance, decreasing = TRUE))[seq_len(rna_hvg_n)]
rna_hvg <- rna_vst[rna_hvg_genes, , drop = FALSE]

rna_summary_table <- tibble(
  metric = c(
    "Genes in raw count table",
    "Duplicated gene symbols collapsed",
    "Numeric count values converted to missing before replacement",
    "Expression filter",
    "Genes retained after low-expression filtering",
    "Highly variable genes retained for PCA/MOFA"
  ),
  value = c(
    format_int(nrow(rna_raw)),
    format_int(rna_duplicate_genes),
    format_int(rna_numeric_na_count),
    paste0("count >= 10 in at least ", min_samples_expressed, " samples (5% of RNA samples)"),
    format_int(nrow(rna_counts_filtered)),
    format_int(nrow(rna_hvg))
  )
)

rna_mean_var <- tibble(
  gene = rownames(rna_vst),
  mean_vst = rowMeans(rna_vst),
  variance_vst = matrixStats::rowVars(rna_vst),
  selected_hvg = gene %in% rna_hvg_genes
)

set.seed(seed_value)
raw_count_vector <- as.numeric(rna_counts)
if (length(raw_count_vector) > 250000) {
  raw_count_vector <- sample(raw_count_vector, 250000)
}
raw_count_distribution <- tibble(log10_count_plus_1 = log10(raw_count_vector + 1))

rna_pca <- run_pca(rna_hvg, scale_features = FALSE)
rna_pca_plot <- make_pca_plot(rna_pca, metadata, "RNA-seq PCA", metadata_color, metadata_shape)
rna_scree_plot <- make_scree_plot(rna_pca, "RNA-seq PCA scree plot")

rna_pca_assoc <- if (isTRUE(rna_pca$ok)) {
  association_tests(rna_pca$scores, metadata, choose_pc_cols(rna_pca, 10))
} else {
  tibble()
}
rna_pca_assoc_compare <- if (isTRUE(rna_pca$ok)) {
  association_method_comparison(rna_pca$scores, metadata, choose_pc_cols(rna_pca, 10))
} else {
  tibble()
}
```

```{r rna-summary-table}
knitr::kable(
  rna_summary_table,
  caption = "RNA-seq preprocessing summary."
)
```

The RNA-seq preprocessing table shows that the raw matrix begins with `r format_int(nrow(rna_counts))` unique gene symbols after duplicate symbols are collapsed. The low-expression filter keeps `r format_int(nrow(rna_counts_filtered))` genes and removes `r format_int(nrow(rna_counts) - nrow(rna_counts_filtered))` genes that are not expressed at count \>= 10 in at least `r min_samples_expressed` samples. The final PCA and MOFA expression view then uses the top `r format_int(nrow(rna_hvg))` highly variable genes, so the multivariate analyses are driven by genes that vary across AML samples rather than by nearly constant or mostly zero genes.

```{r rna-library-size, fig.cap="Distribution of RNA-seq library sizes across AML samples. Large differences in total counts motivate library-size normalisation before multivariate analysis."}
ggplot(rna_library_sizes, aes(x = library_size)) +
  geom_histogram(bins = 35, fill = "#2B6CB0", color = "white") +
  scale_x_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  labs(title = "RNA-seq library size distribution", x = "Total counts per sample", y = "Number of samples") +
  theme_report()
```

The library-size histogram shows substantial sequencing-depth variation across samples. Total counts range from `r format_int(min(rna_library_sizes$library_size))` to `r format_int(max(rna_library_sizes$library_size))`, with a median of `r format_int(median(rna_library_sizes$library_size))`. This spread means that raw counts would partly measure library size rather than biology. DESeq2 size-factor normalisation is therefore required before PCA or MOFA, because otherwise samples with deeper sequencing could appear artificially separated.

```{r rna-count-density, fig.cap="Sampled distribution of raw RNA-seq counts on a log10(count + 1) scale. The spike near zero reflects many genes with very low or absent expression."}
ggplot(raw_count_distribution, aes(x = log10_count_plus_1)) +
  geom_histogram(bins = 60, fill = "#4C78A8", color = "white") +
  labs(title = "Raw count distribution", x = "log10(count + 1)", y = "Number of sampled gene-sample values") +
  theme_report()
```

The raw-count distribution is strongly right-skewed even after plotting on the log10(count + 1) scale. About `r format_pct(mean(rna_counts == 0), accuracy = 0.1)` of all gene-sample entries are zero, and many additional entries are very small. This explains why low-count filtering is important as a large part of the raw matrix contains weak signal that would add noise to distance-based methods. The long tail also shows why a variance-stabilising transformation is needed before treating expression values as approximately continuous.

```{r rna-mean-variance, fig.cap="Mean-variance relationship after DESeq2 variance-stabilising transformation. Highly variable genes selected for PCA and MOFA are highlighted."}
ggplot(rna_mean_var, aes(x = mean_vst, y = variance_vst, color = selected_hvg)) +
  geom_point(alpha = 0.35, size = 0.8) +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "#C2410C"), labels = c("Other genes", "Selected HVGs")) +
  labs(title = "RNA-seq variance after VST", x = "Mean VST expression", y = "Variance", color = NULL) +
  theme_report()
```

After VST, most genes have modest variance across samples: the median gene variance is `r format_num(median(rna_mean_var$variance_vst), 2)`. The selected HVGs have much higher variability, with a median variance of `r format_num(median(rna_mean_var$variance_vst[rna_mean_var$selected_hvg]), 2)` and a maximum of `r format_num(max(rna_mean_var$variance_vst[rna_mean_var$selected_hvg]), 2)`. The orange points are therefore genes that still differ meaningfully between AML samples after normalisation, not genes with the largest raw counts. These genes are the most informative features for identifying expression programmes.

```{r rna-pca-plot, fig.cap="PCA of variance-stabilised RNA-seq data using the 2,000 most variable genes."}
rna_pca_plot
```

The RNA PCA shows the strongest sample-level gradients in the VST-transformed HVG expression matrix. PC1 explains `r format_num(rna_pca$variance$variance_percent[1], 1)`% of the selected-gene variance and PC2 explains `r format_num(rna_pca$variance$variance_percent[2], 1)`%, so the first two axes together capture `r format_num(sum(rna_pca$variance$variance_percent[1:2]), 1)`%. This is enough to reveal broad expression structure, but it is not a complete summary of the transcriptome. Samples far from the centre may represent strong biological expression states or technical extremes worth checking against metadata.

```{r rna-scree-plot, fig.cap="Variance explained by the first ten RNA-seq principal components."}
rna_scree_plot
```

The RNA scree plot decreases gradually from PC1 to PC10. PC1 accounts for `r format_num(rna_pca$variance$variance_percent[1], 1)`% of the variance, while PCs 2-10 each explain between `r format_num(min(rna_pca$variance$variance_percent[2:10]), 1)`% and `r format_num(max(rna_pca$variance$variance_percent[2:10]), 1)`%. This pattern indicates that RNA-seq variation is spread across several expression programmes rather than being driven by one dominant axis. Indicating that no single PC cleanly separates the whole cohort as different PCs likely capture different mixtures of AML biology, cell composition, disease status, genotype, and technical variation.

```{r rna-pca-association-table}
knitr::kable(
  rna_pca_assoc %>% slice_head(n = 8) %>% mutate(across(c(p_value, p_adj), format_p)),
  caption = "Top exploratory ANOVA metadata associations with the first ten RNA-seq PCs."
)
```

The strongest RNA PCA metadata association is: `r format_top_association(rna_pca_assoc)` For categorical variables, the estimate column is eta-squared, meaning the fraction of PC-score variance explained by the group labels in that one-way ANOVA. For continuous variables, the estimate is the fitted slope from the linear model. The p_value is the raw probability under the null model of no association, while p_adj is the Benjamini-Hochberg adjusted p-value across all tested PC-metadata pairs. This table therefore shows which metadata variables best explain the leading RNA expression axes, rather than simply reporting that a test was run.

```{r rna-pca-anova-kruskal-table}
knitr::kable(
  format_comparison_for_table(rna_pca_assoc_compare, n = 8),
  caption = "Side-by-side one-way ANOVA and Kruskal-Wallis comparisons for categorical metadata associations with RNA-seq PCs."
)
```

The Kruskal-Wallis comparison is restricted to categorical metadata because it is a group-comparison test. For RNA-seq PCs, the strongest categorical comparison is: `r format_top_method_comparison(rna_pca_assoc_compare)` When ANOVA and Kruskal-Wallis both give small adjusted p-values, the result is more convincing because the group difference is visible both in the group means and in the rank distribution. If the two tests disagree, I would treat the association more cautiously because the result may depend on outliers, uneven group spread, or small group sizes.

```{r rna-pca-association-heatmap, fig.width=10, fig.height=6, fig.cap="ANOVA p-value heatmap for metadata associations with the first ten RNA-seq PCs. Darker red indicates stronger adjusted evidence that a metadata variable explains a PC score."}
plot_association_heatmap(rna_pca_assoc, "RNA-seq PC metadata association heatmap")
```

The RNA association heatmap is useful because it shows the full PC-by-metadata pattern instead of only the top row of a table. If one metadata variable is significant across several PCs, it suggests a broad effect. If it is restricted to one PC, it suggests a more specific axis. In this dataset the heatmap should be read alongside the scree plot as a low-numbered PC with a strong adjusted p-value is more influential than the same p-value on a PC explaining only a small fraction of expression variance.

```{r rna-pca-center-plot, fig.cap="RNA-seq PCA coloured by sequencing/processing centre. This directly checks whether a centre-associated PC looks like a possible batch axis."}
if ("center_id" %in% names(metadata)) {
  make_pca_plot(rna_pca, metadata, "RNA-seq PCA coloured by centre", "center_id", metadata_shape)
} else {
  empty_plot("RNA-seq PCA coloured by centre", "center_id is not available in the metadata.")
}
```

```{r rna-center-interpretation, results='asis'}
center_hit <- extract_association(rna_pca_assoc, "PC2", "center_id")
if (nrow(center_hit) > 0) {
  cat(
    "For the specific centre/batch question, RNA PC2 has an adjusted ANOVA p-value of **",
    format_p(center_hit$p_adj),
    "** for centre_id, with eta-squared **", format_num(center_hit$estimate, 3),
    "**. If the centre-coloured PCA shows one centre shifted mainly along PC2 while other centres overlap, this is consistent with a centre-associated batch component rather than a clean biological subgroup. If the colours overlap strongly, the p-value should be interpreted more cautiously because a statistical hit can also arise from a few influential samples or confounding between centre and disease composition.\n",
    sep = ""
  )
} else {
  cat("The ANOVA screen did not return a usable PC2-centre test, so there is no direct statistical evidence here that RNA PC2 is centre-driven.\n")
}
```

```{r rna-heatmap, fig.cap="Heatmap of the 50 most variable VST-transformed RNA-seq genes. Rows are scaled genes and columns are samples."}
if (nrow(rna_hvg) >= 2 && ncol(rna_hvg) >= 2) {
  top_heatmap_genes <- head(rna_hvg_genes, 50)
  pheatmap(
    rna_hvg[top_heatmap_genes, , drop = FALSE],
    scale = "row",
    show_colnames = FALSE,
    fontsize_row = 6,
    clustering_method = "ward.D2",
    main = "Top variable RNA-seq genes"
  )
}
```

The RNA heatmap shows the top 50 variable genes directly, with each row scaled to highlight high and low expression across samples. The presence of blocks of samples with coordinated high or low expression supports the PCA result. The cohort contains structured transcriptional variation as the scree plot shows several relevant PCs rather than one dominant component, these heatmap patterns should be interpreted as multiple expression programmes rather than a single binary split of the samples.

## RNA-seq interpretation

After filtering, `r format_int(nrow(rna_counts_filtered))` of `r format_int(nrow(rna_counts))` genes are retained. The filter removes genes with little evidence of expression, while keeping genes with counts of at least 10 in at least `r min_samples_expressed` RNA-seq samples. The top `r format_int(nrow(rna_hvg))` genes by VST variance are then used for PCA and as the expression input to MOFA. This selection keeps the analysis focused on expression features that vary across patients, rather than on genes that are nearly constant after normalisation.

The first two RNA-seq principal components explain `r if (isTRUE(rna_pca$ok)) paste0(format_num(rna_pca$variance$variance_percent[1], 1), "% and ", format_num(rna_pca$variance$variance_percent[2], 1), "%") else "not available"` of the RNA-seq variance among the selected highly variable genes. The most distant samples in the PC1/PC2 plane are `r top_pc_outliers(rna_pca)`. These samples are worth checking as potential biological extremes or technical outliers, although PCA distance alone is not sufficient evidence to exclude them.

The exploratory metadata tests indicate: `r format_top_association(rna_pca_assoc)` This does not prove that the metadata variable causes the RNA-seq separation, but it helps interpret whether a visible axis is plausibly linked to disease stage, sample source, genotype, or another clinical feature. If the PCA plot shows overlapping groups, the correct interpretation is that the dominant RNA expression gradients are either continuous, multi-factorial, or not well captured by the selected metadata annotations.

## Question 1: why not use raw RNA-seq counts directly for PCA or MOFA?

Raw RNA-seq counts are not directly suitable for PCA or Gaussian latent-factor models because they are discrete, highly heteroscedastic, and strongly affected by library size. A sample with more total reads can appear globally different even if its biological expression profile is similar to another sample. In addition, count variance increases with the mean, so highly expressed genes can dominate Euclidean distances and principal components simply because of their measurement scale.

DESeq2 size-factor normalisation addresses differences in sequencing depth, and the variance-stabilising transformation makes the expression values more comparable across the dynamic range. PCA and the Gaussian RNA view in MOFA then operate on values that better approximate continuous, homoscedastic measurements. This does not make the data perfect, but it makes the assumptions of the downstream methods much more defensible than using raw counts.

## Question 2: why filter lowly expressed genes?

Filtering lowly expressed genes reduces noise, multiple-testing burden in downstream analyses, and unnecessary computational complexity. Genes that are zero or near-zero in almost all samples contribute little reliable biological information to PCA or MOFA, but they can add unstable variance from sampling noise, mapping uncertainty, or occasional outlier counts.

In this dataset, filtering removes `r format_int(nrow(rna_counts) - nrow(rna_counts_filtered))` genes before transformation and feature selection. This is a technical choice with biological consequences: very rare cell-state markers could be removed if they are expressed in only a small subset of samples. The chosen threshold is therefore moderate rather than aggressive, retaining genes expressed in at least 5% of RNA-seq samples.

# Task 3: Mutation Preprocessing

```{r mutation-preprocessing}
mutation_matrix <- mutation_raw %>%
  transmute(
    sampleID = as.character(sample_id),
    gene = as.character(gene),
    value = as.integer(as.numeric(value) > 0)
  ) %>%
  filter(!is.na(sampleID), !is.na(gene), gene != "") %>%
  group_by(gene, sampleID) %>%
  summarise(value = max(value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = sampleID, values_from = value, values_fill = 0) %>%
  column_to_rownames("gene") %>%
  as.matrix()
storage.mode(mutation_matrix) <- "numeric"
mutation_matrix[is.na(mutation_matrix)] <- 0
mutation_matrix[mutation_matrix != 0] <- 1

mutation_frequencies <- tibble(
  gene = rownames(mutation_matrix),
  mutated_samples = rowSums(mutation_matrix),
  mutation_frequency = mutated_samples / ncol(mutation_matrix)
) %>%
  arrange(desc(mutated_samples), gene)

mutation_filter_min_samples <- 2
mutation_matrix_filtered <- mutation_matrix[rowSums(mutation_matrix) > mutation_filter_min_samples, , drop = FALSE]
mutation_burden <- tibble(
  sampleID = colnames(mutation_matrix),
  mutation_burden = colSums(mutation_matrix)
) %>%
  left_join(metadata, by = "sampleID")

mutation_summary_table <- tibble(
  metric = c(
    "Genes in mutation panel",
    "Samples with mutation information",
    "Assumption for absent sample-gene entries",
    "Filtering rule",
    "Genes retained after filtering",
    "Overall binary mutation density"
  ),
  value = c(
    format_int(nrow(mutation_matrix)),
    format_int(ncol(mutation_matrix)),
    "Absent sample-gene combinations encoded as 0",
    paste0("gene mutated in > ", mutation_filter_min_samples, " samples"),
    format_int(nrow(mutation_matrix_filtered)),
    format_pct(mean(mutation_matrix == 1), accuracy = 0.1)
  )
)

mutation_pca <- run_pca(mutation_matrix_filtered, scale_features = TRUE)
mutation_pca_plot <- make_pca_plot(mutation_pca, metadata, "Mutation PCA", metadata_color, metadata_shape)
mutation_scree_plot <- make_scree_plot(mutation_pca, "Mutation PCA scree plot")
mutation_pca_assoc <- if (isTRUE(mutation_pca$ok)) {
  association_tests(mutation_pca$scores, metadata, choose_pc_cols(mutation_pca, 10))
} else {
  tibble()
}
mutation_pca_assoc_compare <- if (isTRUE(mutation_pca$ok)) {
  association_method_comparison(mutation_pca$scores, metadata, choose_pc_cols(mutation_pca, 10))
} else {
  tibble()
}
```

```{r mutation-summary-table}
knitr::kable(
  mutation_summary_table,
  caption = "Mutation preprocessing summary."
)
```

The mutation preprocessing table shows that the targeted panel contains `r format_int(nrow(mutation_matrix))` genes across `r format_int(ncol(mutation_matrix))` samples. After converting the long table to a binary gene-by-sample matrix, `r format_int(nrow(mutation_matrix_filtered))` genes remain after removing genes mutated in `r mutation_filter_min_samples` or fewer samples. The overall mutation density is only `r format_pct(mean(mutation_matrix == 1), accuracy = 0.1)`, meaning that most gene-sample combinations are unmutated. This sparsity is the main reason the mutation layer needs binary modelling rather than expression-like Gaussian modelling.

```{r mutation-frequency-table}
knitr::kable(
  mutation_frequencies %>%
    mutate(mutation_frequency = percent(mutation_frequency, accuracy = 0.1)) %>%
    slice_head(n = 15),
  caption = "Mutation frequencies for genes in the targeted panel."
)
```

The frequency table shows that the mutation signal is concentrated in a small number of recurrent AML genes. `r mutation_frequencies$gene[1]` is mutated in `r format_int(mutation_frequencies$mutated_samples[1])` samples (`r format_pct(mutation_frequencies$mutation_frequency[1], accuracy = 0.1)`), followed by `r mutation_frequencies$gene[2]` in `r format_int(mutation_frequencies$mutated_samples[2])` samples (`r format_pct(mutation_frequencies$mutation_frequency[2], accuracy = 0.1)`) and `r mutation_frequencies$gene[3]` in `r format_int(mutation_frequencies$mutated_samples[3])` samples (`r format_pct(mutation_frequencies$mutation_frequency[3], accuracy = 0.1)`). This indicates that mutation PCs and mutation-driven MOFA factors are likely to be strongly influenced by NPM1/FLT3-related genotype structure rather than by many equally common events.

```{r mutation-top-genes-plot, fig.cap="Top mutated genes in the targeted AML mutation panel."}
mutation_frequencies %>%
  slice_head(n = 15) %>%
  mutate(gene = fct_reorder(gene, mutated_samples)) %>%
  ggplot(aes(x = mutated_samples, y = gene)) +
  geom_col(fill = "#7C3AED") +
  labs(title = "Most frequently mutated genes", x = "Number of mutated samples", y = NULL) +
  theme_report()
```

The bar plot makes the imbalance in mutation frequencies visually clear. The top two events, `r mutation_frequencies$gene[1]` and `r mutation_frequencies$gene[2]`, are much more common than most other panel genes, while several genes are present in only a small minority of samples. This means that any broad mutation axis is likely to contrast samples with these recurrent events against samples without them, rather than representing a smooth burden gradient across all genes.

```{r mutation-burden-plot, fig.cap="Sample-level mutation burden across the targeted gene panel."}
ggplot(mutation_burden, aes(x = mutation_burden)) +
  geom_histogram(binwidth = 1, fill = "#7C3AED", color = "white", boundary = -0.5) +
  labs(title = "Mutation burden per sample", x = "Number of mutated panel genes", y = "Number of samples") +
  theme_report()
```

The mutation-burden histogram shows that most samples carry few events in this targeted panel. The median burden is `r format_num(median(mutation_burden$mutation_burden), 0)` mutation, the mean is `r format_num(mean(mutation_burden$mutation_burden), 2)`, and the maximum is `r format_num(max(mutation_burden$mutation_burden), 0)` mutations. There are `r format_int(sum(mutation_burden$mutation_burden == 0))` samples with no mutation among the encoded panel features. This low burden explains why mutation PCA tends to be driven by specific recurrent events and co-mutation patterns rather than by a continuous genome-wide mutational load.

```{r mutation-pca-plot, fig.cap="Exploratory PCA of the filtered binary mutation matrix. PCA is used here as a visual summary, not as a generative model for binary data."}
mutation_pca_plot
```

The mutation PCA summarises the filtered binary mutation matrix as a continuous two-dimensional plot. PC1 explains `r format_num(mutation_pca$variance$variance_percent[1], 1)`% and PC2 explains `r format_num(mutation_pca$variance$variance_percent[2], 1)`% of the filtered mutation variance. These percentages are lower than a clean single-genotype split would produce, which fits the biology of AML as samples are defined by partly overlapping recurrent events such as NPM1, FLT3-ITD, DNMT3A and NRAS rather than by one mutation pattern shared by all patients.

```{r mutation-scree-plot, fig.cap="Variance explained by the first ten mutation principal components."}
mutation_scree_plot
```

The mutation scree plot is fairly flat compared with the drug-response scree plot. PC1 explains `r format_num(mutation_pca$variance$variance_percent[1], 1)`%, and PCs 2-10 each still explain between `r format_num(min(mutation_pca$variance$variance_percent[2:10]), 1)`% and `r format_num(max(mutation_pca$variance$variance_percent[2:10]), 1)`%. As only `r format_int(nrow(mutation_matrix_filtered))` mutation features are retained, this means that several different recurrent mutation events contribute separate axes of variation. It is not appropriate to interpret PC1 as "the" mutation phenotype for the cohort.

```{r mutation-pca-association-table}
knitr::kable(
  mutation_pca_assoc %>% slice_head(n = 8) %>% mutate(across(c(p_value, p_adj), format_p)),
  caption = "Top exploratory ANOVA metadata associations with the first ten mutation PCs."
)
```

The mutation PCA association table confirms that the mutation PCs are linked to known genotype annotations. The top hit is: `r format_top_association(mutation_pca_assoc)` For binary or categorical genotype labels such as NPM1 and FLT3-ITD, the ANOVA p-value tests whether samples in the different mutation-status groups have different PC scores. A small adjusted p-value therefore supports the interpretation that the PC is separating recurrent AML genotype patterns rather than only reflecting random projection of sparse 0/1 data.

```{r mutation-pca-anova-kruskal-table}
knitr::kable(
  format_comparison_for_table(mutation_pca_assoc_compare, n = 8),
  caption = "Side-by-side one-way ANOVA and Kruskal-Wallis comparisons for categorical metadata associations with mutation PCs."
)
```

The side-by-side mutation table checks whether the categorical genotype associations are still visible when scores are ranked rather than compared by group means. The strongest categorical comparison is: `r format_top_method_comparison(mutation_pca_assoc_compare)` Agreement between the two adjusted p-values supports a stable genotype interpretation. A disagreement would suggest that the mutation PC separation is influenced by a few extreme PC scores, unequal group variance, or a distributional shift that is not simply a difference in group means.

```{r mutation-pca-association-heatmap, fig.width=10, fig.height=6, fig.cap="ANOVA p-value heatmap for metadata associations with the first ten mutation PCs."}
plot_association_heatmap(mutation_pca_assoc, "Mutation PC metadata association heatmap")
```

The mutation p-value heatmap helps identify whether genotype annotations are concentrated on one mutation PC or spread across several PCs. This is important because the mutation scree plot is relatively flat: more than one PC can carry meaningful event structure. If NPM1, FLT3-ITD, or related AML variables appear in the darkest cells, that supports a genotype-axis interpretation. However, if centre or specimen variables dominate instead, the mutation PCA would need to be treated more cautiously.

```{r mutation-heatmap, fig.cap="Heatmap of the most frequent binary mutation events. Purple indicates a mutation call and light grey indicates no mutation call."}
if (nrow(mutation_matrix_filtered) >= 2 && ncol(mutation_matrix_filtered) >= 2) {
  top_mutation_genes <- mutation_frequencies %>%
    filter(gene %in% rownames(mutation_matrix_filtered)) %>%
    slice_head(n = 15) %>%
    pull(gene)
  pheatmap(
    mutation_matrix_filtered[top_mutation_genes, , drop = FALSE],
    color = c("#F1F5F9", "#7C3AED"),
    breaks = c(-0.1, 0.5, 1.1),
    show_colnames = FALSE,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    main = "Frequent mutation events"
  )
}
```

The mutation heatmap keeps the binary nature of the data visible showing that most cells are absent calls, with purple blocks marking recurrent events. The rows for `r mutation_frequencies$gene[1]` and `r mutation_frequencies$gene[2]` contain the densest bands, matching the frequency table and the PCA associations. This plot shows that mutation structure is event-based and sparse. The PCA is only a compressed visual summary of these 0/1 patterns.

## Mutation interpretation

The targeted mutation layer contains `r format_int(nrow(mutation_matrix))` genes. The most frequently mutated genes are `r paste(head(mutation_frequencies$gene, 5), collapse = ", ")`, with NPM1 and FLT3-ITD typically expected to be prominent in AML cohorts. After filtering genes mutated in more than `r mutation_filter_min_samples` samples, `r format_int(nrow(mutation_matrix_filtered))` genes are retained for PCA and MOFA. The binary mutation density is only `r format_pct(mean(mutation_matrix == 1), accuracy = 0.1)`, confirming that this view is sparse relative to RNA expression and drug response.

The first two mutation PCs explain `r if (isTRUE(mutation_pca$ok)) paste0(format_num(mutation_pca$variance$variance_percent[1], 1), "% and ", format_num(mutation_pca$variance$variance_percent[2], 1), "%") else "not available"` of the filtered mutation variance. This is a modest proportion, which is expected for sparse binary mutation data where different genes define partially overlapping patient subgroups. The most distant samples in mutation PC space are `r top_pc_outliers(mutation_pca)`, likely reflecting unusual combinations or burdens of panel mutations rather than smooth continuous biology.

PCA on binary mutation data should be interpreted cautiously. It is useful as an exploratory visualisation of co-occurrence and burden, but it does not model the Bernoulli nature of mutation calls. The metadata association screen indicates: `r format_top_association(mutation_pca_assoc)` If visible separation is weak, that is not a failure of the data as mutation profiles in targeted AML panels are often sparse, modular, and better represented by discrete events than by a small number of continuous PCs.

## Question 3: why is a Bernoulli likelihood appropriate for mutation data in MOFA?

A Bernoulli likelihood is appropriate because each mutation feature is encoded as a binary event for each sample, mutated or not mutated. The observations are therefore not continuous measurements with approximately Gaussian noise as they are 0/1 outcomes. Using a Bernoulli likelihood tells MOFA that the mutation view should be modelled as binary probability data, so the latent factors explain changes in mutation probability rather than changes in a continuous expression-like value.

This is especially important for sparse mutation panels. Treating mutation calls as Gaussian values would imply that intermediate values and symmetric continuous residuals are meaningful, which is not biologically correct for the encoded data. The Bernoulli likelihood is still an approximation because it does not fully model mutation calling uncertainty or clonality, but it is better matched to the data type than a Gaussian likelihood.

# Task 4: Drug AUC Preprocessing

```{r drug-preprocessing}
drug_matrix <- drug_raw %>%
  transmute(
    inhibitor = as.character(inhibitor),
    sampleID = as.character(sample_id),
    auc = as.numeric(auc)
  ) %>%
  filter(!is.na(inhibitor), inhibitor != "", !is.na(sampleID)) %>%
  group_by(inhibitor, sampleID) %>%
  summarise(auc = mean(auc, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = sampleID, values_from = auc) %>%
  column_to_rownames("inhibitor") %>%
  as.matrix()
storage.mode(drug_matrix) <- "numeric"

drug_missing_by_drug <- tibble(
  inhibitor = rownames(drug_matrix),
  missing_fraction = rowMeans(is.na(drug_matrix)),
  observed_samples = rowSums(!is.na(drug_matrix))
) %>%
  arrange(desc(missing_fraction), inhibitor)

drug_missing_by_sample <- tibble(
  sampleID = colnames(drug_matrix),
  missing_fraction = colMeans(is.na(drug_matrix)),
  observed_drugs = colSums(!is.na(drug_matrix))
) %>%
  left_join(metadata, by = "sampleID") %>%
  arrange(desc(missing_fraction))

drug_missing_threshold <- 0.20
drug_matrix_filtered <- drug_matrix[drug_missing_by_drug$missing_fraction[match(rownames(drug_matrix), drug_missing_by_drug$inhibitor)] <= drug_missing_threshold, , drop = FALSE]
drug_matrix_imputed <- impute_row_mean(drug_matrix_filtered)
drug_matrix_scaled_for_pca <- row_center_scale(drug_matrix_imputed)

drug_summary_table <- tibble(
  metric = c(
    "Inhibitors in raw drug-response table",
    "Samples with drug-response data",
    "Overall missing AUC fraction",
    "Drug filtering rule",
    "Drugs retained after missingness filtering",
    "AUC range",
    "Median AUC"
  ),
  value = c(
    format_int(nrow(drug_matrix)),
    format_int(ncol(drug_matrix)),
    format_pct(mean(is.na(drug_matrix)), accuracy = 0.1),
    paste0("remove drugs with > ", percent(drug_missing_threshold), " missing AUC values"),
    format_int(nrow(drug_matrix_filtered)),
    paste0(format_num(min(drug_matrix, na.rm = TRUE), 2), " to ", format_num(max(drug_matrix, na.rm = TRUE), 2)),
    format_num(median(drug_matrix, na.rm = TRUE), 2)
  )
)

drug_pca <- run_pca(drug_matrix_scaled_for_pca, scale_features = FALSE)
drug_pca_plot <- make_pca_plot(drug_pca, metadata, "Drug-response PCA", metadata_color, metadata_shape)
drug_scree_plot <- make_scree_plot(drug_pca, "Drug-response PCA scree plot")
drug_pca_assoc <- if (isTRUE(drug_pca$ok)) {
  association_tests(drug_pca$scores, metadata, choose_pc_cols(drug_pca, 10))
} else {
  tibble()
}
drug_pca_assoc_compare <- if (isTRUE(drug_pca$ok)) {
  association_method_comparison(drug_pca$scores, metadata, choose_pc_cols(drug_pca, 10))
} else {
  tibble()
}

drug_mnn <- build_mnn_graph(drug_pca, k = 10, max_pcs = 10)
drug_mnn_cluster_table <- if (isTRUE(drug_mnn$ok)) {
  drug_mnn$clusters %>%
    dplyr::count(drug_mnn_cluster, name = "n_samples") %>%
    arrange(desc(n_samples))
} else {
  tibble()
}

drug_auc_distribution <- as_tibble(as.data.frame(as.table(drug_matrix)), .name_repair = "minimal") %>%
  setNames(c("inhibitor", "sampleID", "auc")) %>%
  filter(!is.na(auc))
```

```{r drug-summary-table}
knitr::kable(
  drug_summary_table,
  caption = "Drug-response preprocessing summary."
)
```

The drug preprocessing table shows that the raw AUC matrix contains `r format_int(nrow(drug_matrix))` inhibitors across `r format_int(ncol(drug_matrix))` samples, with `r format_pct(mean(is.na(drug_matrix)), accuracy = 0.1)` missing values overall. Applying the \>20% missingness filter retains `r format_int(nrow(drug_matrix_filtered))` drugs and removes `r format_int(nrow(drug_matrix) - nrow(drug_matrix_filtered))`. The important point is that imputation (filling in missing values with a reasonable guess so that an analysis can run on a complete matrix) is used only for PCA and heatmap display. The MOFA drug view keeps remaining missing AUC values because MOFA can model incomplete observations directly.

```{r drug-missing-top-table}
knitr::kable(
  drug_missing_by_drug %>%
    mutate(missing_fraction = percent(missing_fraction, accuracy = 0.1)) %>%
    slice_head(n = 12),
  caption = "Drugs with the highest missingness."
)
```

The drug missingness table shows that missingness is concentrated in a subset of inhibitors. The most incomplete drug is `r drug_missing_by_drug$inhibitor[1]`, missing in `r format_pct(drug_missing_by_drug$missing_fraction[1], accuracy = 0.1)` of samples, followed by `r drug_missing_by_drug$inhibitor[2]` at `r format_pct(drug_missing_by_drug$missing_fraction[2], accuracy = 0.1)`. These drugs are excluded because their observed values would represent too small and potentially biased a subset of the cohort.

```{r drug-missing-sample-table}
knitr::kable(
  drug_missing_by_sample %>%
    transmute(
      sampleID,
      missing_fraction = percent(missing_fraction, accuracy = 0.1),
      observed_drugs
    ) %>%
    slice_head(n = 12),
  caption = "Samples with the highest drug-response missingness."
)
```

The sample-level missingness table shows that incomplete drug coverage is also uneven across patients. The most incomplete sample, `r drug_missing_by_sample$sampleID[1]`, has AUC values for only `r format_int(drug_missing_by_sample$observed_drugs[1])` drugs and is missing `r format_pct(drug_missing_by_sample$missing_fraction[1], accuracy = 0.1)` of the panel. Such samples can influence PCA because PCA requires imputation, so the drug PCA is interpreted cautiously. MOFA is less sensitive to individual missing entries because it can leave them as missing after drug filtering.

```{r drug-missingness-plot, fig.cap="Distribution of missing AUC fractions across drugs. Drugs above the dashed line are excluded from PCA and MOFA."}
ggplot(drug_missing_by_drug, aes(x = missing_fraction)) +
  geom_histogram(bins = 30, fill = "#0F766E", color = "white") +
  geom_vline(xintercept = drug_missing_threshold, linetype = "dashed", color = "#B91C1C", linewidth = 0.8) +
  scale_x_continuous(labels = percent_format()) +
  labs(title = "Drug-level missingness", x = "Missing fraction", y = "Number of drugs") +
  theme_report()
```

The missingness histogram shows a main group of drugs with moderate missingness and a long tail of poorly covered drugs. The median drug-level missingness is `r format_pct(median(drug_missing_by_drug$missing_fraction), accuracy = 0.1)`, but the maximum is `r format_pct(max(drug_missing_by_drug$missing_fraction), accuracy = 0.1)`, which explains why the mean missingness is higher. The dashed 20% threshold removes the high-missingness tail while keeping most drugs with reasonably broad sample coverage.

```{r drug-auc-density, fig.cap="Distribution of observed drug AUC values. Lower AUC corresponds to stronger drug sensitivity in a viability-response setting."}
ggplot(drug_auc_distribution, aes(x = auc)) +
  geom_density(fill = "#0F766E", alpha = 0.35, color = "#0F766E") +
  labs(title = "Observed AUC distribution", x = "AUC", y = "Density") +
  theme_report()
```

The AUC distribution spans from `r format_num(min(drug_matrix, na.rm = TRUE), 2)` to `r format_num(max(drug_matrix, na.rm = TRUE), 2)`, with a median of `r format_num(median(drug_matrix, na.rm = TRUE), 2)`. Lower AUC means stronger sensitivity and higher AUC means resistance, the concentration of values around the mid-to-high part of the range suggests that many drug-sample combinations show partial or limited response rather than uniformly strong sensitivity. The distribution is still wide enough to support PCA and MOFA analysis of relative sensitivity and resistance patterns.

```{r drug-missing-bar, fig.cap="Top missing drug-response features before filtering."}
drug_missing_by_drug %>%
  slice_head(n = 15) %>%
  mutate(inhibitor = fct_reorder(inhibitor, missing_fraction)) %>%
  ggplot(aes(x = missing_fraction, y = inhibitor)) +
  geom_col(fill = "#0F766E") +
  scale_x_continuous(labels = percent_format()) +
  labs(title = "Drugs with highest missingness", x = "Missing fraction", y = NULL) +
  theme_report()
```

The missingness bar plot identifies the specific inhibitors removed or most affected by missing data. The top entries, including `r paste(head(drug_missing_by_drug$inhibitor, 3), collapse = ", ")`, have far less complete coverage than the retained majority. Naming these drugs is important because missingness may not be random, it can reflect assay availability, drug panel changes, or measurement failures for particular inhibitors.

```{r drug-pca-plot, fig.cap="PCA of row-mean-imputed and row-scaled drug AUC profiles after excluding highly missing drugs."}
drug_pca_plot
```

The drug PCA shows a stronger leading axis than the RNA or mutation PCA. PC1 explains `r format_num(drug_pca$variance$variance_percent[1], 1)`% of scaled AUC variation and PC2 explains `r format_num(drug_pca$variance$variance_percent[2], 1)`%, so the first two PCs together capture `r format_num(sum(drug_pca$variance$variance_percent[1:2]), 1)`%. This suggests that many drugs share a coordinated response pattern across samples, consistent with a broad sensitivity/resistance gradient. Row-mean imputation is used only so PCA can be computed as the imputed values are not treated as observed measurements.

```{r drug-scree-plot, fig.cap="Variance explained by the first ten drug-response principal components."}
drug_scree_plot
```

The drug scree plot drops sharply after PC1. PC1 explains `r format_num(drug_pca$variance$variance_percent[1], 1)`%, PC2 explains `r format_num(drug_pca$variance$variance_percent[2], 1)`%, and PCs 3-10 each explain only `r format_num(min(drug_pca$variance$variance_percent[3:10]), 1)`% to `r format_num(max(drug_pca$variance$variance_percent[3:10]), 1)`%. This pattern is consistent with one broad drug-response axis plus smaller secondary programmes, rather than many equally strong independent response patterns.

```{r drug-pca-association-table}
knitr::kable(
  drug_pca_assoc %>% slice_head(n = 8) %>% mutate(across(c(p_value, p_adj), format_p)),
  caption = "Top exploratory ANOVA metadata associations with the first ten drug-response PCs."
)
```

The drug PCA association table suggests which clinical or genotype variables align with the strongest AUC gradients. The top hit is: `r format_top_association(drug_pca_assoc)` In this ANOVA screen, categorical estimates are eta-squared values and continuous estimates are fitted slopes. The adjusted p-values indicate that these are not only the smallest nominal p-values from many tests, but associations that remain notable after multiple-testing correction.

```{r drug-pca-anova-kruskal-table}
knitr::kable(
  format_comparison_for_table(drug_pca_assoc_compare, n = 8),
  caption = "Side-by-side one-way ANOVA and Kruskal-Wallis comparisons for categorical metadata associations with drug-response PCs."
)
```

The Kruskal-Wallis comparison asks whether the categorical metadata groups have different ranked drug PC scores. The strongest categorical comparison is: `r format_top_method_comparison(drug_pca_assoc_compare)` If this agrees with the one-way ANOVA result, the drug-response association is less likely to be driven only by a few extreme AUC profiles. If the tests differ, the PCA association should be interpreted cautiously because the group pattern may be non-normal, unevenly spread, or influenced by a small number of samples.

```{r drug-pca-association-heatmap, fig.width=10, fig.height=6, fig.cap="ANOVA p-value heatmap for metadata associations with the first ten drug-response PCs."}
plot_association_heatmap(drug_pca_assoc, "Drug-response PC metadata association heatmap")
```

The drug p-value heatmap shows whether drug-response variation is dominated by one broad clinical/genotype axis or by several weaker axes. A dark cell on PC1 is especially important because drug PC1 explains `r format_num(drug_pca$variance$variance_percent[1], 1)`% of the scaled drug-response variance. If centre appears strongly here, the drug structure may partly reflect assay or processing differences. If the genotype or disease variables dominate, the PCA could be capturing biological drug sensitivity/resistance patterns.

## Drug KNN/MNN neighbourhood analysis

```{r drug-mnn-cluster-table}
if (isTRUE(drug_mnn$ok)) {
  knitr::kable(
    drug_mnn_cluster_table,
    caption = "Drug-response mutual nearest-neighbour cluster sizes."
  )
} else {
  cat("Drug-response MNN graph could not be computed.")
}
```

The KNN/MNN analysis uses the first `r if (isTRUE(drug_mnn$ok)) length(drug_mnn$pc_cols) else 0` drug PCs as a denoised drug-response space. For each sample, the KNN step finds its nearest drug-response neighbours, and the MNN step keeps only pairs that select each other. This is stricter than ordinary KNN because it reduces one-sided neighbour links caused by outlying samples or uneven density. The resulting clusters are drug-response neighbourhoods, not cell types. They indicate groups of AML samples with similar ex vivo AUC profiles.

```{r drug-mnn-pca-plot, fig.cap="Mutual nearest-neighbour graph overlaid on the drug-response PCA. Lines connect samples that are reciprocal nearest neighbours in drug PC space."}
make_mnn_pca_plot(drug_pca, drug_mnn, metadata, "drug_mnn_cluster")
```

```{r drug-mnn-interpretation, results='asis'}
if (isTRUE(drug_mnn$ok) && nrow(drug_mnn_cluster_table) > 0) {
  largest_cluster <- drug_mnn_cluster_table$drug_mnn_cluster[1]
  largest_n <- drug_mnn_cluster_table$n_samples[1]
  n_clusters <- nrow(drug_mnn_cluster_table)
  graph_degree <- drug_mnn$clusters$mnn_degree
  cat(
    "The MNN graph separates the drug-response profiles into **", n_clusters,
    "** reciprocal-neighbour clusters. The largest cluster is **", largest_cluster,
    "** with **", largest_n, "** samples. Samples with high MNN degree sit inside dense drug-response neighbourhoods, whereas low-degree samples lie at the edge of the PCA cloud or have more unusual AUC patterns. In the PCA overlay, a cluster stretched mainly along PC1 should be interpreted as part of the broad sensitivity/resistance gradient, while a compact island separated on both PC1 and PC2 would be stronger evidence for a distinct drug-response subgroup.\n",
    sep = ""
  )
} else {
  cat("The MNN graph did not contain enough reciprocal drug-response neighbours to support cluster-level interpretation.\n")
}
```

```{r drug-heatmap, fig.cap="Heatmap of the 40 most variable filtered drug-response profiles. Values are row-scaled AUC after row-mean imputation for visualisation only."}
if (nrow(drug_matrix_scaled_for_pca) >= 2 && ncol(drug_matrix_scaled_for_pca) >= 2) {
  drug_var <- matrixStats::rowVars(drug_matrix_scaled_for_pca)
  top_drugs <- names(sort(drug_var, decreasing = TRUE))[seq_len(min(40, length(drug_var)))]
  pheatmap(
    drug_matrix_scaled_for_pca[top_drugs, , drop = FALSE],
    show_colnames = FALSE,
    fontsize_row = 6,
    clustering_method = "ward.D2",
    main = "Variable drug-response profiles"
  )
}
```

The drug heatmap shows the most variable filtered drugs after row scaling. Since AUC scales differ between inhibitors, row scaling makes each drug contribute a relative sensitivity/resistance profile rather than letting drugs with larger raw ranges dominate the display. Blocks of samples with consistently high or low scaled AUC across groups of drugs support the PCA result that drug response contains coordinated patterns, not just isolated drug-specific noise.

## Drug-response interpretation

The drug-response matrix contains `r format_int(nrow(drug_matrix))` inhibitors and `r format_int(ncol(drug_matrix))` samples before filtering. Overall missingness is `r format_pct(mean(is.na(drug_matrix)), accuracy = 0.1)`, but missingness is not evenly distributed across inhibitors. Filtering drugs with more than `r percent(drug_missing_threshold)` missing values retains `r format_int(nrow(drug_matrix_filtered))` inhibitors. For PCA, the remaining missing values are imputed by the drug-specific row mean, because PCA cannot operate with missing entries. For MOFA, the filtered AUC matrix is kept with missing values because MOFA can accommodate missing individual observations.

Drug PCA is performed on row-scaled AUC values. Scaling is important because drugs can have different absolute AUC ranges. Without scaling, high-variance drugs would dominate the first PCs even if that variance mainly reflects assay range rather than coordinated sensitivity/resistance structure. The first two drug PCs explain `r if (isTRUE(drug_pca$ok)) paste0(format_num(drug_pca$variance$variance_percent[1], 1), "% and ", format_num(drug_pca$variance$variance_percent[2], 1), "%") else "not available"` of the scaled drug-response variation. The most distant drug-response profiles in PC space are `r top_pc_outliers(drug_pca)`.

The exploratory metadata screen indicates: `r format_top_association(drug_pca_assoc)` If drug-response groups appear in the PCA plot, they may reflect broad sensitivity or resistance programmes, specific pathway dependencies, differences in disease stage, or sample quality effects. Because the drug layer has substantial missingness and ex vivo AUC measurements can be assay-sensitive, these clusters should be interpreted together with RNA-seq, mutation status, and MOFA factors rather than in isolation.

## Question 4: what does a high or low AUC mean?

In a viability-response assay, the AUC summarises the area under the cell viability curve across drug concentrations. A high AUC means that viability remains high across much of the dose range, so the sample is relatively resistant or insensitive to that drug. A low AUC means that viability drops more strongly across the dose range, so the sample is relatively sensitive.

In this dataset, observed AUC values range from `r format_num(min(drug_matrix, na.rm = TRUE), 2)` to `r format_num(max(drug_matrix, na.rm = TRUE), 2)`, with a median of `r format_num(median(drug_matrix, na.rm = TRUE), 2)`. The absolute units depend on how the original dose-response curves were integrated, so the safest interpretation is comparative within this dataset: lower AUC indicates greater sensitivity, and higher AUC indicates greater resistance.

# Task 5: Single-Omics PCA Comparison

```{r pca-rna, fig.width=5, fig.height=4, fig.cap="PCA of RNA-seq data (PC1 vs PC2)."}
rna_pca_plot
```

```{r scree-rna, fig.width=5, fig.height=4, fig.cap="Scree plot for RNA-seq PCA."}
rna_scree_plot
```

```{r pca-mutation, fig.width=5, fig.height=4, fig.cap="PCA of mutation data (PC1 vs PC2)."}
mutation_pca_plot
```

```{r scree-mutation, fig.width=5, fig.height=4, fig.cap="Scree plot for mutation PCA."}
mutation_scree_plot
```

```{r pca-drug, fig.width=5, fig.height=4, fig.cap="PCA of drug-response data (PC1 vs PC2)."}
drug_pca_plot
```

```{r scree-drug, fig.width=5, fig.height=4, fig.cap="Scree plot for drug-response PCA."}
drug_scree_plot
```

The three PCA panels show that each omics layer has a different variance structure. RNA PC1 explains `r format_num(rna_pca$variance$variance_percent[1], 1)`%, with a gradual scree plot, so expression variation is distributed across several transcriptional programmes. Mutation PC1 explains only `r format_num(mutation_pca$variance$variance_percent[1], 1)`%, and the mutation scree is relatively flat, which fits a sparse panel where different genes define different genotype axes. Drug PC1 explains `r format_num(drug_pca$variance$variance_percent[1], 1)`%, the strongest leading PC among the three views, suggesting a broad coordinated drug-response gradient. Comparing the panels shows why integration is useful: no single omics layer gives the full structure of the cohort.

```{r pca-association-combined, fig.width=12, fig.height=7, fig.cap="Combined ANOVA p-value heatmap for metadata associations with the first ten PCs from each omics layer."}
pca_assoc_all <- bind_rows(
  RNAseq = rna_pca_assoc,
  Mutations = mutation_pca_assoc,
  Drug_AUC = drug_pca_assoc,
  .id = "view"
)

if (nrow(pca_assoc_all) > 0) {
  pca_assoc_all %>%
    mutate(
      component = factor(component, levels = paste0("PC", 1:10)),
      metadata_label = str_replace_all(metadata, "_", " "),
      neg_log10_padj = pmin(-log10(p_adj), 12),
      label = case_when(
        p_adj < 0.001 ~ "***",
        p_adj < 0.01 ~ "**",
        p_adj < 0.05 ~ "*",
        TRUE ~ ""
      )
    ) %>%
    ggplot(aes(x = component, y = metadata_label, fill = neg_log10_padj)) +
    geom_tile(color = "white", linewidth = 0.2) +
    geom_text(aes(label = label), size = 2.6) +
    facet_wrap(~view, nrow = 1) +
    scale_fill_gradient(low = "grey96", high = "#B91C1C") +
    labs(title = "PC-metadata association patterns across omics views", x = "Principal component", y = NULL, fill = "-log10 adjusted p") +
    theme_report(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
} else {
  empty_plot("PC-metadata association patterns", "No association tests were available.")
}
```

The combined p-value heatmap is the clearest way to compare what different PCs are capturing. A vertical stripe in one view means a metadata variable is linked to several components, while a single dark tile means a specific PC is carrying that association. This is particularly helpful for the centre/batch question: if centre is mainly significant for RNA PC2 and the centre-coloured PCA shows a PC2 shift, PC2 is plausibly partly batch-related. If centre is weak or scattered across later PCs, it is less likely to be the main driver of the visible PCA pattern.

```{r pca-anova-kruskal-comparison-heatmap, fig.width=12, fig.height=8, fig.cap="Side-by-side one-way ANOVA and Kruskal-Wallis adjusted p-value heatmaps for categorical metadata associations with the first ten PCs."}
pca_assoc_compare_all <- bind_rows(
  RNAseq = rna_pca_assoc_compare,
  Mutations = mutation_pca_assoc_compare,
  Drug_AUC = drug_pca_assoc_compare,
  .id = "view"
)

plot_method_comparison_heatmap(
  pca_assoc_compare_all,
  "One-way ANOVA and Kruskal-Wallis PC association comparison"
)
```

This comparison heatmap shows whether categorical PC associations are consistent between the parametric and rank-based tests. The one-way ANOVA panels focus on group mean differences, while the Kruskal-Wallis panels focus on ranked group separation. When the same metadata variable appears dark in both methods for the same PC, the association is less dependent on the distributional assumptions of ANOVA. When only one method is dark, the result is still useful but should be interpreted as more sensitive to the test choice.

## UMAP visualisation of sample neighbourhoods

```{r umap-preparation}
rna_umap <- run_umap_from_pca(rna_pca, max_pcs = 10, n_neighbors = 15, min_dist = 0.25)
mutation_umap <- run_umap_from_pca(mutation_pca, max_pcs = 10, n_neighbors = 15, min_dist = 0.25)
drug_umap <- run_umap_from_pca(drug_pca, max_pcs = 10, n_neighbors = 15, min_dist = 0.25)

rna_umap_pc1_cor <- umap_pc1_alignment(rna_umap, rna_pca)
mutation_umap_pc1_cor <- umap_pc1_alignment(mutation_umap, mutation_pca)
drug_umap_pc1_cor <- umap_pc1_alignment(drug_umap, drug_pca)
```

```{r rna-umap, fig.width=5, fig.height=4, fig.cap="UMAP embedding from the first ten RNA-seq PCs."}
make_umap_plot(
  rna_umap,
  metadata,
  "RNA-seq UMAP",
  metadata_color,
  metadata_shape
)
```

```{r mutation-umap, fig.width=5, fig.height=4, fig.cap="UMAP embedding from the first ten mutation PCs."}
make_umap_plot(
  mutation_umap,
  metadata,
  "Mutation UMAP",
  metadata_color,
  metadata_shape
)
```

```{r drug-umap, fig.width=5, fig.height=4, fig.cap="UMAP embedding from the first ten drug-response PCs."}
make_umap_plot(
  drug_umap,
  metadata,
  "Drug-response UMAP",
  metadata_color,
  metadata_shape
)
```

PCA first reduces noise, a nearest-neighbour structure is built in the reduced space, and UMAP then places similar observations close together in two dimensions. The important difference is that these observations are bulk AML samples. Therefore, the UMAP can only show whether samples form visual neighbourhoods that may reflect genotype, disease state, sample source, or drug-response phenotype. The RNA UMAP has a UMAP1-PC1 correlation of `r format_num(rna_umap_pc1_cor, 2)`, the mutation UMAP has `r format_num(mutation_umap_pc1_cor, 2)`, and the drug UMAP has `r format_num(drug_umap_pc1_cor, 2)`. A high absolute value means the UMAP is mostly unfolding the first PC gradient, whereas a low value suggests that later PCs or local neighbourhoods also shape the embedding.

```{r drug-umap-mnn, fig.cap="Drug-response UMAP coloured by mutual nearest-neighbour clusters from the drug PCA space."}
if (isTRUE(drug_umap$ok) && isTRUE(drug_mnn$ok)) {
  drug_umap$embedding %>%
    left_join(drug_mnn$clusters, by = "sampleID") %>%
    ggplot(aes(x = UMAP1, y = UMAP2, color = drug_mnn_cluster)) +
    geom_point(size = 2.6, alpha = 0.9) +
    labs(title = "Drug-response UMAP by MNN cluster", x = "UMAP1", y = "UMAP2", color = "MNN cluster") +
    theme_report()
} else {
  empty_plot("Drug-response UMAP by MNN cluster", "Drug UMAP or MNN graph unavailable.")
}
```

The drug UMAP coloured by MNN cluster connects the UMAP visualisation to the KNN/MNN analysis. If the MNN clusters occupy compact UMAP regions, the drug-response neighbourhoods are stable across both methods. If the clusters are stretched along a continuous UMAP arc, the data are better interpreted as a gradient of relative sensitivity/resistance rather than as sharply separated subgroups. This distinction matters because UMAP can make gradual structure look cluster-like. The reciprocal-neighbour graph helps check whether apparent islands correspond to strong local similarity.

## Question 5: limitations of interpreting each omics layer separately

Interpreting each omics layer separately can miss biology that is distributed across layers. For example, a mutation may only be functionally important in samples with a particular expression state, or a drug-response phenotype may reflect both genotype and transcriptional adaptation. Single-omics PCA also makes it difficult to distinguish biological axes from hidden confounding, because technical variation in one assay may look like a real molecular subgroup unless checked against other data types.

The limitations differ by view. RNA-seq is rich but vulnerable to batch, cellular composition, and library-size effects. Mutation data are biologically interpretable but sparse and low-dimensional, so weak PCA separation does not mean mutations are unimportant. Drug AUC data are functional and clinically appealing, but missingness and assay variability can influence structure. Integration is useful because a shared factor model can identify latent axes supported by more than one view, while also allowing view-specific factors when a signal is present in only one data type.

# Task 6: MOFA Input Preparation

```{r mofa-data-preparation}
complete_case_mofa_samples <- Reduce(intersect, list(
  colnames(rna_hvg),
  colnames(mutation_matrix_filtered),
  colnames(drug_matrix_filtered)
))
complete_case_mofa_samples <- sort(complete_case_mofa_samples)

mofa_all_samples <- sort(unique(c(
  colnames(rna_hvg),
  colnames(mutation_matrix_filtered),
  colnames(drug_matrix_filtered)
)))

if (length(mofa_all_samples) < 10) {
  stop("Too few samples for MOFA integration: ", length(mofa_all_samples))
}

align_view_to_samples <- function(mat, all_samples) {
  mat <- as.matrix(mat)
  out <- matrix(
    NA_real_,
    nrow = nrow(mat),
    ncol = length(all_samples),
    dimnames = list(rownames(mat), all_samples)
  )
  present <- intersect(colnames(mat), all_samples)
  out[, present] <- mat[, present, drop = FALSE]
  out
}

rna_mofa <- align_view_to_samples(rna_hvg, mofa_all_samples)
mutation_mofa <- align_view_to_samples(mutation_matrix_filtered, mofa_all_samples)
drug_mofa <- align_view_to_samples(drug_matrix_filtered, mofa_all_samples)

rna_mofa_scaled <- row_center_scale(rna_mofa)
drug_mofa_scaled <- row_center_scale(drug_mofa)

mofa_data <- list(
  RNAseq = rna_mofa_scaled,
  Mutations = mutation_mofa,
  Drug_AUC = drug_mofa_scaled
)

mofa_dimension_table <- tibble(
  view = names(mofa_data),
  features = map_int(mofa_data, nrow),
  samples = map_int(mofa_data, ncol),
  measured_samples = map_int(mofa_data, ~sum(colSums(!is.na(.x)) > 0)),
  missing_values = map_int(mofa_data, ~sum(is.na(.x))),
  missing_fraction = map_dbl(mofa_data, ~mean(is.na(.x))),
  likelihood = c("gaussian", "bernoulli", "gaussian")
)

mofa_sample_coverage <- tibble(sampleID = mofa_all_samples) %>%
  mutate(
    RNAseq = sampleID %in% colnames(rna_hvg),
    Mutations = sampleID %in% colnames(mutation_matrix_filtered),
    Drug_AUC = sampleID %in% colnames(drug_matrix_filtered),
    measured_views = RNAseq + Mutations + Drug_AUC
  )

mofa_input_object <- create_mofa(mofa_data)
```

```{r mofa-dimensions-table}
knitr::kable(
  mofa_dimension_table %>%
    mutate(missing_fraction = percent(missing_fraction, accuracy = 0.1)),
  caption = "Final MOFA view dimensions after feature filtering and union-sample alignment."
)
```

The integrated analysis now uses the union of samples across the selected views, giving `r format_int(length(mofa_all_samples))` samples rather than restricting MOFA to only the `r format_int(length(complete_case_mofa_samples))` complete-case samples. This follows the MOFA practical from the course: if a sample is missing an entire view, that view is represented by `NA` values and MOFA learns from the views that are available for that sample. The final views contain `r format_int(nrow(rna_mofa_scaled))` RNA HVGs, `r format_int(nrow(mutation_mofa))` recurrent mutation features, and `r format_int(nrow(drug_mofa_scaled))` filtered drug AUC features. The missing values are therefore meaningful assay missingness, not imputed measurements.

```{r mofa-sample-coverage-table}
knitr::kable(
  mofa_sample_coverage %>%
    dplyr::count(measured_views, name = "n_samples") %>%
    arrange(desc(measured_views)),
  caption = "Number of MOFA samples with one, two, or three measured omics views."
)
```

This sample-coverage table makes the benefit of the union approach explicit. Samples with all three views are still the most informative for cross-omics factors, but samples with only one or two views do not have to be discarded. Their factor values are inferred from the views that exist for them, while missing whole-view blocks remain missing. This is more appropriate than mean-imputing an absent assay because an unmeasured RNA, mutation, or drug view is not a biological zero.

```{r mofa-data-overview, fig.cap="MOFA input overview showing observed and missing blocks across views before training."}
plot_data_overview(mofa_input_object)
```

The MOFA data overview visualises exactly where each view is observed. White/coloured blocks represent measured values and grey blocks represent missing values. This plot is important because it confirms that the updated MOFA input no longer relies only on complete-case sample overlap. Instead, it preserves incomplete samples and leaves missing drug responses or missing whole views for MOFA to handle during model fitting.

# Task 7: Train a MOFA2 Model

```{r mofa-training}
mofa_model_file <- file.path(work_dir, "MoBi_MultiOmic_4745778_Chua_mofa_model.hdf5")
mofa_train_error <- NULL
mofa_available <- FALSE
mofa_model <- NULL

train_mofa_model <- function() {
  if (file.exists(mofa_model_file)) {
    loaded <- tryCatch(
      MOFA2::load_model(mofa_model_file, load_data = TRUE, verbose = FALSE),
      error = function(e) NULL
    )
    if (inherits(loaded, "MOFA")) {
      loaded_samples <- tryCatch(
        unique(as.character(get_factors(loaded, as.data.frame = TRUE)$sample)),
        error = function(e) character()
      )
      loaded_views <- tryCatch(
        sort(unique(as.character(get_weights(loaded, as.data.frame = TRUE)$view))),
        error = function(e) character()
      )
      if (setequal(loaded_samples, mofa_all_samples) && setequal(loaded_views, names(mofa_data))) {
        return(loaded)
      }
    }
  }

  if (file.exists(mofa_model_file)) {
    file.remove(mofa_model_file)
  }

  model <- mofa_input_object
  data_options <- get_default_data_options(model)
  model_options <- get_default_model_options(model)
  training_options <- get_default_training_options(model)

  model_options$num_factors <- min(8, max(2, length(mofa_all_samples) - 1))
  model_options$likelihoods <- c(
    RNAseq = "gaussian",
    Mutations = "bernoulli",
    Drug_AUC = "gaussian"
  )

  training_options$maxiter <- 600
  training_options$convergence_mode <- "medium"
  training_options$seed <- seed_value
  training_options$verbose <- FALSE

  prepared_model <- prepare_mofa(
    model,
    data_options = data_options,
    model_options = model_options,
    training_options = training_options
  )

  run_mofa(
    prepared_model,
    outfile = mofa_model_file,
    save_data = TRUE,
    use_basilisk = FALSE
  )
}

mofa_model <- tryCatch(
  train_mofa_model(),
  error = function(e) {
    mofa_train_error <<- conditionMessage(e)
    NULL
  }
)
mofa_available <- inherits(mofa_model, "MOFA")

mofa_status <- if (mofa_available) {
  paste(mofa_model@status, collapse = ", ")
} else {
  paste("MOFA model not available:", mofa_train_error)
}

mofa_convergence <- if (mofa_available && "converged" %in% names(mofa_model@training_stats)) {
  as.character(mofa_model@training_stats$converged)
} else if (mofa_available) {
  "Training completed; explicit convergence flag not present in training_stats."
} else {
  "Not available"
}
```

```{r mofa-basic-summary}
cat("MOFA status:", mofa_status, "\n")
cat("Convergence:", mofa_convergence, "\n\n")
if (mofa_available) {
  print(mofa_model)
}
```

The model uses explicit likelihoods matched to each data type. Gaussian for variance-stabilised RNA expression, Bernoulli for binary mutation calls, and Gaussian for continuous drug AUC values. The requested number of latent factors is set to `r min(8, max(2, length(mofa_all_samples) - 1))`, which is enough to capture several independent axes while remaining conservative for `r format_int(length(mofa_all_samples))` union-aligned samples. The model is saved to `MoBi_MultiOmic_4745778_Chua_mofa_model.hdf5`, so repeated rendering can reuse the trained model when the saved sample set matches the current MOFA input. If the input dimensions change, the model is retrained rather than mixing old complete-case results with the updated missing-view design.

If the MOFA backend fails on a local machine, the report records the backend error rather than silently fabricating results. In a successful render, the sections below are populated directly from the trained MOFA model.

## Question 6: MOFA factors versus PCA components

A principal component is a linear axis computed from one data matrix, chosen to explain maximal variance in that matrix under orthogonality constraints. PCA is therefore a single-view dimensionality reduction method: RNA PCA explains RNA variation, mutation PCA explains mutation variation, and drug PCA explains drug-response variation separately.

A MOFA factor is a latent variable inferred jointly across multiple data views with view-specific weights and likelihoods. A factor can be shared across RNA, mutation, and drug response, or it can be mostly specific to one view. This makes MOFA more flexible than PCA for multi-omics data: it can identify coordinated cross-omics biology while still allowing one layer to contain signals that are not present in the others.

# Task 8: Interpret MOFA Factors

```{r mofa-extract-results}
if (mofa_available) {
  ve_raw <- get_variance_explained(mofa_model, as.data.frame = TRUE)
  ve_factor <- ve_raw$r2_per_factor %>%
    as_tibble() %>%
    mutate(
      factor = as.character(factor),
      view = as.character(view),
      r2_percent = value
    )
  ve_total <- ve_raw$r2_total %>%
    as_tibble() %>%
    mutate(
      view = as.character(view),
      r2_percent = value
    )

  factor_values <- get_factors(mofa_model, as.data.frame = TRUE) %>%
    as_tibble() %>%
    mutate(factor = as.character(factor), sampleID = as.character(sample)) %>%
    select(sampleID, group, factor, value)

  factor_scores_wide <- factor_values %>%
    select(sampleID, factor, value) %>%
    pivot_wider(names_from = factor, values_from = value)

  weights_df <- get_weights(mofa_model, as.data.frame = TRUE) %>%
    as_tibble() %>%
    mutate(
      factor = as.character(factor),
      view = as.character(view),
      feature = as.character(feature),
      abs_value = abs(value),
      direction = if_else(value >= 0, "Positive", "Negative")
    )

  mofa_factor_cols <- setdiff(names(factor_scores_wide), "sampleID")
  mofa_assoc <- association_tests(factor_scores_wide, metadata, mofa_factor_cols)
  mofa_assoc_compare <- association_method_comparison(factor_scores_wide, metadata, mofa_factor_cols)

  factor_strength <- ve_factor %>%
    group_by(factor) %>%
    summarise(
      total_r2_percent = sum(r2_percent, na.rm = TRUE),
      max_view_r2_percent = max(r2_percent, na.rm = TRUE),
      strongest_view = view[which.max(r2_percent)],
      n_views_above_1pct = sum(r2_percent >= 1, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(total_r2_percent))

  strongest_overall_factor <- factor_strength$factor[[1]]
  drug_factor <- ve_factor %>%
    filter(view == "Drug_AUC") %>%
    arrange(desc(r2_percent)) %>% dplyr::slice(1) %>%
    pull(factor)
  mutation_factor <- ve_factor %>%
    filter(view == "Mutations") %>%
    arrange(desc(r2_percent)) %>% dplyr::slice(1) %>%
    pull(factor)
  multi_view_factor_candidates <- ve_factor %>%
    group_by(factor) %>%
    summarise(
      total_r2_percent = sum(r2_percent, na.rm = TRUE),
      n_views_above_1pct = sum(r2_percent >= 1, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(n_views_above_1pct >= 2) %>%
    arrange(desc(total_r2_percent))
  biological_factor <- if (nrow(multi_view_factor_candidates) > 0) {
    multi_view_factor_candidates$factor[[1]]
  } else {
    strongest_overall_factor
  }

  top_weights_for <- function(view_name, factor_name, n = 12) {
    weights_df %>%
      filter(view == view_name, factor == factor_name) %>%
      arrange(desc(abs_value)) %>%
      slice_head(n = n)
  }

  top_drug_weights <- top_weights_for("Drug_AUC", drug_factor, 12)
  top_mutation_weights <- top_weights_for("Mutations", mutation_factor, 12)
  top_bio_weights <- weights_df %>%
    filter(factor == biological_factor) %>%
    group_by(view) %>%
    arrange(desc(abs_value), .by_group = TRUE) %>%
    slice_head(n = 8) %>%
    ungroup()
} else {
  ve_factor <- tibble()
  ve_total <- tibble()
  factor_values <- tibble()
  factor_scores_wide <- tibble()
  weights_df <- tibble()
  mofa_assoc <- tibble()
  mofa_assoc_compare <- tibble()
  factor_strength <- tibble()
  strongest_overall_factor <- NA_character_
  drug_factor <- NA_character_
  mutation_factor <- NA_character_
  biological_factor <- NA_character_
  top_drug_weights <- tibble()
  top_mutation_weights <- tibble()
  top_bio_weights <- tibble()
}
```

```{r mofa-factor-strength-table}
if (mofa_available) {
  knitr::kable(
    factor_strength %>%
      mutate(across(ends_with("percent"), ~format_num(.x, 2))),
    caption = "MOFA factors ranked by total variance explained across views."
  )
} else {
  cat("MOFA results are unavailable because model training did not complete: ", mofa_train_error)
}
```

The factor-strength table shows that `r if (mofa_available) strongest_overall_factor else "the strongest factor"` is the largest overall factor, explaining `r if (mofa_available) format_num(factor_strength$total_r2_percent[1], 2) else "NA"`% summed variance across views. Its strongest contribution is from `r if (mofa_available) factor_strength$strongest_view[1] else "NA"`, but it has non-trivial support in `r if (mofa_available) format_int(factor_strength$n_views_above_1pct[1]) else "NA"` views. This is why it is treated as the main integrated biological factor rather than as a purely view-specific axis.

```{r mofa-variance-plot, fig.cap="Variance explained by each MOFA factor in each omics view."}
if (mofa_available) {
  ggplot(ve_factor, aes(x = factor, y = r2_percent, fill = view)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    labs(title = "MOFA variance explained per factor", x = "Factor", y = "Variance explained (%)", fill = "View") +
    theme_report() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
} else {
  empty_plot("MOFA variance explained", "MOFA model unavailable")
}
```

The per-factor variance plot shows which omics layer drives each latent factor. In this model, `r if (mofa_available) strongest_overall_factor else "Factor 1"` explains `r if (mofa_available) format_num(ve_factor$r2_percent[ve_factor$factor == strongest_overall_factor & ve_factor$view == "RNAseq"], 2) else "NA"`% of RNAseq variance, `r if (mofa_available) format_num(ve_factor$r2_percent[ve_factor$factor == strongest_overall_factor & ve_factor$view == "Drug_AUC"], 2) else "NA"`% of Drug_AUC variance, and `r if (mofa_available) format_num(ve_factor$r2_percent[ve_factor$factor == strongest_overall_factor & ve_factor$view == "Mutations"], 2) else "NA"`% of mutation variance. The strongest mutation-specific signal is `r if (mofa_available) mutation_factor else "NA"`, which explains `r if (mofa_available) format_num(max(ve_factor$r2_percent[ve_factor$view == "Mutations"], na.rm = TRUE), 2) else "NA"`% of the mutation view. This separation is useful: some factors are shared across views, while others are mainly capturing genotype or expression-specific structure.

```{r mofa-total-variance-plot, fig.cap="Total variance explained by the complete MOFA model in each view."}
if (mofa_available) {
  ggplot(ve_total, aes(x = fct_reorder(view, r2_percent), y = r2_percent, fill = view)) +
    geom_col(width = 0.7, show.legend = FALSE) +
    coord_flip() +
    labs(title = "Total variance explained by view", x = NULL, y = "Variance explained (%)") +
    theme_report()
} else {
  empty_plot("Total MOFA variance explained", "MOFA model unavailable")
}
```

The total variance plot shows that MOFA captures the largest fraction of the RNAseq view (`r if (mofa_available) format_num(ve_total$r2_percent[ve_total$view == "RNAseq"], 2) else "NA"`%), followed by Drug_AUC (`r if (mofa_available) format_num(ve_total$r2_percent[ve_total$view == "Drug_AUC"], 2) else "NA"`%) and Mutations (`r if (mofa_available) format_num(ve_total$r2_percent[ve_total$view == "Mutations"], 2) else "NA"`%). This is expected because RNA expression has many continuous features and strong covariance structure, while the mutation view is sparse and binary. The lower mutation variance explained does not mean mutations are unimportant. It means the targeted mutation panel contains discrete events that are harder to summarise with a small number of continuous latent factors.

```{r mofa-factor-scatter, fig.cap="MOFA factor values for each union-aligned sample. Points are coloured using the selected metadata variable when available."}
if (mofa_available && all(c("Factor1", "Factor2") %in% names(factor_scores_wide))) {
  mofa_scatter_df <- factor_scores_wide %>% left_join(metadata, by = "sampleID")
  p <- ggplot(mofa_scatter_df, aes(x = Factor1, y = Factor2))
  if (!is.null(metadata_color) && metadata_color %in% names(mofa_scatter_df)) {
    p <- p + geom_point(aes(color = .data[[metadata_color]]), size = 2.6, alpha = 0.9)
  } else {
    p <- p + geom_point(color = "#2B6CB0", size = 2.6, alpha = 0.9)
  }
  p +
    labs(
      title = "MOFA Factor 1 versus Factor 2",
      x = "Factor 1",
      y = "Factor 2",
      color = if (!is.null(metadata_color)) str_replace_all(metadata_color, "_", " ") else NULL
    ) +
    theme_report()
} else {
  empty_plot("MOFA factor scatter", "Factor scores unavailable")
}
```

The Factor 1 versus Factor 2 scatter plot shows sample positions in the integrated latent space. Factor 1 is the strongest multi-view axis in this model, with support from RNAseq and Drug_AUC and a smaller mutation contribution. Factor 2 is less dominant overall but is strongly associated with metadata in the association table below, especially white blood cell count. Therefore, separation in this plot should be read as integrated molecular and functional variation, not as the same object as RNA PC1 or drug PC1.

```{r mofa-factor-metadata-plot, fig.cap="Distribution of the strongest overall MOFA factor across the selected metadata grouping."}
if (mofa_available && !is.na(strongest_overall_factor) && !is.null(metadata_color)) {
  plot_df <- factor_values %>%
    filter(factor == strongest_overall_factor) %>%
    left_join(metadata, by = "sampleID")
  if (metadata_color %in% names(plot_df) && n_distinct(plot_df[[metadata_color]], na.rm = TRUE) <= 8) {
    ggplot(plot_df, aes(x = .data[[metadata_color]], y = value, fill = .data[[metadata_color]])) +
      geom_boxplot(alpha = 0.45, outlier.shape = NA, show.legend = FALSE) +
      geom_jitter(width = 0.15, alpha = 0.75, size = 1.8, show.legend = FALSE) +
      labs(
        title = paste("MOFA", strongest_overall_factor, "by", str_replace_all(metadata_color, "_", " ")),
        x = str_replace_all(metadata_color, "_", " "),
        y = paste(strongest_overall_factor, "value")
      ) +
      theme_report() +
      theme(axis.text.x = element_text(angle = 35, hjust = 1))
  } else {
    empty_plot("MOFA factor values", "Selected metadata variable is not suitable for grouped plotting")
  }
} else {
  empty_plot("MOFA factor values", "MOFA model or metadata grouping unavailable")
}
```

The box and jitter plot checks whether the strongest MOFA factor differs across the selected metadata grouping (`r ifelse(is.null(metadata_color), "not available", str_replace_all(metadata_color, "_", " "))`). The jittered points are important because the MOFA analysis contains `r format_int(length(mofa_all_samples))` union-aligned samples with uneven view coverage, so a group summary alone could hide uneven sample sizes, missing-view structure, or outliers. If the boxes overlap strongly, the factor is better interpreted through its feature weights and continuous metadata associations than as a clean categorical clinical split.

```{r mofa-metadata-associations}
if (mofa_available) {
  knitr::kable(
    mofa_assoc %>%
      slice_head(n = 12) %>%
      mutate(across(c(p_value, p_adj), format_p)),
    caption = "Top exploratory metadata associations with MOFA factors."
  )
} else {
  cat("MOFA factor-metadata associations are unavailable.")
}
```

The MOFA metadata association table gives a more detailed interpretation of the factor scores. The top association is: `r format_top_association(mofa_assoc)` For categorical metadata, the estimate is eta-squared from one-way ANOVA, so larger values mean the group labels explain more of the factor-score spread. For continuous metadata, the estimate is the linear-model slope. The p_value is the raw exploratory p-value, and p_adj is the Benjamini-Hochberg adjusted value across all factor-metadata tests. The strongest adjusted associations suggest which clinical or sample variables help explain a factor, but they should still be treated as hypotheses because MOFA itself is unsupervised.

```{r mofa-anova-kruskal-table}
if (mofa_available) {
  knitr::kable(
    format_comparison_for_table(mofa_assoc_compare, n = 10),
    caption = "Side-by-side one-way ANOVA and Kruskal-Wallis comparisons for categorical metadata associations with MOFA factors."
  )
} else {
  cat("MOFA factor method comparisons are unavailable.")
}
```

The MOFA method comparison table asks whether categorical factor associations are supported by both group-mean testing and rank-based testing. The strongest categorical comparison is: `r if (mofa_available) format_top_method_comparison(mofa_assoc_compare) else "not available"` Agreement between ANOVA and Kruskal-Wallis suggests that the factor separates metadata groups in a way that is not only caused by one or two extreme samples. If the tests disagree, I would interpret the metadata association cautiously and check the factor boxplot before giving it biological meaning.

```{r mofa-association-heatmap, fig.width=10, fig.height=6, fig.cap="ANOVA p-value heatmap for metadata associations with MOFA factor scores."}
if (mofa_available) {
  plot_association_heatmap(mofa_assoc, "MOFA factor metadata association heatmap")
} else {
  empty_plot("MOFA factor metadata association heatmap", "MOFA model unavailable")
}
```

The MOFA association heatmap is interpreted differently from the PCA heatmaps. A dark cell here means that a metadata variable is associated with an integrated latent factor, not only with one omics matrix. If a factor also explains variance in multiple views and has interpretable feature weights, the metadata hit becomes more biologically useful. If a factor is mostly view-specific or explains little variance, even a significant adjusted p-value should be treated as a weaker exploratory signal.

```{r mofa-anova-kruskal-heatmap, fig.width=10, fig.height=7, fig.cap="Side-by-side one-way ANOVA and Kruskal-Wallis adjusted p-value heatmaps for categorical metadata associations with MOFA factors."}
if (mofa_available) {
  plot_method_comparison_heatmap(
    mofa_assoc_compare,
    "One-way ANOVA and Kruskal-Wallis MOFA factor comparison"
  )
} else {
  empty_plot("MOFA factor method comparison", "MOFA model unavailable")
}
```

The side-by-side MOFA heatmap makes the test comparison visible rather than leaving it only in a table. Factors with dark tiles in both rows for the same metadata variable have more stable categorical associations. Factors that are dark only in the ANOVA row may be driven mainly by mean differences, while factors that are dark only in the Kruskal-Wallis row may reflect ranked group separation or distributional differences that are not captured well by a mean comparison.

```{r mofa-weight-plots, fig.width=11, fig.height=7, fig.cap="Top absolute feature weights for selected MOFA factors. Positive and negative weights represent opposite directions along the same latent factor."}
if (mofa_available) {
  factors_to_plot <- unique(c(biological_factor, drug_factor, mutation_factor, strongest_overall_factor))
  factors_to_plot <- factors_to_plot[!is.na(factors_to_plot)]
  factors_to_plot <- head(factors_to_plot, 2)
  plot_weight_panel <- function(factor_name) {
    weights_df %>%
      filter(factor == factor_name) %>%
      group_by(view) %>%
      arrange(desc(abs_value), .by_group = TRUE) %>%
      slice_head(n = 8) %>%
      ungroup() %>%
      mutate(
        feature_label = paste(view, feature, sep = ": "),
        feature_label = fct_reorder(feature_label, value)
      ) %>%
      ggplot(aes(x = value, y = feature_label, fill = direction)) +
      geom_col(width = 0.75) +
      scale_fill_manual(values = c("Positive" = "#B45309", "Negative" = "#2563EB")) +
      labs(title = paste("Top weights for", factor_name), x = "Weight", y = NULL, fill = "Direction") +
      theme_report(base_size = 10)
  }
  if (length(factors_to_plot) == 1) {
    plot_weight_panel(factors_to_plot[[1]])
  } else {
    plot_weight_panel(factors_to_plot[[1]]) | plot_weight_panel(factors_to_plot[[2]])
  }
} else {
  empty_plot("MOFA weights", "MOFA model unavailable")
}
```

The weight plots identify the features that define the selected MOFA factors. For the main biological factor, the largest drug weights include `r if (mofa_available) paste(head(top_bio_weights$feature[top_bio_weights$view == "Drug_AUC"], 4), collapse = ", ") else "not available"`, while the mutation weights include `r if (mofa_available) paste(head(top_bio_weights$feature[top_bio_weights$view == "Mutations"], 4), collapse = ", ") else "not available"`. Features with large absolute weights are the ones most responsible for moving samples along that factor. The sign is interpreted relative to the factor value: high-factor samples tend to have higher scaled values for positive-weight features and lower scaled values for negative-weight features, but the overall sign of a MOFA factor is arbitrary.

```{r mofa-factor-heatmap, fig.cap="Heatmap of MOFA factor values across union-aligned samples."}
if (mofa_available && ncol(factor_scores_wide) > 2) {
  factor_mat <- factor_scores_wide %>%
    column_to_rownames("sampleID") %>%
    as.matrix()
  pheatmap(
    t(factor_mat),
    scale = "row",
    show_colnames = FALSE,
    clustering_method = "ward.D2",
    main = "MOFA factor values"
  )
}
```

The factor heatmap shows the eight learned factors across the `r format_int(length(mofa_all_samples))` union-aligned samples. Distinct row patterns indicate that the model has learned several different axes rather than repeatedly estimating the same signal. This complements the variance plot: factors with low variance explained may still capture specific substructure, but the main interpretation should focus on factors with clear variance support, interpretable weights, and plausible metadata associations.

## MOFA interpretation summary

```{r mofa-interpretation-text, results='asis'}
feature_list <- function(df, n = 6) {
  if (nrow(df) == 0) return("not available")
  df %>%
    slice_head(n = n) %>%
    mutate(txt = paste0(feature, " (", direction, ", weight ", format_num(value, 3), ")")) %>%
    pull(txt) %>%
    paste(collapse = ", ")
}

view_support_sentence <- function(factor_name) {
  if (!mofa_available || is.na(factor_name)) return("not available")
  ve_factor %>%
    filter(factor == factor_name) %>%
    arrange(desc(r2_percent)) %>%
    mutate(txt = paste0(view, " ", format_num(r2_percent, 2), "%")) %>%
    pull(txt) %>%
    paste(collapse = "; ")
}

if (mofa_available) {
  cat(
    "The strongest overall factor is **", strongest_overall_factor, "**, with support across views summarised as ",
    view_support_sentence(strongest_overall_factor), ". The most drug-associated factor is **", drug_factor,
    "**, and the most mutation-associated factor is **", mutation_factor, "**. The selected multi-view biological factor is **",
    biological_factor, "**; this choice prioritises factors with non-trivial variance explained in at least two views when such factors are present.\n\n",
    sep = ""
  )

  cat(
    "For the selected biological factor, the most strongly weighted features across available views are ",
    feature_list(top_bio_weights, 10),
    ". These weights identify features that move most strongly along the factor axis, but the sign of a factor is arbitrary: reversing all factor values and weights would represent the same model.\n\n",
    sep = ""
  )

  cat(
    "The top metadata association with MOFA factors is: ",
    format_top_association(mofa_assoc),
    " This association is exploratory and should be read together with the variance explained and feature weights rather than treated as a standalone clinical result.\n",
    sep = ""
  )
} else {
  cat(
    "MOFA model results are not available in this render because training did not complete. The recorded error was: `",
    mofa_train_error,
    "`. The report code still contains the full MOFA2 workflow and will populate this section when the MOFA backend is available.\n",
    sep = ""
  )
}
```

## Question 7: biological interpretation of one MOFA factor

```{r question-7-answer, results='asis'}
if (mofa_available) {
  bio_support <- ve_factor %>%
    filter(factor == biological_factor) %>%
    arrange(desc(r2_percent))
  bio_views <- bio_support %>%
    filter(r2_percent >= 1 | row_number() == 1) %>%
    mutate(txt = paste0(view, " (", format_num(r2_percent, 2), "% variance explained)")) %>%
    pull(txt) %>%
    paste(collapse = ", ")
  bio_top_by_view <- top_bio_weights %>%
    group_by(view) %>%
    summarise(features = paste(head(feature, 5), collapse = ", "), .groups = "drop") %>%
    mutate(txt = paste0(view, ": ", features)) %>%
    pull(txt) %>%
    paste(collapse = "; ")

  cat(
    "I selected **", biological_factor, "** for biological interpretation because it is the strongest factor with evidence from more than one view when such a factor is available. Its view support is: ",
    bio_views,
    ". The leading weighted features are ",
    bio_top_by_view,
    ".\n\n",
    sep = ""
  )
  cat(
    "A cautious interpretation is that ", biological_factor,
    " represents a coordinated AML axis linking the views listed above. If RNAseq contributes strongly, the factor is likely to capture a transcriptional programme or cellular composition/state gradient. If Drug_AUC contributes strongly, the same axis is also related to ex vivo sensitivity or resistance. If Mutations contributes strongly, discrete genotype events help define the samples at one end of the factor. The exact biological label should be assigned cautiously because MOFA is unsupervised and the assignment would ideally be strengthened by pathway enrichment of RNA weights, validation of mutation groups, and independent drug-response replication.\n",
    sep = ""
  )
} else {
  cat(
    "Question 7 cannot be answered from trained MOFA factors in this render because the MOFA backend did not complete. Once the model trains, the code above selects the strongest multi-view factor and reports its supporting views and top weighted features directly from the model.\n"
  )
}
```

## Question 8: drug-driven MOFA factor and contributing drugs

```{r question-8-answer, results='asis'}
if (mofa_available) {
  drug_view_r2 <- ve_factor %>%
    filter(view == "Drug_AUC", factor == drug_factor) %>%
    pull(r2_percent)
  cat(
    "The factor most strongly driven by drug response is **", drug_factor,
    "**, which explains ", format_num(drug_view_r2, 2),
    "% of the Drug_AUC view variance. The drugs with the largest absolute weights on this factor are ",
    feature_list(top_drug_weights, 10),
    ".\n\n",
    sep = ""
  )
  cat(
    "As lower AUC indicates sensitivity and higher AUC indicates resistance, the sign of each drug weight should be interpreted relative to the sample factor value. Samples with high values of ",
    drug_factor,
    " tend to have higher scaled AUC for positively weighted drugs and lower scaled AUC for negatively weighted drugs, with the reverse pattern for samples with low factor values. The factor therefore describes a coordinated drug-response profile rather than sensitivity to a single inhibitor alone.\n",
    sep = ""
  )
} else {
  cat(
    "Question 8 cannot be answered from trained MOFA factors in this render because the MOFA backend did not complete. When MOFA trains, the code identifies the factor with the highest Drug_AUC variance explained and lists the top weighted drugs.\n"
  )
}
```

## Question 9: mutation-driven MOFA factor and contributing genes

```{r question-9-answer, results='asis'}
if (mofa_available) {
  mutation_view_r2 <- ve_factor %>%
    filter(view == "Mutations", factor == mutation_factor) %>%
    pull(r2_percent)
  cat(
    "The factor most strongly driven by mutations is **", mutation_factor,
    "**, which explains ", format_num(mutation_view_r2, 2),
    "% of the Mutations view variance. The mutation genes with the largest absolute weights on this factor are ",
    feature_list(top_mutation_weights, 10),
    ".\n\n",
    sep = ""
  )
  cat(
    "This factor should be interpreted as a genotype axis: samples at opposite ends differ in the probability of carrying the listed mutation events. Because the mutation panel is sparse and targeted, the factor may reflect one dominant mutation such as NPM1 or FLT3-ITD, a co-mutation pattern, or a contrast between mutually exclusive mutation groups. The strongest interpretation is obtained when the mutation factor also explains variance in RNAseq or Drug_AUC, because that would suggest a functional consequence beyond mutation status alone.\n",
    sep = ""
  )
} else {
  cat(
    "Question 9 cannot be answered from trained MOFA factors in this render because the MOFA backend did not complete. When MOFA trains, the code identifies the factor with the highest Mutations variance explained and lists the top weighted mutation genes.\n"
  )
}
```

# Discussion and limitations

This AML dataset contains complementary molecular and pharmacological information. RNA-seq provides a broad view of transcriptional state, mutation calls provide interpretable genotype events, and drug AUC measurements provide functional response phenotypes. The sample overlap analysis shows that only `r format_int(length(complete_case_mofa_samples))` samples have all three selected views, but the updated MOFA input keeps `r format_int(length(mofa_all_samples))` union-aligned samples by leaving missing views as `NA`. This is preferable to discarding every incomplete sample, although complete-case samples remain the most informative for learning factors shared across all views.

The single-omics analyses show why integration is valuable. RNA-seq contains many features and captures continuous variation, but its PCs can reflect a mixture of biology, specimen composition, and technical factors. Mutation data are sparse and clinically meaningful, but many AML-relevant events are discrete and may not form strong continuous PCA gradients. Drug response is functionally important, yet missingness and assay-specific variability can influence its structure. MOFA addresses these issues by allowing factors to be shared or view-specific and by using likelihoods appropriate to each data type.

The main biological conclusions should remain cautious. If MOFA identifies a drug-dominated factor, it suggests coordinated sensitivity/resistance across groups of inhibitors, but additional pathway annotation and dose-response quality control would be needed to assign mechanism. If a mutation-dominated factor is observed, top weights can identify genotype structure, but targeted mutation panels do not capture all genomic lesions, clonal fractions, copy-number changes, or epigenetic states. If a factor is supported by RNA and drug response, it is a promising candidate for a functional transcriptional programme, but pathway enrichment and external validation would strengthen the claim.

Important technical limitations include moderate complete-case sample size, missing drug-response values, sparse binary mutation features, and preprocessing choices such as expression filtering, highly variable gene selection, drug missingness threshold, and feature scaling. The ANOVA heatmaps help flag possible metadata or centre effects, but they are exploratory and do not replace a full confounder-adjusted model. MOFA is unsupervised, so factors are not guaranteed to align with known clinical variables or causal mechanisms. Future analyses could include RNA pathway enrichment on factor weights, formal differential expression across mutation-defined groups, supervised models of drug response, sensitivity analyses for filtering thresholds, batch/source adjustment where appropriate, and validation in an independent AML cohort.

# Reproducibility information

```{r session-info}
sessionInfo()
```

)----------------------------------------'
# -----------------------------------------------------------------------------
# Write the Quarto file and render HTML
# -----------------------------------------------------------------------------

if (file.exists(qmd_file) && !overwrite_qmd) {
  message("The QMD already exists and overwrite_qmd is FALSE: ", qmd_file)
} else {
  writeLines(qmd_text, qmd_file, useBytes = TRUE)
  message("Wrote Quarto file: ", qmd_file)
}

if (setup_mofa_python) {
  setup_mofa_backend(project_dir)
}

if (render_report) {
  quarto_bin <- find_quarto()
  run_command(quarto_bin, c("render", shQuote(qmd_file), "--to", "html"), "render Quarto HTML")
  if (file.exists(html_file)) {
    message("Rendered HTML report: ", html_file)
  } else {
    warning("Render command completed, but expected HTML file was not found: ", html_file)
  }
}

message("Done.")

