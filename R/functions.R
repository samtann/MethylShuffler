# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'


setwd("/data/projects/punim1157/")

library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
# library(org.Hs.eg.db)


generate_random_CpGs <- function(n_CpGs, array_type = "450k") {

  if(array_type == "EPIC") {
    if(!exists("EPIC_data")) {
      EPIC_data <<- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
      message("Annotation data loaded.")
    } else {
      message("Annotation data already loaded.")
    }
    random_CpGs = sample(EPIC_data$Name, n_CpGs, replace = FALSE)
  }
  else if(array_type == "450k") {
    if(!exists("k450_data")) {
      k450_data <<- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
      message("Annotation data loaded.")
    } else {
      message("Annotation data already loaded.")
    }
  random_CpGs = sample(k450_data$Name, n_CpGs, replace = FALSE)
  }
  else {
    stop("Invalid array type specified. PLease specify 'EPIC' or '450k'.")
  }
  return(random_CpGs)
}

generate_random_CpGs(n_CpGs = 100, array_type = "450k")


generate_random_genes = function(n_genes, array_type) {

  CpGs = generate_random_CpGs(n_genes*2, array_type)

  genes = EPIC_data$UCSC_RefGene_Name[match(CpGs, EPIC_data$Name)]
  genes_list = unlist(strsplit(genes, ";"))
  return(unique(genes_list)[1:n_genes])
}

generate_random_genes(50, "EPIC")




extract_beta_values <- function(beta_matrix, CpG_list) {

  valid_CpGs = CpG_list[CpG_list %in% rownames(beta_matrix)]

  if(length(valid_CpGs < length(CpG_list))) {
    message(paste0("Only ", length(valid_CpGs), " of ", length(CpG_list), " requested CpGs available in this dataset."))
  }

  beta_values = data.frame(t(beta_matrix[valid_CpGs, ]))

  return(beta_values)
}

betas = extract_beta_values(ex, generate_random_CpGs(100, "450k"))



library(ggplot2)
library(progress)

compute_empirical_p_value <- function(beta_matrix, covariates_df, outcome_str, covariate_str, repeats, n_CpGs, p_value_interest, array_type = "450k", plot_distribution = TRUE) {

  pb = progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                        total = repeats,  # Number of total iterations
                        complete = "=",   # Completion bar character
                        incomplete = "-", # Incomplete bar character
                        current = ">", # Current bar character
                        clear = FALSE, # If TRUE, clears the bar when finish
                        width = 100) # Width of the progress bar

  print("Computing principal component 1 across DNA methylation dataset to use as covariate...")
  pc1_beta_matrix = prcomp(t(beta_matrix), scale. = TRUE)$x[, 1]

  print("Generating empirical null distribution...")

  empirical_p_values = replicate(repeats, {

    pb$tick()
    beta_values = suppressMessages(extract_beta_values(beta_matrix, generate_random_CpGs(n_CpGs, array_type)))
    effect_sizes = rnorm(nrow(beta_values))  # random effect sizes from normal distribution
    # score = rowSums(beta_values * effect_sizes)
    pca_result = prcomp(beta_values, scale. = TRUE)  # [,1:150]
    score = pca_result$x[, 1]

    model_formula = as.formula(paste(outcome_str, "~ score + pc1_beta_matrix", covariate_str))
    model = glm(model_formula, data = cbind(covariates_df, score), family = binomial())

    p_value = summary(model)$coefficients["score", "Pr(>|z|)"]

  })

  empirical_p_value = mean(empirical_p_values < p_value_interest)

  if (plot_distribution) {

    # Create Q-Q plot of random p-values
    p_values_df = data.frame(empirical_p_values = empirical_p_values)
    theoretical_quantiles = qunif(ppoints(repeats))
    qq_plot = ggplot(p_values_df, aes(sample = empirical_p_values)) +
      stat_qq(distribution = qunif, dparams = list(min = 0, max = 1)) +
      geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
      labs(title = "Q-Q Plot of Empirical p-values vs. Uniform Distribution",
           x = "Theoretical Quantiles",
           y = "Empirical Quantiles") +
      theme_minimal()
    print(qq_plot) # This line will force the plot to be displayed

    # Create histogram of random p-values with a vertical line for the p-value of interest
    # Create a logical vector for fill color based on p-value interest
    fill_colors = empirical_p_values < p_value_interest
    histogram_plot = ggplot(data.frame(p_value = empirical_p_values, fill_color = fill_colors), aes(x = p_value, fill = fill_color)) +
      geom_histogram(breaks = seq(0, 1, by = 3/repeats),
                     binwidth = 1/repeats,
                     show.legend = FALSE) +
      scale_fill_manual(values = c("TRUE" = "orange", "FALSE" = "grey30")) +
      geom_vline(xintercept = p_value_interest, color = "blue", size = 1) +
      labs(title = "Histogram of Random p-values",
           x = "p-values",
           y = "Frequency") +
      theme_minimal()
      # annotate("text", x = 0.95, y = Inf, label = paste0("Empirical p-value =", round(empirical_p_value, 4)),
      #          hjust = 1, vjust = 1, size = 5, color = "black", parse = TRUE)
    print(histogram_plot)
  }

  return(empirical_p_value)

  # add qq plot
  # allow choice of score or pc1
}



pheno_data = gset@phenoData@data
pheno_data$ASD = ifelse(pheno_data$`disease state:ch1` == "control", 0, 1)
pheno_data$ASD = as.factor(pheno_data$ASD)
pheno_data$sex = as.factor(pheno_data$`Sex:ch1`)
pheno_data$age = as.numeric(pheno_data$`age:ch1`)
pheno_data$brain_region = as.factor(pheno_data$`brain region:ch1`)
compute_empirical_p_value(ex, pheno_data, "ASD", "+ sex + age + brain_region ", 10000, 150, 0.03, "450k", plot_distribution = TRUE)


