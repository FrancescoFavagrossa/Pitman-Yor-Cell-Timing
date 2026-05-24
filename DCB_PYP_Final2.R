
suppressPackageStartupMessages({
  library(data.table)
  library(tickTack) 
  library(dplyr)
})

set.seed(123)



## Parametri modello PYP
MIN_MUT_PER_SEGMENT <- 5      
N_ITER_GIBBS        <- 300   
SIGMA_LIK           <- 0.10   
PURITY              <- 1.0    
ALPHA_PYP           <- 2.0    
DISCOUNT_PYP        <- 0.5    

## Percorsi dati PCAWG
DATA_DIR <- file.path(getwd(), "Data")
FILE_MUT <- file.path(DATA_DIR, "October_2016_whitelist_2583.snv_mnv_indel.maf.xena.nonUS")
FILE_CNV <- file.path(DATA_DIR, "20170119_final_consensus_copynumber_donor")
FILE_METADATA <- file.path(DATA_DIR, "pcawg_specimen_histology_August2016_v9_donor")


load_pcawg_data <- function(file_cnv, file_mut, file_metadata = NULL) {
  cnv <- fread(file_cnv, sep = "\t", showProgress = FALSE)
  setnames(cnv, 
           old = c("chr", "start", "end", "sampleID", "total_cn", "major_cn", "minor_cn"),
           new = c("Chromosome", "Start", "End", "Sample", "Total_CN", "Major_CN", "Minor_CN"))
  cnv[, Chromosome := as.character(gsub("chr", "", Chromosome, ignore.case = TRUE))]
  cnv[, Is_Amplification := (Total_CN > 2)]
  
  mut <- fread(file_mut, showProgress = FALSE)
  setnames(mut,
           old = c("chr", "start", "DNA_VAF"),
           new = c("Chromosome", "Position", "VAF"))
  mut[, Chromosome := as.character(gsub("chr", "", Chromosome, ignore.case = TRUE))]
  mut <- mut[VAF > 0.05]
  
  metadata <- fread(file_metadata, showProgress = FALSE)
  if ("icgc_specimen_id" %in% colnames(metadata)) {
    setnames(metadata, old = "icgc_specimen_id", new = "Sample")
  }
  
  list(cnv = cnv, mut = mut, metadata = metadata)
}

get_patient_info <- function(patient_id, metadata) {
  if (is.null(metadata) || nrow(metadata) == 0) {
    return(list(tumor_type = "Unknown", full_name = patient_id))
  }
  
  patient_row <- metadata[Sample == patient_id]
  if (nrow(patient_row) == 0) {
    return(list(tumor_type = "Unknown", full_name = patient_id))
  }
  
  tumor_type <- "Unknown"
  if ("histology_abbreviation" %in% colnames(patient_row)) {
    tumor_type <- as.character(patient_row$histology_abbreviation[1])
  }
  
  list(tumor_type = tumor_type, full_name = sprintf("%s (%s)", patient_id, tumor_type))
}


determine_cn_type <- function(nA, nB) {
  total <- nA + nB
  if (total == 2 && (nA == 0 || nB == 0)) return("CNLOH")
  if (total == 3) return("Trisomy")
  if (total == 4 && nA == 2 && nB == 2) return("Tetrasomy")
  if (total > 4) return("HighAmp")
  return("Other")
}

estimate_segment_timing <- function(p_cnv, p_mut, purity = PURITY) {
  results_list <- list()
  
  for (i in 1:nrow(p_cnv)) {
    seg <- p_cnv[i]
    in_segment <- p_mut[
      Chromosome == seg$Chromosome & 
        Position >= seg$Start & 
        Position <= seg$End
    ]
    
    if (nrow(in_segment) >= MIN_MUT_PER_SEGMENT) {
      mean_vaf <- mean(in_segment$VAF, na.rm = TRUE)
      nA <- seg$Major_CN
      nB <- seg$Minor_CN
      total_cn <- nA + nB
      
      p1 <- (purity * 1) / (total_cn * purity + 2 * (1 - purity))
      p2 <- (purity * 2) / (total_cn * purity + 2 * (1 - purity))
      
      if (abs(p2 - p1) < 1e-6) next
      
      theta2 <- (mean_vaf - p1) / (p2 - p1)
      theta2 <- pmax(0, pmin(1, theta2))
      theta1 <- 1 - theta2
      
      cn_type <- determine_cn_type(nA, nB)
      
      if (cn_type == "CNLOH") {
        time_est <- (2 * theta2) / (2 * theta2 + theta1)
      } else if (cn_type == "Trisomy") {
        time_est <- (3 * theta2) / (theta1 + 2 * theta2)
      } else if (cn_type == "Tetrasomy") {
        time_est <- (2 * theta2) / (2 * theta2 + theta1)
      } else {
        time_est <- (theta2 * total_cn) / (theta1 + theta2 * total_cn)
      }
      
      time_est <- pmax(0, pmin(1, time_est))
      
      results_list[[length(results_list) + 1]] <- data.frame(
        Chromosome = seg$Chromosome,
        Start = seg$Start,
        End = seg$End,
        Major = nA,
        Minor = nB,
        CN_Type = cn_type,
        Mean_VAF = mean_vaf,
        Time = time_est,
        N_Mutations = nrow(in_segment)
      )
    }
  }
  
  if (length(results_list) == 0) return(NULL)
  do.call(rbind, results_list)
}


fit_pitman_yor_clustering <- function(times, alpha = ALPHA_PYP, discount = DISCOUNT_PYP) {
  times <- as.numeric(times)
  times <- times[!is.na(times)]
  N <- length(times)           
  K_TRUNC <- 15 
  
  z <- sample(1:K_TRUNC, N, replace = TRUE)  
  mu <- runif(K_TRUNC)
  v <- rbeta(K_TRUNC, 1 - discount, alpha) 
  w <- numeric(K_TRUNC)
  
  posterior_mu <- list()
  posterior_z  <- list()
  n_burnin <- floor(N_ITER_GIBBS * 0.5)
  
  for (iter in 1:N_ITER_GIBBS) {
    w[1] <- v[1]
    for (k in 2:K_TRUNC) {
      w[k] <- v[k] * prod(1 - v[1:(k-1)])
    }
    
    for (k in 1:K_TRUNC) {
      idx <- which(z == k)
      if (length(idx) > 0) {
        prec <- length(idx) / SIGMA_LIK^2
        mu[k] <- rnorm(1, mean = mean(times[idx]), sd = 1 / sqrt(prec + 1e-6))
      } else {
        mu[k] <- runif(1)
      }
    }
    
    for (i in 1:N) {
      log_probs <- log(w + 1e-300) + dnorm(times[i], mu, SIGMA_LIK, log = TRUE)
      probs <- exp(log_probs - max(log_probs))
      z[i] <- sample(1:K_TRUNC, 1, prob = probs / sum(probs))
    }
    
    counts <- table(factor(z, levels = 1:K_TRUNC))
    greater <- rev(cumsum(rev(counts))) - counts
    for (k in 1:(K_TRUNC - 1)) {
      v[k] <- rbeta(1, 
                    shape1 = 1 - discount + counts[k],
                    shape2 = alpha + discount * k + greater[k])
    }
    
    if (iter > n_burnin) {
      ord <- order(mu)
      posterior_mu[[iter - n_burnin]] <- mu[ord]
      posterior_z[[iter - n_burnin]]  <- z
    }
  }
  
  avg_mu <- colMeans(do.call(rbind, posterior_mu))
  active <- unique(z)
  mapping <- rank(mu[active])
  
  return(mapping[match(z, active)])
}

analyze_patient_pyp <- function(data, patient_id) {
  patient_info <- get_patient_info(patient_id, data$metadata)
  
  p_cnv <- data$cnv[Sample == patient_id & Is_Amplification == TRUE]
  p_mut <- data$mut[Sample == patient_id]
  
  segments <- estimate_segment_timing(p_cnv, p_mut)
  if (is.null(segments) || nrow(segments) == 0) return(NULL)
  
  segments$Cluster <- fit_pitman_yor_clustering(segments$Time)
  n_clusters <- length(unique(segments$Cluster))
  
  cat(sprintf("🎯 PYP: %d cluster identificati\n", n_clusters))
  
  muts_with_cluster <- list()
  for (i in 1:nrow(segments)) {
    seg <- segments[i, ]
    sub_mut <- p_mut[
      Chromosome == seg$Chromosome & 
        Position >= seg$Start & 
        Position <= seg$End
    ]
    if (nrow(sub_mut) > 0) {
      sub_mut$Cluster <- seg$Cluster
      muts_with_cluster[[i]] <- sub_mut
    }
  }
  
  mutations <- rbindlist(muts_with_cluster, fill = TRUE)
  
  list(
    segments = segments,
    mutations = mutations,
    id = patient_id,
    patient_info = patient_info,
    n_clusters = n_clusters,
    method = "PYP"
  )
}


compute_retrospective_likelihood <- function(result, data_mut, purity = PURITY) {
  cat("\n📊 Likelihood Retrospettiva...\n")
  
  segments <- result$segments
  ll_total <- 0
  ll_per_cluster <- rep(0, result$n_clusters)
  
  for (i in 1:nrow(segments)) {
    seg <- segments[i, ]
    
    muts <- data_mut[
      Sample == result$id &
        Chromosome == seg$Chromosome & 
        Position >= seg$Start & 
        Position <= seg$End
    ]
    
    if (nrow(muts) == 0) next
    
    tau <- seg$Time
    total_cn <- seg$Major + seg$Minor
    
    p1 <- (purity * 1) / (total_cn * purity + 2 * (1 - purity))
    p2 <- (purity * 2) / (total_cn * purity + 2 * (1 - purity))
    
    cn_type <- seg$CN_Type
    if (cn_type == "CNLOH") {
      theta2 <- tau / (2 - tau)
    } else if (cn_type == "Trisomy") {
      theta2 <- tau / (3 - 2*tau)
    } else if (cn_type == "Tetrasomy") {
      theta2 <- tau / (2 - tau)
    } else {
      theta2 <- tau / (total_cn - tau*(total_cn - 2))
    }
    theta2 <- pmax(0, pmin(1, theta2))
    theta1 <- 1 - theta2
    
    depth <- 100
    for (j in 1:nrow(muts)) {
      vaf <- muts$VAF[j]
      n_alt <- round(vaf * depth)
      
      lik <- theta1 * dbinom(n_alt, depth, p1) + 
        theta2 * dbinom(n_alt, depth, p2)
      ll <- log(lik + 1e-300)
      
      ll_total <- ll_total + ll
      ll_per_cluster[seg$Cluster] <- ll_per_cluster[seg$Cluster] + ll
    }
  }
  
  cat(sprintf("✓ Log-likelihood: %.2f (%.2f/seg)\n", 
              ll_total, ll_total / nrow(segments)))
  
  list(ll_total = ll_total, ll_per_cluster = ll_per_cluster)
}



posterior_predictive_vaf <- function(result, data_mut, purity = PURITY) {
  cat("\n🎲 Posterior Predictive Checks...\n")
  
  observed <- list()
  predicted <- list()
  ks_tests <- data.frame()
  
  for (k in 1:result$n_clusters) {
    segs <- result$segments[result$segments$Cluster == k, ]
    obs_vaf <- c()
    pred_vaf <- c()
    
    for (i in 1:nrow(segs)) {
      seg <- segs[i, ]
      
      muts <- data_mut[
        Sample == result$id &
          Chromosome == seg$Chromosome & 
          Position >= seg$Start & 
          Position <= seg$End
      ]
      
      if (nrow(muts) == 0) next
      obs_vaf <- c(obs_vaf, muts$VAF)
      
      tau <- seg$Time
      total_cn <- seg$Major + seg$Minor
      
      p1 <- (purity * 1) / (total_cn * purity + 2 * (1 - purity))
      p2 <- (purity * 2) / (total_cn * purity + 2 * (1 - purity))
      
      cn_type <- seg$CN_Type
      if (cn_type == "CNLOH") {
        theta2 <- tau / (2 - tau)
      } else if (cn_type == "Trisomy") {
        theta2 <- tau / (3 - 2*tau)
      } else if (cn_type == "Tetrasomy") {
        theta2 <- tau / (2 - tau)
      } else {
        theta2 <- tau / (total_cn - tau*(total_cn - 2))
      }
      theta2 <- pmax(0, pmin(1, theta2))
      theta1 <- 1 - theta2
      
      for (j in 1:nrow(muts)) {
        comp <- sample(1:2, 1, prob = c(theta1, theta2))
        p <- ifelse(comp == 1, p1, p2)
        n_alt <- rbinom(1, 100, p)
        pred_vaf <- c(pred_vaf, n_alt / 100)
      }
    }
    
    observed[[k]] <- obs_vaf
    predicted[[k]] <- pred_vaf
    
    if (length(obs_vaf) >= 5 && length(pred_vaf) >= 5) {
      ks <- ks.test(obs_vaf, pred_vaf)
      ks_tests <- rbind(ks_tests, data.frame(
        Cluster = k,
        N = length(obs_vaf),
        KS_stat = ks$statistic,
        KS_pvalue = ks$p.value
      ))
      
      status <- ifelse(ks$p.value > 0.05, "✓", "✗")
      cat(sprintf("  Cluster %d: p=%.4f %s\n", k, ks$p.value, status))
    }
  }
  
  list(observed = observed, predicted = predicted, ks_tests = ks_tests)
}



convert_to_ticktack <- function(segments, mutations, patient_id, purity = PURITY) {
  cat("\n🔄 Conversione formato TickTack...\n")
  
  cn <- segments %>%
    transmute(
      chr = as.numeric(Chromosome),
      from = Start,
      to = End,
      Major = Major,
      minor = Minor,
      CCF = 1
    )
  
  muts <- mutations %>%
    transmute(
      chr = as.numeric(Chromosome),
      from = Position,
      to = Position,
      ref = "G",
      alt = "A",
      DP = 100,
      NV = round(VAF * 100),
      VAF = VAF,
      CCF = 1
    )
  
  cat(sprintf("✓ %d segmenti, %d mutazioni\n", nrow(cn), nrow(muts)))
  
  list(
    cna = cn,
    mutations = muts,
    metadata = data.frame(sample = patient_id, purity = purity)
  )
}


analyze_patient_ticktack <- function(data, patient_id, 
                                     min_mutations = MIN_MUT_PER_SEGMENT,
                                     tolerance = 0.001,
                                     max_attempts = 100) {
  
  cat("\n🎯 TickTack Full: Hierarchical Clustering\n")
  
  patient_info <- get_patient_info(patient_id, data$metadata)
  p_cnv <- data$cnv[Sample == patient_id & Is_Amplification == TRUE]
  p_mut <- data$mut[Sample == patient_id]
  
  segments <- estimate_segment_timing(p_cnv, p_mut)
  if (is.null(segments) || nrow(segments) == 0) return(NULL)
  
  muts_with_segment <- list()
  for (i in 1:nrow(segments)) {
    seg <- segments[i, ]
    sub_mut <- p_mut[
      Chromosome == seg$Chromosome & 
        Position >= seg$Start & 
        Position <= seg$End
    ]
    if (nrow(sub_mut) > 0) {
      sub_mut$Segment <- i
      muts_with_segment[[i]] <- sub_mut
    }
  }
  mutations <- rbindlist(muts_with_segment, fill = TRUE)
  
  tt_data <- convert_to_ticktack(segments, mutations, patient_id)
  
  cat(sprintf("  min_mut=%d, tol=%.4f, attempts=%d\n",
              min_mutations, tolerance, max_attempts))
  cat("  Esecuzione tickTack::fit_h()...\n")
  
  tryCatch({
    results <- tickTack::fit_h(
      tt_data, 
      max_attempts = max_attempts,
      INIT = TRUE,
      tolerance = tolerance,
      min_mutations_number = min_mutations
    )
    
    model_selection <- tickTack::model_selection_h(results$results_timing)
    best_K <- model_selection$best_K
    best_fit <- model_selection$best_fit
    
    cat(sprintf("✓ TickTack: %d cluster (AIC)\n", best_K))
    
    segments$Cluster <- best_fit$summarized_results$clock
    mutations$Cluster <- segments$Cluster[mutations$Segment]
    
    list(
      segments = segments,
      mutations = mutations,
      id = patient_id,
      patient_info = patient_info,
      n_clusters = best_K,
      cluster_means = best_fit$summarized_results$clock_mean,
      method = "TickTack",
      model_selection = model_selection
    )
    
  }, error = function(e) {
    cat(sprintf("✗ Errore: %s\n", e$message))
    NULL
  })
}


compare_methods <- function(result_pyp, result_ticktack) {
  cat("\n📊 COMPARAZIONE: PYP vs TICKTACK\n")
  cat("==================================\n\n")
  
  cat("Numero cluster:\n")
  cat(sprintf("  PYP:      %d\n", result_pyp$n_clusters))
  cat(sprintf("  TickTack: %d\n\n", result_ticktack$n_clusters))
  
  pyp_sizes <- table(result_pyp$segments$Cluster)
  tt_sizes <- table(result_ticktack$segments$Cluster)
  
  cat("Dimensioni cluster:\n")
  cat("  PYP:     ", paste(sprintf("%d", pyp_sizes), collapse=", "), "\n")
  cat("  TickTack:", paste(sprintf("%d", tt_sizes), collapse=", "), "\n\n")
  
  calc_wcv <- function(times, clusters) {
    sizes <- table(clusters)
    n_clust <- length(unique(clusters))
    sum(tapply(times, clusters, var, na.rm=TRUE) * (sizes - 1), na.rm=TRUE) / 
      (length(times) - n_clust)
  }
  
  times <- result_pyp$segments$Time
  pyp_wcv <- calc_wcv(times, result_pyp$segments$Cluster)
  tt_wcv <- calc_wcv(times, result_ticktack$segments$Cluster)
  
  cat("Within-cluster variance:\n")
  cat(sprintf("  PYP:      %.4f\n", pyp_wcv))
  cat(sprintf("  TickTack: %.4f\n", tt_wcv))
  cat(sprintf("  Improvement: %.1f%%\n\n", 100 * (tt_wcv - pyp_wcv) / tt_wcv))
  
  pyp_means <- tapply(times, result_pyp$segments$Cluster, mean)
  tt_means <- tapply(times, result_ticktack$segments$Cluster, mean)
  
  cat("Inter-cluster separation:\n")
  cat(sprintf("  PYP:      %.3f\n", mean(dist(pyp_means))))
  cat(sprintf("  TickTack: %.3f\n\n", mean(dist(tt_means))))
  
  cat("Singleton clusters:\n")
  cat(sprintf("  PYP:      %d\n", sum(pyp_sizes == 1)))
  cat(sprintf("  TickTack: %d\n\n", sum(tt_sizes == 1)))
  
  if (!is.null(result_ticktack$model_selection)) {
    ms <- result_ticktack$model_selection
    cat("Model Selection (TickTack):\n")
    cat(sprintf("  AIC: %.2f\n", ms$best_fit$AIC))
    cat(sprintf("  BIC: %.2f\n\n", ms$best_fit$BIC))
  }
  
  list(
    n_clusters = c(PYP = result_pyp$n_clusters, 
                   TickTack = result_ticktack$n_clusters),
    wcv = c(PYP = pyp_wcv, TickTack = tt_wcv),
    singletons = c(PYP = sum(pyp_sizes == 1), 
                   TickTack = sum(tt_sizes == 1))
  )
}



run_complete_analysis <- function(data, patient_id, output_dir = ".") {
  
  cat("\n")
  cat("╔════════════════════════════════════════╗\n")
  cat("║  ANALISI COMPLETA: PYP + TICKTACK      ║\n")
  cat("╚════════════════════════════════════════╝\n")
  
  prefix <- file.path(output_dir, patient_id)
  
  ## 1. Analisi PYP
  cat("\n▶ STEP 1: PYP Clustering\n")
  result_pyp <- analyze_patient_pyp(data, patient_id)
  if (is.null(result_pyp)) {
    cat("✗ PYP fallito\n")
    return(NULL)
  }
  
  ## 2. Validazione PYP
  cat("\n▶ STEP 2: Validazione PYP\n")
  ll_pyp <- compute_retrospective_likelihood(result_pyp, data$mut)
  ppc_pyp <- posterior_predictive_vaf(result_pyp, data$mut)
  
  ## 3. Analisi TickTack
  cat("\n▶ STEP 3: TickTack Clustering\n")
  result_ticktack <- analyze_patient_ticktack(data, patient_id)
  if (is.null(result_ticktack)) {
    cat("✗ TickTack fallito\n")
    return(NULL)
  }
  
  ## 4. Validazione TickTack
  cat("\n▶ STEP 4: Validazione TickTack\n")
  ll_ticktack <- compute_retrospective_likelihood(result_ticktack, data$mut)
  
  cat(sprintf("\n→ ΔLL (PYP - TickTack): %.2f\n", 
              ll_pyp$ll_total - ll_ticktack$ll_total))
  
  ## 5. Comparazione
  cat("\n▶ STEP 5: Comparazione Metodi\n")
  comparison <- compare_methods(result_pyp, result_ticktack)
  
  ## 6. Salva risultati
  results <- list(
    patient_id = patient_id,
    result_pyp = result_pyp,
    result_ticktack = result_ticktack,
    likelihood_pyp = ll_pyp,
    likelihood_ticktack = ll_ticktack,
    ppc_pyp = ppc_pyp,
    comparison = comparison
  )
  
  saveRDS(results, file = paste0(prefix, "_CompleteResults.rds"))
  
  cat("\n")
  cat("╔════════════════════════════════════════╗\n")
  cat("║  ✓ ANALISI COMPLETATA                  ║\n")
  cat("╚════════════════════════════════════════╝\n")
  cat(sprintf("\nRisultati salvati: %s_CompleteResults.rds\n\n", prefix))
  
  return(results)
}


cat("\n")
cat("╔════════════════════════════════════════╗\n")
cat("║  CARICAMENTO DATI PCAWG                ║\n")
cat("╚════════════════════════════════════════╝\n\n")

data <- load_pcawg_data(FILE_CNV, FILE_MUT, FILE_METADATA)

cat("Selezione paziente con più amplificazioni...\n")
candidates <- data$cnv[Is_Amplification == TRUE, .N, by = Sample][order(-N)]

result <- NULL
chosen_pat <- NULL

for (pat in head(candidates$Sample, 10)) {
  segments_test <- estimate_segment_timing(
    data$cnv[Sample == pat & Is_Amplification == TRUE],
    data$mut[Sample == pat]
  )
  if (!is.null(segments_test) && nrow(segments_test) >= 50) {
    chosen_pat <- pat
    cat(sprintf("✓ Paziente selezionato: %s (%d segmenti)\n\n", 
                chosen_pat, nrow(segments_test)))
    break
  }
}

if (is.null(chosen_pat)) {
  stop("Nessun paziente valido trovato")
}

## Esegui analisi completa
results <- run_complete_analysis(data, chosen_pat, output_dir = ".")

## Stampa sommario finale
cat("\n")
cat("╔════════════════════════════════════════╗\n")
cat("║  SOMMARIO RISULTATI                    ║\n")
cat("╚════════════════════════════════════════╝\n\n")

cat(sprintf("Paziente: %s\n", results$result_pyp$patient_info$full_name))
cat(sprintf("N. segmenti: %d\n\n", nrow(results$result_pyp$segments)))

cat("CLUSTERING:\n")
cat(sprintf("  PYP:      %d cluster\n", results$result_pyp$n_clusters))
cat(sprintf("  TickTack: %d cluster\n\n", results$result_ticktack$n_clusters))

cat("LIKELIHOOD:\n")
cat(sprintf("  PYP:      %.2f\n", results$likelihood_pyp$ll_total))
cat(sprintf("  TickTack: %.2f\n", results$likelihood_ticktack$ll_total))
cat(sprintf("  Δ:        %.2f (PYP migliore)\n\n", 
            results$likelihood_pyp$ll_total - results$likelihood_ticktack$ll_total))

cat("POSTERIOR PREDICTIVE CHECKS (PYP):\n")
passed <- sum(results$ppc_pyp$ks_tests$KS_pvalue > 0.05, na.rm = TRUE)
total <- nrow(results$ppc_pyp$ks_tests)
cat(sprintf("  Cluster con buon fit: %d/%d (%.1f%%)\n\n", 
            passed, total, 100 * passed / total))

cat("✓ Analisi completata con successo!\n\n")
