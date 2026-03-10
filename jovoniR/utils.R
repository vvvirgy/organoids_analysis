
#SCE_PATH = "/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/data/scexp_smad2_karyo_all_organoids.rds"
SCE_PATH = "/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/data/scexp_karyo_all_organoids_filt.rds"
META_PATH = "/orfeo/cephfs/scratch/cdslab/vgazziero/organoids_prj/data/karyotypes_genes_filtered_scrna.rds"

my_fit_devil <- function(
    input_matrix,
    design_matrix,
    overdispersion = "MOM",
    init_overdispersion = NULL,
    init_beta_rough = FALSE,offset = 0,
    size_factors = NULL,
    verbose = FALSE,
    max_iter = 200,
    tolerance = 1e-3,
    CUDA = FALSE,
    batch_size = 1024L,
    parallel.cores = 1
) {
  # - Input parameters ----
  # Read general info about input matrix and design matrix
  gene_names <- rownames(input_matrix)
  ngenes <- nrow(input_matrix)
  nfeatures <- ncol(design_matrix)
  
  # Detect cores to use
  max.cores <- parallel::detectCores()
  if (is.null(parallel.cores)) {
    n.cores <- max.cores
  } else {
    if (parallel.cores > max.cores) {
      message("Requested ", parallel.cores, " cores, but only ", max.cores, " available.")
    }
    n.cores <- min(max.cores, parallel.cores)
  }
  
  # Check if CUDA is available
  CUDA_is_available <- FALSE
  if (CUDA) {
    if (!exists("beta_fit_gpu", mode = "function")) {
      warning(
        "CUDA support was not enabled during package installation. ",
        "The beta_fit_gpu function is not available. ",
        "Reinstall with configure.args='--with-cuda' to enable GPU acceleration: ",
        "devtools::install_github(\"caravagnalab/devil\", force=TRUE, configure.args=\"--with-cuda\"). ",
        "Falling back to CPU computation."
      )
      CUDA <- FALSE
      CUDA_is_available <- FALSE
    } else {
      message("CUDA support detected - using GPU acceleration")
      CUDA_is_available <- TRUE
    }
  }
  
  # - CPU and GPU common part (i.e. size_factors and offset_vectors) ----
  ## - Compute size factors ----
  if (!is.null(size_factors)) {
    
    
    if (class(size_factors) == "character") {
      if (verbose) {
        message("Compute size factors")
      }
      sf <- devil:::calculate_sf(input_matrix, method = size_factors, verbose = verbose)  
    } else {
      if (verbose) {
        message("Using pre-computed size factors")
      }
      sf <- size_factors
    }
  } else {
    sf <- rep(1, nrow(design_matrix))
  }
  
  ## - Compute offset vector ----
  offset_vector <- devil:::compute_offset_vector(offset, input_matrix, sf)
  
  # - Start GPU vs CPU branch ----
  if (CUDA & CUDA_is_available) {
    
    ## - GPU branch ----
    remainder <- ngenes %% batch_size
    extra_genes <- remainder
    genes_batch <- ngenes - extra_genes
    
    message("Fit beta, using CUDA acceleration")
    start_time <- Sys.time()
    res_beta_fit <- beta_fit_gpu(
      input_matrix[seq_len(genes_batch), ],
      design_matrix,
      offset_vector,
      max_iter = max_iter,
      eps = tolerance,
      batch_size = batch_size,
      TEST = FALSE
    )
    
    if (remainder > 0) {
      res_beta_fit_extra <- beta_fit_gpu(
        input_matrix[(genes_batch + 1):ngenes, ],
        design_matrix,
        offset_vector,
        max_iter = max_iter,
        eps = tolerance,
        batch_size = extra_genes,
        TEST = FALSE
      )
    }
    
    end_time <- Sys.time()
    message("[TIMING] Beta fit computing (GPU):", difftime(end_time, start_time, units = "secs"))
    
    # Extract beta and theta from GPU results
    beta <- res_beta_fit$mu_beta
    theta <- res_beta_fit$theta
    
    # gpu_k <- NULL
    # gpu_beta_init <- NULL
    
    if (remainder > 0) {
      beta_extra <- res_beta_fit_extra$mu_beta
      theta_extra <- res_beta_fit_extra$theta
      beta <- rbind(beta, beta_extra)
      theta <- c(theta, theta_extra)
      beta_iters <- c(res_beta_fit$iter, res_beta_fit_extra$iter)
    } else {
      beta_iters <- c(res_beta_fit$iter)
    }
    
    if (is.null(dim(beta))) {
      beta <- matrix(beta, ncol = 1)
    }
    
    rownames(beta) <- gene_names
    
    # Create fit_res structure to match CPU branch
    fit_res <- list(
      beta = beta,
      theta = theta,
      iterations = list(
        beta_iters = beta_iters,
        theta_iters = 0L # GPU uses MOM, no iterative fitting
      )
    )
    
  } else {
    ## - CPU branch ----
    fit_res <- devil:::cpu_fit(
      input_matrix = input_matrix,
      design_matrix = design_matrix,
      offset_vector = offset_vector,
      init_overdispersion = init_overdispersion,
      init_beta_rough = init_beta_rough,
      overdispersion = overdispersion,
      n.cores = n.cores, max_iter = max_iter,
      tolerance = tolerance, verbose = verbose
    )
  }
  
  return(list(
    beta = fit_res$beta,
    overdispersion = fit_res$theta,
    iterations = fit_res$iterations,
    size_factors = sf,
    offset_vector = offset_vector,
    design_matrix = design_matrix,
    input_matrix = input_matrix,
    input_parameters = list(max_iter = max_iter, tolerance = tolerance, parallel.cores = n.cores)
  ))
}

my_test_de = function (devil.fit, contrast, pval_adjust_method = "BH", max_lfc = 10, 
                       clusters = NULL, parallel.cores = 1) {
  max.cores <- parallel::detectCores()
  if (is.null(parallel.cores)) {
    n.cores <- max.cores
  } else {
    if (parallel.cores > max.cores) {
      message("Requested ", parallel.cores, " cores, but only ", max.cores, " available.")
    }
    n.cores <- min(max.cores, parallel.cores)
  }
  
  ngenes   <- nrow(devil.fit$input_matrix)
  nsamples <- nrow(devil.fit$design_matrix)
  contrast <- as.array(contrast)
  
  lfcs <- (devil.fit$beta %*% contrast) %>% unlist() %>% unname() %>% c()
  
  if (!is.null(clusters) & !is.numeric(clusters)) {
    message("Converting clusters to numeric factors")
    clusters <- as.numeric(as.factor(clusters))
  }
  
  out_list <- parallel::mclapply(seq_len(nrow(devil.fit$input_matrix)), function(gene_idx) {
    mu_test <- lfcs[gene_idx]  # natural-log scale
    
    # Cluster-robust (if clusters provided)
    H <- devil:::compute_sandwich(
      devil.fit$design_matrix,
      devil.fit$input_matrix[gene_idx, ],
      devil.fit$beta[gene_idx, ],
      devil.fit$overdispersion[gene_idx],
      devil.fit$size_factors,
      clusters
    )
    var1 <- as.numeric(t(contrast) %*% H %*% contrast)
    se1  <- sqrt(var1)
    p1   <- 2 * stats::pt(abs(mu_test) / se1, df = nsamples - 2, lower.tail = FALSE)
    
    if (!is.null(clusters)) {
      # "Null" (non-clustered) as in your original code
      Hnull <- devil:::compute_sandwich(
        devil.fit$design_matrix,
        devil.fit$input_matrix[gene_idx, ],
        devil.fit$beta[gene_idx, ],
        devil.fit$overdispersion[gene_idx],
        devil.fit$size_factors,
        NULL
      )
      var0 <- as.numeric(t(contrast) %*% Hnull %*% contrast)
      se0  <- sqrt(var0)
      p0   <- 2 * stats::pt(abs(mu_test) / se0, df = nsamples - 2, lower.tail = FALSE)
      
      # Keep your conservative p-value rule, and report matching SE
      if (p0 >= p1) {
        p_final  <- p0
        se_final <- se0
      } else {
        p_final  <- p1
        se_final <- se1
      }
    } else {
      p_final  <- p1
      se_final <- se1
    }
    
    c(pval = p_final, se = se_final)
  }, mc.cores = n.cores)
  
  out_mat  <- do.call(rbind, out_list)
  p_values <- out_mat[, "pval"]
  se_mu    <- out_mat[, "se"]          # SE on natural-log scale
  se_lfc   <- se_mu / log(2)           # SE on log2 fold-change scale
  
  result_df <- dplyr::tibble(
    name     = rownames(devil.fit$beta),
    pval     = as.numeric(p_values),
    adj_pval = stats::p.adjust(as.numeric(p_values), method = pval_adjust_method),
    lfc      = lfcs / log(2),
    se_lfc   = as.numeric(se_lfc)
    # optionally also return se on natural-log scale:
    # se      = as.numeric(se_mu)
  )
  
  result_df <- result_df %>%
    dplyr::mutate(lfc = ifelse(.data$lfc >= max_lfc,  max_lfc,  .data$lfc)) %>%
    dplyr::mutate(lfc = ifelse(.data$lfc <= -max_lfc, -max_lfc, .data$lfc))
  
  df <- nsamples - 2
  tcrit <- stats::qt(0.975, df = df)
  result_df <- result_df %>%
    dplyr::mutate(
      ci_low  = .data$lfc - tcrit * .data$se_lfc,
      ci_high = .data$lfc + tcrit * .data$se_lfc
    )
  
  return(result_df)
}

estimate_mom_dispersion <- function(count_matrix,
                                    design_matrix,
                                    beta_matrix,
                                    sf) {
  
  G <- nrow(count_matrix)   # genes
  n <- ncol(count_matrix)   # cells
  n_design <- nrow(design_matrix)
  p <- ncol(design_matrix)
  
  if (n_design != n)
    stop("design_matrix must have nrow == ncol(count_matrix)")
  
  if (nrow(beta_matrix) != G)
    stop("beta_matrix must have nrow == nrow(count_matrix)")
  
  if (ncol(beta_matrix) != p)
    stop("beta_matrix must have ncol == ncol(design_matrix)")
  
  if (length(sf) != n)
    stop("sf must have length ncol(count_matrix)")
  
  corr <- n / (n - p)
  
  theta <- numeric(G)
  
  for (g in seq_len(G)) {
    num <- 0
    den <- 0
    for (j in seq_len(n)) {
      
      # linear predictor
      eta <- sum(design_matrix[j, ] * beta_matrix[g, ])
      
      mu   <- sf[j] * exp(eta)
      y    <- count_matrix[g, j]
      diff <- y - mu
      
      num <- num + (diff^2 - mu)
      
      if (is.na(num)) {
        print(j)
        stop()
      }
        
      den <- den + mu^2
    }
    
    th <- 0
    if (den > 0) {
      th <- corr * num / den
      if (th < 0) th <- 0  # truncate to non-negative
    }
    
    theta[g] <- th
  }
  
  return(theta)
}
