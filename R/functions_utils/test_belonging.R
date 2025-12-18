test_belonging_alterations = function(x, 
                                      pth = .05, 
                                      n_bootstrap = 10000) {
  x %>% 
    group_by(hgnc_symbol) %>% 
    group_modify(~ {
      
      ref <- .x %>% filter(alteration_classes == "No alteration") %>% pull(observed_expr)
      alt <- .x %>% filter(alteration_classes == "Alteration")
      
      print(ref)
      
      if (length(ref) <= 1 || nrow(alt) == 0) {
        return(tibble())
      }
      
      if (n_bootstrap > 0) {
        
        bootstrap_data <- (mosaic::do(n_bootstrap) * mosaic::resample(ref))
        
        sigma <- apply(bootstrap_data, 1, sd) %>% mean
        mu <- apply(bootstrap_data, 1, mean) %>% mean
        
        ci_int <- abs(pth/2 + c(0, -1))
        
        boot_means <- apply(bootstrap_data, 1, mean)
        
        # ---- FIX: detect failed bootstrap ----
        if (all(is.na(boot_means))) {
          message("⚠️ Bootstrap failed for gene: ", unique(.x$hgnc_symbol))
          return(tibble())
        }
        
        ci <- quantile(boot_means, probs = ci_int, na.rm = TRUE)
        
      } else {
        mu <- mean(ref, na.rm = TRUE)
        sigma <- sd(ref, na.rm = TRUE)
        ci <- NULL
      }
      
      alt_test = alt %>%
        mutate(
          sd_not_alt = sigma,
          mu_not_alt = mu,
          z = (observed_expr - mu) / sigma,
          p_value = 2 * pnorm(-abs(z)),
          significant = p_value <= pth,
          CI_lower = if (!is.null(ci)) min(ci) else NA,
          CI_upper = if (!is.null(ci)) max(ci) else NA,
          CI = pth
        ) 
      
      mut_test = alt_test %>% 
        # dplyr::filter(IMPACT %in% c('HIGH', 'MODERATE')) %>% 
        dplyr::filter(score >= 0.446 | IMPACT %in% c('HIGH', 'MODERATE')) %>% 
        mutate(
          sd_not_alt = sigma,
          mu_not_alt = mu,
          z_h1 = (observed_expr_h1 - mu) / sigma,
          p_value_h1 = 2 * pnorm(-abs(z_h1)),
          significant_h1 = p_value_h1 <= pth,
          CI_lower = if (!is.null(ci)) min(ci) else NA,
          CI_upper = if (!is.null(ci)) max(ci) else NA,
          CI = pth
        )
        
      full_join(alt_test, mut_test)
      
    })
}
