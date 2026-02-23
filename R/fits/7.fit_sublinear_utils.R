# Sub-linear model for VST-normalised expression vs copy number ############################################################
# Model:
#   y_vst = β0 + h * log(K / K_ref) + ε,   ε ~ N(0, σ^2)

# fit sublinear -----------------------------------------------------------

fit_sublinear <- function(data, gene, K_ref = 2) {
  
  data = data %>% 
    dplyr::filter(hgnc_symbol == gene) %>% 
    filter(mutation_status == 'Wild-type') %>% 
    mutate(scaled_cna = log(tot_cna / K_ref))

  fit <- lm(expression ~ scaled_cna, data = data)
  
  out <- list(
    lm     = fit,
    K_ref  = K_ref
  )

  return(out)
}

# prediction sublinear -----------------------------------------------------------
predict_sublinear <- function(fit, K_new, K_ref = 2, interval = c("none","confidence")) {
  
  x_new <- log(K_new / K_ref)
  
  # perform the prediction
  preds <- predict(fit, newdata = data.frame(x = x_new),
                   interval = "confidence")
  as_tibble(preds) %>% 
    mutate(K = K_new)
}

# hypothesis testing linear scaling (h = 1) vs sub-linear (h ≠ 1) -----------------------------------------------------------

test_sublinearity <- function(data, gene, K_ref = 2) {
  
  print(gene)
  
  data = data %>% 
    dplyr::filter(hgnc_symbol == gene) %>% 
    filter(mutation_status == 'Wild-type') %>% 
    mutate(scaled_cna = log(tot_cna / K_ref))
  
  if(nrow(data) > 2){
    
    # Null model: slope fixed to 1 → y = beta0 + 1 * x + ε
    fit_null <- lm(expression ~ 1 + offset(1 * scaled_cna), data = data)
    
    # Alternative model: slope free → y = beta0 + h * x + ε
    fit_alt  <- lm(expression ~ scaled_cna, data = data)
    
    an <- anova(fit_null, fit_alt)
    rownames(an)[2] <- "H1: h != 1"
    
    list(
      'fit_null' = fit_null, 
      'fit_alt' = fit_alt, 
      'anova' = an, 
      'gene' = gene
    )}
  
}

# plotting -----------------------------------------------------------
plot_vst_sublinear <- function(object, y_vst, K) {
  x <- log(K / object$K_ref)
  par(mfrow = c(1,2))
  
  # Residuals vs Fitted
  plot(fitted(object$lm), resid(object$lm),
       xlab = "Fitted (VST)",
       ylab = "Residuals",
       pch = 16, col = "grey40",
       main = "Residuals vs Fitted")
  abline(h = 0, col = "red", lwd = 2)
  
  # QQ-plot of residuals
  qqnorm(resid(object$lm), pch = 16, col = "grey40")
  qqline(resid(object$lm), col = "red", lwd = 2)
  
  par(mfrow = c(1,1))
}


