glm_fit = function(data, 
                   response, 
                   model = as.formula('~ tot_cna + mutation_status'), 
                   alphas = seq(0, 1, 0.1)) {
  
  lapply(unique(data$hgnc_symbol), function(x) {
    
    df = data %>% 
      dplyr::filter(hgnc_symbol == x)
    print(x)
    
    # remove those genes with too little observations
    
    if(nrow(df) > 20) {
      
      mod_mat = model.matrix(model, df)[,-1, drop = FALSE]
      response_variable = df %>%
        pull(response)
      
      # fit cross validation and test different values of alpha
      # alphas = alpha_list
      
      cv = lapply(alphas, function(a) {
        print(a)
        
        tryCatch( { glmnet::cv.glmnet(x = mod_mat, 
                                      y = response_variable, 
                                      # nfolds = 4, 
                                      alpha = as.numeric(a))  }
                  , error = function(e) {NA})
        
      })
      
      names(cv) = alphas
      
      cv = Filter(function(x) !all(is.na(x)), cv)
      if(length(cv) > 0) {
        cvm = sapply(cv, function(x) {
          min(x[['cvm']])
        })
        cvm_best = min(cvm)
        # extract the best alpha
        best_alpha = names(cvm[which(cvm == cvm_best)])
        
        coeffs = coef(cv[[best_alpha]], s = 'lambda.1se')
        
        return(list(
          cv_glm_best_alpha = cv[[best_alpha]], 
          alpha = best_alpha, 
          coefficients = coeffs, 
          gene = x
        ))}
    }
  })
}
