compute_errors = function(glm_fit, 
                       data, 
                       n_obs = 5,
                       response, 
                       model = as.formula('~ tot_cna + mutation_status'), 
                       lambda = "lambda.min") {
  
  lapply(names(glm_fit), function(x) {
    
    df = data %>% 
      dplyr::filter(hgnc_symbol == x) %>% 
      drop_na(all_of(response))
    print(x)
    
    # remove those genes with too little observations
    
    if(nrow(df) < n_obs) {
      print(paste0('Less than ', n_obs, ' observations, discarding gene')) 
    } else {
      
      mod_mat = model.matrix(model, df)[,-1, drop = FALSE]
      response_variable = df %>%
        dplyr::select(all_of(response)) %>% 
        drop_na() %>% 
        as.matrix()
      
      # predict using the model results
      prediction = predict(glm_fit[[x]]$cv_glm_best_alpha, newx = mod_mat, s = lambda)[,,1] 
      rmse_total <- sqrt(mean((response_variable - prediction)^2))
      rmse_per_variable = sqrt(colMeans((response_variable - prediction)^2)) %>% 
        as_tibble_row()
      colnames(rmse_per_variable) = paste0(colnames(rmse_per_variable),'_rmse')
      r2 <- 1 - colSums((response_variable - prediction)^2) / colSums((response_variable - colMeans(response_variable))^2) %>% 
        as_tibble_row()
      colnames(r2) = paste0(colnames(r2), '_R2')
      
      res = tibble(
        gene = x
      ) %>% 
        bind_cols(rmse_per_variable) %>% 
        bind_cols(r2) %>% 
        mutate(
          rmse_total = rmse_total
        )
      
      return(res)
    }
  })
}


get_r2 = function(glm_fit) {
  
  ffit = glm_fit$cv_glm_best_alpha
  
  r2 <- ffit$glmnet.fit$dev.ratio[which(ffit$glmnet.fit$lambda == ffit$lambda.min)]
  return(r2)
}
