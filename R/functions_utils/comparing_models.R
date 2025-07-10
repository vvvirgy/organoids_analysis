compare_models = function(model1, model2) {
  
  cmv_lambda_min_1 = model1$cvm[model1$lambda == model1$lambda.min]
  cmv_lambda_min_2 = model2$cvm[model2$lambda == model2$lambda.min]
  
  cmv_lambda_min_1 - cmv_lambda_min_2
  
}

compute_improvement = function(model1, model2) {
  cmv_lambda_min_1 = model1$cvm[model1$lambda == model1$lambda.min]
  cmv_lambda_min_2 = model2$cvm[model2$lambda == model2$lambda.min]
  
  ((cmv_lambda_min_1-cmv_lambda_min_2)/cmv_lambda_min_1)*100
}

get_coefficients_multiresponse = function(x) {
  lapply(x$coefficients %>% names, function(y) {
    x$coefficients[[y]] %>% 
      as.matrix %>% 
      as.data.frame() %>% 
      tibble::rownames_to_column('predictor') %>% 
      dplyr::rename(coefficient = '1') %>% 
      dplyr::mutate(response = y)
  }) %>% 
    dplyr::bind_rows() %>% 
    dplyr::filter(predictor != '(Intercept)') %>% 
    tidyr::pivot_wider(names_from = c(response, predictor), values_from = coefficient)
}
