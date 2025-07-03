compare_models = function(model1, model2) {
  
  cmv_lambda_min_1 = model1$cvm[model1$lambda == model1$lambda.min]
  cmv_lambda_min_2 = model2$cvm[model2$lambda == model2$lambda.min]
  
  cmv_lambda_min_1 - cmv_lambda_min_2
  
}