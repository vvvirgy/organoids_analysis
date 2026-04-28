survival_analysis = function(data, time, event, what, gene) {
  
  surv_obj <- Surv(time = (data %>% dplyr::pull(time)), event = (data %>% dplyr::pull(event)) )
  
  # Create formula dynamically
  # formula_obj <- as.formula(paste("surv_obj ~", what, '+ Inferred_menopause_status'))
  # formula_obj2 <- as.formula(paste0("Surv(", time, ",", event, ") ~ ", what, '+ Inferred_menopause_status'))
  formula_obj <- as.formula(paste("surv_obj ~", what))
  formula_obj2 <- as.formula(paste0("Surv(", time, ",", event, ") ~ ", what))
  
  # Fit
  fit_km <- survfit(formula_obj, data = data)
  diff <- survdiff(formula_obj, data = data)
  
  # perform also cox fit to have the comparison among each category 
  # data_coxph <- data %>%
  #   mutate(surv_object = Surv(  time  = .data[[time]],
  #                               event = .data[[event]]))
  
  data_df = as.data.frame(data)
  
  fit_coxph = coxph(formula = formula_obj2, ties = 'exact', data = data_df)
  
  hr_plot = tryCatch(expr = {
    ggforest(fit_coxph, data = data_df)}, 
    error = function(e) {
      ggplot()
    })
  
  
  # summary(fit_coxph)
  # coxfit <- coxph(
  #   Surv(os_years, os_status) ~ what,
  #   data = data,
  #   ties = 'exact')
  
  
  # plot
  pp = survfit2(formula_obj2, data = data) %>% 
    ggsurvfit(linewidth = 1.5, ) +
    # add_risktable() %>% 
    labs( x = 'Time (years)', y = 'OS probability')  + 
    # add_confidence_interval() +
    scale_color_brewer(palette = "Set1") +
    add_quantile() + 
    # scale_color_manual(values = wes_palette(3, name = "Darjeeling1", type = "discrete")) +
    # scale_fill_manual(values = wes_palette(3, name = "Darjeeling1", type = "discrete")) + 
    add_pvalue(location = 'annotation') +
    # ggtitle(what) + 
    theme(axis.text = element_text(family = "Arial", size = 13),
          axis.title = element_text(family = "Arial", size = 13),
          legend.text = element_text(family = "Arial", size = 12),
          panel.border = element_rect(linewidth = 1.5), 
          title = element_text(size = 15), 
          legend.direction = 'vertical') + 
    guides(color = guide_legend(ncol = 2))
  
  list('fit' = fit_km, 
       'Diff' = diff, 
       'coxph' = fit_coxph,
       'hr_plot' = hr_plot,
       'km_plot' = pp)
}
