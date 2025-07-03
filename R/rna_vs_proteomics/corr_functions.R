# correlation functions

cor_data.fun<-function(x){
  pvalue<-x$p.value
  cor_value<-x$estimate
  tstat<-x$statistic
  myvector<-c(cor_value,tstat,pvalue)
  return(myvector)
}

Genescorr <- function(subtype_df, 
                      method = 'spearman') {
  
  corrs_gene = lapply(subtype_df$hgnc_symbol, function(x) {
    prot = subtype_df %>% 
      filter(hgnc_symbol == x) %>% 
      pull(mean_intensity)
    
    rna = subtype_df %>% 
      filter(hgnc_symbol == x) %>% 
      pull(value)
    
    corr = cor.test(as.numeric(prot), 
             as.numeric(rna), 
             method = method)
    return(corr)
  })

  corrs_gene.df <- t(sapply(corrs_gene, function(x){
    cor_data.fun(x)
  }))
  colnames(corrs_gene.df) <- c("correlation", "t_stat", "pValue")
  
  corrs.qvalue <- p.adjust(corrs_gene.df[, "pValue"], method = "BH")
  corrs_gene.df <- cbind(corrs_gene.df, qValue = corrs.qvalue)
  
  corrs_gene.df <- as.data.frame(corrs_gene.df)
  return(corrs_gene.df)
}
