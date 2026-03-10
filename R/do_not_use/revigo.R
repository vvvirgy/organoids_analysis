library(httr)
library(stringi)

# args = commandArgs(trailingOnly=TRUE)

# Read user data from a file
# fileName <- 'data/compensation/test_enrichment.tsv'
# userData <- readChar(fileName,file.info(fileName)$size)

# run revigo dimension reduction

data = lapply(enrichment_groups, function(x) {
  x[['GO']]@result %>% 
    dplyr::select(ID, p.adjust) 
  }) 

run_revigo = function(go_data, name) {
  
  userData = go_data %>% 
    mutate(data = paste(ID, p.adjust, sep = '\t')) %>% 
    pull(data) %>% 
    paste(., collapse = '\n') %>% 
    paste0('ID\tp.adjust\n', .)
  
  # Default output file name is result.csv
  fileNameOutput <- paste0("data/compensation/revigo/", name, "_result.html")
  
  # Submit job to Revigo
  httr::POST(
    url = "http://revigo.irb.hr/Revigo",
    body = list(
      cutoff = "0.7",
      valueType = "pvalue",
      speciesTaxon = "9606",
      measure = "SIMREL",
      goList = userData
    ),
    # application/x-www-form-urlencoded
    encode = "form"
  ) -> res
  
  dat <- httr::content(res, as="text", encoding="UTF-8")
  # 
  # # Write results to a file
  # dat <- stri_replace_all_fixed(dat, "\r", "")
  # cat(dat, file=fileNameOutput, fill = FALSE)
  
  txt <- httr::content(res, as = "text", encoding = "UTF-8")
  
  # Extract table block
  table_text <- stringr::str_extract(
    txt,
    "Term ID.*?(?=\\n\\s*\\n)"
  )
  
  df <- readr::read_tsv(table_text, show_col_types = FALSE)
  
  df$group <- name
  
  return(df)
}

revigo_res = lapply((data %>% names)[1:16], function(x) {
  
  gene_group = gsub(':', '_', x) %>% 
    gsub('\\.', '', .) %>% 
    gsub('\\(', '_', .) %>% 
    gsub('\\)', '', .) 
  
  run_revigo(data[[x]], name = gene_group)
})




