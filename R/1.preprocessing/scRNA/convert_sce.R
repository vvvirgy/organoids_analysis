library(zellkonverter)

d_path = '/orfeo/scratch/cdslab/vgazziero/organoids_prj/data/scRNA/reanalysis/splitted_h5/'
d_path = list.files(d_path, full.names = T)

lapply(d_path, function(x) {
  
  dest_path = gsub('splitted_h5', 'sce_new', x)
  dest_path = gsub('h5ad', 'rds', dest_path)
  
  sce = readH5AD(x)
  saveRDS(sce, dest_path)
})
