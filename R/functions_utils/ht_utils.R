# create annotation for the heatmap
create_annotation = function(x, # data 
                             ann_colors, 
                             position) {
  
  ComplexHeatmap::HeatmapAnnotation(df = x, 
                                    col = ann_colors, 
                                    which = position, 
                                    annotation_legend_param = list(nrow = 3, width = 12, by_row = T))
  
}
