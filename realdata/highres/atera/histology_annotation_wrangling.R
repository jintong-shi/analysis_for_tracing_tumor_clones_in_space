
setwd('<DATA_ROOT>/high_res_ST/atera/data/histology')
csv_list <- lapply(list.files(), function(file){
  region <- paste0(strsplit(file,"_")[[1]][1],"_",strsplit(file,"_")[[1]][2])
  csv <- read.csv(file)
  csv$histology <- region
  return(csv)
})%>%bind_rows()
table(csv_list$histology)
write.csv(csv_list,'<DATA_ROOT>/high_res_ST/atera/data/histology/histology_annotation.csv')
