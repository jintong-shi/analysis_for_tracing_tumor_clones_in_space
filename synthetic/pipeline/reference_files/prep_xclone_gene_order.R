# Annotate gene_order.csv with cytoband / chromosome arm -> gene_order_xclone_format.csv
library(dplyr)
setwd('<REF_ROOT>')

cyto <- read.table("cytoBand.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(cyto) <- c("chr","start","end","band","stain")

cyto <- cyto %>%
  mutate(arm = substr(band, 1, 1))  # p or q

gene_order <- read.csv("gene_order.csv", stringsAsFactors = FALSE) %>%
  mutate(mid = (start + end)/2)

cyto <- cyto %>% arrange(chr, start)
cyto_split <- split(cyto, cyto$chr)

get_arm_one_chr <- function(df_gene, df_cyto){
  idx <- findInterval(df_gene$mid, df_cyto$start)
  idx[idx < 1] <- NA
  idx[idx > nrow(df_cyto)] <- NA

  # findInterval only matches on start, so check the gene is inside the band
  ok <- !is.na(idx) & df_gene$mid <= df_cyto$end[idx]
  idx[!ok] <- NA

  df_gene$band <- ifelse(is.na(idx), NA, df_cyto$band[idx])
  df_gene$arm  <- ifelse(is.na(idx), NA, df_cyto$arm[idx])
  df_gene
}

out <- bind_rows(lapply(split(gene_order, gene_order$chr), function(g){
  c <- cyto_split[[ unique(g$chr) ]]
  if (is.null(c)) { g$band <- NA; g$arm <- NA; return(g) }
  get_arm_one_chr(g, c)
}))
out <- out %>% dplyr::mutate(chr_arm = paste0(chr))

out %>% select(gene, chr, start, end, band, arm)

final_df <- out %>%
  transmute(
    GeneName = gene,
    chr      = sub("^chr", "", chr),    # chr1 -> 1
    start    = start,
    stop     = end,
    arm      = arm,
    chr_arm  = paste0(sub("^chr", "", chr), arm),
    band     = band
  )
write.csv(final_df,'gene_order_xclone_format.csv',quote = FALSE,row.names = FALSE)