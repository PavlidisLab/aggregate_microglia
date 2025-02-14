# Microglia ccRE-gene correlations
# Human paper https://www.science.org/doi/10.1126/science.adf7044
# Mouse paper https://www.nature.com/articles/s41586-023-06824-9
# Download home http://catlas.org/catlas_hub/

# Note: human has file with stats (PCC and Pval) while mouse just has positions

# http://catlas.org/catlas_downloads/humanbrain/conns/MGC.pos.stat.txt
# http://catlas.org/catlas_downloads/humanbrain/conns/MGC.pos.bedpe
# http://catlas.org/catlas_downloads/mousebrain/conns/MGL.pos.conns.bedpe
# ------------------------------------------------------------------------------

library(tidyverse)
source("R/00_config.R")


url_hg <- "http://catlas.org/catlas_downloads/humanbrain/conns/MGC.pos.bedpe"
url_mm <- "http://catlas.org/catlas_downloads/mousebrain/conns/MGL.pos.conns.bedpe"


# Download, format, and save microglia ccre - promoter - expression correlation


download_and_format_ccre <- function(path, url) {
  
  if (!file.exists(path)) {
    
    cols <- c(
      "Chromosome_proximal",
      "Start_proximal",	
      "End_proximal",	
      "Chromosome_distal",	
      "Start_distal",	
      "End_distal",
      "posPair"
    )
    
    download.file(url, path)
    
    ccre <- read.delim(path, header = FALSE)
    ccre <- ccre[, 1:7]
    colnames(ccre) <- cols
    
    
    # Formatting chromosome name for consistency with gene table
    ccre <- ccre %>% 
      mutate(
        Chromosome_proximal = str_replace(Chromosome_proximal, "^chr", ""),
        Chromosome_distal = str_replace(Chromosome_distal, "^chr", "")
        )
    
    
    # Seperate symbol from cCRE ID. Human and mouse have flipped columns
    gene_ccre <- str_split(ccre[["posPair"]], pattern = "\\|", simplify = TRUE)
    gene_ccre <- data.frame(gene_ccre)
    colnames(gene_ccre) <- c("Symbol", "cCRE_ID")
    
    if (str_detect(gene_ccre[1, 1], "^cCRE.*")) {
      colnames(gene_ccre) <- c("cCRE_ID", "Symbol")

    }
    
    ccre <- data.frame(Symbol = gene_ccre[["Symbol"]],
                       select(ccre, -posPair),
                       cCRE_ID = gene_ccre["cCRE_ID"])
    
    write.table(ccre, sep = "\t", quote = FALSE, file = path)
    
    return(invisible(NULL))
  }
  
}


download_and_format_ccre(mcg_ccre_path_hg, url_hg)
download_and_format_ccre(mcg_ccre_path_mm, url_mm)
