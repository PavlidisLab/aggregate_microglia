## TODO
## -----------------------------------------------------------------------------

library(tidyverse)
library(aggtools)
source("R/00_config.R")
source("R/utils/functions.R")


mcg_meta <- read.delim(mcg_meta_path)
ids_hg <- unique(filter(mcg_meta, Species == "Human")$ID)
ids_mm <- unique(filter(mcg_meta, Species == "Mouse")$ID)


dat_l <- readRDS(mcg_dat_path)


avg_hg <- do.call(cbind, lapply(dat_l[ids_hg], function(x) rowMeans(x$Mat)))
avg_mm <- do.call(cbind, lapply(dat_l[ids_mm], function(x) rowMeans(x$Mat)))



rowMeans(dat_l$GSE180928$Mat)
