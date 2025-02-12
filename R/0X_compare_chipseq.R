## TODO
## -----------------------------------------------------------------------------

library(tidyverse)
library(data.table)
source("R/00_config.R")
source("R/utils/functions.R")
source("R/utils/vector_comparison_functions.R")

k <- 200  # TopK overlap
set.seed(5)

# Dataset meta and data IDs
mcg_meta <- read.delim(mcg_meta_dedup_path)
ids_hg <- unique(filter(mcg_meta, Species == "Human")$ID)
ids_mm <- unique(filter(mcg_meta, Species == "Mouse")$ID)

# Relevant gene tables
pc_hg <- read.delim(ref_hg_path)
tfs_hg <- read.delim(tfs_hg_path)
pc_mm <- read.delim(ref_mm_path)
tfs_mm <- read.delim(tfs_mm_path)

# List of measurement info to keep filtered genes
count_summ <- readRDS(mcg_count_summ_list_path)

# The final aggregate/global profile used for comparison
agg_mm <- readRDS("/space/scratch/amorin/aggregate_microglia/Cormats/Mm_pcor/aggregate_cormat_FZ_mm.RDS")
agg_hg <- readRDS("/space/scratch/amorin/aggregate_microglia/Cormats/Hg_pcor/aggregate_cormat_FZ_hg.RDS")


# Average bind scores and bind matrix
bind_dat_path <- "/space/scratch/amorin/R_objects/processed_unibind_data.RDS"
bind_summary_path <- "/space/scratch/amorin/R_objects/unibind_Permissive_bindscore_summary.RDS"
bind_summary <- readRDS(bind_summary_path)
bind_dat <- readRDS(bind_dat_path)


# Keep genes measured in microglia
msr_hg <- count_summ$Human$Summ_df
msr_mm <- count_summ$Mouse$Summ_df

keep_hg <- count_summ$Human$Filter_genes
keep_mm <- count_summ$Mouse$Filter_genes


gene_hg <- "SPI1"
gene_mm <- "Spi1"


# Human SPI1 no microglia data
all_bind_ids_hg <- bind_dat$Permissive_hg$Meta %>% 
  filter(Symbol == gene_hg) %>% 
  pull(ID)


mcg_bind_ids_mm <- bind_dat$Permissive_mm$Meta %>% 
  filter(str_to_title(Symbol) == gene_mm & str_detect(ID, ".*microglia.*")) %>% 
  pull(ID)


all_bind_ids_mm <- bind_dat$Permissive_mm$Meta %>% 
  filter(str_to_title(Symbol) == gene_mm) %>% 
  pull(ID)


non_mcg_bind_ids_mm <- setdiff(all_bind_ids_mm, mcg_bind_ids_mm)

# TODO: show agreement of SPI1 coexpr aggregation between:
# 1) SPI1 microglia bind profiles
# 2) SPI1 all bind profiles (all + shuffle same n)


all_bind_hg <- bind_dat$Permissive_hg$Mat_qnl[, all_bind_ids_hg]
all_bind_mm <- bind_dat$Permissive_mm$Mat_qnl[, all_bind_ids_mm]
mcg_bind_mm <- bind_dat$Permissive_mm$Mat_qnl[, mcg_bind_ids_mm]
nonmcg_bind_mm <- bind_dat$Permissive_mm$Mat_qnl[, non_mcg_bind_ids_mm]


identical(pc_hg$Symbol, rownames(agg_hg$Agg_mat))
identical(pc_hg$Symbol, rownames(bind_summary$Human_TF))

identical(pc_mm$Symbol, rownames(agg_mm$Agg_mat))
identical(pc_mm$Symbol, rownames(bind_summary$Mouse_TF))

identical(rowMeans(all_bind_hg), bind_summary$Human_TF[, gene_hg])


bind_df_hg <- data.frame(
  Symbol = pc_hg$Symbol,
  Coexpr = agg_hg$Agg_mat[, gene_hg],
  Bind_all = bind_summary$Human_TF[, gene_hg],
  N_msr = msr_hg$N_msr
)

bind_df_hg <- filter(bind_df_hg, Symbol != gene_hg)



bind_df_mm <- data.frame(
  Symbol = pc_mm$Symbol,
  Coexpr = agg_mm$Agg_mat[, gene_mm],
  Bind_all = bind_summary$Mouse_TF[, gene_mm],
  Bind_mcg = rowMeans(mcg_bind_mm),
  Bind_nonmcg = rowMeans(nonmcg_bind_mm),
  N_msr = msr_mm$N_msr
)
bind_df_mm <- filter(bind_df_mm, Symbol != gene_mm)


# TODO: inspecct high bind low/no msr
bind_df_hg <- filter(bind_df_hg, Symbol %in% keep_hg)
bind_df_mm <- filter(bind_df_mm, Symbol %in% keep_mm)



plot(bind_df_hg$Coexpr, bind_df_hg$Bind_all)
plot(bind_df_mm$Coexpr, bind_df_mm$Bind_all)
plot(bind_df_mm$Coexpr, bind_df_mm$Bind_mcg)
plot(bind_df_mm$Coexpr, bind_df_mm$Bind_nonmcg)

cor(select_if(bind_df_hg, is.numeric), method = "spearman", use = "pairwise.complete.obs")
cor(select_if(bind_df_mm, is.numeric), method = "spearman", use = "pairwise.complete.obs")

colwise_topk_intersect(select_if(bind_df_hg, is.numeric), k = 200, check_k_arg = FALSE)
colwise_topk_intersect(select_if(bind_df_mm, is.numeric), k = 200, check_k_arg = FALSE)


# data.frame(
#   Coexpr = c(bind_df_mm$Coexpr, bind_df_mm$Coexpr),
#   Bind = c(bind_df_mm$Bind_mcg, bind_df_mm$Bind_nonmcg),
#   Group = c(rep("Microglia", nrow(bind_df_mm)), rep("Non-microglia", nrow(bind_df_mm)))
# )


p1a <- ggplot(bind_df_mm, aes(x = Bind_mcg, y = Coexpr)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic()

p1b <- ggplot(bind_df_mm, aes(x = Bind_nonmcg, y = Coexpr)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic()

egg::ggarrange(p1a, p1b, nrow = 1)


p2 <- ggplot(bind_df_mm, aes(x = Bind_all, y = Coexpr)) +
  geom_point() +
  geom_smooth(aes(x = Bind_mcg, y = Coexpr), method = "lm", colour = "red") +
  geom_smooth(aes(x = Bind_nonmcg, y = Coexpr), method = "lm", colour = "blue") +
  theme_classic()


fit1 <- lm(Coexpr ~ Bind_mcg + Bind_nonmcg, data = bind_df_mm)
fit2 <- lm(Coexpr ~ Bind_mcg * Bind_nonmcg, data = bind_df_mm)
anova(fit1, fit2)


topk_mat <- cbind(bind_summary$Mouse_TF[keep_mm, ], Coexpr = agg_mm$Agg_mat[keep_mm, gene_mm])
topk_mat <- colwise_topk_intersect(topk_mat, k = k)
topk_df <- mat_to_df(topk_mat, value = "Topk")


topk_mat <- cbind(all_bind_mm[keep_mm, ], Coexpr = agg_mm$Agg_mat[keep_mm, gene_mm])
topk_mat <- colwise_topk_intersect(topk_mat, k = k)
topk_df <- mat_to_df(topk_mat, value = "Topk")


bind_mat_mm <- bind_dat$Permissive_mm$Mat_qnl[keep_mm, ]
# bind_mat_mm <- bind_dat$Permissive_mm$Mat_qnl
agg_vec <- topk_sort(agg_mm$Agg_mat[keep_mm, "Spi1"], k = k)


topk_l <- mclapply(1:ncol(bind_mat_mm), function(x) {
  topk_intersect(agg_vec, topk_sort(bind_mat_mm[, x], k = k))
}, mc.cores = ncore)


stopifnot(identical(colnames(bind_mat_mm), bind_dat$Permissive_mm$Meta$ID))


topk_df <- data.frame(
  ID = bind_dat$Permissive_mm$Meta$ID,
  Topk = unlist(topk_l),
  Symbol = bind_dat$Permissive_mm$Meta$Symbol,
  Mcg = bind_dat$Permissive_mm$Meta$ID %in% mcg_bind_ids_mm
)

topk_df <- topk_df %>%
  mutate(
    Group = case_when(
      Symbol == "SPI1" & Mcg ~ "SPI mcg",
      Symbol == "SPI1" & !Mcg ~ "SPI Non-mcg",
      Symbol != "SPI1" ~ "Non-SPI1"
    ),
    Group = factor(Group, levels = c("Non-SPI1", "SPI Non-mcg", "SPI mcg"))
  )




ggplot(topk_df, aes(x = Group, y = Topk)) +
  geom_violin(fill = "slategrey") +
  # geom_violin(fill = "slategrey", alpha = 0.1) +
  geom_jitter(data = filter(topk_df, Group != "Non-SPI1"),
              shape = 21, size = 2.4, width = 0.1) +
  geom_boxplot(width = 0.05) +
  
  ylab("Top200") +
  xlab(NULL) +
  theme_classic() +
  theme(text = element_text(size = 25))



     