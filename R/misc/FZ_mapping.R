library(tidyverse)
library(cowplot)

pcor_raw <- readRDS("~/sc_output/Microglia/Mm_pcor/Spi1_rawcor.tsv")
pcor_fz <- DescTools::FisherZ(pcor_raw)
pcor_fz_inv <- DescTools::FisherZInv(pcor_fz)

max_pcor <- max(pcor_raw[setdiff(rownames(pcor_raw), "Spi1"), ], na.rm = TRUE)
min_pcor <- min(pcor_raw[setdiff(rownames(pcor_raw), "Spi1"), ], na.rm = TRUE)


# seq_pcor <- seq(-1, 1, length.out = 1000)
# seq_fz <- DescTools::FisherZ(seq(-1, 1, length.out = 1000))
# seq_fzinv <- DescTools::FisherZInv(seq_fz)


plot_df1 <- data.frame(
  Pcor = seq(-1, 1, length.out = 1000),
  FZ = DescTools::FisherZ(seq(-1, 1, length.out = 1000))
)


p1 <- 
  ggplot(plot_df1, aes(x = Pcor, y = FZ)) +
  geom_point(shape = 21) +
  geom_vline(xintercept = max_pcor, col = "red") +
  geom_vline(xintercept = min_pcor, col = "blue") +
  xlab("Pearson's cor") +
  ylab("Fisher's Z") +
  theme_classic() +
  theme(text = element_text(size = 25))



p2a <- 
  data.frame(
  X = rowMeans(pcor_raw, na.rm = TRUE),
  Y = rowMeans(pcor_fz, na.rm = TRUE)) %>% 
  ggplot(., aes(x = X, y = Y)) +
  geom_point() +
  xlim(c(-0.2, 0.2)) +
  ylim(c(-0.2, 0.2)) +
  xlab("Average Pearson's cor") +
  ylab("Average Fisher's Z") +
  ggtitle("Spi1 aggregate") +
  theme_classic() +
  theme(text = element_text(size = 25))




p2b <- 
  data.frame(
  X = rowMeans(pcor_raw, na.rm = TRUE),
  Y = rowMeans(pcor_fz_inv, na.rm = TRUE)) %>% 
  ggplot(., aes(x = X, y = Y)) +
  geom_point() +
  xlim(c(-0.2, 0.2)) +
  ylim(c(-0.2, 0.2)) +
  xlab("Average Pearson's cor") +
  ylab("Average Inverse Fisher's Z") +
  ggtitle("Spi1 aggregate") +
  theme_classic() +
  theme(text = element_text(size = 25))




p2c <- 
  data.frame(
  X = rowMeans(pcor_fz, na.rm = TRUE),
  Y = rowMeans(pcor_fz_inv, na.rm = TRUE)) %>% 
  ggplot(., aes(x = X, y = Y)) +
  geom_point() +
  xlim(c(-0.2, 0.2)) +
  ylim(c(-0.2, 0.2)) +
  xlab("Average Fisher's Z") +
  ylab("Average Inverse Fisher's Z") +
  ggtitle("Spi1 aggregate") +
  theme_classic() +
  theme(text = element_text(size = 25))


plot_grid(p2a, p2b, p2c, nrow = 1)
