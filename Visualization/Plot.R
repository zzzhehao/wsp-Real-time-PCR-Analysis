library(ggpubr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("AES.R")
source("Pool seq infection strain distr.R")
source("Pool seq infection level.R")
source("Tetracycline.R")
source("Reproduction.R")
source("Offspring infection status.R")

ggsave(filename = "f55.png", plot = f5.5, width = 8, height = 5)
ggsave(filename = "f56.png", plot = f5.6, width = 5, height = 5)

f5.7 <- ggarrange(f5.7A.bar.infection_categories, f5.7B.quantity_group.quantified, nrow = 1, align = "h", labels = c("A", "B"))
ggsave(filename = "f57.png", plot = f5.7, width = 11, height = 6)

f5.8 <- ggarrange(f.58A, f5.8b.f1in, f5.8C.dsxp, nrow = 1, labels = c("A", "B", "C"), widths = c(2.5,2,2))
ggsave(filename = "f58.png", plot = f5.8, width = 11, height = 6)
