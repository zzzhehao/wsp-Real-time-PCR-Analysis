library(readxl)
library(ggpubr)
library(dplyr) 
library(stringr)
library(rstatix)
library(tidyverse)
library(broom)
library(ggplot2)
library(ggsignif)

# Import data
raw <- read_xlsx("/Users/hu_zhehao/Library/Mobile Documents/com~apple~CloudDocs/UHH/B.Sc. Biologie/Bachelorarbeit/DNA Samples/Zhehao_Hu_Bachelorthesis_Data.xlsx", sheet = "Pool Screen", range = "A4:P113")
raw <- raw %>% mutate(Haplotype = case_when(
  Haplotype == "HT1st" ~ "HT1*",
  Haplotype == "HT2/2st" ~ "HT2/2*",
  .default = Haplotype
))

# Test

raw <- raw %>% ungroup() %>% filter(Haplotype != "HT3" & Haplotype != "2023 HT1*")

shapiro <- raw %>% 
  filter(!is.na(Haplotype)) %>%
  group_by(Haplotype) %>%
  nest() %>%
  mutate(Shapiro = map(data, ~shapiro.test(.x$`Initial Target Copies`)))
shapiro_g <- shapiro %>% 
  mutate(shapiro_g = Shapiro %>% map(glance)) %>%
  unnest(shapiro_g)
shapiro_g

stat.result <- raw %>% ungroup() %>% filter(!is.na(`Initial Target Copies`)) %>% wilcox_test(`Initial Target Copies` ~ Haplotype)
stat.result

# Show significance
stat.result <- stat.result %>% add_xy_position(x="Haplotype", dodge = 0.8)
stat.result <- mutate(stat.result, y.position = log(y.position, base = 10))

#### Plot ####
n_fun.pl <- function(x){
  return(data.frame(y = 10,
                    label = length(x)))
} # Add sample size

f5.6 <- raw %>% filter(!is.na(`Ct Mean`)) %>% filter(!is.na(Haplotype)) %>% filter(!is.na(`Suggest Result`)) %>% filter(!is.na(`Initial Target Copies`) & Haplotype != "2023 HT1*") %>%
  group_by(Haplotype, `Suggest Result`) %>%
  ggplot(aes(x = Haplotype, y = `Initial Target Copies`)) +
  geom_boxplot(aes(fill = Haplotype)) +
  stat_summary(fun.data = n_fun.pl, geom = "text", vjust = -0.2) +
  theme_USGS_box() +
  theme(axis.line = element_blank(),
        legend.position = c(0.5,-0.2), #0.5, -0.2
        legend.direction = "horizontal",
        axis.text = element_text(size = 12),
        aspect.ratio = 0.8) +
  scale_fill_manual(values = cbp4) +
  scale_color_manual(values = c("#FBD173","#FF904D","#DA5858")) +
  ylab("") +
  xlab("")+
  labs(title = "") +
  stat_pvalue_manual(stat.result, step.increase = 0.05, size = 3) +
  scale_y_continuous(trans = scales::log10_trans(), 
                     breaks = c(1e5,1e7,1e9),
                     expand = expansion(mult = c(0.07,0.12)))
f5.6

#### Histo ####

histo.poolseq.lvl <-raw %>% group_by(Haplotype) %>% filter(!is.na(Haplotype)) %>%
  ggplot()+
  geom_histogram(aes(x=`Initial Target Copies`, fill = Haplotype)) +
  facet_wrap(vars(Haplotype),scales = "free_x") +
  theme_USGS_box() +
  theme(strip.text.x.top = element_text(face = "bold"),
        legend.position = c(-2.3,0.9),
        legend.direction = "horizontal",
        aspect.ratio = 1,
        axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_manual(values = cbp4) +
  scale_colour_manual(values = cbp4)

