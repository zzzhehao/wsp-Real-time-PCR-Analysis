library(readxl)
library(dplyr) 
library(stringr)
library(ggplot2)
library(ggpubr)
library(plyr)
library(reshape2)
library(gridExtra)
library(ggsignif)
library(rstatix)
library(tidyverse)
library(broom)

# Import data
raw <- read_xlsx("/Users/hu_zhehao/Library/Mobile Documents/com~apple~CloudDocs/UHH/B.Sc. Biologie/Bachelorarbeit/DNA Samples/Zhehao_Hu_Bachelorthesis_Data.xlsx", sheet = "Tetracycline", range = "A3:U61") %>% filter(!is.na(Group)) %>% filter(!is.na(Gl_ID))

# Data manipulation
raw <- raw %>% select(c(1,4,5,14,15,16,20,21))
colnames(raw)[8] <- "quantity"
raw <- raw %>% 
  mutate(Interpretation = case_when(
    `Suggest result` == "Negative" ~ "negative",
    `quantity` == "Positive" ~ "positive",
  )) %>%
  mutate(Interpretation = replace(Interpretation, is.na(Interpretation), "quantified"))
raw$Interpretation <- as.factor(raw$Interpretation)
r.quantified <- filter(raw, Interpretation == "quantified") %>% filter(!is.na(Haplotype) & !is.na(Group))
r.quantified$`quantity` <- as.numeric(r.quantified$`quantity`) 

#### Test ####
tc.test <- r.quantified

shapiro <- tc.test %>% 
  filter(!is.na(Haplotype)) %>%
  group_by(Group, Haplotype) %>%
  nest() %>%
  mutate(Shapiro = map(data, ~shapiro.test(.x$quantity)))
shapiro_g <- shapiro %>% 
  mutate(shapiro_g = Shapiro %>% map(glance)) %>%
  unnest(shapiro_g)
shapiro_g

stat.result <- tc.test %>% ungroup() %>% filter(!is.na(quantity)) %>% group_by(Haplotype) %>% t_test(quantity ~ Group) %>% add_significance()
stat.result

stat.result <- stat.result %>% add_xy_position(x = "Group", dodge = 0.8)
stat.result <- mutate(stat.result, y.position = log(y.position, base = 10))

#### Fig 5.7B Quantified: quantity_groups vars Haplotype ####
n_fun.p1 <- function(x){
  return(data.frame(y = 8.3,
                    label = length(x)))
}
n.haplotypes.p1 <- c(
  'HT1' = "HT1",
  'HT1st' = "HT1*",
  'HT2/2st' = "HT2/2*"
)
colnames(r.quantified)[8] <- "quantity"
f5.7B.quantity_group.quantified <- 
  r.quantified %>% filter(!is.na(Haplotype)) %>%
  ggplot(aes(x=Group,y=quantity)) + 
  geom_boxplot(aes(fill=Group)) +
  scale_y_log10() +
  facet_grid(.~Haplotype,labeller=as_labeller(n.haplotypes.p1),scales='free_x') +
  stat_summary(fun.data = n_fun.p1, geom = "text", hjust = 0.5) +
  theme_USGS_box() +
  theme(strip.text.x.top = element_text(face = "bold", size = 10),
        aspect.ratio = 3,
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10)
  ) +
  scale_fill_manual(values = cbp5) +
  xlab("") +
  ylab("Initial wsp-Gene Copies") +
  stat_pvalue_manual(stat.result, step.increase = 0.01, size = 4, tip.length = 0.01)+
  scale_y_continuous(trans = scales::log10_trans(), 
                     breaks = c(1e5,1e7,1e8),
                     expand = expansion(mult = c(0.07,0.10)))
f5.7B.quantity_group.quantified

#### Fig 5.7A Bar: categories of infection status ####
n_fun.pA <- function(x){
  return(data.frame(y = -18,
                    label = length(x)))
}
n.haplotypes.pA <- c(
  'HT1' = "HT1",
  'HT1st' = "HT1*",
  'HT2/2st' = "HT2/2*"
)
f5.7A.bar.infection_categories <- 
  raw %>% filter(!is.na(Haplotype)) %>%
  ggplot(aes(x = Interpretation, fill = Group))+
  geom_histogram(stat = "count", binwidth = 1500) +
  stat_count(binwidth = 1500, 
             geom = "text", 
             color = "white", 
             aes(label = after_stat(count),
                 group = Group),
             size = 3.5,
             position = position_stack(vjust = 0.55)) +
  facet_wrap(vars(Haplotype),labeller=as_labeller(n.haplotypes.pA),scales='free_x') +
  theme_USGS_box() +
  theme(strip.text.x.top = element_text(face = "bold", size = 10),
        legend.position = "none",
        aspect.ratio = 3,
        axis.text.x = element_text(angle = 45, hjust = 1, size = 11)
  ) +
  scale_fill_manual(values = cbp5) +
  scale_color_manual(values = c(NA, NA, "red"))+
  xlab("") +
  ylab("Number of Observations")
f5.7A.bar.infection_categories

#### All group test ####
raw$quantity <- as.numeric(raw$quantity)
rawtest <- raw %>% filter(!is.na(Haplotype) & !is.na(Group)) %>% ungroup() %>% filter(!is.na(quantity)) %>% mutate(testg = as.factor(paste0(Haplotype, Group))) %>% t_test(quantity ~ testg)
rawtest