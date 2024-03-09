# Calling required packages
library(readxl)
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(purrr)

# Import Data
raw <- read_xlsx("Zhehao_Hu_Bachelorthesis_Data.xlsx", sheet = "Reproduction", range = "A3:P61") %>% filter(!is.na(Group))
colnames(raw)[5:16]  <- c("sum_layed","sum_hatched", "sum_grown", "6/27_layed", "6/27_hatched", "6/27_grown","7/4_layed", "7/4_hatched", "7/4_grown","7/11_layed", "7/11_hatched", "7/11_grown")
raw <- raw %>% mutate(across(c(3,4), as.factor))

# calculate hatch rate
raw <- raw %>% mutate(sum_hatch_rate = sum_hatched / sum_layed)

library(rstatix)
library(broom)

rp.test <- raw # Create submatrix for statistic test

shapiro <- rp.test %>%
  filter(!is.na(Haplotype)) %>% # Drop NAs
  group_by(Group, Haplotype) %>% # Samples with same treatment (Group) and haplotype were grouped together
  nest() %>%
  mutate(Shapiro = map(data, ~shapiro.test(.x$sum_hatch_rate))) # Store the test result
shapiro_g <- shapiro %>% 
  mutate(shapiro_g = Shapiro %>% map(glance)) %>%
  unnest(shapiro_g)

shapiro_g # Display the test result
stat.result <- 
  rp.test %>% 
  ungroup() %>% 
  filter(!is.na(sum_hatch_rate)) %>% # Drop NAs
  group_by(Haplotype) %>% # Apply seperate tests by grouping samples according to haplotype
  t_test(sum_hatch_rate ~ Group) %>% # T-test comparing hatch rate between different treatment
  add_significance()
stat.result # Display the dif test result

# Show significance
stat.result <- stat.result %>% add_xy_position(x="Group", dodge = 0.8)

# test all groups
stat.result.g <- 
  rp.test %>% 
  ungroup() %>% 
  filter(!is.na(sum_hatch_rate)) %>% # Drop NAs
  group_by(Group) %>%
  t_test(sum_hatch_rate ~ Haplotype) %>% # T-test comparing hatch rate between different treatment
  add_significance()

stat.result.g # Display the t test result -> control HT1 vs control HT2 interesting
stat.result.g <- stat.result.g[2,]
stat.result.g <- stat.result.g %>% add_xy_position(x="Haplotype", dodge = 0.8)
stat.result.g # extract

# Labelling haplotypes
n.hatchrate <- c(
  'HT1' = "HT1",
  'HT1st' = "HT1*",
  'HT2/2st' = "HT2/2*"
)

# Labelling sample size
n_fun.hr <- function(x){
  return(data.frame(y = 1.2,
                    label = length(x)))
}

f.58A <- 
  ggplot(raw,aes(x=Group,y=sum_hatch_rate)) + 
  geom_boxplot(aes(fill=Group)) + 
  facet_grid(.~Haplotype,labeller=as_labeller(n.hatchrate),scales='free_x') +
  stat_summary(fun.data = n_fun.hr, geom = "text", vjust = -0.2) +
  theme_USGS_box() +
  theme(strip.text.x.top = element_text(face = "bold"),
        legend.position = "none",
        aspect.ratio = 3,
        axis.text.x = element_text(angle = 45, hjust = 1),
  ) +
  scale_fill_manual(values = cbp5) +
  scale_colour_manual(values = cbp5) +
  xlab("") +
  ylab("Egg hatch rate") +
  stat_pvalue_manual(stat.result, step.increase = 0.01, size = 4, tip.length = 0.01)+
  scale_y_continuous(labels = scales::percent,
                     breaks = c(0,0.2,0.4,0.6,0.8,1),
                     expand = expansion(mult = c(0.07,0.15)))
f.58A
