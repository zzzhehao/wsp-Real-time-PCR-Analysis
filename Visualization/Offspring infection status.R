library(readxl)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(reshape2)
library(scales)
library(rstatix)
library(ggh4x)
library(ggsignif)
library(ggpubr)

# Import data
datapath <- "/Zhehao_Hu_Bachelorthesis_Data.xlsx"
# Parents data
mamas <- read_xlsx(datapath, sheet = "Tetracycline", range = "A3:V61")
# F1 data
kids <- read_xlsx(datapath, sheet = "Tetracycline Offspring", range = "A3:Y51")
dsx <- kids

# Data preparation
mamas <- mamas %>% filter(!is.na(Gl_ID)) %>% mutate(Partner_M = NA) %>% filter(`Initial Target Copies` != "Positive")
colnames(mamas)[2:3] <- c("PID","OID")
mamas$OID <- as.character(mamas$OID)
mamas$Gl_ID <- as.character(mamas$Gl_ID)

kids <- kids %>% filter(!is.na(Gl_ID)) %>% select(-c(6:14,20,22)) %>% filter(`Initial Target Copies` != "Positive")
kids$`Initial Target Copies` <- as.numeric(kids$`Initial Target Copies`)
mamas$`Initial Target Copies` <- as.numeric(mamas$`Initial Target Copies`)
colnames(kids)[2:3] <- c("PID","OID")

# Merge parents and F1 data
df <- bind_rows(mamas, kids) %>% arrange(PID) %>% 
  mutate(Name = case_when( # assign unique ID
    is.na(OID) ~ paste0(PID, " Parent"),
    .default = paste0(PID, " ", OID)), 
    .before = PID)

df <- df %>% filter(PID %in% select(filter(df,duplicated(df$PID)), PID)$PID)
df$`Initial Target Copies` <- as.numeric(df$`Initial Target Copies`)

# calculate percentage of F1 infection level divided by parents (Not shown)
df <- df %>% group_by(PID) %>% mutate(per = `Initial Target Copies`/max(`Initial Target Copies`))
df <- df %>% mutate(disp = paste0(round(per*100, digit = 3),"%"))

#### Fig 5.8C dsx review ####
dsx <- dsx %>% filter(!is.na(`dsx Sex`))
f5.8C.dsxp <- dsx %>% group_by(Group) %>%
  ggplot() +
  geom_histogram(aes(x = factor(interaction(Haplotype, Group), levels = c("HT1.Control", "HT1st.Tetracycline", "HT2/2st.Control", "HT2/2st.Tetracycline")), fill = `dsx Sex`),stat = "count", position = position_stack()) +
  scale_x_discrete(labels = c("HT1 Control", "HT1* Tetracycline", "HT2/2* Control", "HT2/2* Tetracycline")) +
  theme_USGS_box() +
  theme(strip.text.x.top = element_text(face = "bold", size = 10),
        legend.direction = "horizontal",
        legend.position = c(0.3,0.78),
        aspect.ratio = 1.2,
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10)
  ) +
  scale_fill_manual(values = c("#4249C1", "#DC75A8"), name = "dsx Result", labels = c("female", "male")) +
  ylab("Count") +
  xlab("")
f5.8C.dsxp

#### Fig 5.8B F1 offspring infection level ####
df <- df %>% mutate(generation = case_when(
  is.na(OID) ~ "Parent",
  !is.na(OID) ~ "F1"
))
# test
library(rstatix)
library(tidyverse)
library(broom)

f1in.test <- df %>% filter(!is.na(`Initial Target Copies`) & !is.na(Haplotype))

shapiro <- f1in.test %>% 
  filter(!is.na(Haplotype)) %>%
  group_by(generation, Group) %>%
  nest() %>%
  mutate(Shapiro = map(data, ~shapiro.test(.x$`Initial Target Copies`)))
shapiro_g <- shapiro %>% 
  mutate(shapiro_g = Shapiro %>% map(glance)) %>%
  unnest(shapiro_g)
shapiro_g

stat.result <- f1in.test %>% ungroup() %>% group_by(Group) %>% t_test(`Initial Target Copies` ~ generation) %>% add_significance()
stat.result

stat.result <- stat.result %>% add_xy_position(x = "generation", dodge = 0.8)
stat.result <- mutate(stat.result, y.position = (log(y.position, base = 10)+0.5) - 0.05)

# Fig 5.8B Visualization
n_fun.fin <- function(x){
  return(data.frame(y = 9,
                    label = length(x)))
}
df <- df %>% mutate(generation = case_when(
  .default = generation,
  generation == "Parent" ~ "Mother"
  )
)
pfin.ord <- c("Mother", "F1")
f5.8b.f1in <- 
  df %>% filter(!is.na(Group)) %>% filter(!is.na(Haplotype)) %>%
  ggplot(aes(x = factor(generation, level = pfin.ord), y = `Initial Target Copies`)) +
  geom_boxplot(aes(fill = Group)) +
  facet_grid(.~Group,scales='free_x') +
  stat_summary(fun.data = n_fun.fin, geom = "text", vjust = -0.2) +
  theme_USGS_box() +
  theme(strip.text.x.top = element_text(face = "bold", size = 10),
        aspect.ratio = 3,
        legend.direction = "horizontal",
        legend.position = c(0.25,0.2),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.line = element_blank()) +
  scale_fill_manual(values = cbp5) +
  xlab("") +
  ylab("wsp gene copies/µg DNA") +
  stat_pvalue_manual(stat.result, step.increase = 0.05, size = 3) +
  scale_y_continuous(trans = scales::log10_trans(), 
                     breaks = c(1e3,1e5,1e7, 1e9),
                     expand = expansion(mult = c(0.07,0.12)))
f5.8b.f1in

#### NOT SHOWN - Pbar family overview new ####
pbar <-
  df %>% filter(!is.na(Group)) %>% filter(!is.na(Haplotype)) %>% arrange(OID) %>%
  ggplot(aes(x=Name, y=`Initial Target Copies`, fill=PID))+
  geom_bar(stat = "identity") +
  geom_text(aes(label = disp), size = 2, vjust = 3) +
  geom_text(aes(label = Haplotype), size = 2, vjust = 4.5) +
  scale_y_log10(limits = c(1e0,1e9)) +
  theme_USGS_box() +
  theme(strip.text.x.top = element_text(face = "bold", size = 10),
        aspect.ratio = 0.2,
        legend.position = "none",
        legend.direction = "horizontal",
        axis.text.x = element_text(angle = 45),
        axis.line = element_blank()) +
  facet_wrap(.~factor(generation, level = c("Parent", "F1"))~Group, nrow=2, scales = "free_x", labeller = function(x) {x[2]}) +
  xlab("") +
  ylab("") 
pbar
