library(readxl)
library(dplyr)
library(ggplot2)
library(ggpattern)

# Import data
raw <- read_xlsx("Zhehao_Hu_Bachelorthesis_Data.xlsx", sheet = "Pool Screen", range = "A4:N113") %>% filter(!is.na(`Suggest Result`)) %>% filter(`Suggest Result`!="Anomaly") %>% filter(!is.na(`Haplotype`)) 

# Data manipulation
raw <- raw %>% mutate(Haplotype = case_when(
  Haplotype == "HT1st" ~ "HT1*",
  Haplotype == "HT2/2st" ~ "HT2/2*",
  .default = Haplotype
)) 
raw <- raw %>% mutate(`Suggest Result` = case_when(
  `Suggest Result` == "wA2/wB" ~ "wLytA2 / wLytB",
  `Suggest Result` == "wA1" ~ "wLytA1",
  `Suggest Result` == "wA2 + wB" ~ "wLytA2 + wLytB",
  .default = `Suggest Result`
))

# since ggplot do not have native pie chart geom
# summarize each groups
raw2 <- raw %>% group_by(Haplotype, `Suggest Result`) %>% tally()
# filter out unrelated data
raw2 <- raw2 %>% filter(Haplotype != "HT3" & Haplotype != "2023 HT1*")
# make result as factor instead of string
raw2$`Suggest Result` <- as.factor(raw2$`Suggest Result`)
# calculate total n from each group
raw2 <- merge(raw2, aggregate(n~Haplotype, raw2, sum), by = "Haplotype")
# calculate each group's portion in concerned haplotype
raw2 <- raw2 %>% ungroup() %>% group_by(Haplotype) %>% mutate(labels = `n.x`/`n.y`)
raw2

# Define theme
theme_USGS_box <- function(base_family = "serif", ...){
  theme_bw(base_family = base_family, ...) +
    theme(
      plot.title = element_text(size = 12),
      legend.background = element_rect(color = "black", fill = "white"),
      panel.background = element_rect(fill = "grey96", colour = "grey20")
    )
}

#### Plot ####
f5.5 <- 
  raw2 %>% filter(Haplotype != "HT3" & Haplotype != "2023 HT1*") %>% group_by(Haplotype) %>%
  ggplot(aes(x="", y=labels, group = `Suggest Result`))+
  geom_bar(aes(fill=`Suggest Result`), 
           stat = "identity") + # make a barplot of the result frequency in each haplotype (facet)
  ylim(0,1) + # unite y axis range since y is frequencies
  coord_polar(theta = "y") + # make barplot a pie chart
  geom_text(aes(label = round(`n.x`, digits = 2), x = 1.2),
            position = position_stack(vjust = 0.5),
            size = 2.5,
            vjust = -1,
            hjust = 0.4
  ) +
  facet_wrap(vars(Haplotype)) + 
  theme_USGS_box() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = c(0.5,-0.11),
        legend.direction = "horizontal",
        strip.text.x.top = element_text(face = "bold"),
        axis.line = element_blank(),
        panel.grid.major  = element_blank(),
        panel.grid.minor  = element_blank(),) +
  scale_fill_manual(values = cbp3) +
  xlab("") +
  ylab("") +
  labs(fill = "Wolbachia infection")+
  labs(title = "")
f5.5
