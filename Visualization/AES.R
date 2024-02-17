theme_USGS_box <- function(base_family = "serif", ...){
  theme_bw(base_family = base_family, ...) +
    theme(
      plot.title = element_text(size = 12),
      axis.ticks.length = unit(-0.05, "in"),
      axis.text.y = element_text(margin=unit(c(0.3,0.3,0.3,0.3), "cm")),
      axis.text.x = element_text(margin=unit(c(0.3,0.3,0.3,0.3), "cm")),
      axis.ticks.x = element_blank(),
      aspect.ratio = 1,
      legend.background = element_rect(color = "black", fill = "white"),
      panel.background = element_rect(fill = "grey96", colour = "grey20")
    )
}

#### Colorblindness-friendly palette ####
# https://davidmathlogic.com/colorblind/#%23648FFF-%23785EF0-%23DC267F-%23FE6100-%23FFB000

cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

cbp2 <- (c("#4249C1","#8F7CE4","#DC75A8","#FF904D","#FBD173"))

cbp3 <- (c("#999999","#FBD173","#FF904D","#DA5858")) # Wolbachia strains

cbp4 <- (c("#4249C1","#8F7CE4","#DC75A8")) # Haplotype

cbp5 <- (c("#888888","#35CCC9")) # Antibiotic treatment