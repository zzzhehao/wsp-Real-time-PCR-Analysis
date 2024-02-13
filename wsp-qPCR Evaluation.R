# Setup
# Manually create a new folder called "Evaluation", all the report files will be exported there

wsp_analysis <- function() {
  
  library(readxl)
  library(plyr)
  library(purrr)
  library(rlist)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(stringr)
  library(ggpubr)
  library(gridExtra)
  library(knitr)
  library(log4r) #citation
  
  # Set working directory to current location
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  
  # Set path to xls. report file with all datasheets generated from StepOne Real-time PCR System
  file.xls.report <- paste0(xls_filePath, plate_ID,".xls")
  
  # Standard curve parameters (wA1)
  slope.a1 <- -3.595
  y_intcpt.a1 <- 39.274
  eff.a1 <- 89.732
  # Standard curve parameter (Average)
  slope.a2b <- -3.531
  y_intcpt.a2b <- 39.951
  eff.a2b <- 91.994
  
  # Log Setup
  logFile <- "wsp qPCR Evaluation.log"
  fileLogger <- logger(threshold = "INFO",appenders = file_appender(logFile))
  info(fileLogger, paste0("\n",
                          "================================================================= \n",
                          "## wsp qPCR Evaluation by Zhehao Hu (2024) ## \n",
                          "Plate ID: ", plate_ID, "\n",
                          "Date: ", format(Sys.time(), "%a %b %e %H:%M:%S %Y"),"\n",
                          "GitHub: https://github.com/zzzhehao/BA-Analysis", "\n",
                          "================================================================="))
  
  # IC parameter
  ic4 <- 19.72
  ic4SD <- 0.56
  
  #### Import ####
  # Melt norm Data import
  Mc.norm.raw <- as_tibble(readxl::read_xls(file.xls.report, sheet = "Melt Region Normalized Data", range = "A8:DQ104")) %>% select(-c(3,4)) %>% drop_na()
  print(file.xls.report)
  colnames(Mc.norm.raw)[c(1,2)] <- c("Well","Well_name")
  
  # Melt derv Data import
  Mc.drv.raw <- as_tibble(readxl::read_xls(file.xls.report, sheet = "Melt Region Derivative Data", range = "A8:DQ104")) %>% select(-c(3,4)) %>% drop_na()
  colnames(Mc.drv.raw)[c(1,2)] <- c("Well","Well_name")
  
  # Melt temp Data import
  Mc.temp.raw <- as_tibble(readxl::read_xls(file.xls.report, sheet = "Melt Region Temperature Data", range = "A8:DQ104")) %>% select(-c(3,4)) %>% drop_na()
  colnames(Mc.temp.raw)[c(1,2)] <- c("Well","Well_name")
  
  # Import result table
  output <- as_tibble(readxl::read_xls(file.xls.report, sheet = "Results", range = "A8:S104"))
  colnames(output)[c(1,3,7,8)] <- c("Well_name","Sample","Ct","Ct_Mean")
  
  # set well no. in output
  output <- left_join(Mc.temp.raw[,c(1,2)],output, by = "Well_name")
  
  # Import amplification data
  ampl <- as_tibble(readxl::read_xls(file.xls.report, sheet = "Amplification Data", range = "A8:E3848")) %>% drop_na()
  
  # aes setups
  theme_USGS_box <- function(base_family = "serif", ...){
    theme_bw(base_family = base_family, ...) +
      theme(
        plot.title = element_text(size = 12),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text = element_text(size = 10),
        plot.caption = element_text(size = 10),
        axis.ticks.length = unit(-0.05, "in"),
        axis.text.y = element_text(margin=unit(c(0.3,0.3,0.3,0.3), "cm")),
        axis.text.x = element_text(margin=unit(c(0.3,0.3,0.3,0.3), "cm")),
        axis.ticks.x = element_blank(),
        legend.background = element_rect(color = "black", fill = "white"),
        panel.background = element_rect(fill = "grey96", colour = "grey20"),
      )
  }
  
  #### Data preperation ####
  # Find NTC
  output <- output %>% mutate(`Sample` = case_when(
    Task == "NTC" ~ "NTC",
    .default = `Sample`
  ))
  
  # Assign global sample ID
  output <- output %>%
    mutate(Gl_ID = case_when(
      str_detect(Sample, "Sample ") ~ paste0(plate_ID, sprintf('%02d', as.numeric(gsub("Sample ", "", Sample)))),
      .default = paste0(plate_ID, "_", Sample))
    ) %>%
    relocate(Gl_ID, .before = Sample)
  
  # Initial melt curve signal
  output <- left_join(output,Mc.norm.raw[,c(2,3)], by = "Well_name")
  colnames(output)[22] <- "iM60"
  
  # 74-peak initial signal (Reading 34: 70.1ºC)
  output <- left_join(output,Mc.norm.raw[,c(2,36)], by = "Well_name")
  colnames(output)[23] <- "iM70"
  
  # wa1 involving DB-drop signal (Reading 40: 72ºC)
  output <- left_join(output,Mc.norm.raw[,c(2,42)], by = "Well_name")
  colnames(output)[24] <- "iM72"
  
  # 79-peak initial signal (Reading 58: 77.3ºC)
  output <- left_join(output,Mc.norm.raw[,c(2,60)], by = "Well_name")
  colnames(output)[25] <- "iM77"
  
  # 79-drop percentage
  output <- output %>% mutate(M77_70 = iM77/iM70)
  
  # Find DB drops
  neg.readings <- rowSums(Mc.drv.raw[c(3:43)] < 0)
  DB <- neg.readings > 0
  output <- left_join(output,as_tibble(data.frame(Mc.drv.raw[,2],DB,neg.readings)), by = "Well_name")
  
  # 79-peak percentage of initial signal
  output <- output %>% mutate(M77_70 = iM77/iM70)
  
  # detecting wa1-involving anomalies
  wa1DB <- output$iM72 < output$iM77
  output <- left_join(output,as_tibble(data.frame(wa1DB, Mc.drv.raw[,2])), by = "Well_name")
  
  # Generate tabular report
  output <- output %>% mutate(Tm2 = case_when(is.na(Tm2)~0,.default = Tm2))
  sreport <- output %>%
    group_by(Gl_ID) %>%
    reframe(
      Sample = Sample, 
      Ct_Mean = median(Ct), 
      Tm1 = median(Tm1), 
      Tm2 = median(Tm2), 
      iM60 = median(`iM60`),
      iM70 = median(`iM70`),
      iM72 = median(`iM72`),
      iM77 = median(`iM77`),
      M77_70 = median(`M77_70`),
      DB = DB,
      wa1DB = wa1DB,
      neg.readings = median(neg.readings)
    ) %>% unique()
  output <- output %>% mutate(Tm2 = case_when(Tm2 == 0 ~ NA,.default = Tm2))
  sreport[3] <- as.numeric(unlist(sreport[3]))
  sreport[c(3:6)] <- round(sreport[c(3:6)], digits = 2)
  sreport <- sreport %>% mutate(Tm2 = case_when(Tm2 == 0 ~ NA,.default = Tm2))
  
  #### Evaluation ####
  # 1.algorithm: positive or negative?
  sreport <- sreport %>% mutate(agr1 =  0 + 
                                  ifelse(!is.na(Tm2), 4, 0) + 
                                  ifelse(Tm1 > 77.5, 2, 0) + 
                                  ifelse(M77_70 > 0.6, 1, 0)
  )
  eva1 <- as_tibble(data.frame(agr1=c(0:7),eva1=c("Negative", rep("Inspect", 2), "Quantify", "Negative", rep("Inspect", 2), "Positive")))
  sreport <- left_join(sreport, eva1, by = "agr1")
  
  # 2.algorithm: infection type
  sreport <- sreport %>% mutate(agr2 =  0+
                                  ifelse(Tm1 > 79.5, 4, 0) + 
                                  ifelse(DB == TRUE, 2, 0) + 
                                  ifelse(iM72 < iM77, 1, 0)
  )
  eva2 <- as_tibble(data.frame(agr2=c(0:7),eva2=c("wA2/wB", "Inspect", "wA2 + wB", "Inspect", "wA1", rep("Inspect",3))))
  sreport <- left_join(sreport, eva2, by = "agr2")
  
  # 3.algorithm: negative double check
  sreport <- sreport %>% mutate(agr3 =  0+
                                  ifelse(Ct_Mean < 30, 2, 0) + 
                                  ifelse(iM70>5.3, 1, 0) 
  )
  eva3 <- as_tibble(data.frame(agr3=c(0:3),eva3=c("Negative", rep("Inspect",3))))
  sreport <- left_join(sreport, eva3, by = "agr3")
  
  # Summarize result
  sreport <- sreport %>% mutate(`Suggest result` = case_when(
    eva1 == "Negative" & eva3 == "Negative" ~ eva1,
    eva1 == "Negative" & eva3 != "Negative" ~ eva3,
    eva1 == "Quantify" | eva1 == "Positive" ~ eva2,
    eva1 == "Inspect" ~ eva1,
    .default = "ERROR"
  ))
  sreport <- sreport %>% mutate(Warning = case_when(
    `Suggest result` == "Inspect" ~ "Inspect",
    eva1 == "Positive" ~ "Contamination",
    .default = ""
  )) %>%
    mutate(EvaCode = paste0("E-",as.character(agr1), as.character(agr2), as.character(agr3)))
  
  # Calculate initial target copies
  icfactor <- as.numeric((sreport %>% filter(Sample == "IC4") %>% select(Ct_Mean) / ic4 )) ^ -1 # Calculate factor for calibrating IC
  icfactor <- c(rep(icfactor, length(sreport$Sample)))
  sreport <- sreport %>% mutate(Ct_cal = Ct_Mean * icfactor) #Calibrate IC
  
  sreport <- sreport %>% mutate(`Initial target copies` = case_when(
    eva1 == "Quantify" ~ case_when(
      eva2 == "wA1" ~ formatC(round(((10^(-1/slope.a1))^(y_intcpt.a1-Ct_Mean))*100, digit = 0), format = "e", digits = 2),
      .default = formatC(round(((10^(-1/slope.a2b))^(y_intcpt.a2b-Ct_Mean))*100, digit = 0), format = "e", digits = 2)
    )
  ))
  
  #### Log ####
  # Overview
  info(fileLogger, paste0(
    "  No. of Samples: ", nrow(sreport)-2, "\n",
    "  Double Infection: ", nrow((sreport %>% filter(DB == TRUE))), "\n",
    "  ⚠️ Inspection required: ", nrow((sreport %>% filter(`Suggest result` == "Inspect")))
  ))
  # Warning NTC
  CtNTC <- (sreport %>% filter(Sample == "NTC"))$Ct_Mean
  if (CtNTC < 30) {
    warn(fileLogger, paste0("⚠️ Irregular NTC, Ct low   Ct: ",
                            CtNTC))
  } else {
    info(fileLogger, paste0("✅ NTC in range   Ct: ",
                            CtNTC))
  }
  
  # Clean up evaluation ram
  sreport <- sreport %>% select(-c("agr1", "agr2", "agr3", "eva1", "eva2", "eva3"))
  
  #### Plots N Mc ####
  # Reshape melt data: Make a list of every well with it's melt curve raw data
  gthr <- function (df) {
    gather(df)[-1,]
  }
  l.temp <- dlply(Mc.temp.raw[,-c(2)],.(Well), gthr)
  l.norm <- dlply(Mc.norm.raw[,-c(2)],.(Well), gthr)
  l <- map2(l.temp, l.norm, left_join, by = "key")
  
  # add result summarise to each list
  linfo <- select(
    left_join(select(output, c("Gl_ID", "Well")), select(sreport, c("Gl_ID", "Ct_Mean","Tm1", "Tm2", "iM70", "DB", "neg.readings", "M77_70", "Suggest result","Initial target copies")), by = "Gl_ID"), 
    c("Well","Gl_ID","Ct_Mean", "Tm1", "Tm2", "iM70", "DB", "neg.readings", "M77_70", "Suggest result","Initial target copies"))
  linfo <- dlply(linfo,.(Well))
  l <- map2(l,linfo,list)
  
  # Group list of wells by samples
  merge <- function(list) {
    list.append(list, Target = as.character(list[[2]][1,2]))
  }
  l <- map(l, merge)
  l <- list.group(l, Target)
  
  # Define plot functions
  plot <- ggplot()
  drawplot <- function (list) {
    line <- geom_line(data = list[[1]], aes(x = value.x, y = value.y))
  }
  ap.layer <- function (list) {
    plot$layers <- list.append(plot$layers, list)
  }
  all.plot <- function (list) {
    plot$layers <- list()
    lines <- map(list, drawplot)
    plot$layers <- unlist(map(lines, ap.layer))
    
    plot <- plot +
      ylim(-1, 13) +
      xlab("Temperature (ºC)") +
      ylab("Normalized Reporter (Rn)") +
      labs(title = paste0(list[[1]][[2]][1,2], " Melt Curve"),
           caption = paste0(
             "Gl_ID: ", list[[1]][[2]][1,2], "\n\n",
             "Tm1: ", round(list[[1]][[2]][1,4], digit = 2), "\n",
             "Tm2: ", round(list[[1]][[2]][1,5], digit = 2), "\n",
             "70ºC Signal (Rn): ", round(list[[1]][[2]][1,6], digit = 2), "\n",
             "M77/70: ", round(list[[1]][[2]][1,9], digit = 4)*100, "%\n", "\n",
             "Double Infection: ", list[[1]][[2]][1,7], "\n",
             "Inhibited Reading Counts: ", list[[1]][[2]][1,8], "\n\n",
             "Suggest Result: ", list[[1]][[2]][1,10], "\n"
           ),
      ) +
      theme_USGS_box() +
      theme(plot.caption.position = "panel",
            aspect.ratio = 1,
            plot.caption = element_text(hjust = 0))
    
    return(plot)
    plot$layers <- list()
  }
  
  # Draw plots
  mcplots <- map(l, all.plot)
  
  #### Plots Amp ####
  # Reshape melt data
  colnames(ampl)[1] <- "Well_name"
  ampl <- left_join(ampl, select(output, c(1,2,4)))
  
  ampl <- dlply(ampl,.(Well), list)
  ampl <- map2(ampl,linfo,list)
  
  # Group list of wells by samples
  merge <- function(list) {
    list.append(list, Target = as.character(list[[1]][[1]][1,7]))
  }
  ampl <- map(ampl, merge)
  ampl <- list.group(ampl, Target)
  
  # Define plot functions
  plota <- ggplot()
  drawplot.a <- function (list) {
    linea <- geom_line(data = list[[1]][[1]], aes(x = Cycle, y = ΔRn))
  }
  ap.layer.a <- function (list) {
    plota$layers <- list.append(plota$layers, list)
  }
  all.plot.a <- function (list) {
    plota$layers <- list()
    linesa <- map(list, drawplot.a)
    plota$layers <- unlist(map(linesa, ap.layer.a))
    
    plota <- plota +
      scale_y_log10(limits = c(1e-5, 1e2)) +
      xlab("Cycle") +
      ylab("Normalized Reporter (Rn)") +
      labs(title = paste0("Amplification Plot"),
           caption = paste0(
             "Ct: ", round(list[[1]][[2]][3], digit = 2), "\n",
             "Initial Target Copies: ", list[[1]][[2]][11], "\n"
           )) +
      theme_USGS_box() +
      theme(plot.caption.position = "panel",
            aspect.ratio = 1,
            plot.caption = element_text(hjust = 0))
    
    return(plota)
    plota$layers <- list()
  }
  
  # Draw plots
  amplots <- map(ampl, all.plot.a)
  
  
  
  #### Plots info ####
  lgd <- ggplot()+
    labs(title = paste0("Evaluation Code Inspection Infos"),
         caption = paste0("\n","1.Digit","\n",
                          "1,5: M77/70 too high; 2,6: M77/70 too low","\n\n",
                          "2. Digit","\n",
                          "1,3,5,6,7,9,11,13,14,15: wA1 involving DB","\n\n",
                          "3. Digit","\n",
                          "1: iM70 too high; 2: Ct too low; 3: iM70 too high, Ct too low"),
    ) +
    theme_USGS_box() +
    theme(plot.caption.position = "panel",
          aspect.ratio = 0.01,
          plot.caption = element_text(hjust = 0))
  
  #### Export ####
  # Export Report
  sreport[6:10] <- round(sreport[6:10], digit = 2)
  sreport <- sreport %>% mutate(Status = "raw")
  write.csv(sreport, file = as.character(paste0("Evaluation/",plate_ID, " Evaluation Results.csv")))
  info(fileLogger, paste0("Tabular results are successfully saved at: " ,"Evaluation/",plate_ID, " Evaluation Results.csv"))
  
  # Export plots
  combinelist <- function (list1, list2) {
    list(list1, list2)
  }
  allplots <- map2(mcplots, amplots, combinelist)
  combineplots <- function (list) {
    cp <- ggarrange(list[[1]], list[[2]], align = "h")
    return(cp)
  }
  allplots <- map(allplots, combineplots)
  
  pagealn <- function(plt) {
    plt <- plt + theme(aspect.ratio = 0.8)
  }
  allplots <- map(allplots, pagealn)
  allplots <- list.append(allplots, lgd)
  
  ggsave(file=as.character(paste0("Evaluation/",plate_ID, " Amplification and Melt Curve Plots.pdf")), marrangeGrob(grobs = allplots, nrow = 2, ncol = 1, top = NULL), width = 210, height = 297, units = "mm", device = "pdf")
  info(fileLogger, paste0("Graphic results are successfully saved at: " ,"Evaluation/",plate_ID, " Amplification and Melt Curve Plots.pdf"))
  
  print(paste0("Evaluation of Plate No. ", plate_ID, " has been succesfully completed."))
  info(fileLogger, paste0("✅ Evaluation of Plate ID. ", plate_ID, " has been succesfully completed. No fatal events occured."))
  info(fileLogger, paste0("Please be aware, the automatic analysis cannot replace manual inspection. Evaluation code will provide further detailed analysis results. \n", "# \n", "# \n", "# End of evaluation \n", "# \n", "#"))

}