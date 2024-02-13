# Set working directory to current location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Set the Plate List (will be used for Gl_ID and file names)
# ❗️Please ensure all .xls files contains all datasheets and named with the format "Plate_ID.xls" and saved under the xls_filePath
xls_filePath <- "xls report/"
plate_IDs <- c(
  20230406,
  20230502,
  20230503,
  20230504,
  20230505,
  20230522,
  20230714,
  20230717,
  20230816,
  20230823,
  20230828,
  20230831,
  20231020
)

source("wsp-qPCR Evaluation.R") # Calling wsp Evaluation Script

logFile <- "wsp qPCR Evaluation.log"
fileLogger <- logger(threshold = "INFO",appenders = file_appender(logFile)) # Logger set up

info(fileLogger, paste0("\n",
                        "================================================================= \n",
                        "## wsp qPCR Evaluation MasterScript ## \n",
                        "Initialized", "\n",
                        "Date: ", format(Sys.time(), "%a %b %e %H:%M:%S %Y"),"\n",
                        "GitHub: https://github.com/zzzhehao/BA-Analysis", "\n",
                        "================================================================= \n",
                        "# \n", "Evaluation begins \n", "#"))

for (plate_ID in plate_IDs) {
  wsp_analysis()
} # Applying Evaluation

info(fileLogger, paste0("# \n", "# \n", "# End of Evaluation MasterScript \n", "\n", "\n"))