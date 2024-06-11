# _wsp_ qPCR Analysis

R functions for evaluation of StepOne Realtime PCR System result using _wsp_qpcr_ Primer for detection of _wsp_ (_Wolbachia_ Surface Protein) fragments from _Altica lythri_.

This is a part of the supplements from my bachelor thesis **How’s _Wolbachia_ doing? Various _Wolbachia_ communities in the female-biased flea beetle _Altica lythri_.**

## 1 Map of content

```
|	
+-- Evaluation MasterScript.R # this script is outdated no longer useful
+-- wsp-qPCR Evaluation.R # Code of the function
|
+-- Visualization/ # Scripts for visualization of the thesis data
	+-- AES.R # Global aesthetics setting
	...
```

## 2 Setup

All analysis were performed using R Statistical Software (v4.2.0; R Core Team, 2022) and RStudio (v2023.9.0.463).

## 3 Features after thesis Submission

The original version corresponding to thesis submission is at the [commit 01edbdf](https://github.com/zzzhehao/wsp-Real-time-PCR-Analysis/tree/01ebdbf8a724d36446cb1e920add2505835d4866) 

- **feat:** add wsp_analysis_batch() as bulk analysis option
- **feat:** add arguments in function to customise IC and standard curve parameter
- **perf:** no need to specify plate_ID anymore, it also accepts characters  
- **feat:** no need to manually check for required packages  
- **perf:** improved log formatting and layout  
- **docs:** removed evaluation code information in plot document, expected to be added back as independent document in the future
- **fix**: removed customization of standard curve parameters, this now should be manually altered in **wsp-qPCR Evaluation.R** at line 42-49

## 4 Analyse Data

### 4.1 Raw file

The _.xls_ file exported by StepOne Real time PCR system must contain following sheets in the same file, named as the plate ID, which will be assigned to sample name and used in reports later.

- Sample Setup
- Results
- Amplification Data
- Melt Region Temperature Data
- Melt Region Normalized Data
- Melt Region Derivative Data

### 4.2 Analyse data with functions

#### 4.2.1 Description

The wsp_analysis functions perform the initial data preparation and interpretation regarding *Wolbachia* infection status.
- `wsp_analysis` analyse single report file.
- `wsp_analysis_batch` analyse multiple files in the same folder as a batch and multiple runs belonging to the same project.

#### 4.2.2 Usage

```R
wsp_analysis(xls_filePath, Ct_ic = NULL, SD_ic = NULL, in_batch = FALSE )

wsp_analysis_batch(xls_folderPath, Ct_ic = NULL, SD_ic = NULL)
```

#### 4.2.3 Arguments

`xls_filePath`: A path to the analysed _.xls_ file.

`xls_folderPath`: A path to the folder containing analysed _.xls_ files, all files must be placed directly under the folder.

`Ct_ic`: Optional. Customizable input for average Ct from interplate calibrator.

`SD_ic`: Optional. Customizable input for standard deviation of Ct from interplate calibrator. Plate with IC Ct outside this range will trigger warning in report. No direct effect on analyses.

`in_batch`: Internal. Used by `wsp_analysis_batch` to process all file together.

#### 4.2.4 Outputs

Operating the function will export two reports to output folder. 
- When analyzing single file, the output folder is named as _Evaluation_.
- When analyzing multiple files, the output folder is named as _Batch_Evaluation\__\[current-date]

**Evaluation Results.csv** presents tabular results of each sample with following critical value:
- *Gl_ID*: unique sample ID
- *Sample*: Sample name inherited from StepOne repots
- *Ct_Mean*: Mean value of Ct from all technical replicates
- *Tm1*: Temperature of the first melting curve peak
- *Tm2*: Temperature of the second melting curve peak
- *iM60*, *iM70*, *iM72*, *iM77*: Nomalized melting curve signal at each temperature
- *DB*: Double infection
- *wa1DB*: Double infection involving strain *w*LytA1
- *Suggest result*: Suggested result, yet not fully reliable
- *Warning*: If *Inspect* shows, manual interpretation is mandatory
- *EvaCode*: Traceback of the evaluation logic
- *Ct_cal*: Adjusted Ct regarding to interplate calibrator
- *Initial target copies*: (if available) Fragment copies per µg DNA

**Amplification and Melt Curve Plots.pdf** presents graphical result, i.e. amplification plot and melting curve plot of each sample. All technical replicates are plotted together, with short report underneath.

#### 4.2.5 Examples

```R
# Call the functions
source("wsp-qPCR Evaluation.R")

# Analyse single file
wsp_analysis("Sample Data/20230406.xls")

# Analyse multiple files
wsp_analysis_batch("Sample Data")

# Customize parameters
wsp_analysis("Sample Data/20230512.xls", Ct_ic = 18.28, SD_ic = 0.43)
```

