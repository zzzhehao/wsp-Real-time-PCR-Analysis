# _wsp_ qPCR Analysis

R functions for evaluation of StepOne Realtime PCR System result using _wsp_qpcr_ Primer for detection of _wsp_ (_Wolbachia_ Surface Protein) fragments from _Altica lythri_.

This is a part of the supplements from my bachelor thesis **How’s _Wolbachia_ doing? Various _Wolbachia_ communities in the female-biased flea beetle _Altica lythri_.**

## 1 Map of content

```
|	
+-- Evaluation MasterScript.R # this script is outdated and no longer useful
+-- wsp-qPCR Evaluation.R # Code of the function
|
+-- Visualization/ # Scripts for visualization of the thesis data
	+-- AES.R # Global aesthetics setting
	...
```

## 2 Setup

All analyses were performed using R Statistical Software (v4.2.0; R Core Team, 2022) and RStudio (v2023.9.0.463).

## 3 Features after thesis submission

The original version corresponding to thesis submission is at the [commit 01edbdf](https://github.com/zzzhehao/wsp-Real-time-PCR-Analysis/tree/01ebdbf8a724d36446cb1e920add2505835d4866) 

- **feat:** add wsp_analysis_batch() as bulk analysis option
- **feat:** add arguments in function to customise IC and standard curve parameter
- **perf:** no need to specify plate_ID anymore, it also accepts characters  
- **feat:** no need to manually check for required packages  
- **perf:** improve log formatting and layout  
- **docs:** remove evaluation code information in plot document, expected to be added back as independent document in the future
- **fix**: removed customization of standard curve parameters, this now should be manually altered in **wsp-qPCR Evaluation.R** at line 42-49
- **docs**: add EvaCode documentation in README
- **feat**: adapt new workflow: interplate calibrator now has new default: "IPC", can also be specified manually; old workflow can be switched by specifying `ancient = True`
- **feat**: adapt new workflow: sample name now should be assigned to 'sample name' in StepOne software in default

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
- `wsp_analysis` analyses single report file.
- `wsp_analysis_batch` analyses multiple files in the same folder as a batch and multiple runs belonging to the same project.

#### 4.2.2 Usage

```R
wsp_analysis(xls_filePath, Ct_ic = NULL, SD_ic = NULL, in_batch = FALSE, ancient = FALSE, IC_name = "IPC")

wsp_analysis_batch(xls_folderPath, Ct_ic = NULL, SD_ic = NULL, ancient = FALSE, IC_name = "IPC")
```

#### 4.2.3 Arguments

`xls_filePath`: A path to the analysed _.xls_ file.

`xls_folderPath`: A path to the folder containing analysed _.xls_ files, all files must be placed directly under the folder.

`Ct_ic`: Optional. Customizable input for average Ct from interplate calibrator.

`SD_ic`: Optional. Customizable input for standard deviation of Ct from interplate calibrator. Plate with IC Ct outside this range will trigger warning in report. No direct effect on analyses.

`in_batch`: Internal. Used by `wsp_analysis_batch` to process all file together.

`ancient`: Optional. If true, assuming original workflow.

`IC_name`: Optional. If interplate calibrator is not named as default "IPC", the name should be specified here. It is important that this remains consistent for batch analysis.

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
- *EvaCode*: Traceback of the evaluation logic, see chapter 5
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

## 5 EvaCode

**EvaCode** is a three digit code assigned by `wsp_analysis`, it performs initial interpretation of the sample result, but is not fully reliable, taken rare occurring events into account. It provides a first view of the sample results and help for final interpretation.

Each code digit corresponds to a category after a set of examination. The category, i.e. the decimal code which embed in EvaCode is one-hot encoded answer for a set of binary examination. 

E.g., an evaluation code of “E-343” suggests “Positive sample; No Tm2, Tm1 in *wsp* range (positive), confirmed by melting curve; Clear for quantification; *w*LytA1 infection”, the last code (3) can be ignored since the first code already clearly assigned the sample as positive.

### 5.1 Digit 1: Negative or Positive?

| Digit of binary code | Test              | Binary options                                                                                                                                                                                               |
| -------------------- | ----------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| 1                    | `!is.na(Tm2)`     | 0: No Tm2; 1: Tm2 exists.                                                                                                                                                                                    |
| 2                    | `Tm1 > 77.5`      | 0: Tm1 < 77.5ºC, amplification is not in _wsp_ range, thus suggest as negative;<br>1: Tm1 > 77.5ºC, amplification is in _wsp_ range, thus suggest as positive.                                               |
| 3                    | `iM77/iM70 > 0.6` | 0: Main drop of melting curve is within 70ºC – 77ºC, likely 74ºC, thus suggest as negative;<br>1: Main drop of melting curve is not within 70ºC – 77ºC, likely within 77ºC – 79ºC, thus suggest as positive. |

| Binary Code | Decimal Code (1. EvaCode) | Result   | Description                                                                                                              | Warning message |
| ----------- | ------------------------- | -------- | ------------------------------------------------------------------------------------------------------------------------ | --------------- |
| 000         | 0                         | Negative | No Tm2, Tm1 not in _wsp_ range (negative), confirmed by melting curve.                                                   |                 |
| 001         | 1                         | Inspect  | No Tm2, Tm1 not in _wsp_ range (negative), inconsistent to melting curve.                                                | Inspect         |
| 010         | 2                         | Inspect  | No Tm2, Tm1 in _wsp_ range (positive) but inconsistent to melting curve.                                                 | Inspect         |
| 011         | 3                         | Quantify | No Tm2, Tm1 in _wsp_ range (positive), confirmed by melting curve, clear for quantification                              |                 |
| 100         | 4                         | Negative | Tm2 exists, Tm1 not in _wsp_ range (negative), confirmed by melting curve.                                               |                 |
| 101         | 5                         | Inspect  | Tm2 exists, Tm1 not in _wsp_ range (negative) but inconsistent to melting curve.                                         | Inspect         |
| 110         | 6                         | Inspect  | Tm2 exists, Tm1 in _wsp_ range (positive) but inconsistent to melting curve.                                             | Inspect         |
| 111         | 7                         | Positive | Tm2 exists, Tm1 in _wsp_ range (positive), confirmed by melting curve, presence of Tm2 suggests potential contamination. | Contamination   |

### 5.2 Digit 2: Infection Type

| Digit of binary code | Test          | Binary options                                                                                                                               |
| -------------------- | ------------- | -------------------------------------------------------------------------------------------------------------------------------------------- |
| 1                    | `Tm1 > 79.5`  | 0: Tm1 < 79.5ºC, non-*w*LytA1 infection;<br>1: Tm1 > 79.5ºC, *w*LytA1-infection.                                                             |
| 2                    | `is.na(DB)`   | 0: No double infection;<br>1: Double infected.                                                                                               |
| 3                    | `iM72 < iM77` | 0: No peak within range 72ºC-77ºC, thus suggest as *w*LytA2 + *w*LytB DB;<br>1: Peak within range 72ºC-77ºC, suggests *w*LytA1 involving DB. |

| Binary Code | Decimal Code (2. EvaCode) | Result             | Description                                                                                                                      | Warning message |
| ----------- | ------------------------- | ------------------ | -------------------------------------------------------------------------------------------------------------------------------- | --------------- |
| 000         | 0                         | *w*LytA2/*w*LytB   | *w*LytA2/*w*LytB infection.                                                                                                      |                 |
| 001         | 1                         | Inspect            | *w*LytA2/*w*LytB infection but inconsistent to melting curve.                                                                    | irregular DB    |
| 010         | 2                         | *w*LytA2 + *w*LytB | *w*LytA2 + *w*LytB double infection.                                                                                             |                 |
| 011         | 3                         | Inspect            | *w*LytA2 + *w*LytB double infection but inconsistent to melting curve.                                                           | irregular DB    |
| 100         | 4                         | *w*LytA1           | *w*LytA1 infection.                                                                                                              |                 |
| 101         | 5                         | Inspect            | *w*LytA1 involving double infection, but no DB detected.                                                                         | irregular DB    |
| 110         | 6                         | Inspect            | Tm1 suggests *w*LytA1 involving double infection, but inconsistent to melting curve.                                             | irregular DB    |
| 111         | 7                         | Inspect            | *w*LytA1 involving double infection, confirmed by melting curve, matches expectation but not observed in natural sample to date. | irregular DB    |

### 5.3 Digit 3: Negative Check

| Digit of binary code | Test         | Binary options                                                                                                                             |
| -------------------- | ------------ | ------------------------------------------------------------------------------------------------------------------------------------------ |
| 1                    | `Ct < 29`    | 0: Ct in range;<br>1: Ct too low.                                                                                                          |
| 2                    | `iM70 < 5.3` | 0: Melting curve matches NTC profile;<br>1: Melting curve suggest large amount of fluorescence signal at 70ºC, likely something amplified. |

| Binary Code | Decimal Code (3. EvaCode) | Result   | Description                               | Warning message |
| ----------- | ------------------------- | -------- | ----------------------------------------- | --------------- |
| 00          | 0                         | Negative |                                           |                 |
| 01          | 1                         | Inspect  | High fluorescence signal at 70ºC.         | Inspect         |
| 10          | 2                         | Inspect  | Ct low                                    | Inspect         |
| 11          | 3                         | Inspect  | Ct low, high fluorescence signal at 70ºC. | Inspect         |
