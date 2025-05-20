# Ptarmigan Telomere Analysis
Python and R code for calculating the relative telomere lengths of Icelandic rock ptarmigans (*Lagopus muta*).

Authors: Jasmine Baclig, Molly J. Hansen

## Content
- indiv_effic_files: folder containing output of linregpcr.py in CSV format.
- lm_info_files: folder containing files needed to label and identify each ptarmigan sample.
- rdml_files: folder containing RDML files from the LightCycler software.
- txt_files: folder containing TXT files of the qPCR data from the LightCycler software.
- data_aggregration.R: R code that combines all ptarmigan-related data.
- linregpcr.py: Python code for running LinRegPCR on all the files included in rdml_files and their corresponding files in txt_files.
- main_script.R: R code for running all other R scripts included here.
- make_all_samples.R: R code that creates a dataframe about the year, sex, and age of all the ptarmigan samples (written by MJH).
- output.txt: TXT file containing the console output of linregpcr.py, which can be used to track the progress of the code.
- rtl_calculation.R: R code that analyzes the qPCR data and calculates the relative telomere lengths of each ptarmigan sample.
- sex_age.R: R code that extracts sex and age information from the input files in lm_info_files (written by MJH).

## How To Use
- If different data will be used, do the following for each LightCycler file to be analyzed:
  - Exclude samples that appears multiple times (not including replicate reactions) so that each sample has a passing result on each assay only once. This may include samples that were repeated for reasons other than a bad Cq error or a bad efficiency. If unsure, continue anyway, and the rtl_calculation.R code will indicate which samples have more than one data point in the "duplicates" dataframe.
  - Currently, the rtl_calculation.R code recognizes negative control with the sample name "Negative" and positive control with the sample name "GB". Either change the sample names in the LightCycler files or edit the R code so that the naming remains consistent. Standards are not used at all throughout the data and code.
  - Go to "File" > "Save As", and in the "Save as type" drop-down menu, select "RDML file". Save this in rdml_files.
  - In the "Analysis" tab, find the "Results Table" and click the top-leftmost cell to select all cells. Right-click anywhere in the middle, and click "Export to File". Save this TXT file in txt_files, and make sure that the name matches with the corresponding RDML file in rdml_files.
- Run linregpcr.py (Note: this can be run in RStudio). Afterward, check if there are the same number of files in indiv_effic_files as there are in rdml_files.
- If ptarmigan data already included in this repository will be used, run main_script.R.
- If different data will be used only for calculating relative telomere lengths, run rtl_calculation.R.