# Ptarmigan Telomere Analysis
Python and R code for calculating the relative telomere lengths of Icelandic rock ptarmigans (*Lagopus muta*).

Authors: Jasmine Baclig, Molly J. Hansen

## Content
- indiv_effic_files: folder containing output of linregpcr.py in CSV format.
- lm_info_files: folder containing files needed to label and identify each ptarmigan sample.
- rdml_files: folder containing RDML files from the LightCycler software.
- txt_files: folder containing TXT files of the qPCR data from the LightCycler software.
- linregpcr.py: Python code for running LinRegPCR on all the files included in rdml_files and their corresponding files in txt_files.
- main_script.R: R code for running all other R scripts included here.
- make_all_samples.R: R code that creates a dataframe about the year, sex, and age of all the ptarmigan samples (written by MJH).
- output.txt ; TXT file containing the console output of linregpcr.py, which can be used to track the progress of the code.
- rtl_calculation.R: R code that analyzes the qPCR data and calculates the relative telomere lengths of each ptarmigan sample.
- sex_age.R: R code that extracts sex and age information from the input files in lm_info_files (written by MJH).

## How To Use
- Do the following for each LightCycler file to be analyzed:
  - Optional: exclude
  - Go to "File" > "Save As", and in the "Save as type" drop-down menu, select "RDML file".