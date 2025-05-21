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
  - Additionally, the rtl_calculation.R code recognizes two gene names that was used in this study: "TOX" (reference gene) and "TELO" (gene of interest). Either change the gene names in the LightCycler files or edit the R code so that the naming remains consistent. Note that the positive and negative controls should also have a gene label.
  - Go to "File" > "Save As", and in the "Save as type" drop-down menu, select "RDML file". Save this in rdml_files.
  - In the "Analysis" tab, find the "Results Table" and click the top-leftmost cell to select all cells. Right-click anywhere in the middle, and click "Export to File". Save this TXT file in txt_files, and make sure that the name matches with the corresponding RDML file in rdml_files (or at least they are in the same order/position in their own folder after the files are sorted alphabetically by rtl_calculation.R).
- Run linregpcr.py (Note: this can be run in RStudio). Afterward, check if there are the same number of files in indiv_effic_files as there are in rdml_files.
- If ptarmigan data already included in this repository will be used, run main_script.R.
- If different data will be used only for calculating relative telomere lengths, run rtl_calculation.R. Here are some important dataframes in this code:
  - bad_effic_data: indicates which samples have a bad efficiency that was calculated by linregpcr.py.
    - If a GB sample appears here, the whole plate should be repeated.
    - Note: there is no dataframe included that indicates which samples have a bad Cq error.
  - duplicates: indicates which samples have more than one instance where it passed in one assay.
    - The dataframe sample_data can be used to find out which run the duplicate data are in.
  - missing_data: indicates which samples are not yet done based on a given list of all samples available
    - If working with a different set of samples than what is included, line 71 of the code must be edited by changing what list of all samples should be used by the script.
  - calculated_data: includes the relative telomere lengths for all complete samples.
    - An additional line of code should be written to export this dataframe into a CSV file.