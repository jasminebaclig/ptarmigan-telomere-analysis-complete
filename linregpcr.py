# Calculates individual efficiencies for qPCR assays
# Author: Jasmine Baclig

import rdmlpython as rdml
import os
import sys

# Stores path for input and output files
input_folder = ".//rdml_files//"
output_folder = ".//indiv_effic_files//"

# Gets all RDML files in specified directory
file_list = list(filter(lambda name: ".rdml" in name, os.listdir(input_folder)))

# Prints output to file instead of console
sys.stdout = open("output.txt", "w")

for file in file_list:
    print("========== WORKING ON " + file + " ==========")

    # Validates file
    cli_validate = rdml.Rdml()
    cli_resValidate = cli_validate.validate(filename = input_folder + file)
    print(cli_resValidate)

    # Runs LinRegPCR
    cli_linRegPCR = rdml.Rdml(input_folder + file)
    cli_expList = cli_linRegPCR.experiments()
    cli_exp = cli_expList[0]
    cli_runList = cli_exp.runs()
    cli_run = cli_runList[0]
    cli_result = cli_run.linRegPCR(updateRDML = True, saveResultsCSV = True, timeRun = True, verbose = True)

    # Saves output to csv
    with open(output_folder + file[:-4] + "csv", "w") as cli_f:
        cli_f.write(cli_result["resultsCSV"])

    print("\n\n")

# Closes output file
sys.stdout.close()
