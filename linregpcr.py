import rdmlpython as rdml
import os

# Gets all RDML files in specified directory
path = "C://Users//angelajasminebac//Documents//Python//rdml_files"
file_list = list(filter(lambda name: ".rdml" in name, os.listdir(path)))

for file in file_list:
    # Validates file
    cli_validate = rdml.Rdml()
    cli_resValidate = cli_validate.validate(filename = file)
    print(cli_resValidate)

    # Runs LinRegPCR
    cli_linRegPCR = rdml.Rdml(file)
    cli_expList = cli_linRegPCR.experiments()
    cli_exp = cli_expList[0]
    cli_runList = cli_exp.runs()
    cli_run = cli_runList[0]
    cli_result = cli_run.linRegPCR(updateRDML = True, saveResultsCSV = True, timeRun = True, verbose = True)

    with open("//indiv_effic_files" + file[:-4] + "csv", "w") as cli_f:
        cli_f.write(cli_result["resultsCSV"])