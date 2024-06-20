import rdmlpython as rdml
import os

path = "C://Users//angelajasminebac//Documents//Python"
print(list(filter(lambda name: ".rdml" in name, os.listdir(path))))


"""
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

with open("result.csv", "w") as cli_f:
    cli_f.write(cli_result["resultsCSV"])
"""