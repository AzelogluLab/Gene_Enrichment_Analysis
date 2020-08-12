def fix_go_names(directory, file):
    ##Description
    ##Reads a .csv file with the GO-Ontology Data, then fixes the term names to not show the ID, then writes back to the file with the appended name. 
    ##Arguments
    ##directory: The directory where the csv files are
    ##file: The name of the file that you want to fix the names for
    ##Returns:
    ##None

    import os 
    import csv
    import pandas as pd
    #Change directory
    os.chdir(directory)
    #Read in the file
    data = pd.read_csv(file)
    #Loop through the terms and remove the tags 
    for i in range(0, len(data["Terms"])):
        data.at[i, "Terms"] = data["Terms"][i][0:-13]
    #Save to csv
    data.to_csv(file)
    

def fix_mgi_names(directory, file):
    ##Description
    ##Reads a .csv file with the MGI Data, then fixes the term names to not show the ID, then writes back to the file with the appended name. 
    ##Arguments
    ##directory: The directory where the csv files are
    ##file: The name of the file that you want to fix the names for
    ##Returns:
    ##None

    import os 
    import csv
    import pandas as pd 

    #Change directory
    os.chdir(directory)
    #Read in the file
    data = pd.read_csv(file)
    #Loop through the terms and remove the tags
    for i in range(0, len(data["Terms"])):
        data.at[i, "Terms"] = data["Terms"][i][11:]
    #Save to csv
    data.to_csv(file)
    

def fix_wiki_names(directory, file):
    ##Description
    ##Reads a .csv file with the WIKI Pathways Data, then fixes the term names to not show the ID, then writes back to the file with the appended name. 
    ##Arguments
    ##directory: The directory where the csv files are
    ##file: The name of the file that you want to fix the names for
    ##Returns:
    ##None

    import os 
    import csv
    import pandas as pd  

    #Change directory
    os.chdir(directory)
    #Read in the file
    data = pd.read_csv(file)
    #Loop through the terms and remove the tags
    for i in range(0, len(data["Terms"])):
        fix = data["Terms"][i].split()[0:-1]
        fix = ' '.join(fix)
        data.at[i, "Terms"] = fix
    #Save to csv 
    data.to_csv(file)
