

#Setup functions
def getUserID(directory, file, title, number):
    ##getUserID
    ##Description
    ##Takes Gene IDs and sends them to the EnrichR API, returns an ID associated to the call
    ##Arguments:
    ##directory: The directory of your input gene file
    ##file: The filename of the gene file
    ##title: The name of the title you want associated with the API call (Does not affect the statistical outcome)
    ##number: The number of top genes you want as your input sorting by p-value
    ##Returns:
    ##user_list_id: The id number associated with the EnrichR API

    import json
    import requests
    import csv
    import os 
    import pandas as pd

    ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/addList'
    #Change directory
    os.chdir(directory)
    #Load in the file
    data = pd.read_csv(file)
    #Sort the data by the p-value
    data = data.sort_values(by = ['PValue'])
    #Take the top number of genese specified in the parameters
    genes = data.head(number)
    #Format the values for EnrichR
    genes = genes[genes.columns[1]].tolist()
    genes = '\n'.join(genes)

    #Organize API Call
    description = title
    payload = {
        'list': (None, genes),
        'description': (None, description)
    }

    response = requests.post(ENRICHR_URL, files = payload)
    if not response.ok:
        raise Exception('Error analyzing gene list')

    data = json.loads(response.text)
    user_list_id = data['userListId']
    return(user_list_id)


def getAllGenesID(directory, file1, file2, title, number):
    ##getAllGenesData
    ##Description
    ##Takes Gene IDs from both the upregulated and downregulated gene lists, makes and API call and retrieves
    ##the associated id.
    ##Arguments:
    ##directory: The directory of your input gene files
    ##file1: The filename of the gene file (upregulated)
    ##file2: The filename of the gene file (downregulated)
    ##title: The name of the title you want associated with the API call (Does not affect the statistical outcome)
    ##number: The number of top genes you want as your input sorting by p-value
    ##Returns:
    ##user_list_id: The id number associated with the EnrichR API
    import json
    import requests
    import csv
    import os 
    import pandas as pd

    ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/addList'
    #Change directory
    os.chdir(directory)
    #Load in the data from file1
    data = pd.read_csv(file1)
    #Sort the data by p-value
    data = data.sort_values(by = ['PValue'])
    #Take the top number of gene speicified in the parameters
    genes = data.head(number)
    #Format the values for EnrichR
    genes1 = genes[genes.columns[1]].tolist()
    #Load in the data from file2
    data = pd.read_csv(file2)
    #Sort the data by p-value
    data = data.sort_values(by = ['PValue'])
    #Take the top number of genes specified in the parameters
    genes = data.head(number)
    #Format the values for EnrichR
    genes2 = genes[genes.columns[1]].tolist()
    genes = genes1 + genes2
    genes = '\n'.join(genes)

    #Organize API Call
    description = title
    payload = {
        'list': (None, genes),
        'description': (None, description)
    }

    response = requests.post(ENRICHR_URL, files = payload)
    if not response.ok:
        raise Exception('Error analyzing gene list')

    data = json.loads(response.text)
    print(data)
    user_list_id = data['userListId']
    return(user_list_id)


def getEnrichResults(ID, library):
    ##getEnrichResults
    ##Description
    ##Fetches the results from the api call
    ##Arguments:
    ##directory: The directory of your input gene files
    ##ID: The ID associated with API call
    ##library: The library you want enrich the genes with
    ##Returns:
    ##data: The data from the api call 
    import json
    import requests

    ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/enrich'
    query_string = '?userListId=%s&backgroundType=%s'
    user_list_id = ID
    gene_set_library = library
    response = requests.get(
        ENRICHR_URL + query_string % (user_list_id, gene_set_library)
    )
    if not response.ok:
        raise Exception('Error fetching enrichment results')
    data = json.loads(response.text)
    return(data)

def makeDataframe(data, library, directory, file):
    ##makeDataframe
    ##Description
    ##Takes the data from the API and converts it into a dataframe
    ##Arguments:
    ##data: The data from the api call
    ##library: The library you want enrich the genes with 
    ##directory: The directory of your input gene files
    ##file: The file you want the dataframe to save to 
    ##Returns:
    ##df: The dataframe created from the data retrieved from the EnrichR API 
    import os 
    import pandas as pd

    #Extract data
    data  = data[library]
    #Create empty vectors
    term_list = []
    p_val_list = []
    genes_list = []
    p_val_adj_list = []    
    #Parse data and fill the vectors 
    for i in range(0, len(data) - 1):
        output = data[i]
        term = output[1]
        p_val = output[2]
        genes = output[5]
        p_val_adj = output[6]
        term_list.insert(i, term)
        p_val_list.insert(i, p_val)
        genes_list.insert(i, genes)
        p_val_adj_list.insert(i, p_val_adj) 
    #Convert vectors into a dataframe
    data_list =  {'Terms':term_list, 'Genes':genes_list, 'p-value':p_val_list, 'Adjusted p-value':p_val_adj_list}
    df = pd.DataFrame(data_list)
    #Change Directory
    os.chdir(directory)
    #Save dataframe to a file 
    df.to_csv(file, header = True)
    return(df)

def EnrichR_pipeline(directory, output_directory, up_file, down_file, title, number, library):
    ##Description
    ##Take the gene files and run them through EnrichR and return the results in a .csv format.
    ##Arguments
    ##directory: The directory of the gene files
    ##output_directory: The directory where you want to save the data from EnrichR
    ##up_file: The name of the file that contains the upregulated genes
    ##down_file: The name of the file that contains the downregulated genes
    ##title: The name of the title you want associated with the API call (Does not affect the statistical outcome)
    ##number: The number of top genes you want as your input sorting by p-value
    ##library: The library you want enrich the genes with
    ##Returns:
    ##None

    #Obtain the IDs aorresponding to the API calls
    upID = getUserID(directory = directory, file = up_file, title = title, number = number)
    downID = getUserID(directory = directory, file = down_file, title = title, number = number)
    up_down_ID = getAllGenesID(directory = directory, file1 = up_file, file2 = down_file, title = title, number = number)
    #Extract data from the API call using the IDs 
    up_data = getEnrichResults(ID = upID, library = library)
    down_data = getEnrichResults(ID = downID, library = library)
    up_down_data = getEnrichResults(ID = up_down_ID, library = library)
    #Turn data into dataframes
    up_data_df = makeDataframe(data = up_data, library = library, directory = output_directory, file = str("Upreg_" + library + ".csv"))
    down_data_df = makeDataframe(data = down_data, library = library, directory = output_directory, file = str("Downreg_" + library + ".csv"))
    up_down_df = makeDataframe(data = up_down_data, library = library, directory = output_directory, file = str("Upreg_Downreg_" + library + ".csv"))



#Run Script
directory = 
output_directory = 
up_file = 
down_file = 
number_of_genes = 
lib1 = 
EnrichR_pipeline(directory = directory, output_directory = output_directory, up_file = up_file, down_file = down_file, title = 'YAP', number = number_of_genes, library = lib)

