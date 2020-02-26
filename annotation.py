import csv
import mygene
import pandas as pd

# This function will annotate the table of results with the
# HUGO gene name and the gene description.
def annotation(file):

    mg = mygene.MyGeneInfo()
    f = open(file, 'r')
    csv_reader = csv.reader(f, delimiter=",")

    # Getting a list of the gene_IDs available in the first column
    # of the csv available
    gene_ID_list = []
    for row in csv_reader:
        gene_ID = row[0]
        gene_ID_list.append(gene_ID)

    # Removing the title and putting our list through mygene, to get more information
    # such as the HUGO name (symbol) and the summary. This will all be in dictionary form in
    # results
    gene_ID_list.remove('gene_ID')
    results = mg.querymany(gene_ID_list,scopes='ensembl.gene', fields='symbol,summary')

    symbol = []
    summary = []
    # If the key of the dictionary (results) contains symbol or summary,
    # we can append the value of that key to their respective lists
    # If the keys doesn't exist, this loop will just append "NA" for that index
    for hit in results:

        if "symbol" in hit:
            symbol.append(hit['symbol'])
        else:
            symbol.append('NA')

        if "summary" in hit:
            summary.append(hit['summary'])
        else:
            summary.append('NA')

    list_symbol = symbol
    list_summary = summary

    # Now we're using PANDAS to output our lists (list_symbol, list_summary) into
    # 2 additional columns to our original file. We are giving them the
    # title "HUGO Gene Name" and "Gene Description" respectively. We're going to
    # create a new file that contains the columns of the original file and our 2
    # additional columns. We will name this file "results_L-vs-NL_full_annotation.csv"

    df = pd.read_csv(file)
    df['HUGO Gene Name'] = pd.DataFrame({'HUGO Gene Name': list_symbol})
    df['Gene Description'] = pd.DataFrame({'Gene Description': list_summary})
    df.to_csv("results_L-vs-NL_full_annotation.csv")


    f.close()

