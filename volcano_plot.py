import matplotlib.pyplot as plt
import csv
import numpy as np

# This function will initialize the x and y values we need for our volcano
# plot into a list form. We need the log2FoldChain values, padj values
# and the gene names
def form_lists(file):
    f = open(file, 'r')
    csv_reader = csv.reader(f, delimiter=",")

    log2FC_list = []
    padj_list = []
    gene_name_list = []
    # In this loop we are appending the log2FoldChain values,
    # padj values and gene names found in column 4, 8, and 9
    # of the annotated file
    f.readline()
    for row in csv_reader:

        log2FC = float(row[3])
        log2FC_list.append(log2FC)

        padj = row[7]
        padj_list.append(padj)

        gene_name = row[8]
        gene_name_list.append(gene_name)

    return log2FC_list, padj_list, gene_name_list


# Now we need to clean our lists. The padj list has several
# empty values,because there is no padj value for that gene.
# We want to remove them. If we remove the padj value, we also
# need to remove the log_2FC value with the same index, and similarly
# the gene name with the same index. This way we can do the scatterplot
# correctly.
def clean_lists():

    log2FC_list, padj_list, gene_name_list = form_lists("results_L-vs-NL_full_annotation.csv")


    # This for loop is making an index list, for all the indexes of
    # all the empty values in the padj list.
    index_list = []
    for index in range(len(padj_list)):
        if padj_list[index] == '':
            index_list.append(index)

    # This loop is appending all the empty values (in padj_list) as noted down by the
    # index list, into empty_padj_list. List comprehension is used to say
    # that for values in padj_list, if they are not in empty_padj_list
    # (if they arent empty values and are therefore actual values), they can
    # stay and be assigned to the variable, cleaned_padj_list

    # The next lines are saying that for the all indexes in index list,
    # the index value in log2FC_list will be changed to ''. The value will be
    # be changed empty for that index. List comprehension is used to say that
    # all values in log2FC_list, they can stay if they are not empty and will be
    # assigned to the variable removed_log2FC_list. The same logic applies to
    # gene_name_list
    empty_padj_list = []
    for index in index_list:
        empty_padj_list.append(padj_list[index])
        log2FC_list[index] = ''
        gene_name_list[index] = ''

    cleaned_padj_list = ([i for i in padj_list if i not in empty_padj_list])
    cleaned_log2FC_list = ([i for i in log2FC_list if i != ''])
    cleaned_gene_name_list = ([i for i in gene_name_list if i != ''])

    # This loop is taking all the values in cleaned_padj_list, and
    # getting the -log10 value of it. It is appending the -log10 values
    # into log10pvaluelist
    log10pvaluelist = []
    for value in cleaned_padj_list:
        float_value = float(value)
        new_value = abs(np.log10(float_value))
        log10pvaluelist.append(new_value)

    return cleaned_log2FC_list, log10pvaluelist, cleaned_gene_name_list

# This function will plot an annotated volcano plot, with a line
# drawn at FDR = 0.05 and Log2(fold change) = 2
def volcano_plot():

    cleaned_log2FC_list, log10pvaluelist, cleaned_gene_name_list = clean_lists()

    plt.scatter(cleaned_log2FC_list, log10pvaluelist, color = "red")
    plt.ylim([0, 30])
    plt.xlim([-5, 4])
    plt.xlabel("log2 Fold Change")
    plt.ylabel("-log10(P-value)")
    plt.title("Volcano Plot")
    # The lines below (plt.axvline and plt.axhline) enables us to draw a straight line
    # that crosses the x and y line in the values listed below
    plt.axvline(x=2, color='k')
    plt.axhline(y=1.3010299956639813, color='k', linestyle='-')

    # We are now going to annotate the outliers. The for loop below
    # is getting the index for x and y in a specific range. It is then
    # appending the indexes for those x and y values into their respective lists.
    # We need our x and y indexes to match because that will be the coordinate
    # for our candidate gene. We make the lists into sets and get the intersection.
    # We assign the intersection into the variable, upperleft_index_gene_names.

    #UPPERLEFT OUTLIERS
    upperleft_index_x_list = []
    upperleft_index_y_list = []
    for x,y in zip(cleaned_log2FC_list,log10pvaluelist):
        if x < -3 and y > 17:
            index_x = cleaned_log2FC_list.index(x)
            upperleft_index_x_list.append(index_x)

            index_y = log10pvaluelist.index(y)
            upperleft_index_y_list.append(index_y)

    upperleft_index_gene_names = (list(set(upperleft_index_y_list).intersection(upperleft_index_x_list)))

    # The indexes in upperleft_index_gene_names can be used to return values
    # such as gene names we want to annotate and their coordinates. The x and
    # y coordinates are cleaned_log2FC_list and log10pvaluelist

    for index in upperleft_index_gene_names:
        plt.text(cleaned_log2FC_list[index] - 0.1, log10pvaluelist[index] - 0.1, cleaned_gene_name_list[index], horizontalalignment='right', size='small', color='black',
                weight='semibold')

    # This follows the same logic as above
    #UPPERRIGHT OUTLIERS
    upperright_index_x_list = []
    upperright_index_y_list = []
    for x,y in zip(cleaned_log2FC_list,log10pvaluelist):
        if x > 2.2 and y > 8.5:
            index_y = log10pvaluelist.index(y)
            upperright_index_y_list.append(index_y)

            index_x = cleaned_log2FC_list.index(x)
            upperright_index_x_list.append(index_x)


    upperright_index_gene_names = (list(set(upperright_index_y_list).intersection(upperright_index_x_list)))

    for names in upperright_index_gene_names:
        plt.text(cleaned_log2FC_list[names] - 0.1, log10pvaluelist[names] - 0.1, cleaned_gene_name_list[names], horizontalalignment='right', size='small', color='black',
                weight='semibold')


    plt.show()



