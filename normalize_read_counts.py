# These functions will take as input the gene_counts.tsv file
# and will normalize all the read counts per gene.
def normalize_values(file):
    f = open(file, 'r')

    gene_ID = []
    sum_SC1 = 0
    sum_SC2 = 0
    sum_SC3 = 0
    sum_SC4 = 0
    sum_SC5 = 0
    sum_SC6 = 0

    # Here we are getting the sum of the all the
    # read counts in the sample column (SC1, SC2 etc.).
    # line[i] specifies the read counts in a specific column.
    # line[2] is the read counts in column SC2 for example,
    # because the 2nd index of line is column SC2.
    # This sum is the library size.

    # We are also appending all the gene IDs
    # found in the first column into the list, gene_ID
    f.readline()
    for line in f:

        line = line.replace('\n', "")
        line = line.split("\t")

        gene_ID.append(line[0])

        sum_SC1 += int(line[1])
        sum_SC2 += int(line[2])
        sum_SC3 += int(line[3])
        sum_SC4 += int(line[4])
        sum_SC5 += int(line[5])
        sum_SC6 += int(line[6])

    normalized_SC1 = []
    normalized_SC2 = []
    normalized_SC3 = []
    normalized_SC4 = []
    normalized_SC5 = []
    normalized_SC6 = []

    # Here we are getting the normalized value, by dividing
    # each read count with the library size for that sample
    # column. We then append the normalized value into a new list
    # specific to that column

    f = open(file, 'r')
    f.readline()
    for line in f:

        line = line.replace('\n', "")
        line = line.split("\t")

        norm_value_SC1 = int(line[1]) / sum_SC1
        normalized_SC1.append(norm_value_SC1)

        norm_value_SC2 = int(line[2]) / sum_SC2
        normalized_SC2.append(norm_value_SC2)

        norm_value_SC3 = int(line[3]) / sum_SC3
        normalized_SC3.append(norm_value_SC3)

        norm_value_SC4 = int(line[4]) / sum_SC4
        normalized_SC4.append(norm_value_SC4)

        norm_value_SC5 = int(line[5]) / sum_SC5
        normalized_SC5.append(norm_value_SC5)

        norm_value_SC6 = int(line[6]) / sum_SC6
        normalized_SC6.append(norm_value_SC6)

    return gene_ID, normalized_SC1, normalized_SC2, normalized_SC3, normalized_SC4, normalized_SC5, normalized_SC6

