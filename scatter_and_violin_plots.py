import normalize_read_counts as nrc
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# This function will return as a list the -log2 values of
# the normalized read counts for each sample of a
# candidate gene.
def log_list_outliers():
    # Here are are getting the gene_ID and normalized sample
    # lists as already calculated by normalize_values
    gene_ID, normalized_SC1, normalized_SC2, normalized_SC3, normalized_SC4, normalized_SC5, normalized_SC6 = nrc.normalize_values("gene_counts.tsv")

    # We know that for the gene name MCAM, the ID is ENSG00000076706.
    # We can get the index of that ID as it falls in the gene ID list
    # This same index will return the normalized read count value
    # for each of the 6 different samples. This is because they the gene
    # ID and normalized values for each sample are on the same row and have
    # the same index. We are putting the normalized value for each
    # sample into a list. For all the values of that list we are going to
    # get the -log2 value of it and will append it into the log2_MCAM list
    # The log 2 lists will be the y value of our scatterplot

    log2_MCAM = []
    index_MCAM = gene_ID.index("ENSG00000076706")
    list_MCAM = [normalized_SC1[index_MCAM], normalized_SC2[index_MCAM], normalized_SC3[index_MCAM]
                ,normalized_SC4[index_MCAM], normalized_SC5[index_MCAM], normalized_SC6[index_MCAM]]
    for number in list_MCAM:
        log2_MCAM.append(abs(np.log2(number)))

    # The code below follows the same logic as above except with
    # different genes

    log2_KCNIP3 = []
    index_KCNIP3 = gene_ID.index("ENSG00000115041")
    list_KCNIP3 = [normalized_SC1[index_KCNIP3], normalized_SC2[index_KCNIP3], normalized_SC3[index_KCNIP3],
                   normalized_SC4[index_KCNIP3], normalized_SC5[index_KCNIP3], normalized_SC6[index_KCNIP3]]
    for number in list_KCNIP3:
        log2_KCNIP3.append(abs(np.log2(number)))

    log2_KRT75 = []
    index_KRT75 = gene_ID.index("ENSG00000170454")
    list_KRT75 = [normalized_SC1[index_KRT75], normalized_SC2[index_KRT75], normalized_SC3[index_KRT75],
                  normalized_SC4[index_KRT75], normalized_SC5[index_KRT75], normalized_SC6[index_KRT75]]
    for number in list_KRT75:
        log2_KRT75.append(abs(np.log2(number)))

    log2_CPA4 = []
    index_CPA4 = gene_ID.index("ENSG00000128510")
    list_CPA4 = [normalized_SC1[index_CPA4], normalized_SC2[index_CPA4], normalized_SC3[index_CPA4],
                 normalized_SC4[index_CPA4], normalized_SC5[index_CPA4], normalized_SC6[index_CPA4]]
    for number in list_CPA4:
        log2_CPA4.append(abs(np.log2(number)))

    log2_LMX1B = []
    index_LMX1B = gene_ID.index("ENSG00000136944")
    list_LMX1B = [normalized_SC1[index_LMX1B], normalized_SC2[index_LMX1B], normalized_SC3[index_LMX1B],
                  normalized_SC4[index_LMX1B], normalized_SC5[index_LMX1B], normalized_SC6[index_LMX1B]]
    for number in list_LMX1B:
        log2_LMX1B.append(abs(np.log2(number)))

    log2_ADAM8 = []
    index_ADAM8 = gene_ID.index("ENSG00000151651")
    list_ADAM8 = [normalized_SC1[index_ADAM8], normalized_SC2[index_ADAM8], normalized_SC3[index_ADAM8],
                  normalized_SC4[index_ADAM8], normalized_SC5[index_ADAM8], normalized_SC6[index_ADAM8]]
    for number in list_ADAM8:
        log2_ADAM8.append(abs(np.log2(number)))

    log2_CDT1 = []
    index_CDT1 = gene_ID.index("ENSG00000167513")
    list_CDT1 = [normalized_SC1[index_CDT1], normalized_SC2[index_CDT1], normalized_SC3[index_CDT1],
                 normalized_SC4[index_CDT1], normalized_SC5[index_CDT1], normalized_SC6[index_CDT1]]
    for number in list_CDT1:
        log2_CDT1.append(abs(np.log2(number)))

    log2_RINL = []
    index_RINL = gene_ID.index("ENSG00000187994")
    list_RINL = [normalized_SC1[index_RINL], normalized_SC2[index_RINL], normalized_SC3[index_RINL],
                 normalized_SC4[index_RINL], normalized_SC5[index_RINL], normalized_SC6[index_RINL]]
    for number in list_RINL:
        log2_RINL.append(abs(np.log2(number)))

    log2_TFPI2 = []
    index_TFPI2 = gene_ID.index("ENSG00000105825")
    list_TFPI2 = [normalized_SC1[index_TFPI2], normalized_SC2[index_TFPI2], normalized_SC3[index_TFPI2],
                  normalized_SC4[index_TFPI2], normalized_SC5[index_TFPI2], normalized_SC6[index_TFPI2]]
    for number in list_TFPI2:
        log2_TFPI2.append(abs(np.log2(number)))

    return log2_MCAM, log2_KCNIP3, log2_KRT75, log2_CPA4, log2_LMX1B, log2_ADAM8, log2_CDT1, log2_RINL, log2_TFPI2


# This function will read a file (samples_data.tsv) and will
# put the 3rd column of the file (type) into list form.
# This will serve as the x value in our scatterplot and histogram
def type_list(file):

    f = open(file, 'r')

    f.readline()
    type = []
    for line in f:
        line = line.replace('\n', "")
        line = line.split("\t")
        type.append(line[2])

    return type

# This function will give us a scatterplot of
# the normalized gene counts per sample type for each
# candidate gene
def scatterplot():

    log2_MCAM, log2_KCNIP3, log2_KRT75, log2_CPA4, log2_LMX1B, log2_ADAM8, log2_CDT1, log2_RINL, log2_TFPI2 = log_list_outliers()
    type = type_list("samples_data.tsv")


    plt.subplot(3, 3, 1)
    plt.title("MCAM")
    plt.xlabel("Sample Type")
    plt.ylabel("-log2(Normalized Gene Counts)", size=8)
    plt.scatter(type, log2_MCAM, color="red")

    plt.subplot(3, 3, 2)
    plt.title("KCNIP3")
    plt.xlabel("Sample Type")
    plt.ylabel("-log2(Normalized Gene Counts)", size=8)
    plt.scatter(type, log2_KCNIP3, color="blue")

    plt.subplot(3, 3, 3)
    plt.title("KRT75")
    plt.xlabel("Sample Type")
    plt.ylabel("-log2(Normalized Gene Counts)", size=8)
    plt.scatter(type, log2_KRT75, color="black")

    plt.subplot(3, 3, 4)
    plt.title("CPA4")
    plt.xlabel("Sample Type")
    plt.ylabel("-log2(Normalized Gene Counts)", size=8)
    plt.scatter(type, log2_CPA4, color="yellow")

    plt.subplot(3, 3, 5)
    plt.title("LMX1B")
    plt.xlabel("Sample Type")
    plt.ylabel("-log2(Normalized Gene Counts)", size=8)
    plt.scatter(type, log2_LMX1B, color="purple")

    plt.subplot(3, 3, 6)
    plt.title("ADAM8")
    plt.xlabel("Sample Type")
    plt.ylabel("-log2(Normalized Gene Counts)", size=8)
    plt.scatter(type, log2_ADAM8, color="orange")

    plt.subplot(3, 3, 7)
    plt.title("CDT1")
    plt.xlabel("Sample Type")
    plt.ylabel("-log2(Normalized Gene Counts)", size=8)
    plt.scatter(type, log2_CDT1, color="green")

    plt.subplot(3, 3, 8)
    plt.title("RINL")
    plt.xlabel("Sample Type")
    plt.ylabel("-log2(Normalized Gene Counts)", size=8)
    plt.scatter(type, log2_RINL, color="pink")

    plt.subplot(3, 3, 9)
    plt.title("TFPI2")
    plt.xlabel("Sample Type")
    plt.ylabel("-log2(Normalized Gene Counts)", size=8)
    plt.scatter(type, log2_TFPI2, color="grey")

    plt.show()

# This function will give us a violin plot of
# the normalized gene counts per sample type for each
# candidate gene
def violin_plot():

    log2_MCAM, log2_KCNIP3, log2_KRT75, log2_CPA4, log2_LMX1B, log2_ADAM8, log2_CDT1, log2_RINL, log2_TFPI2 = log_list_outliers()
    type = type_list("samples_data.tsv")


    plt.subplot(3, 3, 1)
    plt.xlabel("Sample Type")
    plt.ylabel("-log2(Normalized Gene Counts)", size= 8)
    plt.title("MCAM")
    ax = sns.violinplot(x=type, y=log2_MCAM, palette = "BuGn_r")

    plt.subplot(3, 3, 2)
    plt.xlabel("Sample Type")
    plt.ylabel("-log2(Normalized Gene Counts)", size= 8)
    plt.title("KCNIP3")
    ax = sns.violinplot(x=type, y=log2_KCNIP3, palette = "BuGn_r")

    plt.subplot(3, 3, 3)
    plt.xlabel("Sample Type")
    plt.ylabel("-log2(Normalized Gene Counts)", size= 8)
    plt.title("KRT75")
    ax = sns.violinplot(x=type, y=log2_KRT75, palette = "BuGn_r")

    plt.subplot(3, 3, 4)
    plt.xlabel("Sample Type")
    plt.ylabel("-log2(Normalized Gene Counts)", size= 8)
    plt.title("CPA4")
    ax = sns.violinplot(x=type, y=log2_CPA4, palette = "BuGn_r")

    plt.subplot(3, 3, 5)
    plt.xlabel("Sample Type")
    plt.ylabel("-log2(Normalized Gene Counts)", size= 8)
    plt.title("LMX1B")
    ax = sns.violinplot(x=type, y=log2_LMX1B, palette = "BuGn_r")

    plt.subplot(3, 3, 6)
    plt.xlabel("Sample Type")
    plt.ylabel("-log2(Normalized Gene Counts)", size= 8)
    plt.title("ADAM8")
    ax = sns.violinplot(x=type, y=log2_ADAM8, palette = "BuGn_r")

    plt.subplot(3, 3, 7)
    plt.xlabel("Sample Type")
    plt.ylabel("-log2(Normalized Gene Counts)", size= 8)
    plt.title("CDT1")
    ax = sns.violinplot(x=type, y=log2_CDT1, palette = "BuGn_r")

    plt.subplot(3, 3, 8)
    plt.xlabel("Sample Type")
    plt.ylabel("-log2(Normalized Gene Counts)", size= 8)
    plt.title("RINL")
    ax = sns.violinplot(x=type, y=log2_RINL, palette = "BuGn_r")

    plt.subplot(3, 3, 9)
    plt.xlabel("Sample Type")
    plt.ylabel("-log2(Normalized Gene Counts)", size= 8)
    plt.title("TFIP2")
    ax = sns.violinplot(x=type, y=log2_TFPI2, palette = "BuGn_r")

    plt.show()

