import matplotlib.pyplot as plt
import seaborn as sns; sns.set()
from scipy.spatial import distance
import numpy as np
import normalize_read_counts as nrc


# This function will output a heatmap given the 6 lists
# of normalized values obtained in normalize_read_counts.py
# To get the heatmap values we need to calculate the
# euclidean distance between the 6 lists

def heatmap():

    gene_ID, normalized_SC1, normalized_SC2, normalized_SC3, normalized_SC4, normalized_SC5, normalized_SC6 = nrc.normalize_values("gene_counts.tsv")

    normalized_lists = (normalized_SC1, normalized_SC2, normalized_SC3, normalized_SC4, normalized_SC5, normalized_SC6)
    heatmap_data = []

    # This function will take calculate the euclidean distance
    # between any 2 lists and append them into the list uniform_data
    def euclidean_distance(listx, listy):
        heatmap_data.append(distance.euclidean(listx, listy))
        return heatmap_data

    # This function calls has 2 loops. The outer loop goes
    # through the numbers 0-5, and the inner loop does the same.
    # The innerloop calls the euclidean_distance function.
    # The euclidean_distance function will be performed on
    # all 2 paired combinations of the 6 lists. For example,
    # normalized_lists[0] = normalized_SC1, normalized_lists[0] =
    # normalized_SC1 (because that is the 0th index of normalized_lists).
    # So the first combination will be (normalized_SC1, normalized_SC1)
    # The inner loop will iterate through the entire list, so the next value will be
    # (normalized_SC1, normalized_SC2), (normalized_SC1, normalized_SC3)
    # and so on. When the inner loop finishes iterating through the list, the outer loop
    # loop will iterate  and go to the 1st index from the 0th. The inner
    # loop will reset and iterate again through the list again.
    # We will get (normalized_SC2, normalized_SC1), (normalized_SC2, normalized_SC2)
    # (normalized_SC2, normalized_SC3) and so on.

    # All these values will be stored in heatmap_data

    def data_for_heatmap():
        for x in range(6):
            for y in range(6):
                euclidean_distance(normalized_lists[x], normalized_lists[y])

    # We are calling on data_for_heatmap(), to give us the complete
    # version of heatmap_data, with all 36 values appended. In order
    # for this to be mapped in a heatmap, it needs to be in a 6 by 6 array
    # so we are reshaping the data
    data_for_heatmap()
    heatmap_data = np.reshape(heatmap_data, (6, 6))

    # Now we are going to start plotting the heatmap using seaborn

    # Setting a mask because the upper right triangle of the
    # heatmap is redundant. This will block out the upper right triangle
    sns.set(style="white")
    mask = np.zeros_like(heatmap_data, dtype=np.bool)
    mask[np.triu_indices_from(mask)] = True

    ax = sns.heatmap(heatmap_data, mask=mask)

    row = ["SC1", "SC2", "SC3", "SC4", "SC5", "SC6"]
    column = ["SC1", "SC2", "SC3", "SC4", "SC5", "SC6"]

    ax.set_xticklabels(row)
    ax.set_yticklabels(column)
    ax.set_title("Heatmap of the samples")

    plt.xlabel("Samples")
    plt.ylabel("Samples")
    plt.show()

