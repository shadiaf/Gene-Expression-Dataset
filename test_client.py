import heatmap
import annotation
import volcano_plot
import scatter_and_violin_plots

# The outputs are already in the output_file. Running this client file will
# give you the same thing as the outputs in that file

# Will output a heatmap
print (heatmap.heatmap())

# Will output a new csv file: results_L-vs-NL_full_annotated
# It will show up within the Final Project folder
# This will take some time as it needs to query
print (annotation.annotation("results_L-vs-NL_full.csv"))

# Will output a volcano plot
print (volcano_plot.volcano_plot())

# Will output a scatter plot
print (scatter_and_violin_plots.scatterplot())

# Will output a violin plot
print (scatter_and_violin_plots.violin_plot())