######################################################################################################
import argparse
import os
from datetime import datetime
import numpy as np
import scipy.sparse as sp
import pandas as pd
import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
######################################################################################################

######################################################################################################
input_directory=''
os.chdir(input_directory)
file_prefix = ''
file_suffix = '.h5ad.gz' # The input file suffix to search for
sample_names = [] # The sample name strings to use for the plots
sample_indices = [] # Indexing information for the samples
# dataset_types = ['m6A_minus_strand', 'm6A_plus_strand', '5mC_minus_strand', '5mC_plus_strand', 'combined_minus_strand', 'combined_plus_strand']
# datasets_to_analyze = [0, 1]
filter_on_coordinates, lower_bound, upper_bound = [True, 0, 3000]
filter_columns_on_nan, position_nan_threshold = [True, 0.6] # Keep only positions that occur with less than the given nan frequency
filter_rows_on_nan, read_nan_threshold = [False, 0.025] # Keep only reads that have less than the given nan frequency
nan0_0minus1= False
fill_nans = False
filter_min_methylation, min_read_methylation = [False, 0] # Minimal methylation quantity in read to keep the read
show_basic_plot, save_basic_plot = [False, False]
n_kmeans_clusters,show_kmeans_plot, save_kmeans_plot = [5, False, False] # Number of k-means clusters to use, whether to show the plot, whether to save the plot
show_cluster_plot = False
show_hierchical_plot, save_hierarchical_plot = [False, False]
pairwise_vmax_value = 1.5 # Maximal pairwise value to use in the pairwise heatmap colorbar
show_pearson_correlation_plot, save_pearson_correlation_plot = [False, False]
show_pairwise_combinations_plot, save_pairwise_combinations_plot = [True, False]
plot_pairwise_single_position, pairwise_positions_to_plot, save_pairwise_single_position = [True, ['2787','2311', '2187'], True]
interaction= ['1-1','1-0','0-1','0-0'] # Pairwise combination interactions
pairwise_combinations = 0 # Index of the pairwise combination to analyze
show_scanpy_plots, save_scanpy_plots = [False, False]
umap_obs_to_plot = ['leiden', 'Sample', 'kmeans_labels_reordered']
n_neighbors = 10 # Number of neighbors for the neighborhood graph
leiden_resolution = 0.01 # Resolution value for leiden clustering. Higher value generates more clusters
######################################################################################################

######################################################################################################
# Get the current date
current_date = datetime.now()
# Format the date as a string
date_string = current_date.strftime("%Y%m%d")
date_string = date_string[2:]
def time_string(sep=':'):
    current_time = datetime.now()
    return current_time.strftime(f"%H{sep}%M{sep}%S")

def adata_to_df(adata):
    """
    Takes in an anndata object and returns a dataframe
    """
    # Extract the data matrix
    data_matrix = adata.X
    # Extract observation (read) annotations
    obs_df = adata.obs
    # Extract variable (position) annotations
    var_df = adata.var
    # Convert data matrix and annotations to pandas DataFrames
    df = pd.DataFrame(data_matrix, index=obs_df.index, columns=var_df.index)
    return df

def filter_adata_by_nan_proportion(adata, threshold, axis='obs'):
    if axis == 'obs':
        # Calculate the proportion of NaN values in each read
        nan_proportion = np.isnan(adata.X).mean(axis=1)
        # Filter reads to keep reads with less than a certain NaN proportion
        filtered_indices = np.where(nan_proportion <= threshold)[0]
        filtered_adata = adata[filtered_indices, :]
    elif axis == 'var':
        # Calculate the proportion of NaN values at a given position
        nan_proportion = np.isnan(adata.X).mean(axis=0)
        # Filter positions to keep positions with less than a certain NaN proportion
        filtered_indices = np.where(nan_proportion <= threshold)[0]
        filtered_adata = adata[:, filtered_indices]
    else:
        raise ValueError("Axis must be either 'obs' or 'var'")
    return filtered_adata

def calculate_correlation_matrix(df):
    # Initialize an empty matrix to store the correlation coefficients
    n_cols = len(df.columns)
    correlation_matrix = np.zeros((n_cols, n_cols), dtype=float)
    print('{0}: Calculating correlation matrix'.format(time_string()))
    # Iterate through each pair of columns
    for i in range(n_cols):
        print('{0}: Progress {1}%'.format(time_string(), 100*i/n_cols))
        for j in range(n_cols):
            # Calculate the correlation coefficient for the pair (i, j)
            correlation_matrix[i, j] = df.iloc[:, i].corr(df.iloc[:, j], method='pearson', min_periods=1)
    return correlation_matrix

def pair_heatmap_plot(matrix, title, xlabel, ylabel, save, save_name, cmap='coolwarm', colorbar_label=False, vmax=False):
    fig, axes = plt.subplots(figsize=(12, 6))
    if save:
        dpi = 3000
    else:
        dpi = 300
    fig.set_dpi(dpi)
    if vmax:
        sns.heatmap(matrix, annot=False, cmap=cmap, fmt='.2f', linewidths=0, cbar=True, ax=axes, vmax=vmax)
    else:
        sns.heatmap(matrix, annot=False, cmap=cmap, fmt='.2f', linewidths=0, cbar=True, ax=axes)
    axes.set_title('')
    axes.set_xlabel(xlabel)
    axes.set_ylabel(ylabel)
    axes.tick_params(axis='x', rotation=70)
    fig.suptitle(title)
    if colorbar_label:
        colorbar = axes.collections[0].colorbar
        colorbar.set_label(colorbar_label)
    else:
        pass
    # axes.set_xticks(np.arange(0.5, len(df_sorted.columns)))
    # axes.set_xticklabels(df_sorted.columns)
    # axes.set_yticks(np.arange(0.5, len(df_sorted.columns)))
    # axes.set_yticklabels(df_sorted.columns)
    plt.tight_layout()
    if save:
        fig.savefig(save_name, bbox_inches='tight', pad_inches=0.1)
    else:
        plt.show()
    return

def calculate_combinations_matrix_numpy(df):
    matrix = df.to_numpy()
    print(matrix)
    matrix = np.nan_to_num(matrix,nan=2)
    n_cols = matrix.shape[1]
    combinations_matrix = np.zeros((n_cols, n_cols, 4))
    print('{0}: Calculating combination matrix'.format(time_string()))
    # Iterate through each pair of columns
    for i in range(n_cols):
        for j in range(n_cols):
            # Count the occurrences of '1-1', '1-0', '0-1', and '0-0' for the pair (i, j)
            count_11 = np.nansum((matrix[:, i] == 1) & (matrix[:, j] == 1))
            count_10 = np.nansum((matrix[:, i] == 1) & (matrix[:, j] == 0))
            count_01 = np.nansum((matrix[:, i] == 0) & (matrix[:, j] == 1))
            count_00 = np.nansum((matrix[:, i] == 0) & (matrix[:, j] == 0))
            
            # Store the counts in the combinations matrix
            combinations_matrix[i, j, 0] = count_11/matrix.shape[0]
            combinations_matrix[i, j, 1] = count_10/matrix.shape[0]
            combinations_matrix[i, j, 2] = count_01/matrix.shape[0]
            combinations_matrix[i, j, 3] = count_00/matrix.shape[0]
    
    return combinations_matrix

######################################################################################################

###################################################
### Get input h5ad file name and load it into adata###
# List all files in the directory
files = os.listdir(os.getcwd())
# get current working directory
cwd = os.getcwd()
# Filter file names that contain the search string in their filename and keep them in a list
hdfs = [hdf for hdf in files if file_suffix in hdf and file_prefix in hdf]
# Sort file list by names and print the list of file names
hdfs.sort()
print('{0}: Sample files found: {1}'.format(time_string(), hdfs))
print('{0}: Reading in {1} h5ad'.format(time_string(), hdfs[0]))
adata = ad.read_h5ad(hdfs[0])
print('{0}: successfully loaded h5ad into memory as adata'.format(time_string()))
###################################################

# Iterate over Datasets
dataset_types = ['6B6_5mC_top', '6BALB_cJ_5mC_top']
datasets_to_analyze = [0, 1]
print('{0}: Begin iterating over datasets: {1}'.format(time_string(), [dataset_types[i] for i in datasets_to_analyze]))
for i, dataset in enumerate(dataset_types):
    # Analyze the Datasets of interest one at a time
    if i in datasets_to_analyze:
        print('{0}: Iterating over {1} Dataset'.format(time_string(), dataset))
        # Subset and keep only the Dataset of interest
        category_of_interest = 'Reference'
        value_to_keep = dataset
        print('{0}: Initial adata {1}'.format(time_string(), adata))
        print('{0}: Subsetting adata to keep {1} Dataset'.format(time_string(), dataset))
        adata_subset = adata[adata.obs[category_of_interest] == value_to_keep]
        print('{0}: adata_subset {1}'.format(time_string(), adata_subset))

        if filter_on_coordinates:
            # Subset and keep positions within the specified coordinate range
            position_list = list(range(lower_bound, upper_bound + 1))
            position_list = [str(pos) for pos in position_list]
            position_set = set(position_list)
            print('{0}: Subsetting adata to keep data between coordinates {1} and {2}'.format(time_string(), lower_bound, upper_bound))
            adata_subset = adata_subset[:, adata_subset.var_names.isin(position_set)]
            print('{0}: adata_subset {1}'.format(time_string(), adata_subset))

        if filter_columns_on_nan:
            # Remove positions that have over a certain frequency of nan values
            print('{0}: Subsetting adata to keep positions with atleast {1} percent of data coverage for coordinates in the defined window'.format(time_string(), position_nan_threshold*100))
            adata_subset = filter_adata_by_nan_proportion(adata_subset, position_nan_threshold, axis='var')
            print('{0}: adata_subset {1}'.format(time_string(), adata_subset))

        if filter_rows_on_nan:
            # Remove reads that have over a certain frequency of nan values
            print('{0}: Subsetting adata to keep reads with atleast {1} percent of data for coordinates in the defined window'.format(time_string(), read_nan_threshold*100))
            adata_subset = filter_adata_by_nan_proportion(adata_subset, read_nan_threshold, axis='obs')
            print('{0}: adata_subset {1}'.format(time_string(), adata_subset))

        # Save a copy of the anndata object without filling
        adata_unfilled = adata_subset.copy()
        print(adata_unfilled.var_names)

        # Forward and back fill each read to replace nan values with nearest positional methylation state value
        if fill_nans:
            print('{0}: Filling nan values in adata_subset with nearest positional methylation state values'.format(time_string()))
            df = pd.DataFrame(adata_subset.X)
            df = df.ffill(axis=1).bfill(axis=1)
            adata_subset.X = df.values
        
        if nan0_0minus1:
            print('{0}: Changing non methlylated values (ie 0) to -1 and filling nan values in adata_subset with 0 and'.format(time_string()))
            old_value, new_value = [0, -1]
            df = adata_to_df(adata_subset)
            df = df.replace(old_value, new_value)
            old_value, new_value = [np.nan, 0]
            df = df.replace(old_value, new_value)
            adata_subset.X = df.values
        
        # Calculate the average methylation density within each read
        print('{0}: Calculating the mean methlyation state of each read'.format(time_string()))
        read_methylation_mean_values = np.nanmean(adata_subset.X, axis=1)
        adata_subset.obs['mean_methylation'] = read_methylation_mean_values

        if filter_min_methylation:
            # Filter reads to keep reads with at least a certain methylation frequency
            filtered_indices = np.where(read_methylation_mean_values >= min_read_methylation)[0]
            print('{0}: Subsetting adata_subset to keep reads with at least {1} percent methylation'.format(time_string(), min_read_methylation*100))
            adata_subset = adata_subset[filtered_indices, :]
            print('{0}: adata_subset {1}'.format(time_string(), adata_subset))       

        sample_set = set(sample for sample in list(adata_subset.obs.Sample))
        sample_list = list(sample_set)
        sample_list.sort()

        if show_basic_plot:
        # Plotting basic heatmap
            for j, sample in enumerate(sample_list):
                if sample in sample_set:
                    print('{0}: sorting {1} data for sample {2}'.format(time_string(), dataset, sample))
                    temp_adata_subset = adata_subset[adata_subset.obs['Sample'] == str(sample)]
                    df = adata_to_df(temp_adata_subset)
                    df['row_sum'] = df.sum(axis=1)
                    # Sort the DataFrame by the row sums in descending order
                    df_sorted = df.sort_values(by='row_sum', ascending=False)
                    # Drop the temporary 'row_sum' column
                    df_sorted = df_sorted.drop(columns=['row_sum'])
                    df_sorted = df_sorted.reset_index(drop=True)
                    # Create heatmap
                    title = '{0} binarized methylation data for {1} from positions {2} to {3}'.format(sample_names[j], dataset, lower_bound, upper_bound)
                    cmap = 'Blues'
                    xlabel = 'Position'
                    ylabel = 'Read'
                    colorbar_label = 'Methylation'
                    save_name = cwd + '/{0}_{1} sorted binarized methylation data for {2} {3}_{4}'.format(date_string, sample_names[j], dataset, lower_bound, upper_bound)
                    print('{0}: Plotting the sorted SMF data'.format(time_string()))
                    pair_heatmap_plot(df_sorted, title, xlabel, ylabel, save_basic_plot, save_name, cmap, colorbar_label)

        if show_kmeans_plot:
            from sklearn.cluster import KMeans
            # Create a KMeans object with the desired number of clusters
            print('{0}: Calculating kmeans clustering of adata_subset for {1} clusters'.format(time_string(), n_kmeans_clusters))
            kmeans = KMeans(n_clusters=n_kmeans_clusters)
            # Fit the KMeans model to the data
            kmeans.fit(adata_subset.X)
            # Get the cluster labels for each data point
            cluster_labels = kmeans.labels_
            # Add the kmeans cluster data as an observation to the anndata object
            adata_subset.obs['kmeans_labels'] = cluster_labels.astype(str)
            # Calculate the mean of each observation categoty of each cluster
            cluster_means = adata_subset.obs.groupby('kmeans_labels').mean()
            print(cluster_means)
            # Sort the cluster indices by mean methylation value
            sorted_clusters = cluster_means.sort_values(by='mean_methylation', ascending=False).index
            # Create a mapping of the old cluster values to the new cluster values
            sorted_cluster_mapping = {old: new for new, old in enumerate(sorted_clusters)}
            # Apply the mapping to create a new observation value: kmeans_labels_reordered
            adata_subset.obs['kmeans_labels_reordered'] = adata_subset.obs['kmeans_labels'].map(sorted_cluster_mapping)
            adata_subset.obs['kmeans_labels_reordered'] = adata_subset.obs['kmeans_labels_reordered'].astype(int)
            # Sort the adata_subset object by the kmeans_labels_reordered observation
            print('{0}: Sorting data on kmeans_labels_reordered observation'.format(time_string()))
            sorted_indices = adata_subset.obs.sort_values(by='kmeans_labels_reordered').index
            adata_subset = adata_subset[sorted_indices].copy()
            #adata_subset.obs = adata_subset.obs.sort_values(by='kmeans_labels_reordered')
            for j, sample in enumerate(sample_list):
                temp_adata_subset = adata_subset[adata_subset.obs['Sample'] == str(sample)]
                df = adata_to_df(temp_adata_subset)
                df = df.reset_index(drop=True)
                # Create heatmap
                title = '{0} Kmeans Clustered (n = {1}) reads {2} from positions {3} to {4}'.format(sample_names[j], n_kmeans_clusters ,dataset, lower_bound, upper_bound)
                cmap = 'coolwarm'  #'Blues'
                xlabel = 'Position'
                ylabel = 'Read'
                colorbar_label = 'Methylation state: Blue=Unmethylated, Red=Methylated'
                save_name = cwd + '/{0}_{1}_{2} clustered binarized methylation data for {3} {4}_{5}'.format(date_string, time_string(sep=''), sample_names[j], dataset, lower_bound, upper_bound)
                print('{0}: Plotting the k-means cluster sorted SMF data'.format(time_string()))
                pair_heatmap_plot(df, title, xlabel, ylabel, save_kmeans_plot, save_name, cmap, colorbar_label)

        if show_pearson_correlation_plot:
            for j, sample in enumerate(sample_list):
                temp_adata_subset = adata_subset[adata_subset.obs['Sample'] == str(sample)]
                df = adata_to_df(temp_adata_subset)
                df = df.reset_index(drop=True)
                #Calculate the correlation matrix
                result_matrix = calculate_correlation_matrix(df)
                result_df = pd.DataFrame(result_matrix)
                result_df.columns = df.columns
                result_df = result_df.set_index(df.columns)
                xlabel = 'Position'
                ylabel = 'Position'
                colorbar_label = 'Pearson correlation'
                title = 'Pearson correlation for {0} {1}'.format(sample_names[j], dataset)
                save_name = cwd + '/{0}_{1} Pearson correlation for {2} {3}_{4}'.format(date_string, sample_names[j],  dataset, lower_bound, upper_bound)
                pair_heatmap_plot(result_df, title, xlabel, ylabel, save_pearson_correlation_plot, save_name, colorbar_label=colorbar_label)

        if show_pairwise_combinations_plot:
            for j, sample in enumerate(sample_list):
                temp_adata_subset = adata_unfilled[adata_unfilled.obs['Sample'] == str(sample)]
                df = adata_to_df(temp_adata_subset)
                df = df.reset_index(drop=True)
                result_matrix = calculate_combinations_matrix_numpy(df)
                result_matrix_single = result_matrix[:,:,pairwise_combinations]
                result_df = pd.DataFrame(result_matrix_single)
                result_df.columns = df.columns
                result_df = result_df.set_index(df.columns)
                interaction_II = ['Methylated_Methylated', 'Methylated_Unmethylated', 'Unmethylated_Methylated', 'Unmethylated_Unmethylated']
                xlabel = 'Position'
                ylabel = 'Position'
                colorbar_label = 'Pairwise Probabilities' # On diagnal is the probability of a square being in the first state. Off the diagnal is pairwise probability
                title = '{0} Pairwise {1} linkage {2}'.format(sample_names[j], interaction_II[pairwise_combinations], dataset)
                save_name = cwd + '/{0}_{1} Pairwise {2} linkage {3} {4}_{5}'.format(date_string, sample_names[j], interaction_II[pairwise_combinations], dataset, lower_bound, upper_bound)
                pair_heatmap_plot(result_df, title, xlabel, ylabel, save_pairwise_combinations_plot, save_name, cmap='viridis', colorbar_label=colorbar_label)
                if plot_pairwise_single_position:
                    plt.figure(figsize=(10, 6))
                    plt.title(f'{sample_names[j]} Pairwise {interaction_II[pairwise_combinations]} probabilities for {dataset}')
                    #plt.ylim(-0.05, 1)  # Adjust y-axis limits
                    for position in pairwise_positions_to_plot:
                        if position in result_df.index:
                            row_to_plot = result_df.loc[position]
                            # Plot the selected row
                            row_df = pd.DataFrame({'Positions': row_to_plot.index, 'Pairwise probability': row_to_plot.values})
                            ax = sns.lineplot(data=row_df, x='Positions', y='Pairwise probability', marker='o', markersize=5, label=f'{position}')
                            num_ticks = 20  # Number of x-ticks you want to display
                            tick_positions = range(0, len(row_df), len(row_df) // num_ticks)  # Calculate tick positions
                            tick_labels = row_df['Positions'][::len(row_df) // num_ticks]  # Get labels for tick positions
                            ax.set_xticks(tick_positions)
                            ax.set_xticklabels(tick_labels, rotation=45)  # Rotate labels if necessary
                            ax.spines['right'].set_visible(False)
                            ax.spines['top'].set_visible(False)
                            ax.legend(loc='upper left')
                    if save_pairwise_single_position:
                        save_name = cwd + '/{0}_{1} Pairwise {2} probability plot for {3} {4}_{5}'.format(date_string, sample_names[j], interaction_II[pairwise_combinations], dataset, lower_bound, upper_bound)
                        plt.savefig(save_name)
                    else:
                        plt.show()


        if show_scanpy_plots:
            import scanpy as sc
            # set scanpy parameters
            sc.settings.verbosity = 0             # verbosity: errors (0), warnings (1), info (2), hints (3)
            sc.logging.print_header()
            sc.settings.set_figure_params(dpi=300, facecolor='white')

            # Calculate QC metrics     
            sc.pp.calculate_qc_metrics(adata_subset, percent_top=None, log1p=False, inplace=True) # Calculates QC metrics

            # Logarithmize the data and save raw counts
            adata_subset.layers["counts"] = adata_subset.X.copy() # preserve raw counts
            sc.pp.log1p(adata_subset) # logarithmize the data
            adata_subset.raw = adata_subset # used for downstream visualization

            ## Identify highly variable positions and save the name of the positions##
            sc.pp.highly_variable_genes(adata_subset, flavor='seurat', layer='counts', n_top_genes=50)
            highly_variable_gene_names = adata_subset.var_names[adata_subset.var['highly_variable']]
            #sc.pl.highly_variable_genes(adata_subset) # Visualize highly variable positions

            ## Scale the data and perform PCA ##
            sc.pp.scale(adata_subset)
            sc.tl.pca(adata_subset, svd_solver='arpack')
            #sc.pl.pca_scatter(adata_subset)

            ## compute the neighborhood graph of cells using the PCA representation of the data matrix ##
            n_pcs = 20
            if len(adata_subset.var) < n_pcs:
                n_pcs = len(adata_subset.var) // 2
            if len(adata_subset.obs) < n_neighbors*2:
                n_neighbors = len(adata_subset.obs) // 2
            sc.pp.neighbors(adata_subset, n_neighbors=n_neighbors, n_pcs=n_pcs, use_rep='X_pca')

            ## Leiden graph-clustering method ##
            sc.tl.leiden(adata_subset, resolution=leiden_resolution)

            ## Compute UMAP ##
            sc.tl.umap(adata_subset)
            if save_scanpy_plots:
                save_name = cwd + '/{0}_{1} UMAP {2} {3} {4}_{5}'.format(date_string, sample_names[j], interaction_II[pairwise_combinations], dataset, lower_bound, upper_bound)
                sc.pl.umap(adata_subset, color=umap_obs_to_plot, wspace=0.5, save=save_name)
            else: 
                sc.pl.umap(adata_subset, color=umap_obs_to_plot, wspace=0.5)
            #sc.tl.dendrogram(adata_subset, groupby='Dataset')
            #sc.pl.dendrogram(adata_subset, groupby='Dataset') 
            #sc.pl.matrixplot(adata_subset, var_names=highly_variable_gene_names ,groupby='Dataset', dendrogram=True)

        if show_hierchical_plot:
            from scipy.cluster.hierarchy import linkage, dendrogram
            for j, sample in enumerate(sample_list):
                temp_adata_subset = adata_subset[adata_subset.obs['Sample'] == str(sample)]
                df = adata_to_df(temp_adata_subset)
                df = df.reset_index(drop=True)
                clustergrid = sns.clustermap(df, cmap='coolwarm', method='average', metric='euclidean', row_cluster=True, col_cluster=False, figsize=(10, 6), cbar=(False))
                clustergrid.fig.suptitle('{0} Hierarchical clustering of methylated reads for {1} from positions {2} to {3}'.format(sample_names[j], dataset, lower_bound, upper_bound), y=0.9)
                clustergrid.ax_heatmap.set_xlabel('Positions')
                clustergrid.ax_heatmap.set_ylabel('Reads')
                plt.subplots_adjust(left=0.1, right=0.9, top=0.8, bottom=0.2)   
                if clustergrid.cax is not None:
                    clustergrid.cax.set_visible(False)         
                if save_hierarchical_plot:
                    save_name=cwd + '/{0}_{1} Hierarchical heatmap for {2} {3}_{4}'.format(date_string, sample_names[j],  dataset, lower_bound, upper_bound)
                    clustergrid.savefig(save_name, bbox_inches='tight', pad_inches=0.1)
                else:
                    plt.show()
               

######################################################################################################
