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
    #matrix = np.nan_to_num(matrix,nan=2)
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
