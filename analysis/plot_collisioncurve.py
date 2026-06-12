import argparse
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde
import matplotlib.patches as mpatches

# Define subplot target dimensions
w_subplot = 11
h_subplot = 7


# Configurations
mutation_model = {0:"MUTATION_MODEL_SIMPLE_SNP_ONLY",1:"MUTATION_MODEL_GEOMETRIC_MUTATOR"}
mutation_expressions = {0:"BAL",1:"SUB-ONLY",2:"DEL-LITE",3:"INS-LITE",4:"SUB-LITE"}

# distance_metric = {0:"None",1:"Hamming Similarity",2:"Jaccard Similarity",3:"Cosine Similarity",4:"Angular Similarity",5:"Edit Similarity"}

def read_collision_data_complete(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

    sections = {}
    data_rows = []
    similarity_values = []
    snp_rate = []
    delRate = []
    insmean = []
    insrate = []
    stayRate = []
    rand_base_params = []   # This is equal to snp_rate for sub-only model and equal to g_mean from geo-ins model

    current_AND_param = ""
    current_OR_param = ""
    # Base metadata that applies to all rows
    base_metadata = {}
    hash_parameters = {}  # Dictionary to store dynamic hash parameters

    for line in lines:
        line = line.strip()
        if not line:
            continue

        if line.startswith(':1:'):
            base_metadata['test_name'] = line[3:].strip()
        elif line.startswith(':2:'):
            base_metadata['column_headers'] = line[3:].strip().split(',')
        elif line.startswith(':3:'):
            line_content = line[3:].strip()
            line_parts = line_content.split(',')
            base_metadata['Hashname'] = line_parts[0].strip()
            base_metadata['SequenceLength'] = int(line_parts[1].strip())
            base_metadata['DistanceMetric'] = str(line_parts[2].strip())
            base_metadata['MutationModel'] = int(line_parts[3].strip())
            base_metadata['MutationExpression'] = int(line_parts[4].strip())
        elif line.startswith(':4.1:'):  # Parameter names
            param_names = line[5:].strip().split(',')
            hash_parameters['param_names'] = [name.strip() for name in param_names]
        elif line.startswith(':4.2:'):  # Parameter descriptions
            param_descs = line[5:].strip().split(',')
            hash_parameters['param_descs'] = [desc.strip() for desc in param_descs]
        elif line.startswith(':4.3:'):  # Parameter values
            param_values = line[5:].strip().split(',')
            # Handle missing or empty parameter values
            hash_parameters['param_values'] = [
                float(value.strip()) if value.strip() else 0.0 for value in param_values
            ]
        elif line.startswith(':5:'):
            line_content = line[3:].strip()
            similarity_values = [float(x.strip()) for x in line_content.split(',')]
        elif line.startswith(':6:'):
            line_content = line[3:].strip()
            snp_rate = [float(x.strip()) for x in line_content.split(',')]
        elif line.startswith(':7:'):
            line_content = line[3:].strip()
            delRate = [float(x.strip()) for x in line_content.split(',')]
        elif line.startswith(':8:'):
            line_content = line[3:].strip()
            insmean = [float(x.strip()) for x in line_content.split(',')]
        elif line.startswith(':9:'):
            line_content = line[3:].strip()
            stayRate = [float(x.strip()) for x in line_content.split(',')]
        elif line.startswith(':10:'):
            # Skip or store if needed
            pass
        elif line.startswith(':11:'):
            line_content = line[4:].strip()
            # Store current AND/OR params for the upcoming :12: line
            amplification_params = line_content.split(',')
            current_AND_param = amplification_params[0]
            current_OR_param = amplification_params[1]
        elif line.startswith(":13:"):
            line_content = line[4:].strip()
            rand_base_params = [float(x.strip()) for x in line_content.split(',')]
        elif line.startswith(":14:"):
            line_content = line[4:].strip()
            insrate = [float(x.strip()) for x in line_content.split(',')]
        elif line.startswith(":12:"):
            line_content = line[4:].strip()
            collision_rates = [float(x.strip()) for x in line_content.split(',')]
            # Create a new data row for this set of collision rates
            row_data = {
                'test_name': base_metadata.get('test_name', ''),
                'hashname': base_metadata.get('Hashname', ''),
                'sequencelength': base_metadata.get('SequenceLength', 0),
                'tokenlength': base_metadata.get('TokenLength', 0),
                'DistanceMetric': base_metadata.get('DistanceMetric', ''),
                'MutationModel': base_metadata.get('MutationModel', 0),
                'MutationExpression': base_metadata.get('MutationExpression', 0),
                'AND_param': current_AND_param,
                'OR_param': current_OR_param,
                'similarity_values': similarity_values.copy(),
                'snp_rate': snp_rate.copy(),
                'delRate': delRate.copy(),
                'insmean': insmean.copy(),
                'insrate': insrate.copy(),
                'stayRate': stayRate.copy(),
                'collision_rates': collision_rates.copy(),
                'rand_base_params': rand_base_params,
            }
            # Add hash parameters dynamically
            if 'param_names' in hash_parameters and 'param_values' in hash_parameters:
                for name, value in zip(hash_parameters['param_names'], hash_parameters['param_values']):
                    row_data[name] = value
                    # print(value)

            # Add a new column to store parameter names
            row_data['parameter_columns'] = hash_parameters.get('param_names', [])

            data_rows.append(row_data)

    return sections, pd.DataFrame(data_rows)


# For verification
def plot_histogram(values, title="Histogram", xlabel="Value", ylabel="Count"):
    fig, ax = plt.subplots()

    ax.hist(
        values,
        bins=26,          # required
        alpha=0.85,
        edgecolor="black", # keep if your current style uses outlined bars
        linewidth=0.4
    )

    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    # keep grid settings if your current style uses them
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    return fig, ax

def plot_distribution_curve_fixed_range_100(
        values,
        ax=None,
        label=None,
        color=None,
        alpha=0.22,      # translucent fill
        bw_method=0.2,   # smoothing control
):
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.figure

    # keep the fixed domain behavior (0..100)
    v = np.asarray(values, dtype=float)
    v = v[(v >= 0) & (v <= 100)]
    if len(v) < 2:
        return fig, ax  # not enough data for KDE

    x = np.linspace(0, 100, 500)   # fixed x-range, never changes
    kde = gaussian_kde(v, bw_method=bw_method)
    y = kde(x)

    ax.plot(x, y, color=color, linewidth=2, label=label)
    ax.fill_between(x, 0, y, color=color, alpha=alpha)

    ax.set_xlim(0, 100)
    ax.set_xlabel("Value")
    ax.set_ylabel("Density")
    ax.grid(True, alpha=0.25)
    return fig, ax

def plot_distribution_curve(
        values,
        ax=None,
        label=None,
        color=None,
        alpha=0.22,
        bw_method=0.1,
        linestyle='-',   # Default: solid line
        dashes=None,     # Optional: Tuple for custom spacing (e.g., (5, 5))
        linewidth=2      # Allow changing thickness
):
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.figure

    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if len(v) < 2:
        return fig, ax

    # Dynamic x-domain logic
    xmin, xmax = v.min(), v.max()
    if xmin == xmax:
        pad = 1e-6
        xmin -= pad
        xmax += pad
    else:
        pad = 0.05 * (xmax - xmin)
        xmin -= pad
        xmax += pad

    x = np.linspace(xmin, xmax, 500)
    kde = gaussian_kde(v, bw_method=bw_method)
    y = kde(x)

    # 1. Create the line object
    (line,) = ax.plot(x, y, color=color, linewidth=linewidth, linestyle=linestyle, label=label)

    # 2. Apply custom dashes if provided (overrides linestyle)
    # Format: (offset, (on_len, off_len, ...))
    if dashes:
        line.set_dashes(dashes)

    ax.fill_between(x, 0, y, color=color, alpha=alpha)

    ax.set_xlabel("Value")
    ax.set_ylabel("Density")
    ax.grid(True, alpha=0.25)

    return fig, ax



# Average per bin curve
def plot_binned_average_curve(
        row,
        ax=None,
        color=None,
        markersize=6,
        linewidth=1,
        alpha=1,
):
    """Add a single bin-averaged collision curve to ax. Returns (fig, ax)."""
    if ax is None:
        # fig, ax = plt.subplots(figsize=(11, 7))
        fig, ax = plt.subplots()
    else:
        fig = ax.figure


    x_values = np.array(row['similarity_values'])
    y_values = np.array(row['collision_rates'])

    # Define bins
    bin_edges = np.arange(0, 1.01, 0.01)  # 0.00, 0.01, 0.02, ..., 1.00
    num_bins = len(bin_edges) - 1  # 100 bins
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2  # Proper bin centers


    # Compute mean for each bin
    bin_means = []
    for i in range(num_bins):
        mask = (x_values >= bin_edges[i]) & (x_values < bin_edges[i + 1])
        if mask.sum() > 0:
            bin_means.append(np.mean(y_values[mask]))
        else:
            bin_means.append(np.nan)

    bin_means = np.array(bin_means)
    valid = ~np.isnan(bin_means)

    mutation_model_used = ""
    if row['MutationModel'] == 0:
        mutation_model_used = "SUB-ONLY"

    if row['MutationModel'] == 1:
        mutation_model_used = "SubIndel-" + mutation_expressions[row["MutationExpression"]]

    # Build label
    if row['AND_param'] == '1' and row['OR_param'] == '1':
        param_labels = ", ".join([f"{name}={row[name]}" for name in row['parameter_columns'] if row[name] != 0.0])
        if param_labels:
            label = f"{row['hashname']} (L={row['sequencelength']}, {param_labels})[{mutation_model_used}]"
        else:
            label = f"{row['hashname']} (L={row['sequencelength']})[{mutation_model_used}]"
    else:
        param_labels = ", ".join([f"{name}={row[name]}" for name in row['parameter_columns'] if row[name] != 0.0])
        if param_labels:
            label = f"{row['hashname']} (L={row['sequencelength']}, {param_labels}, b={row['AND_param']}, r={row['OR_param']})[{mutation_model_used}]"
        else:
            label = f"{row['hashname']} (L={row['sequencelength']}, b={row['AND_param']}, r={row['OR_param']})[{mutation_model_used}]"
            
    ax.plot(
        bin_centers[valid],
        bin_means[valid],
        'o',
        linestyle='solid',
        color=color,
        alpha=alpha,
        markersize=markersize,
        linewidth=linewidth,
        label=label,
    )

    return fig, ax

def finalize_binned_average_plot(
        ax,
        similarity_label="Similarity",
        s1=None,
        s2=None,
        savename=None,
):
    # """Add ideal line, vertical markers, grid, legend, and optionally save."""
    # ax.plot([0, 1], [0, 1], color="black", linewidth=1, alpha=1, label="Theoretical similarity estimator")

    # if s1 is not None:
    #     ax.axvline(x=s1, color='black', label=f'$s_1 = {s1}$', linewidth=0.7, linestyle='dashed')
    # if s2 is not None:
    #     ax.axvline(x=s2, color='black', label=f'$s_2 = {s2}$', linewidth=0.7, linestyle='dashed')

    ax.set_xlim(-0.01, 1.01)
    ax.set_ylim(-0.01, 1.01)
    ax.minorticks_on()
    ax.grid(visible=True, which="both", axis="both", alpha=0.5)
    ax.set_xlabel(similarity_label, fontsize=22)
    ax.set_ylabel('Average Collision Rate per bin', fontsize=22)
    ax.tick_params(axis='both', labelsize=20)  # Add this line
    ax.legend(fontsize=18)
    ax.figure.tight_layout()

    if savename:
        ax.figure.savefig(savename, dpi=600, bbox_inches='tight')

# # Color Scatter Plot for Subs-Only Model
# def plot_color_scatter(
#         df,
#         n_cols=2,
#         xlabel="Similarity",
#         save_filename=None,
#         show_plot=True
# ):
#     """
#     Generates a grid of scatter plots for hash collision rates vs similarity.

#     Parameters:
#     - df: DataFrame containing the data. Must contain columns:
#           'similarity_values', 'collision_rates', 'snp_rate', 'hashname',
#           'sequencelength', 'tokenlength', and optionally 'SubseqHash_k', 'SubseqHash_d'.
#     - n_cols: Number of columns in the subplot grid.
#     - xlabel: Label for the X-axis (variable 'similarity' in original snippet).
#     - save_filename: If provided, saves the figure to this path (e.g., "my_plot").
#                      Suffixes like _snpcolor_multiplot.png will be appended automatically
#                      if you want to match the original style, or provide the full name.
#     """

#     n_rows = len(df)
#     # Calculate required grid rows
#     n_grid_rows = (n_rows + n_cols - 1) // n_cols

#     # Create figure
#     # Note: Adjust figsize calculation if you want a fixed size rather than dynamic
#     # fig, axes = plt.subplots(nrows=n_grid_rows, ncols=n_cols, figsize=(28, 7 * n_grid_rows))
#     fig, axes = plt.subplots(
#         nrows=n_grid_rows, 
#         ncols=n_cols, 
#         figsize=(w_subplot * n_cols, h_subplot * n_grid_rows),
#         constrained_layout=True 
#     )

#     # Flatten axes for easy 1D iteration.
#     # Handle edge case where n_rows=1 (subplots returns a single Axes object, not array)
#     if n_rows == 1:
#         axes = np.array([axes])
#     axes = axes.flatten()

#     # Iterate using enumerate to get a safe 0-based index 'i' for axes
#     for i, (idx, row) in enumerate(df.iterrows()):
#         ax = axes[i]

#         # Extract data (assuming cells contain lists/arrays)
#         x_values = np.array(row['similarity_values'])
#         y_values = np.array(row['collision_rates'])
#         label_values = np.array(row['snp_rate'])

#         # Determine Label based on hashname
#         hash_name = row.get('hashname', 'Unknown')
#         seq_len = row.get('sequencelength', '?')


#         mutation_model_used = ""
#         if row['MutationModel'] == 0:
#             mutation_model_used = "SUB-ONLY"

#         if row['MutationModel'] == 1:
#             mutation_model_used = "SubIndel-" + mutation_expressions[row["MutationExpression"]]

        # # Build label
        # if row['AND_param'] == '1' and row['OR_param'] == '1':
        #     param_labels = ", ".join([f"{name}={row[name]}" for name in row['parameter_columns'] if row[name] != 0.0])
        #     if param_labels:
        #         label = f"{row['hashname']} (L={row['sequencelength']}, {param_labels})[{mutation_model_used}]"
        #     else:
        #         label = f"{row['hashname']} (L={row['sequencelength']})[{mutation_model_used}]"
        # else:
        #     param_labels = ", ".join([f"{name}={row[name]}" for name in row['parameter_columns'] if row[name] != 0.0])
        #     if param_labels:
        #         label = f"{row['hashname']} (L={row['sequencelength']}, {param_labels}, b={row['AND_param']}, r={row['OR_param']})[{mutation_model_used}]"
        #     else:
        #         label = f"{row['hashname']} (L={row['sequencelength']}, b={row['AND_param']}, r={row['OR_param']})[{mutation_model_used}]"

#         # if hash_name == "SubseqHash-64":
#         #     k_val = row.get('SubseqHash_k', '?')
#         #     d_val = row.get('SubseqHash_d', '?')
#         #     label = f"SubseqHash (L={seq_len}, k={k_val}, d={d_val})[{mutation_model_used}]"
#         # else:
#         #     token_len = row.get('tokenlength', '?')
#         #     label = f"{hash_name} (L={seq_len}, k={token_len})[{mutation_model_used}]"

#         # Scatter plot with colorbar
#         scatter = ax.scatter(x_values, y_values, c=label_values, cmap='viridis',
#                              alpha=0.3, edgecolors='black', linewidth=0.3, s=25)

#         # Add colorbar for each subplot
#         cbar = plt.colorbar(scatter, ax=ax)
#         cbar.set_label('SNP Rate', fontsize=14)

#         # Axis settings
#         ax.set_xlim(0, 1.01)
#         ax.set_ylim(0, 1.01)
#         ax.set_xlabel(xlabel, fontsize=22)
#         ax.set_ylabel('Collision Rates', fontsize=22)
#         ax.grid(True, alpha=0.5)

#         # Legend (using the generated label)
#         # We pass [scatter] to legend to ensure the label is associated with the points,
#         # though ax.legend([label]) works, passing the handle is cleaner.
#         ax.legend(handles=[scatter], labels=[label], loc='upper left', fontsize=18)
#         ax.tick_params(axis='both', labelsize=14)

#     # Hide unused subplots (if n_rows < total grid slots)
#     for j in range(n_rows, len(axes)):
#         axes[j].set_visible(False)

#     # plt.tight_layout()

#     # name = f"{row['hashname']}"
#     # plt.savefig(f"{name}_snpcolor_multiplot.png", dpi=600)

#     # Save logic
#     if save_filename:
#         # If user didn't provide an extension, add the style from the original snippet
#         if not save_filename.endswith('.png'):
#             final_name = f"{save_filename}_snpcolor_multiplot.png"
#         else:
#             final_name = save_filename

#         plt.savefig(final_name, dpi=600)
#         print(f"Plot saved to {final_name}")

#     if show_plot:
#         plt.show()

# No Scatter Plot for Subs-Only Model
def plot_monocolor_scatter(
        df,
        n_cols=2,
        xlabel="Similarity",
        save_filename=None,
        show_plot=True
):
    """
    Generates a grid of scatter plots for hash collision rates vs similarity.

    Parameters:
    - df: DataFrame containing the data. Must contain columns:
          'similarity_values', 'collision_rates', 'snp_rate', 'hashname',
          'sequencelength', 'tokenlength', and optionally 'SubseqHash_k', 'SubseqHash_d'.
    - n_cols: Number of columns in the subplot grid.
    - xlabel: Label for the X-axis (variable 'similarity' in original snippet).
    - save_filename: If provided, saves the figure to this path (e.g., "my_plot").
                     Suffixes like _snpcolor_multiplot.png will be appended automatically
                     if you want to match the original style, or provide the full name.
    """

    n_rows = len(df)
    # Calculate required grid rows
    n_grid_rows = (n_rows + n_cols - 1) // n_cols

    # Create figure
    # Note: Adjust figsize calculation if you want a fixed size rather than dynamic
    # fig, axes = plt.subplots(nrows=n_grid_rows, ncols=n_cols, figsize=(28, 7 * n_grid_rows))
    fig, axes = plt.subplots(
        nrows=n_grid_rows, 
        ncols=n_cols, 
        figsize=(w_subplot * n_cols, h_subplot * n_grid_rows),
        constrained_layout=True 
    )

    # Flatten axes for easy 1D iteration.
    # Handle edge case where n_rows=1 (subplots returns a single Axes object, not array)
    if n_rows == 1:
        axes = np.array([axes])
    axes = axes.flatten()

    # Iterate using enumerate to get a safe 0-based index 'i' for axes
    for i, (idx, row) in enumerate(df.iterrows()):
        ax = axes[i]

        # Extract data (assuming cells contain lists/arrays)
        x_values = np.array(row['similarity_values'])
        y_values = np.array(row['collision_rates'])
        # label_values = np.array(row['snp_rate'])

        # Determine Label based on hashname
        hash_name = row.get('hashname', 'Unknown')
        seq_len = row.get('sequencelength', '?')


        mutation_model_used = ""
        if row['MutationModel'] == 0:
            mutation_model_used = "SUB-ONLY"

        if row['MutationModel'] == 1:
            mutation_model_used = "SubIndel-" + mutation_expressions[row["MutationExpression"]]

        # Build label
        if row['AND_param'] == '1' and row['OR_param'] == '1':
            param_labels = ", ".join([f"{name}={row[name]}" for name in row['parameter_columns'] if row[name] != 0.0])
            if param_labels:
                label = f"{row['hashname']} (L={row['sequencelength']}, {param_labels})[{mutation_model_used}]"
            else:
                label = f"{row['hashname']} (L={row['sequencelength']})[{mutation_model_used}]"
        else:
            param_labels = ", ".join([f"{name}={row[name]}" for name in row['parameter_columns'] if row[name] != 0.0])
            if param_labels:
                label = f"{row['hashname']} (L={row['sequencelength']}, {param_labels}, b={row['AND_param']}, r={row['OR_param']})[{mutation_model_used}]"
            else:
                label = f"{row['hashname']} (L={row['sequencelength']}, b={row['AND_param']}, r={row['OR_param']})[{mutation_model_used}]"

        # Scatter plot with colorbar
        scatter = ax.scatter(x_values, y_values,
                             alpha=0.3, edgecolors='black', linewidth=0.3, s=25)

        # Add colorbar for each subplot
        # cbar = plt.colorbar(scatter, ax=ax)
        # cbar.set_label('SNP Rate', fontsize=14)

        # Axis settings
        ax.set_xlim(0, 1.01)
        ax.set_ylim(0, 1.01)
        ax.set_xlabel(xlabel, fontsize=22)
        ax.set_ylabel('Collision Rates', fontsize=22)
        ax.grid(True, alpha=0.5)

        ax.legend( labels=[label], loc='upper left', fontsize=18)
        ax.tick_params(axis='both', labelsize=20)

    # Hide unused subplots (if n_rows < total grid slots)
    for j in range(n_rows, len(axes)):
        axes[j].set_visible(False)

    # plt.tight_layout()

    # name = f"{row['hashname']}"
    # plt.savefig(f"{name}_snpcolor_multiplot.png", dpi=600)

    # Save logic
    if save_filename:
        # If user didn't provide an extension, add the style from the original snippet
        if not save_filename.endswith('.png'):
            final_name = f"{save_filename}_monocolor_multiplot.png"
        else:
            final_name = save_filename

        plt.savefig(final_name, dpi=600, bbox_inches='tight')
        print(f"Plot saved to {final_name}")

    if show_plot:
        plt.show()



# BOX PLOT per bin
def plot_hash_collision_boxplots(
        df,
        n_cols=2,
        xlabel="Similarity",
        save_filename=None,
        num_bins=26,
        show_plot=True
):
    """
    Generates a grid of box plots for hash collision rates binned by similarity.

    Parameters:
    - df: DataFrame containing 'similarity_values', 'collision_rates', and metadata columns.
    - n_cols: Number of columns in the subplot grid.
    - xlabel: Label for the X-axis (variable 'similarity' in original snippet).
    - save_filename: If provided, saves the figure.
    - num_bins: Number of bins to divide the similarity score (0 to 1).
    """

    n_rows = len(df)
    n_grid_rows = (n_rows + n_cols - 1) // n_cols

    # Create figure with dynamic height
    # fig, axes = plt.subplots(nrows=n_grid_rows, ncols=n_cols, figsize=(28, 7 * n_grid_rows))
    fig, axes = plt.subplots(
        nrows=n_grid_rows, 
        ncols=n_cols, 
        figsize=(w_subplot * n_cols, h_subplot * n_grid_rows),
        constrained_layout=True 
    )

    # Handle edge case where n_rows=1 (subplots returns single Axes, not array)
    if n_rows == 1:
        axes = np.array([axes])
    axes = axes.flatten()

    # Define bins
    # Note: bin_edges has num_bins+1 points.
    # The original logic used `int(sim_value * num_bins)` effectively mapping [0,1] to indices 0..num_bins-1
    bin_edges = np.linspace(0, 1.05, num_bins + 1)

    for i, (idx, row) in enumerate(df.iterrows()):
        ax = axes[i]

        # Extract data
        x_values = np.array(row['similarity_values'])
        y_values = np.array(row['collision_rates'])

        # --- Label Generation Logic ---
        hash_name = row.get('hashname', 'Unknown')
        seq_len = row.get('sequencelength', '?')
        and_param = str(row.get('AND_param', '?'))
        or_param = str(row.get('OR_param', '?'))

        # Determine if we show AND/OR params
        show_params = not (and_param == '1' and or_param == '1')

        mutation_model_used = ""
        if row['MutationModel'] == 0:
            mutation_model_used = "SUB-ONLY"

        if row['MutationModel'] == 1:
            mutation_model_used = "SubIndel-" + mutation_expressions[row["MutationExpression"]]

        param_labels = ", ".join([f"{name}={row[name]}" for name in row['parameter_columns'] if row[name] != 0.0])
        base_label = f"{hash_name} (L={seq_len}, {param_labels}"

        if show_params:
            label = f"{base_label}, b={and_param}, r={or_param})[{mutation_model_used}]"
        else:
            label = f"{base_label})[{mutation_model_used}]"

        # --- Binning Logic ---
        # Create empty lists for each bin
        binned_collision_rates = [[] for _ in range(num_bins)]

        for j in range(len(x_values)):
            sim_value = x_values[j]
            collision_value = y_values[j]

            # Map similarity 0.0-1.0 to bin index
            bin_idx = int(sim_value * num_bins)
            if bin_idx >= num_bins:
                bin_idx = num_bins - 1

            binned_collision_rates[bin_idx].append(collision_value)

        # Filter empty bins for plotting
        plot_data = []
        valid_positions = []
        valid_tick_labels = [] # To store labels for x-axis

        for k, bin_data in enumerate(binned_collision_rates):
            if len(bin_data) > 0:
                plot_data.append(bin_data)
                # Position is 1-based index corresponding to bin index
                valid_positions.append(k + 1)
                valid_tick_labels.append(f"{bin_edges[k]:.1f}")

        # --- Plotting ---
        bp = ax.boxplot(
            plot_data,
            positions=valid_positions, # Place boxes at specific x-coordinates
            widths=0.5,
            patch_artist=True,
            showfliers=False
        )

        # Style the boxes
        for box in bp['boxes']:
            box.set(facecolor='#C4C4C4', alpha=1, edgecolor='black', linewidth=1)
        for median in bp['medians']:
            median.set(color='black', linewidth=1.5)
        for whisker in bp['whiskers']:
            whisker.set(color='black', linewidth=1)
        for cap in bp['caps']:
            cap.set(color='black', linewidth=1)

        # Create custom legend patch
        box_patch = mpatches.Patch(facecolor='#C4C4C4', edgecolor='black', alpha=1, label=label)

        # Axis settings
        ax.set_xlim(0.5, num_bins + 0.5)
        ax.set_ylim(0, 1.1)
        ax.set_xlabel(xlabel, fontsize=22)
        ax.set_ylabel('Collision Rates', fontsize=22)
        ax.grid(True, alpha=0.5, axis='y')

        # Legend
        ax.legend(handles=[box_patch], loc='upper left', fontsize=18)

        # --- X-Ticks Logic ---
        # The original code set ticks at specific valid_positions if (pos-1) % 5 == 0.
        # This aligns ticks to the bin edges 0.0, 0.25, 0.5 etc depending on num_bins.
        tick_pos_to_show = [p for p in valid_positions if (p - 1) % 5 == 0]

        # We need to map these positions back to the bin_edges value
        # valid_positions contains 'k+1', so index is 'p-1'
        tick_labels_to_show = [f"{bin_edges[p-1]:.1f}" for p in tick_pos_to_show]

        ax.set_xticks(tick_pos_to_show)
        ax.set_xticklabels(tick_labels_to_show)
        ax.tick_params(axis='both', labelsize=20)

    # Hide unused subplots
    for j in range(n_rows, len(axes)):
        axes[j].set_visible(False)

    # plt.tight_layout()

    # Save logic
    if save_filename:
        if not save_filename.endswith('.png'):
            final_name = f"{save_filename}_boxplot_multiplot.png"
        else:
            final_name = save_filename

        plt.savefig(final_name, dpi=600, bbox_inches='tight')
        print(f"Plot saved to {final_name}")

    if show_plot:
        plt.show()

# def plot_overlapping_collision_boxplots(
#         df,
#         xlabel="Similarity",
#         save_filename=None,
#         num_bins=26,
#         show_plot=True
# ):
#     """
#     Generates a single plot with overlapping boxplots for all hash implementations.
#     Each row in df gets a different color, and boxes are offset within each bin.

#     Parameters:
#     - df: DataFrame containing 'similarity_values', 'collision_rates', and metadata.
#     - xlabel: Label for the X-axis.
#     - save_filename: If provided, saves the figure.
#     - num_bins: Number of bins to divide similarity (0 to 1).
#     - show_plot: Whether to display the plot.
#     """

#     fig, ax = plt.subplots(figsize=(16, 8))

#     # Color palette for different hash implementations
#     colors = plt.cm.Set3(np.linspace(0, 1, len(df)))
    
#     # Define bins
#     bin_edges = np.linspace(0, 1.05, num_bins + 1)
    
#     # Store legend handles
#     legend_handles = []

#     for row_idx, (idx, row) in enumerate(df.iterrows()):
#         # Extract data
#         x_values = np.array(row['similarity_values'])
#         y_values = np.array(row['collision_rates'])

#         # Build label
#         hash_name = row.get('hashname', 'Unknown')
#         seq_len = row.get('sequencelength', '?')
#         and_param = str(row.get('AND_param', '1'))
#         or_param = str(row.get('OR_param', '1'))

#         mutation_model_used = ""
#         if row['MutationModel'] == 0:
#             mutation_model_used = "SUB-ONLY"
#         elif row['MutationModel'] == 1:
#             mutation_model_used = "SubIndel-" + mutation_expressions[row["MutationExpression"]]

#         param_labels = ", ".join([f"{name}={row[name]}" for name in row['parameter_columns'] if row[name] != 0.0])
#         label = f"{hash_name} (L={seq_len}, {param_labels})[{mutation_model_used}]"
        
#         if not (and_param == '1' and or_param == '1'):
#             label += f" (b={and_param}, r={or_param})"

#         # Bin the data
#         binned_collision_rates = [[] for _ in range(num_bins)]

#         for j in range(len(x_values)):
#             sim_value = x_values[j]
#             collision_value = y_values[j]
#             bin_idx = int(sim_value * num_bins)
#             if bin_idx >= num_bins:
#                 bin_idx = num_bins - 1
#             binned_collision_rates[bin_idx].append(collision_value)

#         # Filter empty bins
#         plot_data = []
#         valid_positions = []

#         for k, bin_data in enumerate(binned_collision_rates):
#             if len(bin_data) > 0:
#                 plot_data.append(bin_data)
#                 valid_positions.append(k)

#         # Offset positions for overlapping boxplots
#         # Each row gets a different offset within the bin
#         offset = (row_idx - len(df) / 2) * 0.15  # Adjust width as needed
#         offset_positions = [pos + 1 + offset for pos in valid_positions]

#         # Plot boxplots with offset
#         bp = ax.boxplot(
#             plot_data,
#             positions=offset_positions,
#             widths=0.12,
#             patch_artist=True,
#             showfliers=False,
#             label=label
#         )

#         # Style the boxes with specific color
#         for box in bp['boxes']:
#             box.set(facecolor=colors[row_idx], alpha=0.7, edgecolor='black', linewidth=0.8)
#         for median in bp['medians']:
#             median.set(color='darkred', linewidth=1.5)
#         for whisker in bp['whiskers']:
#             whisker.set(color='black', linewidth=0.8)
#         for cap in bp['caps']:
#             cap.set(color='black', linewidth=0.8)

#         # Create legend patch
#         legend_handles.append(
#             mpatches.Patch(facecolor=colors[row_idx], edgecolor='black', alpha=0.7, label=label)
#         )

#     # Axis settings
#     ax.set_xlim(0.5, num_bins + 0.5)
#     ax.set_ylim(0, 1.1)
#     ax.set_xlabel(xlabel, fontsize=16)
#     ax.set_ylabel('Collision Rates', fontsize=16)
#     ax.grid(True, alpha=0.3, axis='y')

#     # X-axis ticks (every 10th bin)
#     tick_positions = [i + 1 for i in range(0, num_bins, 10)]
#     tick_labels = [f"{bin_edges[i]:.2f}" for i in range(0, num_bins, 10)]
#     ax.set_xticks(tick_positions)
#     ax.set_xticklabels(tick_labels)

#     # Legend
#     ax.legend(handles=legend_handles, loc='upper left', fontsize=11, ncol=1)
#     ax.tick_params(axis='both', labelsize=12)

#     fig.tight_layout()

#     # Save logic
#     if save_filename:
#         # Create directory if it doesn't exist
#         save_dir = os.path.dirname(save_filename)
#         if save_dir:
#             os.makedirs(save_dir, exist_ok=True)
        
#         if not save_filename.endswith('.png'):
#             final_name = f"{save_filename}_overlapping_boxplots.png"
#         else:
#             final_name = save_filename
        
#         plt.savefig(final_name, dpi=600, bbox_inches='tight')
#         print(f"Plot saved to {final_name}")

#     if show_plot:
#         plt.show()

def plot_verification_curves(
        df,
        n_cols=2,
        save_filename=None,
        show_plot=True
):
    """
    Generates a grid of verification distribution curves for all data rows.
    Each subplot shows KDE curves for rand_base_params, insmean, delRate, snp_rate, and insrate.

    Parameters:
    - df: DataFrame from read_collision_data_complete().
    - n_cols: Number of columns in the subplot grid.
    - save_filename: If provided, saves the figure.
    - show_plot: Whether to call plt.show().
    """

    n_rows = len(df)
    n_grid_rows = (n_rows + n_cols - 1) // n_cols

    # fig, axes = plt.subplots(nrows=n_grid_rows, ncols=n_cols, figsize=(14 * n_cols, 6 * n_grid_rows))
    fig, axes = plt.subplots(
        nrows=n_grid_rows, 
        ncols=n_cols, 
        figsize=(w_subplot * n_cols, h_subplot * n_grid_rows),
        constrained_layout=True 
    )

    # Handle edge cases
    if n_rows == 1:
        axes = np.array([axes])
    axes = axes.flatten()

    # Curve definitions: (column_name, label, color, linestyle, dashes, linewidth)
    curve_specs = [
        ('rand_base_params', 'rand_base_params', '#1f77b4', '-',    None,                    3),
        ('insmean',          'insmean',           '#2ca02c', '--',   None,                    2),
        ('delRate',          'delRate',            '#d62728', ':',    None,                    2),
        ('snp_rate',         'snp_rate',           '#ff7f0e', '-',    (5, 10),                 2),
        ('insrate',          'ins_rate',           '#1cd328', '-',    (3, 3, 1, 3, 1, 3),     2),
    ]

    for i, (idx, row) in enumerate(df.iterrows()):
        ax = axes[i]

        # Build subplot title/label
        hash_name = row.get('hashname', 'Unknown')
        seq_len = row.get('sequencelength', '?')

        mutation_model_used = ""
        if row['MutationModel'] == 0:
            mutation_model_used = "SUB-ONLY"
        elif row['MutationModel'] == 1:
            mutation_model_used = "SubIndel-" + mutation_expressions[row["MutationExpression"]]

        param_labels = ", ".join([f"{name}={row[name]}" for name in row['parameter_columns'] if row[name] != 0.0])
        title = f"{hash_name} (L={seq_len}, {param_labels})[{mutation_model_used}]"

        and_param = str(row.get('AND_param', '1'))
        or_param = str(row.get('OR_param', '1'))
        if not (and_param == '1' and or_param == '1'):
            title += f" AND={and_param}, OR={or_param}"

        # Plot each curve
        for col_name, label, color, linestyle, dashes, lw in curve_specs:
            values = row.get(col_name, [])
            if values is not None and len(values) > 1:
                plot_distribution_curve(
                    values,
                    ax=ax,
                    label=label,
                    color=color,
                    linestyle=linestyle,
                    dashes=dashes,
                    linewidth=lw,
                )

        ax.set_title(title, fontsize=24)
        ax.legend(fontsize=18, loc='upper right')
        ax.set_xlabel("Parameter Value", fontsize=22)
        ax.set_ylabel("Density", fontsize=22)
        ax.tick_params(axis='both', labelsize=20)

    # Hide unused subplots
    for j in range(n_rows, len(axes)):
        axes[j].set_visible(False)

    # plt.tight_layout()

    # Save logic
    if save_filename:
        if not save_filename.endswith('.png'):
            final_name = f"{save_filename}_verificationCurves.png"
        else:
            final_name = save_filename
        plt.savefig(final_name, dpi=600, bbox_inches='tight')
        print(f"Plot saved to {final_name}")

    if show_plot:
        plt.show()


def compute_rho(df, s1, s2, bin_width=0.05, eps=1e-6):
    """
    Compute the LSH rho value for hash function performance evaluation.
    """
    # Extract similarity values and collision rates from dataframe
    sim_values = np.asarray(df['similarity_values'], dtype=float)
    collision_probs = np.asarray(df['collision_rates'], dtype=float)
    
    # Find values in bin around r [r - bin_width, r + bin_width]
    mask_s1 = (sim_values >= s1 - bin_width) & (sim_values <= s1 + bin_width)
    collision_probs_s1 = collision_probs[mask_s1]
    
    # Find values in bin around cr [cr - bin_width, cr + bin_width]
    mask_s2 = (sim_values >= s2 - bin_width) & (sim_values <= s2 + bin_width)
    collision_probs_s2 = collision_probs[mask_s2]
    
    # Check that both bins have data
    if len(collision_probs_s1) == 0:
        raise ValueError(
            f"No similarity values found in bin [s1 - bin_width, s1 + bin_width] = "
            f"[{s1 - bin_width}, {s1 + bin_width}]. Consider increasing bin_width."
        )
    
    if len(collision_probs_s2) == 0:
        raise ValueError(
            f"No similarity values found in bin [s2 - bin_width, s2 + bin_width] = "
            f"[{s2 - bin_width}, {s2 + bin_width}]. Consider increasing bin_width."
        )
    
    # Compute average collision probabilities
    p1 = np.mean(collision_probs_s1)
    p2 = np.mean(collision_probs_s2)
    
    # Clamp for numerical stability
    p1 = min(p1, 1 - eps)
    p2 = max(p2, eps)
    
    # Compute rho
    rho = np.log(p1) / np.log(p2)
    
    # Build result dictionary
    result = {
        'p1': p1,
        'p2': p2,
        'rho': rho
    }
    
    # Add warning if rho >= 1
    if p1 <= p2:
        result['warning'] = "rho >= 1: hash family is not effective at these thresholds"
    
    return result


def main():
    parser = argparse.ArgumentParser(
        description="BioLSHasher collision data plotter - reads the BioLSHasher output data files and generates analysis plots."
    )
    parser.add_argument(
        "outdatafile",
        help="Path to the collision results output data file (e.g. results/collisionResults_SubseqHash-64.outdata)"
    )
    # parser.add_argument(
    #     "--similarity", "-s",
    #     default=None,
    #     help="Similarity label for axes (auto-detected from CSV if omitted)"
    # )
    parser.add_argument(
        "--s1", type=float, default=0.94,
        help="Upper similarity threshold for vertical line (default: 0.94)"
    )
    parser.add_argument(
        "--s2", type=float, default=0.62,
        help="Lower similarity threshold for vertical line (default: 0.62)"
    )
    parser.add_argument(
        "--no-verification", action="store_true",
        help="Skip the verification distribution curves"
    )
    parser.add_argument(
        "--no-scatter", action="store_true",
        help="Skip the color scatter plot"
    )
    parser.add_argument(
        "--no-boxplot", action="store_true",
        help="Skip the box plot"
    )
    args = parser.parse_args()

    # Validate file exists
    if not os.path.isfile(args.outdatafile):
        print(f"Error: file not found: {args.outdatafile}", file=sys.stderr)
        sys.exit(1)

    # Read data
    print(f"Reading: {args.outdatafile}")
    sections, df = read_collision_data_complete(args.outdatafile)

    if df.empty:
        print("Error: no data rows parsed from the CSV.", file=sys.stderr)
        sys.exit(1)

    print(f"Loaded {len(df)} data row(s).")

    # Auto-detect similarity label from the distance metric in the first row
    # if args.similarity is not None:
    #     similarity = args.similarity
    # else:
    metric_id = df.iloc[0].get('DistanceMetric', '')
    similarity = metric_id + " " + "Similarity";
    print(f"Auto-detected similarity label: {similarity}")

    similarity_name = df.iloc[0].get('DistanceMetric', 'results')

    # Create output directory for plots based on CSV path and similarity_name
    csv_dir = os.path.dirname(os.path.abspath(args.outdatafile))
    output_dir = os.path.join(csv_dir, similarity_name)
    os.makedirs(output_dir, exist_ok=True)
    print(f"Plots will be saved in: {output_dir}")

    # if not args.no_verification:
    #     print("Plotting verification distribution curves...")
    #     plot_verification_curves(df, save_filename=os.path.join(output_dir, similarity_name), show_plot=False)
    #     # print(f"Saved")

    # --- Average per bin ---
    print("Plotting binned average curve...")
    fig, ax = plt.subplots(figsize=(w_subplot, h_subplot))
    """Add ideal line, vertical markers, grid, legend, and optionally save."""
    # ax.plot([0, 1], [0, 1], color="black", linewidth=1, alpha=1, label="Theoretical similarity estimator")

    # if args.s1 is not None:
    #     ax.axvline(x=args.s1, color='black', linewidth=0.7, linestyle='dashed', label='_nolegend_')
    #     ax.text(args.s1 + 0.01, 0.95, f'$s_1 = {args.s1}$', color='black', fontsize=14, rotation=0)

    # if args.s2 is not None:
    #     ax.axvline(x=args.s2, color='black', linewidth=0.7, linestyle='dashed', label='_nolegend_')
    #     ax.text(args.s2 + 0.01, 0.95, f'$s_2 = {args.s2}$', color='black', fontsize=14, rotation=0)

    for idx, row in df.iterrows():
        plot_binned_average_curve(row, ax=ax)
    savename = os.path.join(output_dir, f"{similarity_name}_binaveraged.png")
    finalize_binned_average_plot(ax, similarity_label=similarity, s1=args.s1, s2=args.s2, savename=savename)
    print(f"Saved: {savename}")
    # plt.show()

    # # --- Color Scatter Plot ---
    # if not args.no_scatter:
    #     print("Plotting color scatter...")
    #     savename = os.path.join(output_dir, f"{similarity_name}_multicolor_scatter_plot.png")
    #     plot_color_scatter(df, xlabel=similarity, n_cols=2,  save_filename = savename, show_plot=False)
    #     # print(f"Saved: {savename}")

    # --- MonoColor Scatter Plot ---
    if not args.no_scatter:
        print("Plotting monocolor scatter...")
        # Dynamically set n_cols: use 1 column if only 1 row, else use 2
        n_cols = 1 if len(df) == 1 else 2
        savename = os.path.join(output_dir, f"{similarity_name}_monocolor_scatter_plot.png")
        plot_monocolor_scatter(df, xlabel=similarity, n_cols=n_cols, save_filename=savename, show_plot=False)

    # --- Box Plot ---
    if not args.no_boxplot:
        print("Plotting box plots...")
        n_cols = 1 if len(df) == 1 else 2
        savename = os.path.join(output_dir, f"{similarity_name}_box_plot.png")
        plot_hash_collision_boxplots(df, xlabel=similarity, n_cols=n_cols, save_filename=savename, show_plot=False)

    # Process all rows and compute rho for each
    rho_results = []
    s1 = 0.95
    s2 = 0.70
    bin_width = 0.05

    for idx, row in df.iterrows():
        try:
            result = compute_rho(row, s1, s2, bin_width=bin_width)
            
            # Extract hash parameters
            param_labels = ", ".join([f"{name}={row[name]}" for name in row['parameter_columns'] if row[name] != 0.0])
            
            rho_results.append({
                'hashname': row['hashname'],
                'sequencelength': row['sequencelength'],
                'hash_params': param_labels,
                'AND_param': row['AND_param'],
                'OR_param': row['OR_param'],
                's1': s1,
                's2': s2,
                'bin_width': bin_width,
                'p1': result['p1'],
                'p2': result['p2'],
                'rho': result['rho'],
                'warning': result.get('warning', None)
            })
            print(f"{row['hashname']} (L={row['sequencelength']}, {param_labels}): s1={s1}, s2={s2}, bin_width={bin_width} → p1={result['p1']:.6f}, p2={result['p2']:.6f}, rho={result['rho']:.4f}")
        except ValueError as e:
            print(f"Error for {row['hashname']}: {e}")

    # Convert to dataframe for analysis
    rho_df = pd.DataFrame(rho_results)
    print("\n" + "="*140)
    print(rho_df.to_string(index=False))

    # Save rho results to CSV with parameters
    rho_output_path = os.path.join(output_dir, f"{similarity_name}_rho_analysis.csv")
    rho_df.to_csv(rho_output_path, index=False)
    print(f"\nRho analysis saved to: {rho_output_path}")

    # plot_overlapping_collision_boxplots(
    #     df,
    #     xlabel=similarity,
    #     save_filename=os.path.join(output_dir, f"{similarity_name}_overlapping_boxplots"),
    #     num_bins=26,
    #     show_plot=False
    # )
    
    # Output processed dataframe
    processed_csv = os.path.join(output_dir, f"{similarity_name}_processed_coll.csv")
    df.to_csv(processed_csv, index=False)
    print(f"Saved processed dataframe to: {processed_csv}")

    print("Done.")


if __name__ == "__main__":
    main()
