import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import itertools
from adjustText import adjust_text
import pandas as pd

# Configurations
mutation_model = {0:"MUTATION_MODEL_SIMPLE_SNP_ONLY",1:"MUTATION_MODEL_GEOMETRIC_MUTATOR"}
mutation_expressions = {0:"BALANCED",1:"SUB_ONLY",2:"DEL_LITE",3:"INS_LITE",4:"SUB_LITE"}


def read_ann_data_complete(filepath):
    """
    Parse the custom ANN results CSV format with multiple experiment blocks.
 
    Format recap:
        :1:   experiment title
        :2:   metadata column headers (Hashname, SequenceLength, etc)
        :3:   metadata values
        :4.1: hash param names
        :4.2: hash param descriptions
        :4.3: hash param values
        :5:   data column headers  (b, r, c, Avg_Recall, etc)
        :6:   one data row per line
 
    Each :1: starts a new experiment block.  All metadata and hash params
    from that block are broadcast onto every data row from :6: lines, so
    the result is a flat DataFrame : one row per (b, r)
    configuration per experiment.
 
    Returns
    -------
    df : pd.DataFrame
        Columns:
          - experiment_id        (int, 1-based, unique per block)
          - test_name            (str, from :1:)
          - Hashname, SequenceLength, DistanceMetric,
            MutationModel, MutationExpression   (from :3:)
          - param_<name>         (float, one column per :4.x: parameter)
          - param_names          (list[str], parameter names : same for all
                                  rows in a block, useful for display)
          - b, r, c              (float to cast to int where appropriate)
          - Avg_Recall, Avg_Precision, Avg_FPR, Avg_Query_Time_ms, Avg_Index_Size_MB  (float, from :6:)
          - any other :5: columns present in the file
    """
    data_rows      = []
    current_exp_id = 0
 
    # Per-block state (reset on every :1:)
    test_name      = ''
    base_meta      = {}
    hash_parameters = {}   # {param_names, param_descs, param_values}
    col_headers    = []    # from :5:
 
    with open(filepath, 'r') as f:
        lines = f.readlines()
 
    for raw_line in lines:
        line = raw_line.strip()
        if not line:
            continue
 
        # ── Identify prefix ───────────────────────────────────────────────────
        if line.startswith(':1:'):
            current_exp_id += 1
            test_name       = line[3:].strip()
            base_meta       = {}
            hash_parameters = {}
            col_headers     = []
 
        elif line.startswith(':2:'):
            base_meta['_meta_keys'] = [x.strip() for x in line[3:].strip().split(',')]
 
        elif line.startswith(':3:'):
            values = [v.strip() for v in line[3:].strip().split(',')]
            # print(values)
            keys   = base_meta.pop('_meta_keys', [])
            # print(keys)
            # Parse known fields with proper types; keep unknowns as strings
            _type_map = {
                'SequenceLength':    int,
                'MutationModel':     int,
                'MutationExpression':int,
            }
            for k, v in zip(keys, values):
                try:
                    base_meta[k] = _type_map.get(k, str)(v)
                except (ValueError, TypeError):
                    base_meta[k] = v
 
        elif line.startswith(':4.1:'):
            hash_parameters['param_names'] = [x.strip() for x in line[5:].strip().split(',')]
 
        elif line.startswith(':4.2:'):
            hash_parameters['param_descs'] = [x.strip() for x in line[5:].strip().split(',')]
 
        elif line.startswith(':4.3:'):
            raw_vals = line[5:].strip().split(',')
            hash_parameters['param_values'] = [
                float(v.strip()) if v.strip() else 0.0 for v in raw_vals
            ]
 
        elif line.startswith(':5:'):
            col_headers   = [h.strip() for h in line[3:].strip().split(',')]
 
        elif line.startswith(':6:'):
 
            values = [v.strip() for v in line[3:].strip().split(',')]
            if len(values) != len(col_headers):
                continue   # malformed row so, skip silently
 
            # ── Build one flat row ────────────────────────────────────────────
            row = {
                'experiment_id':     current_exp_id,
                'test_name':         test_name,
                # Metadata fields (with safe .get fallbacks)
                'Hashname':          base_meta.get('Hashname', ''),
                'SequenceLength':    base_meta.get('SequenceLength', 0),
                'DistanceMetric':    base_meta.get('Distance Metric', ''),
                'MutationModel':     base_meta.get('Mutation Model', 0),
                'MutationExpression':base_meta.get('Mutation Expression', 0),
                # Hash param names list (same for every row in this block)
                'param_names':       hash_parameters.get('param_names', []),
            }
 
            # Hash parameters as individual columns: param_k, param_d, ...
            p_names  = hash_parameters.get('param_names',  [])
            p_values = hash_parameters.get('param_values', [])
            for pname, pval in zip(p_names, p_values):
                row[f'param_{pname}'] = pval
 
            # Data columns from :6: line
            for col, val in zip(col_headers, values):
                try:
                    row[col] = float(val)
                except ValueError:
                    row[col] = val
 
            data_rows.append(row)
 
    # ── Report dropped old-format blocks ──────────────────────────────────────
    n_new   = len({r['experiment_id'] for r in data_rows})
 
    if not data_rows:
        print("[warn] No data rows found. Returning empty DataFrame.")
        return pd.DataFrame()
 
    df = pd.DataFrame(data_rows)
 
    # Tidy up integer columns that arrived as floats
    for col in ('b', 'r', 'c', 'SequenceLength', 'MutationModel', 'MutationExpression'):
        if col in df.columns:
            df[col] = df[col].astype(int)
 
    print(f"[info] Parsed {len(df)} rows across {n_new} experiment(s).")
    return df
 


# ──────────────────────────────────────────────────────────────────────────────
# Label helpers  (operate on a single-experiment sub-DataFrame)
# ──────────────────────────────────────────────────────────────────────────────
 
def experiment_label(edf, exp_id):
    """
    Build a human-readable title for one experiment sub-DataFrame.
 
    Format:  Exp<id>: <Hashname>  <param_k>=<v>, ...
    """
    row = edf.iloc[0]
    parts = []
 
    hash_name = row.get('Hashname', '')
    if hash_name:
        parts.append(str(hash_name))
 
    # Dynamic hash params: columns named param_*
    param_cols = [c for c in edf.columns if c.startswith('param_') and c != 'param_names']
    param_str  = ', '.join(f"{c[len('param_'):]!s}={row[c]}" for c in param_cols)
    if param_str:
        parts.append(param_str)
 
    label = '  '.join(parts) if parts else f'experiment {exp_id}'
    return f"Exp{exp_id}: {label}"
 
 
def short_label(edf, exp_id):
    """Compact single-line label for legends."""
    row       = edf.iloc[0]
    hash_name = row.get('Hashname', f'exp{exp_id}')
    param_cols = [c for c in edf.columns if c.startswith('param_') and c != 'param_names']
    param_str  = ','.join(f"{c[len('param_'):]!s}={row[c]}" for c in param_cols)
    base = f"{hash_name} ({param_str})" if param_str else str(hash_name)
    return f"Exp{exp_id}: {base}"

# ──────────────────────────────────────────────────────────────────────────────
# Shared style helpers
# ──────────────────────────────────────────────────────────────────────────────
 
MARKERS = ['o', 's', '^', 'D', 'v', 'P', '*', 'X', 'h', '<']
 
# def exp_style(exp_idx, n_exps):
#     cmap   = plt.cm.tab10
#     color  = cmap(exp_idx / max(n_exps - 1, 1))
#     marker = MARKERS[exp_idx % len(MARKERS)]
#     return color, marker
 
COLORS = [
'#e41a1c', '#377eb8', '#4daf4a', '#ff7f00', '#984ea3',
'#a65628', '#f781bf', '#999999', '#66c2a5', '#fc8d62'
]  # ColorBrewer Set1 + Set2 mix — perceptually distinct

def exp_style(exp_idx, n_exps):
    color  = COLORS[exp_idx % len(COLORS)]
    marker = MARKERS[exp_idx % len(MARKERS)]
    return color, marker


 
# def b_color_map(all_b_vals):
#     cmap = plt.cm.viridis
#     return {b: cmap(i / max(len(all_b_vals) - 1, 1)) for i, b in enumerate(sorted(all_b_vals))}
 

def b_color_map(all_b_vals):
    palette = ['#1b9e77', '#d95f02', '#7570b3', '#e7298a',
               '#66a61e', '#e6ab02', '#a6761d', '#666666']  # ColorBrewer Dark2
    sorted_b = sorted(all_b_vals)
    return {b: palette[i % len(palette)] for i, b in enumerate(sorted_b)}


def finalize(ax, fig, output_path, texts=None):
    if texts:
        adjust_text(
            texts,
            ax=ax,
            expand=(1.5, 1.5),
            arrowprops=dict(arrowstyle='-', color='gray', lw=0.9, alpha=0.7, shrinkA=12, shrinkB=5),
            # verbose=False,
        )
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved : {output_path}")
    # plt.show()
    plt.close(fig)
    
def iter_experiments(df):
    """
    Yield (exp_id, sub_df) for each unique experiment_id in df,
    in ascending order.  exp_id is the 1-based integer from the file.
    """
    for exp_id in sorted(df['experiment_id'].unique()):
        yield exp_id, df[df['experiment_id'] == exp_id].copy()





## ──────────────────────────────────────────────────────────────────────────────
## ──────────────────────────────────────────────────────────────────────────────

# Plot 1: Query Time vs Recall
def plot_query_time_vs_recall(df, output_path):
    """
    Avg_Query_Time_ms vs Avg_Recall.  Exposes the recall–speed tradeoff.
    Color = experiment, points annotated with (b,r).
    """
    if 'Avg_Query_Time_ms' not in df.columns:
        print("[skip] plot_query_time_vs_recall: column 'Avg_Query_Time_ms' not present.")
        return
 
    fig, ax = plt.subplots(figsize=(11, 7))
    texts   = []
    n       = df['experiment_id'].nunique()
 
    for plot_idx, (exp_id, edf) in enumerate(iter_experiments(df)):
        color, marker = exp_style(plot_idx, n)
        label         = short_label(edf, exp_id)
 
        ax.scatter(edf['Avg_Recall'], edf['Avg_Query_Time_ms'],
                   color=color, marker=marker, s=70, alpha=0.75,
                   edgecolors='black', linewidth=0.4,
                   label=label, zorder=5)
 
        for _, row in edf.iterrows():
            texts.append(ax.text(row['Avg_Recall'], row['Avg_Query_Time_ms'],
                                 f"({row['b']},{row['r']})",
                                 fontsize=10, fontstyle='italic',
                                 alpha=0.9,
                                 bbox=dict(boxstyle='round,pad=0.15',
                                           facecolor='white',
                                           edgecolor='none',
                                           alpha=0.6)))
 
    ax.set_xlabel('Average Recall', fontsize=13)
    ax.set_ylabel('Average Query Time (ms)', fontsize=13)
    ax.legend(loc='best', fontsize=10, title='Experiment')
    finalize(ax, fig, output_path, texts)

# Plot 2: Index Size vs Recall
def plot_index_size_vs_recall(df, output_path):
    """
    Avg_Index_Size_MB vs Avg_Recall.  Exposes the recall vs memory tradeoff.
    """
    if 'Avg_Index_Size_MB' not in df.columns:
        print("[skip] plot_index_size_vs_recall: column 'Avg_Index_Size_MB' not present.")
        return
 
    fig, ax = plt.subplots(figsize=(11, 7))
    texts   = []
    n       = df['experiment_id'].nunique()
 
    for plot_idx, (exp_id, edf) in enumerate(iter_experiments(df)):
        color, marker = exp_style(plot_idx, n)
        label         = short_label(edf, exp_id)
 
        ax.scatter(edf['Avg_Recall'], edf['Avg_Index_Size_MB'],
                   color=color, marker=marker, s=70, alpha=0.75,
                   edgecolors='black', linewidth=0.4,
                   label=label, zorder=5)
 
        for _, row in edf.iterrows():
            texts.append(ax.text(row['Avg_Recall'], row['Avg_Index_Size_MB'],
                                 f"({row['b']},{row['r']})",
                                 fontsize=10,fontstyle='italic',
                                 alpha=0.9,
                                 bbox=dict(boxstyle='round,pad=0.15',
                                           facecolor='white',
                                           edgecolor='none',
                                           alpha=0.6)))
            
 
    ax.set_xlabel('Average Recall', fontsize=13)
    ax.set_ylabel('Average Index Size (MB)', fontsize=13)
    ax.legend(loc='best', fontsize=10, title='Experiment')
    finalize(ax, fig, output_path, texts)
 
# Plot 3: All (b,r) scatter  : experiments separated by color
def plot_recall_vs_fpr_scatter(df, output_path, use_log=False):
    """
    Scatter of every (b,r) point.  Color = experiment, label = (b,r).
    Gives a quick global overview of the full parameter space.
    """
    fig, ax = plt.subplots(figsize=(11, 7))
    texts   = []
    n       = df['experiment_id'].nunique()
 
    for plot_idx, (exp_id, edf) in enumerate(iter_experiments(df)):
        color, marker = exp_style(plot_idx, n)
        label         = short_label(edf, exp_id)
        
        # ax.scatter(edf['Avg_FPR'].clip(lower=1e-8),edf['Avg_Recall'],
        #            color=color, marker=marker, s=70, alpha=0.75,
        #            edgecolors='black', linewidth=0.4,
        #            label=label, zorder=5)

        ax.scatter(edf['Avg_FPR'],edf['Avg_Recall'],
                   color=color, marker=marker, s=70, alpha=0.75,
                   edgecolors='black', linewidth=0.4,
                   label=label, zorder=5)
 
        for _, row in edf.iterrows():
            texts.append(ax.text(row['Avg_FPR'], row['Avg_Recall'],
                                 f"({row['b']},{row['r']})",
                                 fontsize=10, fontstyle='italic',
                                 alpha=0.9,
                                 bbox=dict(boxstyle='round,pad=0.15',
                                           facecolor='white',
                                           edgecolor='none',
                                           alpha=0.6)))
                                 
    scale_lbl = ' (log)' if use_log else ''
    ax.set_xlabel(f'Average FPR{scale_lbl}', fontsize=13)
    ax.set_ylabel('Average Recall', fontsize=13)
    if use_log:
        ax.set_xscale('symlog', linthresh=1e-7)
    
    ax.legend(loc='best', fontsize=10, title='Experiment')
    finalize(ax, fig, output_path, texts)
 
# Plot 4: Recall vs FPR  (all (b, r) points, lines grouped by b) Lined 
def plot_recall_vs_fpr(df, output_path, use_log=False):
    """
        ROC-style: FPR on x, Recall on y.  All experiments on one axes.
        Color = b value (shared colormap across experiments).
        Marker = experiment (so overlapping b values stay distinguishable).
        Each point labelled r=<r>.
    """
    fig, ax = plt.subplots(figsize=(11, 7))
    texts   = []

    # if not use_log:
    #     ax.plot([0, 1], [0, 1], color='gray', linestyle=':', linewidth=1.0,
    #             alpha=0.6, label='random baseline', zorder=1)

    all_b    = sorted(df['b'].unique())
    b_colors = b_color_map(all_b)
    n        = df['experiment_id'].nunique()

    for plot_idx, (exp_id, edf) in enumerate(iter_experiments(df)):
        _, marker = exp_style(plot_idx, n)
        exp_lbl   = short_label(edf, exp_id)

        for b_val in all_b:
            bdf = edf[edf['b'] == b_val].sort_values('Avg_FPR')
            if bdf.empty:
                continue

            ax.plot(bdf['Avg_FPR'], bdf['Avg_Recall'],
                    color=b_colors[b_val], marker=marker, markersize=10,
                    markeredgecolor='black', markeredgewidth=0.5,
                    linestyle='--', linewidth=1.2, alpha=0.8,
                    label=f'b={b_val} ({exp_lbl})', zorder=5)

            for _, row in bdf.iterrows():
                texts.append(ax.text(row['Avg_FPR'], row['Avg_Recall'],
                                    f"r={row['r']}", fontsize=10,
                                    fontstyle='italic',alpha=0.9,
                                    bbox=dict(boxstyle='round,pad=0.15',
                                           facecolor='white',
                                           edgecolor='none',
                                           alpha=0.6)))

    scale_lbl = ' (log)' if use_log else ''
    ax.set_xlabel(f'Average FPR{scale_lbl}', fontsize=13)
    ax.set_ylabel('Average Recall', fontsize=13)
    if use_log:
        ax.set_xscale('symlog', linthresh=1e-7)
    ax.legend(loc='best', fontsize=9, title='b  (experiment)')
    finalize(ax, fig, output_path, texts)
 
# Plot 5: Precision Vs Recall (All (b,r) scatter)
def plot_precision_vs_recall(df, output_path):
    """
    Precision-Recall scatter. All experiments on one axes.
    Color = experiment. Each point labelled (b,r).
    """
    fig, ax = plt.subplots(figsize=(11, 7))   # wider canvas gives adjust_text more room
    texts   = []
    n       = df['experiment_id'].nunique()

    for plot_idx, (exp_id, edf) in enumerate(iter_experiments(df)):
        color, marker = exp_style(plot_idx, n)
        label         = short_label(edf, exp_id)

        ax.scatter(edf['Avg_Recall'], edf['Avg_Precision'],
                   color=color, marker=marker, s=70, alpha=0.75,
                   edgecolors='black', linewidth=0.4,
                   label=label, zorder=5)

        for _, row in edf.iterrows():
            texts.append(ax.text(row['Avg_Recall'], row['Avg_Precision'],
                                 f"({int(row['b'])},{int(row['r'])})",
                                 fontsize=10,
                                 fontstyle='italic',
                                 alpha=0.9,
                                 bbox=dict(boxstyle='round,pad=0.15',
                                           facecolor='white',
                                           edgecolor='none',
                                           alpha=0.6)))

    ax.set_xlabel('Average Recall', fontsize=13)
    ax.set_ylabel('Average Precision', fontsize=13)
    ax.legend(loc='best', fontsize=10, title='Experiment')
    finalize(ax, fig, output_path, texts)

# Plot 6: Best FPR per recall bin
def plot_best_fpr_per_recall_bin(df, output_path, num_bins=15, use_log=False):
    """
        For each experiment, bin all (b,r) points by Avg_Recall.
        Pick the minimum FPR in each bin
        Each point is annotated with the (b,r) pair that achieved the minimum.
    """
    fig, ax = plt.subplots(figsize=(11, 7))
    texts   = []
 
    bin_edges = np.linspace(df['Avg_Recall'].min() - 1e-9,
                            df['Avg_Recall'].max() + 1e-9, num_bins + 1)
 
    exp_ids = sorted(df['experiment_id'].unique())
    n       = len(exp_ids)
 
    for plot_idx, (exp_id, edf) in enumerate(iter_experiments(df)):
        color, marker = exp_style(plot_idx, n)
        label         = short_label(edf, exp_id)
 
        bin_centers, min_fprs, annotations = [], [], []
 
        for i in range(num_bins):
            lo, hi  = bin_edges[i], bin_edges[i + 1]
            mask    = (edf['Avg_Recall'] >= lo) & (edf['Avg_Recall'] < hi)
            sub     = edf[mask]
            if sub.empty:
                continue
            best = sub.loc[sub['Avg_FPR'].idxmin()]
            bin_centers.append((lo + hi) / 2)
            min_fprs.append(best['Avg_FPR'])
            annotations.append(f"({best['b']},{best['r']})")
 
        order        = np.argsort(bin_centers)
        bin_centers  = [bin_centers[i] for i in order]
        min_fprs     = [min_fprs[i]    for i in order]
        annotations  = [annotations[i] for i in order]
 
        ax.plot(bin_centers, min_fprs,
                color=color, marker=marker, markersize=10,
                markeredgecolor='black', markeredgewidth=0.5,
                linestyle='--', linewidth=1.1, alpha=0.5,
                label=label, zorder=5)
 
        for x, y, ann in zip(bin_centers, min_fprs, annotations):
            texts.append(ax.text(x, y, ann, fontsize=10, fontstyle='italic',
                                 alpha=0.9,
                                 bbox=dict(boxstyle='round,pad=0.15',
                                           facecolor='white',
                                           edgecolor='none',
                                           alpha=0.6)))
            
    scale_lbl = ' (log)' if use_log else ''
    ax.set_xlabel('Recall Bin Center', fontsize=13)
    ax.set_ylabel(f'Best (Min) FPR in Bin{scale_lbl}', fontsize=13)
    if use_log:
        ax.set_yscale('symlog', linthresh=1e-7)
    ax.legend(loc='best', fontsize=10, title='Experiment')
    finalize(ax, fig, output_path, texts)


# Entry point
def run_all_plots(filepath, out_dir=None):
    df = read_ann_data_complete(filepath)
 
    if df.empty:
        print("No new-format experiments found. Nothing to plot.")
        return
    
    # --- Directory creation ---
    # Use similarity name for subdirectory
    similarity_name = df.iloc[0].get('DistanceMetric', 'results')
    csv_dir = os.path.dirname(os.path.abspath(filepath))
    output_dir = os.path.join(csv_dir, similarity_name)
    os.makedirs(output_dir, exist_ok=True)
    print(f"Plots and CSV will be saved in: {output_dir}")

    # # Summary
    # for exp_id, edf in iter_experiments(df):
    #     print(f"  [{exp_id}] {experiment_label(edf, exp_id)}  -->  {len(edf)} rows present")
    
    
    # stem = Path(filepath).stem
    # out_dir = Path(out_dir) if out_dir else Path(filepath).parent
    # out_dir.mkdir(parents=True, exist_ok=True)
 
    # --- Plot saving ---
    def out(suffix):
        return os.path.join(output_dir, f"{similarity_name}_{suffix}.png")
    
    plot_query_time_vs_recall(  df, out('query_time_vs_recall'))
    plot_index_size_vs_recall(  df, out('index_size_vs_recall'))

    plot_recall_vs_fpr_scatter( df, out('recall_vs_fpr_scatter'))
    plot_recall_vs_fpr_scatter( df, out('recall_vs_fpr_scatter_log'),    use_log=True)
    
    plot_recall_vs_fpr(         df, out('recall_vs_fpr'))
    plot_recall_vs_fpr(         df, out('recall_vs_fpr_log'),    use_log=True)

    plot_precision_vs_recall(   df, out('precision_vs_recall'))
    
    plot_best_fpr_per_recall_bin(df, out('best_fpr_per_bin'))
    plot_best_fpr_per_recall_bin(df, out('best_fpr_per_bin_log'), use_log=True)
    
    # --- Save processed DataFrame ---
    processed_csv = os.path.join(output_dir, f"{similarity_name}_processed_ann.csv")
    df.to_csv(processed_csv, index=False)
    print(f"Saved processed dataframe to: {processed_csv}")
    
    return df

 
if __name__ == "__main__":
    import sys
    import os
    from pathlib import Path
 
    if len(sys.argv) < 2:
        raise ValueError("Missing required argument: 'filepath'.\nUsage: python script.py <filepath> [out_dir]")
        
    filepath = sys.argv[1]
    out_dir  = sys.argv[2] if len(sys.argv) > 2 else None
 
    df = run_all_plots(filepath, out_dir)
    # if df is not None:
    #     print("\nDataFrame shape:", df.shape)
    #     print(df.to_string())







#Extra codes
# def plot_recall_vs_fpr(df, output_path, use_log=False):
#     """
#     All (b, r) points on a single plot.
#     Lines grouped by b within each experiment.
#     Color indicates b-value, line style/marker indicates experiment.
#     """
#     fig, ax = plt.subplots(figsize=(12, 8))
#     texts = []
 
#     all_b    = sorted(df['b'].unique())
#     b_colors = b_color_map(all_b)
    
#     # Define styles to distinguish different experiments
#     markers = ['o', 's', '^', 'D', 'v', 'p', '*']
#     linestyles = ['-', '--', '-.', ':']
 
#     for plot_idx, (exp_id, edf) in enumerate(iter_experiments(df)):
#         marker    = markers[plot_idx % len(markers)]
#         linestyle = linestyles[plot_idx % len(linestyles)]
#         exp_lbl   = experiment_label(edf, exp_id)
 
#         for b_val in all_b:
#             bdf = edf[edf['b'] == b_val].sort_values('Avg_Recall')
#             if bdf.empty:
#                 continue
 
#             ax.plot(bdf['Avg_Recall'], bdf['Avg_FPR'],
#                     color=b_colors[b_val], marker=marker, markersize=6,
#                     markeredgecolor='black', markeredgewidth=0.5,
#                     linestyle=linestyle, linewidth=1.2, alpha=0.8,
#                     label=f'{exp_lbl} (b={b_val})', zorder=5)
 
#             for _, row in bdf.iterrows():
#                 texts.append(ax.text(row['Avg_Recall'], row['Avg_FPR'],
#                                      f"r={row['r']}", fontsize=8,
#                                      fontweight='bold', alpha=0.85))
 
#     scale_lbl = ' (log)' if use_log else ''
#     ax.set_xlabel('Average Recall', fontsize=12)
#     ax.set_ylabel(f'Average FPR{scale_lbl}', fontsize=12)
#     ax.set_title('Recall vs FPR (All Experiments)', fontsize=13, fontweight='bold')
    
#     if use_log:
#         # Applying the symlog fix from the scatter plot to the Y-axis here
#         ax.set_yscale('symlog', linthresh=1e-5)
        
#     # Move legend outside the plot to avoid covering data points
#     ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.5), 
#               fontsize=9, title='Experiment & b-value')
#     ax.grid(True, alpha=0.3)
    
#     adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', color='gray', lw=0.5))
 
#     plt.tight_layout()
#     plt.savefig(output_path, dpi=300, bbox_inches='tight')
#     print(f"Saved : {output_path}")
#     plt.show()
#     plt.close(fig)



# def plot_recall_vs_fpr(df, output_path, use_log=False):
#     """
#     One subplot per experiment.
#     Within each subplot: color = b value, line connects points with same b
#     sorted by recall, each point labelled r=<r>.
#     """
#     exp_ids = sorted(df['experiment_id'].unique())
#     n       = len(exp_ids)
#     ncols   = min(n, 2)
#     nrows   = (n + ncols - 1) // ncols
#     fig, axes = plt.subplots(nrows, ncols, figsize=(7 * ncols, 5.5 * nrows), squeeze=False)
 
#     all_b    = sorted(df['b'].unique())
#     b_colors = b_color_map(all_b)
 
#     for plot_idx, (exp_id, edf) in enumerate(iter_experiments(df)):
#         ax    = axes[plot_idx // ncols][plot_idx % ncols]
#         texts = []
 
#         for b_val in all_b:
#             bdf = edf[edf['b'] == b_val].sort_values('Avg_Recall')
#             if bdf.empty:
#                 continue
 
#             ax.plot(bdf['Avg_Recall'], bdf['Avg_FPR'],
#                     color=b_colors[b_val], marker='o', markersize=6,
#                     markeredgecolor='black', markeredgewidth=0.5,
#                     linestyle='--', linewidth=1.2, alpha=0.8,
#                     label=f'b={b_val}', zorder=5)
 
#             for _, row in bdf.iterrows():
#                 texts.append(ax.text(row['Avg_Recall'], row['Avg_FPR'],
#                                      f"r={row['r']}", fontsize=8,
#                                      fontweight='bold', alpha=0.85))
 
#         scale_lbl = ' (log)' if use_log else ''
#         ax.set_xlabel('Average Recall', fontsize=12)
#         ax.set_ylabel(f'Average FPR{scale_lbl}', fontsize=12)
#         ax.set_title(experiment_label(edf, exp_id), fontsize=11, fontweight='bold')
#         if use_log:
#             ax.set_yscale('log')
#         ax.legend(loc='best', fontsize=9, title='b value')
#         ax.grid(True, alpha=0.3)
#         adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', color='gray', lw=0.5))
 
#     for i in range(n, nrows * ncols):
#         axes[i // ncols][i % ncols].set_visible(False)
 
#     plt.tight_layout()
#     plt.savefig(output_path, dpi=300, bbox_inches='tight')
#     print(f"Saved : {output_path}")
#     plt.show()
#     plt.close(fig)
