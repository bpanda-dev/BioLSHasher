import argparse
import os
import sys
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
import ast

mutation_expressions = {0: "BALANCED", 1: "SUB_ONLY", 2: "DEL_LITE", 3: "INS_LITE", 4: "SUB_LITE"}


def safe_literal_eval(val):
    if isinstance(val, str):
        try:
            return ast.literal_eval(val)
        except:
            return val
    return val


def read_processed_dataframe(csv_path):
    df = pd.read_csv(csv_path)
    # Reconstruct array columns from strings
    list_cols = [
        'similarity_values', 'snp_rate', 'delRate', 'insmean', 'insrate',
        'stayRate', 'collision_rates', 'rand_base_params', 'parameter_columns'
    ]
    for col in list_cols:
        if col in df.columns:
            df[col] = df[col].apply(safe_literal_eval)
    return df



def _param_str(row, key):
    v = row.get(key, None)
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return '1'
    return str(int(v)) if isinstance(v, (float, np.floating)) and v == int(v) else str(v)


def generate_row_label(row):
    mutation_model_used = ""
    if row['MutationModel'] == 0:
        mutation_model_used = "SUB-ONLY"
    elif row['MutationModel'] == 1:
        mutation_model_used = "SubIndel-" + mutation_expressions.get(row["MutationExpression"], "")

    param_labels = ""
    if 'parameter_columns' in row and isinstance(row['parameter_columns'], list):
        param_labels = ", ".join([f"{name}={row[name]}" for name in row['parameter_columns'] if name in row])

    param_part = f", {param_labels}" if param_labels else ""

    and_p = _param_str(row, 'AND_param')
    or_p  = _param_str(row, 'OR_param')

    and_or_part = f", AND={and_p}, OR={or_p}" if (and_p != '1' or or_p != '1') else ""

    return f"{row['hashname']} (L={row['sequencelength']}{param_part}{and_or_part})[{mutation_model_used}]"


def plot_scatter(df, sim_metric=""):
    fig = go.Figure()

    for idx, row in df.iterrows():
        label = generate_row_label(row)
        x_vals = row['similarity_values']
        y_vals = row['collision_rates']

        fig.add_trace(go.Scatter(
            x=x_vals, y=y_vals,
            mode='markers',
            name=label,
            marker=dict(size=6, opacity=0.6),
            hovertemplate='Similarity: %{x:.4f}<br>Collision Rate: %{y:.4f}'
        ))

    fig.update_layout(
        title="Collision Rates vs Similarity (Scatter)",
        xaxis_title=f"{sim_metric} Similarity" if sim_metric else "Similarity",
        yaxis_title="Collision Rates",
        xaxis=dict(range=[0, 1.01]),
        yaxis=dict(range=[0, 1.01]),
        template="ggplot2",
        colorway=px.colors.qualitative.Safe,
        hoverlabel=dict(namelength=-1),
        legend=dict(title="Hash Configurations", bgcolor='rgba(255,255,255,0.9)')
    )

    return fig


def plot_binned_average(df, sim_metric=""):
    fig = go.Figure()

    bin_edges = np.arange(0, 1.04, 0.04)
    num_bins = len(bin_edges) - 1
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    for idx, row in df.iterrows():
        label = generate_row_label(row)
        print(">>>>", label)
        x_vals = np.array(row['similarity_values'])
        y_vals = np.array(row['collision_rates'])

        bin_means = []
        for i in range(len(bin_edges) - 1):
            mask = (x_vals >= bin_edges[i]) & (x_vals < bin_edges[i + 1])
            if mask.sum() > 0:
                bin_means.append(np.mean(y_vals[mask]))
            else:
                bin_means.append(np.nan)

        bin_means = np.array(bin_means)
        valid = ~np.isnan(bin_means)

        fig.add_trace(go.Scatter(
            x=bin_centers[valid],
            y=bin_means[valid],
            mode='lines+markers',
            name=label,
            marker=dict(size=8),
            line=dict(width=2)
        ))

    fig.update_layout(
        title="Average Collision Rate per Bin",
        xaxis_title=f"{sim_metric} Similarity" if sim_metric else "Similarity",
        yaxis_title="Average Collision Rate",
        xaxis=dict(range=[0, 1.01]),
        yaxis=dict(range=[0, 1.01]),
        template="ggplot2",
        colorway=px.colors.qualitative.Safe,
        hoverlabel=dict(namelength=-1),
        legend=dict(title="Hash Configurations", bgcolor='rgba(255,255,255,0.9)')
    )

    return fig

def plot_box_plot(df, num_bins=51, sim_metric=""):
    fig = go.Figure()

    bin_edges = np.linspace(0, 1.0, num_bins + 1)

    for loop_idx, (idx, row) in enumerate(df.iterrows()):
        label = generate_row_label(row)
        x_vals = np.array(row['similarity_values'])
        y_vals = np.array(row['collision_rates'])

        binned_y = []
        binned_x = []
        for j in range(len(x_vals)):
            sim_value = x_vals[j]
            bin_idx = int(sim_value * num_bins)
            if bin_idx >= num_bins:
                bin_idx = num_bins - 1

            binned_x.append(f"{bin_edges[bin_idx]:.2f}")
            binned_y.append(y_vals[j])

        fig.add_trace(go.Box(
            x=binned_x,
            y=binned_y,
            name=label,
            boxpoints=False,
            marker_color=px.colors.qualitative.Safe[loop_idx % len(px.colors.qualitative.Safe)]
        ))

    fig.update_layout(
        title="Collision Rates Boxplot per Bin",
        xaxis_title=f"{sim_metric} Similarity Bin Start" if sim_metric else "Similarity Bin Start",
        yaxis_title="Collision Rates",
        boxmode='group',
        template="ggplot2",
        colorway=px.colors.qualitative.Safe,
        yaxis=dict(range=[0, 1.01]),
        hoverlabel=dict(namelength=-1),
        legend=dict(title="Hash Configurations", bgcolor='rgba(255,255,255,0.9)')
    )

    fig.update_xaxes(categoryorder='array', categoryarray=[f"{v:.2f}" for v in bin_edges])

    return fig


def get_ann_groups(df):
    if 'experiment_id' in df.columns:
        group_cols = ['experiment_id']
    else:
        group_cols = ['Hashname', 'SequenceLength', 'DistanceMetric', 'MutationModel', 'MutationExpression', 'param_k', 'param_d', 'c']
        group_cols = [col for col in group_cols if col in df.columns]

    groups = []
    for keys, group in df.groupby(group_cols, dropna=False):
        row = group.iloc[0]
        hash_name = row.get('Hashname', f'Exp{keys}')

        param_cols = [
            c for c in group.columns
            if c.startswith('param_') and c != 'param_names'
            and pd.notna(row[c]) and row[c] != 0
        ]
        param_str = ','.join(
            f"{c[len('param_'):]}="
            f"{int(row[c]) if isinstance(row[c], (float, np.floating)) and row[c] == int(row[c]) else row[c]}"
            for c in param_cols
        )

        if 'experiment_id' in df.columns:
            exp_id = keys if not isinstance(keys, tuple) else keys[0]
            label = f"Exp{int(exp_id)}: {hash_name}" + (f" ({param_str})" if param_str else "")
        else:
            label = str(hash_name) + (f" ({param_str})" if param_str else "")

        group = group.sort_values('Avg_Recall')
        groups.append((label, group))
    return groups


def plot_ann_fpr_vs_recall(df):
    fig = go.Figure()
    for label, group in get_ann_groups(df):
        text = [f"({b},{r})" for b, r in zip(group['b'], group['r'])] if 'b' in group.columns and 'r' in group.columns else []

        clipped_fpr = [max(v, 1e-6) for v in group['Avg_FPR']]

        fig.add_trace(go.Scatter(
            x=clipped_fpr, y=group['Avg_Recall'],
            mode='markers', name=label, text=text,
            marker=dict(size=10, opacity=0.75, line=dict(width=1, color='DarkSlateGrey')),
            hovertemplate='FPR: %{x:.4g}<br>Recall: %{y:.4f}<br>%{text}' if text else 'FPR: %{x:.4g}<br>Recall: %{y:.4f}'
        ))

    fig.update_layout(
        title="Avg Recall vs Avg FPR", xaxis_title="Average FPR (log)", yaxis_title="Average Recall",
        template="ggplot2", colorway=px.colors.qualitative.Safe, hoverlabel=dict(namelength=-1),
        legend=dict(title="Experiment", bgcolor='rgba(255,255,255,0.9)')
    )
    fig.update_xaxes(type="log")
    return fig


def plot_ann_precision_vs_recall(df):
    fig = go.Figure()
    for label, group in get_ann_groups(df):
        text = [f"(b,r):({b},{r})" for b, r in zip(group['b'], group['r'])] if 'b' in group.columns and 'r' in group.columns else None
        fig.add_trace(go.Scatter(
            x=group['Avg_Recall'], y=group['Avg_Precision'],
            mode='markers', name=label, text=text,
            marker=dict(size=10, opacity=0.75, line=dict(width=1, color='DarkSlateGrey')),
            hovertemplate='Recall: %{x:.4f}<br>Precision: %{y:.4f}<br>%{text}' if text else 'Recall: %{x:.4f}<br>Precision: %{y:.4f}'
        ))
    fig.update_layout(
        title="Avg Precision vs Avg Recall", xaxis_title="Average Recall", yaxis_title="Average Precision",
        template="ggplot2", colorway=px.colors.qualitative.Safe, hoverlabel=dict(namelength=-1),
        legend=dict(title="Experiment", bgcolor='rgba(255,255,255,0.9)')
    )
    return fig


def plot_ann_min_fpr_per_bin(df, num_bins=50):
    fig = go.Figure()
    bin_edges = np.linspace(df['Avg_Recall'].min() - 1e-9, df['Avg_Recall'].max() + 1e-9, num_bins + 1)

    for label, group in get_ann_groups(df):
        bin_centers, min_fprs, annotations = [], [], []
        for i in range(num_bins):
            lo, hi = bin_edges[i], bin_edges[i + 1]
            mask = (group['Avg_Recall'] >= lo) & (group['Avg_Recall'] < hi)
            sub = group[mask]
            if sub.empty:
                continue
            best_idx = sub['Avg_FPR'].idxmin()
            best = sub.loc[best_idx]
            bin_centers.append((lo + hi) / 2)

            fpr_val = best['Avg_FPR']
            min_fprs.append(fpr_val if fpr_val >= 1e-6 else 1e-6)
            annotations.append(f"({best['b']},{best['r']})")

        order = np.argsort(bin_centers)
        bin_centers  = [bin_centers[i]  for i in order]
        min_fprs     = [min_fprs[i]     for i in order]
        annotations  = [annotations[i]  for i in order]

        fig.add_trace(go.Scatter(
            x=bin_centers, y=min_fprs,
            mode='lines+markers', name=label, text=annotations,
            marker=dict(size=10, line=dict(width=1, color='black')),
            line=dict(dash='dash', width=1.5),
            hovertemplate='Recall Bin Center: %{x:.4f}<br>Min FPR: %{y:.4g}<br>%{text}'
        ))

    fig.update_layout(
        title="Min FPR per Recall Bin", xaxis_title="Recall Bin Center", yaxis_title="Best (Min) FPR in Bin (log)",
        template="ggplot2", colorway=px.colors.qualitative.Safe, hoverlabel=dict(namelength=-1),
        legend=dict(title="Experiment", bgcolor='rgba(255,255,255,0.9)')
    )
    fig.update_yaxes(type="log")
    return fig


def plot_ann_query_time_vs_recall(df):
    fig = go.Figure()
    for label, group in get_ann_groups(df):
        text = [f"({b},{r})" for b, r in zip(group['b'], group['r'])] if 'b' in group.columns and 'r' in group.columns else None

        y_col = 'Avg_Query_Time_ms' if 'Avg_Query_Time_ms' in group.columns else 'query_time'
        if y_col not in group.columns:
            continue

        fig.add_trace(go.Scatter(
            x=group['Avg_Recall'], y=group[y_col],
            mode='markers', name=label, text=text,
            marker=dict(size=10, opacity=0.75, line=dict(width=1, color='DarkSlateGrey')),
            hovertemplate='Recall: %{x:.4f}<br>Query Time: %{y:.4f}<br>%{text}' if text else 'Recall: %{x:.4f}<br>Query Time: %{y:.4f}'
        ))
    fig.update_layout(
        title="Avg Query Time vs Avg Recall", xaxis_title="Average Recall", yaxis_title="Average Query Time (ms)",
        template="ggplot2", colorway=px.colors.qualitative.Safe, hoverlabel=dict(namelength=-1),
        legend=dict(title="Experiment", bgcolor='rgba(255,255,255,0.9)')
    )
    return fig


def main():
    parser = argparse.ArgumentParser(
        description="Generate interactive Plotly visualizations for BioLSHasher processed metrics.")
    parser.add_argument("outdatafile", nargs='+', help="Paths to *_processed-*.outdata files")
    parser.add_argument("--show", action="store_true", help="Launch browser with plots immediately")
    args = parser.parse_args()

    collision_df = None
    ann_df = None
    basename = "BioLSHasher"
    out_dir = ""

    for f in args.outdatafile:
        if not os.path.isfile(f):
            print(f"File not found: {f}")
            continue
        df = read_processed_dataframe(f)
        if df.empty:
            continue
        test_name = df['test_name'].iloc[0] if 'test_name' in df.columns else ""

        if "Collision" in test_name:
            collision_df = df
            basename = os.path.basename(f).replace("_processed_coll.csv", "")
            out_dir = os.path.dirname(os.path.abspath(f))
        elif "Neighbour" in test_name:
            ann_df = df
            if out_dir == "":
                basename = os.path.basename(f).replace("_processed_ann.csv", "").replace("_processed.csv", "")
                out_dir = os.path.dirname(os.path.abspath(f))

    replacements = {
        "show_collision":           "display: none;",
        "scatter_div":              "",
        "binned_div":               "",
        "box_div":                  "",
        "show_ann":                 "display: none;",
        "ann_fpr_recall_div":       "",
        "ann_precision_recall_div": "",
        "ann_min_fpr_div":          "",
        "ann_query_time_div":       ""
    }

    _plotlyjs_injected = False

    def _plotlyjs():
        nonlocal _plotlyjs_injected
        if not _plotlyjs_injected:
            _plotlyjs_injected = True
            return True   # embeds Plotly.js; guaranteed to be present before any chart JS
        return False

    if collision_df is not None:
        sim_metric = collision_df['DistanceMetric'].iloc[0] if 'DistanceMetric' in collision_df.columns else ""
        replacements["show_collision"] = "display: block;"
        replacements["binned_div"]     = plot_binned_average(collision_df, sim_metric).to_html(full_html=False, include_plotlyjs=_plotlyjs())
        replacements["scatter_div"]    = plot_scatter(collision_df, sim_metric).to_html(full_html=False, include_plotlyjs=_plotlyjs())
        replacements["box_div"]        = plot_box_plot(collision_df, sim_metric=sim_metric).to_html(full_html=False, include_plotlyjs=_plotlyjs())

    if ann_df is not None:
        replacements["show_ann"] = "display: block;"
        replacements["ann_fpr_recall_div"]       = plot_ann_fpr_vs_recall(ann_df).to_html(full_html=False, include_plotlyjs=_plotlyjs())
        replacements["ann_precision_recall_div"] = plot_ann_precision_vs_recall(ann_df).to_html(full_html=False, include_plotlyjs=_plotlyjs())
        replacements["ann_min_fpr_div"]          = plot_ann_min_fpr_per_bin(ann_df).to_html(full_html=False, include_plotlyjs=_plotlyjs())
        replacements["ann_query_time_div"]       = plot_ann_query_time_vs_recall(ann_df).to_html(full_html=False, include_plotlyjs=_plotlyjs())

    if collision_df is None and ann_df is None:
        print("No valid CSV files provided or files did not match expected test_names.")
        sys.exit(1)

    import string
    template_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "dashboard_template.html")
    with open(template_path, "r", encoding="utf-8") as f:
        html_template_raw = f.read()

    html_template = string.Template(html_template_raw).safe_substitute(**replacements)

    dashboard_path = os.path.join(out_dir, f"{basename}_dashboard.html")
    with open(dashboard_path, "w", encoding="utf-8") as f:
        f.write(html_template)

    print(f"Dashboard generated successfully at: {dashboard_path}")

    if args.show:
        import webbrowser
        url = 'file://' + os.path.realpath(dashboard_path)
        try:
            if not webbrowser.open(url):
                raise RuntimeError("webbrowser.open returned False")
        except Exception as e:
            print(f"[warn] Could not launch browser automatically ({e}).")
            print(f"       Open manually: {url}")


if __name__ == "__main__":
    main()