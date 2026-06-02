# report_utils.py
from pathlib import Path
from dataclasses import dataclass
import pandas as pd
import matplotlib.pyplot as plt
import io, base64
import numpy as np
import ast
import seaborn as sns
from matplotlib.ticker import FormatStrFormatter
from matplotlib.lines import Line2D
import zipfile 

# ----------------------------------------------------------
@dataclass
class Config:
    AMPlify_TSV : Path 
    tAMPer_CSV  : Path
    PDB_DIR     : Path 
    Report_HTML : Path
    LOGO_PNG    : Path
    TEMPLATE_PATH : Path 

    def __init__(self, amplify_tsv, tamper_csv, pdb_dir, report_html, logo_png, template_path):
        self.AMPlify_TSV = Path(amplify_tsv)
        self.tAMPer_CSV = Path(tamper_csv)
        self.PDB_DIR = Path(pdb_dir)
        self.Report_HTML = Path(report_html)
        self.LOGO_PNG = logo_png
        self.TEMPLATE_PATH = template_path

# ----------------------------------------------------------
def load_AMPlify_results(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str)

def load_tAMPer_results(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, dtype=str)

# ----------------------------------------------------------
def build_summary_table(amplify_df: pd.DataFrame,
                        tamper_df:  pd.DataFrame) -> pd.DataFrame:
    merged = (
        amplify_df
        .merge(tamper_df, left_on="Sequence_ID", right_on="id",
               how="inner", suffixes=("", "_t"))
        .rename(columns={
            "Probability_score": "AMPlify_Overall_Probability",
            "Prediction":        "AMPlify_Prediction",
            "score":             "tAMPer_Toxicity_Probability",
            "prediction":        "tAMPer_Prediction",
        })
    )

    # map 0/1 → Non-toxic / Toxic
    merged["tAMPer_Prediction"] = (
        merged["tAMPer_Prediction"]
        .replace({"0.0": "Non-toxic", "1.0": "Toxic",
                  0.0: "Non-toxic",  1.0: "Toxic"})
    )

    for col in ["AMPlify_Overall_Probability",
                "tAMPer_Toxicity_Probability"]:
        merged[col] = (
            merged[col]
            .astype(float)
            .round(3)
            .map(lambda x: f"{x:.3f}")
        )

    return merged[[
        "Sequence_ID", "Sequence", "Length", "Charge",
        "AMPlify_Overall_Probability", "tAMPer_Toxicity_Probability",
        "AMPlify_Prediction", "tAMPer_Prediction"
    ]]

# ----------------------------------------------------------
def plot_probability_scatter_individual(summary_df: pd.DataFrame,
                                        pad            = 0.02,
                                        dpi            = 500,
                                        figsize        = (6, 6),
                                        label_fontsize = 14,
                                        title_fontsize = 16,
                                        tick_fontsize  = 12):
    """
    Return {Sequence_ID: {'png': b64, 'svg': b64}}
    A PNG/SVG pair where the current sequence is in red and
    all others in black.
    """
    out   = {}
    x_all = summary_df["AMPlify_Overall_Probability"].astype(float)
    y_all = summary_df["tAMPer_Toxicity_Probability"].astype(float)

    for _, row in summary_df.iterrows():
        sid   = row["Sequence_ID"]
        x_cur = float(row["AMPlify_Overall_Probability"])
        y_cur = float(row["tAMPer_Toxicity_Probability"])

        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

        ax.scatter(x_all, y_all,
                   c="black", alpha=1,
                   edgecolors="w", linewidth=0.5,
                   s=100)
        ax.scatter([x_cur], [y_cur],
                   c="red", edgecolors="w",
                   linewidth=0.75, s=100, zorder=5)

        ax.set_xlim(-pad, 1 + pad)
        ax.set_ylim(-pad, 1 + pad)

        ax.set_xlabel("AMPlify Score",
                      fontsize=label_fontsize, fontweight="bold")
        ax.set_ylabel("tAMPer Score",
                      fontsize=label_fontsize, fontweight="bold")
        ax.set_title(sid,
                     fontsize=title_fontsize, fontweight="bold")

        ax.tick_params(axis='both', which='major',
                       labelsize=tick_fontsize)

        ax.plot([0, 1], [0, 1],
                c="gray", ls="--", alpha=0.2)

        handles, labs = plt.gca().get_legend_handles_labels()
        point_of_interest = Line2D([0], [0], label=sid, marker='o', markersize=10, 
                        markeredgecolor='r', markerfacecolor='r', linestyle='')
        all_points = Line2D([0], [0], label="Other Peptides", marker='o', markersize=10, 
                        markeredgecolor='black', markerfacecolor='black', linestyle='')

        # add manual symbols to auto legend
        handles.extend([point_of_interest, all_points])
        ax.legend(handles=handles)

        # encode PNG
        buf_png = io.BytesIO()
        fig.savefig(buf_png, format="png", bbox_inches="tight")
        b64_png = base64.b64encode(buf_png.getvalue()).decode("ascii")

        # encode SVG
        buf_svg = io.BytesIO()
        fig.savefig(buf_svg, format="svg", bbox_inches="tight")
        b64_svg = base64.b64encode(buf_svg.getvalue()).decode("ascii")

        plt.close(fig)
        out[sid] = {"png": b64_png, "svg": b64_svg}

    return out

# ----------------------------------------------------------
def plot_attention_heatmap_single(sequence_id,
                                  sequence,
                                  attention_str,
                                  palette_name    = 'mako',
                                  dpi             = 500,
                                  scale           = 2.0,
                                  title_fontsize  = 8,
                                  tick_fontsize   = 8,
                                  legend_labelsize= 6):
    """
    Return {'png': b64, 'svg': b64} with a 1 x L heat-map rod.
    """
    attn = np.asarray(ast.literal_eval(attention_str), float)
    L = len(attn)

    cmap = sns.color_palette(palette_name, as_cmap=True)

    base_w = L * 0.15
    base_h = 1.2
    fig = plt.figure(figsize=(base_w * scale, base_h * scale), dpi=dpi)
    ax  = fig.add_axes([0.1, 0.35, 0.8, 0.5])

    im = ax.imshow(
        attn[np.newaxis, :],
        aspect='auto',
        cmap=cmap,
        vmin=0, vmax=attn.max(),
        extent=(-0.5, L - 0.5, 0, 1)
    )

    step = max(1, L // 20)
    ax.set_xticks(np.arange(0, L, step))
    ax.set_xticklabels(
        [sequence[i] for i in range(0, L, step)],
        fontweight='bold',
        fontsize=tick_fontsize * scale
    )
    ax.set_yticks([])

    ax.set_title(sequence_id,
                 fontsize=title_fontsize * scale)

    cax = fig.add_axes([0.25, 0.02, 0.5, 0.06])
    cbar = fig.colorbar(im, cax=cax, orientation='horizontal')

    cbar.set_label('AMPlify Attention Score',
                   fontsize=legend_labelsize * scale)
    cbar.ax.tick_params(labelsize=legend_labelsize * scale,
                        pad=2 * scale)
    cbar.ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))

    plt.tight_layout(pad=1.0 * scale)

    buf_png = io.BytesIO()
    fig.savefig(buf_png, format='png', bbox_inches='tight')
    b64_png = base64.b64encode(buf_png.getvalue()).decode('ascii')

    buf_svg = io.BytesIO()
    fig.savefig(buf_svg, format='svg', bbox_inches='tight')
    b64_svg = base64.b64encode(buf_svg.getvalue()).decode('ascii')

    plt.close(fig)
    return {"png": b64_png, "svg": b64_svg}

# ----------------------------------------------------------
def load_pdb_map(summary_df: pd.DataFrame, pdb_dir: Path) -> dict:
    """
    Return {Sequence_ID: base64-encoded PDB text or None},
    searching for '*_relaxed_rank_001_*.pdb' inside
    '{pdb_dir}/{Sequence_ID}.result.zip'.
    """
    out = {}
    pdb_dir = Path(pdb_dir)
    
    for sid in summary_df["Sequence_ID"]:
        zip_path = pdb_dir / f"{sid}.result.zip"
        
        if zip_path.exists():
            try:
                with zipfile.ZipFile(zip_path, 'r') as zf:
                    matching_files = [
                        name for name in zf.namelist() 
                        if name.endswith('.pdb') and '_relaxed_rank_001_' in name
                    ]
                    
                    if matching_files:
                        pdb_bytes = zf.read(matching_files[0])
                        pdb_b64 = base64.b64encode(pdb_bytes).decode("ascii")
                        out[sid] = pdb_b64
                        continue
                        
            except zipfile.BadZipFile:
                print(f"Warning: {zip_path.name} is corrupted or not a valid zip archive.")
                
        out[sid] = None
        
    return out

# ----------------------------------------------------------
def class_breakdown_svg(summary_df: pd.DataFrame,
                        figsize=(8, 4), dpi=500) -> str:
    """
    Two donut charts (AMP vs non-AMP, Toxic vs Non-toxic) with rounded
    ends and the chart title centred in the hole.
    Returns a base64-encoded SVG (transparent background).
    """
    import io, base64, math
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import pandas as pd

    # -- colours ---------------------------------------------------
    green = "#00C04B"
    red   = "#ff2424"
    gray  = "#E5E5EA"

    # -- build binary class counts ---------------------------------
    amp_pred = summary_df["AMPlify_Prediction"].str.strip().str.lower()
    tox_pred = summary_df["tAMPer_Prediction"].str.strip().str.lower()

    amp_counts = pd.Series({
        "AMP":      (amp_pred == "amp").sum(),
        "non-AMP":  (amp_pred != "amp").sum()
    })
    tox_counts = pd.Series({
        "Toxic":      (tox_pred == "toxic").sum(),
        "Non-toxic":  (tox_pred != "toxic").sum()
    })

    # colour order must match counts order
    palette_amp = [green, gray]
    palette_tox = [red,   gray]

    # -- figure / axes ---------------------------------------------
    fig, axes = plt.subplots(
        ncols=2, figsize=figsize, dpi=dpi,
        subplot_kw=dict(aspect="equal")
    )
    fig.patch.set_alpha(0)
    for ax in axes:
        ax.set_facecolor('none')

    ring_width = 0.15          # thickness of the donut ring
    outer_r    = 1.0           # outer radius
    inner_r    = outer_r - ring_width
    mid_r      = inner_r + ring_width / 2      # centre-line of ring
    cap_r      = ring_width / 2 * 0.91              # radius of the cap circles

    # -- helper to add a circular cap at an angle θ (deg) ----------
    def add_cap(ax, theta_deg, color):
        theta = math.radians(theta_deg)
        x = mid_r * math.cos(theta)
        y = mid_r * math.sin(theta)
        cap = mpatches.Circle((x, y), cap_r, color=color, ec="none")
        ax.add_patch(cap)

    # -- helper to draw one donut with rounded ends ----------------
    def donut(ax, counts, colors, title):
        wedges, _ = ax.pie(
            counts.values,
            startangle=90,
            radius=outer_r,
            colors=colors,
            wedgeprops=dict(width=ring_width, edgecolor="white")
        )

        # Add caps only for the coloured slice (first wedge)
        total = counts.sum()
        frac  = counts.iloc[0] / total
        start = 90                      # starting at 12 o’clock
        end   = 90 + frac * 360

        add_cap(ax, start, colors[0])   # cap at start
        add_cap(ax, end,   colors[0])   # cap at end

        # Title in the hole
        ax.text(0, 0, title,
                ha="center", va="center",
                fontsize=14, weight="bold")

        ax.legend(
            wedges,
            [f"{lbl} ({cnt})" for lbl, cnt in counts.items()],
            loc="center left",
            bbox_to_anchor=(1.05, 0, 0.4, 1),
            frameon=False
        )

    donut(axes[0], amp_counts, palette_amp, "AMPlify")
    donut(axes[1], tox_counts, palette_tox, "tAMPer")

    plt.tight_layout()

    # -- export as base64-SVG --------------------------------------
    buf_svg = io.BytesIO()
    fig.savefig(buf_svg, format="svg",
                bbox_inches="tight", transparent=True)
    b64_svg = base64.b64encode(buf_svg.getvalue()).decode("ascii")
    plt.close(fig)
    return b64_svg
