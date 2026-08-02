import argparse
import scanpy as sc
import squidpy as sq
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.cm as cm
from skimage import exposure
from scipy.optimize import nnls
from scipy.sparse import issparse
import math
from scipy.stats import mannwhitneyu
from collections import defaultdict

#/research/groups/plummgrp/projects/Visium_855118_Steve_Mack_Siri_09232025/common/Tuyu_Apr2026/IF_images #IF


THRESHOLD = 6

def marker_unique_excel(exl_file,sheet_name,cluster = "supertype.markers.combo",class_label = "class_label"):
    
    markers_df = pd.read_excel(exl_file,sheet_name=sheet_name)

    #interested in supertype.markers.combo and class_label and class_id_label
    markers_df = markers_df[[cluster,class_label]]
    markers_df = markers_df.dropna(subset=[cluster])
    #markers_df.to_csv("test.csv", index=False)

    markers_df['gene_list'] = markers_df[cluster].str.split(',').apply(
        lambda x: [g.strip() for g in x]
    )

    genes_long = (
        markers_df[[cluster, class_label]]
        .assign(gene=markers_df[cluster].str.split(','))
        .explode('gene')
    )

    genes_long['gene'] = genes_long['gene'].str.strip()
    genes_long = genes_long[genes_long['gene'] != '']  # remove empty strings

    gene_class_count = (
        genes_long.groupby('gene')[class_label]
        .nunique()
        .rename('n_classes')
    )

    unique_genes = gene_class_count[gene_class_count == 1].index

    unique_table = (
        genes_long[genes_long['gene'].isin(unique_genes)]
        [[class_label, 'gene']]
        .drop_duplicates()
        .sort_values([class_label, 'gene'])
        .reset_index(drop=True)
    )


    unique_marker_dict = (
        unique_table.groupby(class_label)['gene']
        .apply(list)
        .to_dict()
    )

    return unique_table, unique_marker_dict

def marker_unique_excel_subclass(exl_file,sheet_name,cluster = "supertype.markers.combo",class_label = "class_label",subclass_label = "subclass_label"):

    """
        split marker unique genes into class and subclass.  Based on Allen Mouse Brain.
        Note: preserving unique genes between supertype but NOT subclass.  May need to change...

    """

    markers_df = pd.read_excel(exl_file,sheet_name=sheet_name)

    markers_df = markers_df[[cluster,class_label,subclass_label]]
    markers_df = markers_df.dropna(subset=[cluster])
    

    markers_df['gene_list'] = markers_df[cluster].str.split(',').apply(
        lambda x: [g.strip() for g in x]
    )

    genes_long = (
        markers_df[[cluster, class_label,subclass_label]]
        .assign(gene=markers_df[cluster].str.split(','))
        .explode('gene')
    )

    genes_long['gene'] = genes_long['gene'].str.strip()
    genes_long = genes_long[genes_long['gene'] != '']  # remove empty strings

    gene_class_count = (
        genes_long.groupby('gene')[class_label]
        .nunique()
        .rename('n_classes')
    )

    unique_genes = gene_class_count[gene_class_count == 1].index

    unique_table = (
        genes_long[genes_long['gene'].isin(unique_genes)]
        [[class_label, subclass_label,'gene']]
        .drop_duplicates()
        .sort_values([class_label, 'gene'])
        .reset_index(drop=True)
    )

  
    unique_marker_dict = defaultdict(lambda: defaultdict(list))

    for _, row in unique_table.iterrows():
        unique_marker_dict[row[class_label]][row[subclass_label]].append(row["gene"])

    unique_marker_dict = {k: dict(v) for k, v in unique_marker_dict.items()}
        
    return unique_table, unique_marker_dict



def score_geneset(adata, genes, name):
    """Mean log-normalised expression of detected genes in the set."""
    present = [g for g in genes if g in adata.var_names]

    if len(present) == 0:
        adata.obs[name] = 0.0
        return 0
    expr = adata[:, present].X
    if hasattr(expr, "toarray"):
        expr = expr.toarray()
    adata.obs[name] = expr.mean(axis=1)
    return len(present)

def minimum_threshold_gene_filter(adata,markers,min_markers=3):

    filtered_markers = {}
    for ct, genes in markers.items():
        present  = [g for g in genes if g in adata.var_names]
        detected = [g for g in present
                if (adata[:, g].X > 0).sum() > adata.n_obs * 0.01]  # >1% spots
        filtered_markers[ct] = detected
        print(f"  {ct:30s}  {len(detected):3d}/{len(genes):3d} detected")

    # Drop cell types with too few markers
    valid_cts = {ct: g for ct, g in filtered_markers.items() if len(g) >= min_markers}
    print(f"\n{len(valid_cts)} cell types with >= {min_markers} detected markers")

    return valid_cts



def plot_celltype_overlay(adata, celltypes, img_key="hires", 
                          spot_size=1.4, alpha=1.0, min_threshold=0.3,
                          figsize=(10, 10),border_color="white"):
    """
    Overlay multiple cell types on the same spatial image.
    Each cell type gets a distinct color; spots are only shown
    if their score exceeds min_threshold (as a percentile).
    """
    library_id = list(adata.uns["spatial"].keys())[0]
    img = adata.uns["spatial"][library_id]["images"][img_key]
    scale = adata.uns["spatial"][library_id]["scalefactors"][f"tissue_{img_key}_scalef"]
    
    # Get spot coordinates (in image pixel space)
    coords = adata.obsm["spatial"] * scale  # shape: (n_spots, 2)
    
    # Assign a color per cell type
    colors = cm.tab20.colors  # 20 distinct colors
    
    fig, ax = plt.subplots(figsize=figsize)
    ax.imshow(img,alpha = 0.4)
    ax.axis("off")
    
    legend_handles = []
    
    for i, ct in enumerate(celltypes):
        if ct not in adata.obs.columns:
            print(f"Warning: {ct} not found in adata.obs, skipping")
            continue
        
        scores = adata.obs[ct].values
        
        # Only plot spots above threshold percentile
        threshold = np.nanpercentile(scores, min_threshold * 100)
        mask = scores > threshold
        
        # Normalize scores for alpha/size scaling
        s_norm = (scores[mask] - scores[mask].min()) / (scores[mask].max() - scores[mask].min() + 1e-8)
        
        color = colors[i % len(colors)]
        
        # ── Draw border layer first (slightly larger, neon color) ─────────────────
        scatter = ax.scatter(
            coords[mask, 0],
            coords[mask, 1],
            s=spot_size * (1 + s_norm * 3), ##spot_size * 1.6, #spot_size * (1 + s_norm * 3),        # slightly bigger than fill
            c=border_color,
            linewidths=1,
            alpha=alpha,
            zorder=2,
        )

        legend_handles.append(scatter)
        
        # ── Draw filled spots on top ───────────────────────────────────────────────
        sc_fill = ax.scatter(
            coords[mask, 0], 
            coords[mask, 1],
            s=spot_size * (1 + s_norm * 3),
            c=[color],
            linewidths=0,
            alpha=alpha,
            zorder=3,
            label=ct, 
        )
        legend_handles.append(sc_fill)
    
    ax.legend(
        loc="upper right",
        bbox_to_anchor=(1.3, 1),
        fontsize=7,
        markerscale=2,
        framealpha=0.8,
    )
    plt.title("Cell type overlay", fontsize=12)
    plt.tight_layout()
    return fig



def plot_spatial_with_border(adata, ct, img_key="hires", 
                              cmap="magma", spot_size=8,
                              border_color="lime", border_width=0.8,
                              min_percentile=50, figsize=(8, 8)):
    """
    Plot spatial cell type scores with a neon border around each spot.
    """
    library_id = list(adata.uns["spatial"].keys())[0]
    img = adata.uns["spatial"][library_id]["images"][img_key]
    scale = adata.uns["spatial"][library_id]["scalefactors"][f"tissue_{img_key}_scalef"]
    
    coords = adata.obsm["spatial"] * scale  # pixel coordinates
    scores = adata.obs[ct].values

    # Only show spots above threshold
    threshold = np.nanpercentile(scores, min_percentile)
    mask = scores > threshold

    # Normalize scores for colormap
    vmin, vmax = scores[mask].min(), np.nanpercentile(scores[mask], min_percentile)
    norm = plt.Normalize(vmin=vmin, vmax=vmax)
    colors = plt.get_cmap(cmap)(norm(scores[mask]))

    fig, ax = plt.subplots(figsize=figsize)
    ax.imshow(img,alpha = 0.4)

    # ── Draw border layer first (slightly larger, neon color) ─────────────────
    ax.scatter(
        coords[mask, 0], coords[mask, 1],
        s=spot_size * 1.6,        # slightly bigger than fill
        c=border_color,
        linewidths=0,
        zorder=2,
    )

    # ── Draw filled spots on top ───────────────────────────────────────────────
    sc_fill = ax.scatter(
        coords[mask, 0], coords[mask, 1],
        s=spot_size,
        c=colors,
        linewidths=0,
        zorder=3,
    )

    # Colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    plt.colorbar(sm, ax=ax, shrink=0.4, label=ct)

    ax.set_title(ct, fontsize=12)
    ax.axis("off")
    plt.tight_layout()

    return fig


def plot_spatial_with_border_v2(adata, ct, ax, fig, img_key="hires",
                              cmap="magma", spot_size=8,
                              border_color="lime", border_width=0.8,
                              min_percentile=50):
    """
    Plot spatial cell type scores with a neon border around each spot.
    """
    library_id = list(adata.uns["spatial"].keys())[0]
    img = adata.uns["spatial"][library_id]["images"][img_key]
    scale = adata.uns["spatial"][library_id]["scalefactors"][f"tissue_{img_key}_scalef"]

    coords = adata.obsm["spatial"] * scale  # pixel coordinates
    scores = adata.obs[ct].values

    #safe guard empy scores
    if np.all(scores == 0) or np.all(np.isnan(scores)):
        ax.imshow(img,alpha = 0.4)
        ax.set_title(f"{ct}\n(no expression)", fontsize=8)
        ax.axis("off")
        return


    # Only show spots above threshold
    threshold = np.nanpercentile(scores, min_percentile)
    mask = scores > threshold

    # Normalize scores for colormap
    vmin, vmax = scores[mask].min(), np.nanpercentile(scores[mask], min_percentile)
    norm = plt.Normalize(vmin=vmin, vmax=vmax)
    colors = plt.get_cmap(cmap)(norm(scores[mask]))

    #fig, ax = plt.subplots(figsize=figsize)
    ax.imshow(img,alpha = 0.4)

    # ── Draw border layer first (slightly larger, neon color) ─────────────────
    ax.scatter(
        coords[mask, 0], coords[mask, 1],
        s=spot_size * 1.6,        # slightly bigger than fill
        c=border_color,
        linewidths=0,
        zorder=2,
    )

    # ── Draw filled spots on top ───────────────────────────────────────────────
    sc_fill = ax.scatter(
        coords[mask, 0], coords[mask, 1],
        s=spot_size,
        c=colors,
        linewidths=0,
        zorder=3,
    )

    # Colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    plt.colorbar(sm, ax=ax, shrink=0.4, label=ct)

    ax.set_title(ct, fontsize=12)
    ax.axis("off")
    plt.tight_layout()

    return fig


def weighted_score(adata, genes, min_detection=0.01):
    """
    Weight each gene by its detection rate across spots.
    Downweights genes that are rarely detected (unreliable at 8um).
    """
    present = [g for g in genes if g in adata.var_names]
    if not present:
        return np.zeros(adata.n_obs)

    expr = adata[:, present].X
    if hasattr(expr, 'toarray'):
        expr = expr.toarray()

    # Detection rate per gene
    detection = (expr > 0).mean(axis=0)

    # Only use genes detected in >1% of spots
    mask = detection >= min_detection
    if mask.sum() == 0:
        return np.zeros(adata.n_obs)

    expr     = expr[:, mask]
    detection = detection[mask]

    # Weight by detection rate
    weights  = detection / detection.sum()

    return expr @ weights


def nnls_deconvolve(adata, markers_dict, n_jobs=-1):
    """
    Non-negative least squares deconvolution.
    Fits cell type proportions to observed expression.
    """
    # Build signature matrix — mean expression per cell type
    # Using only genes in your marker dict
    all_genes = list({g for genes in markers_dict.values()
                      for g in genes if g in adata.var_names})
    ct_names  = list(markers_dict.keys())

    # Signature: average normalized expression of marker genes per cell type
    # Here we use detection-weighted mean as proxy for signature
    sig = np.zeros((len(all_genes), len(ct_names)))
    gene_idx = {g: i for i, g in enumerate(all_genes)}

    for j, (ct, genes) in enumerate(markers_dict.items()):
        present = [g for g in genes if g in adata.var_names]
        for g in present:
            expr = adata[:, g].X
            if issparse(expr):
                expr = expr.toarray().flatten()
            sig[gene_idx[g], j] = expr.mean()

    # Get spot expression for those genes
    X = adata[:, all_genes].X
    if issparse(X):
        X = X.toarray()

    # NNLS per spot
    props = np.zeros((adata.n_obs, len(ct_names)))
    for i in range(adata.n_obs):
        props[i], _ = nnls(sig, X[i])

    # Normalize to sum to 1
    row_sums = props.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1
    props = props / row_sums

    prop_df = pd.DataFrame(props, index=adata.obs_names,
                           columns=[f"nnls_{ct}" for ct in ct_names])
    return prop_df


def high_confidence_assignment(adata, valid_cts,
                                score_pct=75,    # percentile threshold
                                min_methods=2):  # must pass N methods
    """
    Assign cell type only if multiple scoring methods agree.
    """
    ct_names = list(valid_cts.keys())
    n_spots  = adata.n_obs
    votes    = np.zeros((n_spots, len(ct_names)), dtype=int)

    for j, ct in enumerate(ct_names):
        # Method A — scanpy score
        if f"score_{ct}" in adata.obs:
            s = adata.obs[f"score_{ct}"].values
            thresh = np.nanpercentile(s, score_pct)
            votes[:, j] += (s >= thresh).astype(int)

        # Method B — weighted score
        if f"wscore_{ct}" in adata.obs:
            s = adata.obs[f"wscore_{ct}"].values
            thresh = np.nanpercentile(s, score_pct)
            votes[:, j] += (s >= thresh).astype(int)

        # Method C — NNLS proportion
        if f"nnls_{ct}" in adata.obs:
            s = adata.obs[f"nnls_{ct}"].values
            thresh = np.nanpercentile(s, score_pct)
            votes[:, j] += (s >= thresh).astype(int)

    # Assign cell type if >= min_methods agree
    max_votes   = votes.max(axis=1)
    best_ct_idx = votes.argmax(axis=1)

    assignments = np.where(
        max_votes >= min_methods,
        np.array(ct_names)[best_ct_idx],
        "unassigned"
    )

    # Confidence score = fraction of methods agreeing / max possible
    confidence = max_votes / 3.0

    adata.obs['celltype_assignment']  = assignments
    adata.obs['assignment_confidence'] = confidence
    adata.obs['assignment_votes']      = max_votes

    # Summary
    counts = pd.Series(assignments).value_counts()
    print(f"\nAssignment summary (min {min_methods}/3 methods):")
    print(counts.to_string())
    print(f"\nUnassigned: {(assignments=='unassigned').sum()} "
          f"({(assignments=='unassigned').mean()*100:.1f}%)")
    print(f"High confidence (all 3 agree): {(max_votes==3).sum()} spots")

    return adata


def main():

    '''
    continuation after spatial_pipeline.py
    automation version of integrative_analysis_5_23_2026.ipynb
    '''

    parser = argparse.ArgumentParser()
    parser.add_argument('adata_name', type=str, help="path to adata output from spatial_pipeline.py")
    parser.add_argument('label',type=str,help="output label")
    #parser.add_argument('min_pct',type=float,help="score must exceed minimum pct threshold to be visible")
    args = parser.parse_args()

    
    adata = sc.read_h5ad(args.adata_name)

    marker_table_allen_subclass,marker_allen_subclass = marker_unique_excel_subclass('/research/groups/mackgrp/home/demerson/AIB/yao_et_al_Nature_2023.xlsx','supertype_annotation')
    marker_table_allen,marker_allen = marker_unique_excel('/research/groups/mackgrp/home/demerson/AIB/yao_et_al_Nature_2023.xlsx','supertype_annotation')
    
    #remove ct with < 3 genes in adata
    marker_allen = minimum_threshold_gene_filter(adata,marker_allen)

    marker_allen['gfp_probe'] = ['Gfp', 'GFP', 'EGFP', 'eGFP', 'gfp', 'egfp','EGFP-L10a', 'Gfp-L10a', 'td-Tomato', 'mCherry']
    tumor = {}
    instance = pd.read_excel("/research/groups/mackgrp/home/demerson/spatial_rabies/93_core_gene_set_commonly_upregulated_in_human_and_mouse_C11orf95_RELAFUS1_driven_ependymoma.xlsx",sheet_name="Sheet1",header=None)
    tumor['C11orf95_relafus1'] = instance[0].str.capitalize().tolist()

    #method 1 scoring
    for ct, genes in marker_allen.items():
        #n = score_geneset(adata, genes, ct)
        sc.tl.score_genes(
            adata,
            gene_list=genes,
            score_name= f"tlscore_{ct}",
            use_raw=False
        )

    for ct, genes in tumor.items():
        sc.tl.score_genes(
            adata,
            gene_list=genes,
            score_name= f"tltumscore_{ct}",
            use_raw=False
        )

    #method 2 scoring
    #for ct, genes in marker_allen.items():
        #adata.obs[f"wscore_{ct}"] = weighted_score(adata, genes)

    #method 3 scoring
    #prop_df = nnls_deconvolve(adata, marker_allen)
    #adata.obs = adata.obs.join(prop_df)

    #combine all 3 method 4
    #adata = high_confidence_assignment(adata, marker_allen,
    #                               score_pct=75, min_methods=2)        


    score_cols = [c for c in adata.obs.columns if c.startswith("tlscore_")]
    scores = adata.obs[score_cols].copy()
    #scores.columns = [c.replace("score_", "") for c in score_cols]

    top2 = np.sort(scores.values, axis=1)[:, -2:]
    margin = top2[:, 1] - top2[:, 0]   # best minus second-best
    adata.obs["tlcelltype_assigned"] = scores.idxmax(axis=1)
    adata.obs["tlcelltype_score"] = scores.max(axis=1)
    adata.obs["tlcelltype_margin"] = margin

    # Require both absolute score AND separation from runner-up
    low_conf = (adata.obs["tlcelltype_score"] < 0.1) | (adata.obs["tlcelltype_margin"] < 0.05)
    adata.obs.loc[low_conf, "tlcelltype_assigned"] = "Unknown"

    gfp_adata = adata[adata.obs['tlcelltype_assigned'] == "tlscore_gfp_probe"]

    score_cols = [c for c in gfp_adata.obs.columns if c.startswith("tlscore_")]
    scores_df = gfp_adata.obs[score_cols].copy()
    #gfp removed from consideration
    scores_df = scores_df.drop(columns=["tlscore_gfp_probe"])

    top2 = np.sort(scores_df.values, axis=1)[:, -2:]
    margin = top2[:, 1] - top2[:, 0]   # best minus second-best
    gfp_adata.obs["tlcelltype_assigned"] = scores_df.idxmax(axis=1)
    gfp_adata.obs["tlcelltype_score"] = scores_df.max(axis=1)
    gfp_adata.obs["tlcelltype_margin"] = margin

    low_conf = (gfp_adata.obs["tlcelltype_score"] < 0.1) | (gfp_adata.obs["tlcelltype_margin"] < 0.05)
    gfp_adata.obs.loc[low_conf, "tlcelltype_assigned"] = "Unknown"
    gfp_adata_mod =  gfp_adata[gfp_adata.obs["tlcelltype_assigned"] != "Unknown"]
    counts = gfp_adata_mod.obs["tlcelltype_assigned"].value_counts()
    total = counts.sum()
    pcts = counts / total * 100

    fig, axes = plt.subplots(1, 2, figsize=(18, 6))

    counts.plot(
        ax=axes[0],
        kind="bar",
        color="steelblue"
        )

    axes[0].set_xlabel("predicted cell type highest score")
    axes[0].set_ylabel("Count")
    axes[0].set_title("Most frequent max-scoring cell type at GFP probe spot")

    cumulative = 0
    colors = plt.cm.tab20.colors

    for i, (celltype, pct) in enumerate(pcts.items()):
        axes[1].bar(0, pct, bottom=cumulative, color=colors[i % len(colors)],
           width=0.5, label=f"{celltype} ({pct:.1f}%)")
        # Only label segments large enough to read
        if pct > 2:
            axes[1].text(
                0, cumulative + pct / 2,
                f"{celltype}\n{pct:.1f}%",
                ha="center", va="center",
                fontsize=7, color="white", fontweight="bold",
            )
        cumulative += pct

    axes[1].set_xlim(-0.5, 0.5)
    axes[1].set_ylim(0, 100)
    axes[1].set_xticks([])
    axes[1].set_ylabel("Percentage of GFP+ spots (%)")
    axes[1].set_title("GFP+ spots assigned to\nmax-scoring cell type")
    axes[1].legend(
        loc="upper right",
        bbox_to_anchor=(2.0, 1),
        fontsize=8,
        title="Cell type",
    )

    plt.tight_layout()
    plt.savefig("spatial_gfp_spots_top_score_tlscore_" + args.label + "_bar_conf.png",dpi=300,bbox_inches="tight")
    plt.clf()


    # Further subset gfp + immune into subclass


    gfp_immune_adata = gfp_adata[gfp_adata.obs["tlcelltype_assigned"] == "tlscore_Immune"]
    
    #marker_allen_subclass["Immune"]

    print("gfp_immune_adata.var_names",gfp_immune_adata.var_names)
    print(marker_allen_subclass["Immune"])

    for ct, genes in marker_allen_subclass["Immune"].items():
        #n = score_geneset(adata, genes, ct)
        sc.tl.score_genes(
            gfp_immune_adata,
            gene_list=genes,
            score_name= f"tlscore_GFP_Immune_{ct}",
            use_raw=False
        )

    score_cols = [c for c in gfp_immune_adata.obs.columns if c.startswith("tlscore_GFP_Immune_")]
    scores_df = gfp_immune_adata.obs[score_cols].copy()

    top2 = np.sort(scores_df.values, axis=1)[:, -2:]
    margin = top2[:, 1] - top2[:, 0]   # best minus second-best
    gfp_immune_adata.obs["tlcelltype_subclass_assigned"] = scores_df.idxmax(axis=1)
    gfp_immune_adata.obs["tlcelltype_subclass_score"] = scores_df.max(axis=1)
    gfp_immune_adata.obs["tlcelltype_subclass_margin"] = margin

    low_conf = (gfp_immune_adata.obs["tlcelltype_subclass_score"] < 0.1) | (gfp_immune_adata.obs["tlcelltype_subclass_margin"] < 0.05)
    gfp_immune_adata.obs.loc[low_conf, "tlcelltype_subclass_assigned"] = "Unknown"
    gfp_immune_adata_mod =  gfp_immune_adata[gfp_immune_adata.obs["tlcelltype_subclass_assigned"] != "Unknown"]
    counts = gfp_immune_adata.obs["tlcelltype_subclass_assigned"].value_counts()

    ax = counts.plot(
        kind="bar",
        figsize=(12, 6))

    ax.set_xlabel("subclass type highest score")
    ax.set_ylabel("Frequency")
    ax.set_title("Most frequent max-scoring immune subclass at GFP probe spot (conf separation)")
    plt.tight_layout()
    plt.savefig("spatial_gfp_spots_top_score_tlscore_" + args.label + "_bar_conf_immune_subclass.png",dpi=300,bbox_inches="tight")
    plt.clf()
    
  
    gfp_immune_adata_tlscore_conf = gfp_immune_adata.copy()
 
    # Reassign gfp probe to max of remaining scores and plot count (no conf separation)

    celltype_scores = list(marker_allen.keys())

    #Inspecting all GFP probed spots (score > 0) and analyzing which class has the highest score excluding GFP_probe

    gfp_adata = adata[adata.obs['tlcelltype_assigned'] == "tlscore_gfp_probe"].copy()

    score_cols = [c for c in gfp_adata.obs.columns if c.startswith("tlscore_")]
    scores_df = gfp_adata.obs[score_cols].copy()
    #gfp removed from consideration
    scores_df = scores_df.drop(columns=["tlscore_gfp_probe"])

    top2 = np.sort(scores_df.values, axis=1)[:, -2:] 
    margin = top2[:, 1] - top2[:, 0]   # best minus second-best
    gfp_adata.obs["tlcelltype_assigned"] = scores_df.idxmax(axis=1)
    gfp_adata.obs["tlcelltype_score"] = scores_df.max(axis=1)
    gfp_adata.obs["tlcelltype_margin"] = margin

    #low_conf = (gfp_adata.obs["tlcelltype_score"] < 0.1) | (gfp_adata.obs["tlcelltype_margin"] < 0.05)
    #gfp_adata.obs.loc[low_conf, "tlcelltype_assigned"] = "Unknown"
    gfp_adata_mod =  gfp_adata[gfp_adata.obs["tlcelltype_assigned"] != "Unknown"]
    counts = gfp_adata_mod.obs["tlcelltype_assigned"].value_counts()
    total = counts.sum()
    pcts = counts / total * 100

    fig, axes = plt.subplots(1, 2, figsize=(18, 6))

    counts.plot(
        ax=axes[0],
        kind="bar",
        color="steelblue"
        )    

    axes[0].set_xlabel("predicted cell type highest score")
    axes[0].set_ylabel("Count")
    axes[0].set_title("Most frequent max-scoring cell type at GFP probe spot")


    cumulative = 0
    colors = plt.cm.tab20.colors

    for i, (celltype, pct) in enumerate(pcts.items()):
        axes[1].bar(0, pct, bottom=cumulative, color=colors[i % len(colors)],
           width=0.5, label=f"{celltype} ({pct:.1f}%)")
        # Only label segments large enough to read
        if pct > 2:
            axes[1].text(
                0, cumulative + pct / 2,
                f"{celltype}\n{pct:.1f}%",
                ha="center", va="center",
                fontsize=7, color="white", fontweight="bold",
            )
        cumulative += pct

    axes[1].set_xlim(-0.5, 0.5)
    axes[1].set_ylim(0, 100)
    axes[1].set_xticks([])
    axes[1].set_ylabel("Percentage of GFP+ spots (%)")
    axes[1].set_title("GFP+ spots assigned to\nmax-scoring cell type")
    axes[1].legend(
        loc="upper right",
        bbox_to_anchor=(2.0, 1),
        fontsize=8,
        title="Cell type",
    )

    plt.tight_layout()
    plt.savefig("spatial_gfp_spots_top_score_tlscore_" + args.label + "_bar.png",dpi=300,bbox_inches="tight")

    # Further subset gfp + immune into subclass

    
    gfp_immune_adata = gfp_adata[gfp_adata.obs["tlcelltype_assigned"] == "tlscore_Immune"]
    

    #marker_allen_subclass["Immune"]
    for ct, genes in marker_allen_subclass["Immune"].items():
        #n = score_geneset(adata, genes, ct)
        sc.tl.score_genes(
            gfp_immune_adata,
            gene_list=genes,
            score_name= f"tlscore_GFP_Immune_{ct}",
            use_raw=False
        )

    score_cols = [c for c in gfp_immune_adata.obs.columns if c.startswith("tlscore_GFP_Immune_")]
    scores_df = gfp_immune_adata.obs[score_cols].copy()

    top2 = np.sort(scores_df.values, axis=1)[:, -2:]
    margin = top2[:, 1] - top2[:, 0]   # best minus second-best
    gfp_immune_adata.obs["tlcelltype_subclass_assigned"] = scores_df.idxmax(axis=1)
    gfp_immune_adata.obs["tlcelltype_subclass_score"] = scores_df.max(axis=1)
    gfp_immune_adata.obs["tlcelltype_subclass_margin"] = margin

    #low_conf = (gfp_immune_adata.obs["tlcelltype_subclass_score"] < 0.1) | (gfp_immune_adata.obs["tlcelltype_subclass_margin"] < 0.05)
    #gfp_immune_adata.obs.loc[low_conf, "tlcelltype_subclass_assigned"] = "Unknown"
    #gfp_immune_adata_mod =  gfp_immune_adata[gfp_immune_adata.obs["tlcelltype_subclass_assigned"] != "Unknown"]
    counts = gfp_immune_adata.obs["tlcelltype_subclass_assigned"].value_counts()

    ax = counts.plot(
        kind="bar",
        figsize=(12, 6))

    ax.set_xlabel("subclass type highest score")
    ax.set_ylabel("Frequency")
    ax.set_title("Most frequent max-scoring subclass at GFP probe spot")
    plt.tight_layout()
    plt.savefig("spatial_gfp_spots_top_score_tlscore_" + args.label + "_bar_immune_subclass.png",dpi=300,bbox_inches="tight")
    

    gfp_immune_adata_tlscore = gfp_immune_adata.copy()

    adata.uns["spatial"].keys()
    library_id = list(adata.uns["spatial"].keys())[0]
    adata.uns["spatial"][library_id].keys()
    adata.uns["spatial"][library_id]["scalefactors"]


    scores = adata.obs['tlscore_gfp_probe'].values
    threshold = np.nanpercentile(scores, THRESHOLD)

    mask = scores > threshold

    coords = adata.obsm["spatial"]

    mean_x = coords[mask, 0].mean()
    mean_y = coords[mask, 1].mean()

    #plot around GFP probe spots

    sq.pl.spatial_scatter(
        adata,
        img_res_key="hires",
        img_alpha = 0.4,
        color='tlscore_gfp_probe',
        size=5,
        alpha=1,
        dpi=300,
        title='tlscore gfp probe',
    )
  
    plt.savefig("celltype_gfp_tlscore_" + args.label + ".png",dpi=300,bbox_inches="tight",facecolor="white")

    # plot fusion scored spot

    sq.pl.spatial_scatter(
        adata,
        img_res_key="hires",
        img_alpha = 0.4,
        color='tltumscore_C11orf95_relafus1',
        size=5,
        alpha=1,
        dpi=300,
        title='tlscore C11orf95_relafus1',
    )

    plt.savefig("celltype_fusion_c11orf95_relafus1_tlscore_" + args.label + ".png",dpi=300,bbox_inches="tight",facecolor="white")

    #normalized plots for tlscore

    score_cols = ['tlscore_gfp_probe', 'tlscore_Astro-Epen', 'tlscore_Immune', 
              'tlscore_CNU-LGE GABA', 'tltumscore_C11orf95_relafus1']
    score_cols_norm = []

    norm_threshold = THRESHOLD/100 #0.07  # normalized score threshold


    for col in score_cols:
        if col not in adata.obs.columns:
            print(f"Skipping '{col}': not found in adata.obs")
            continue

        norm_col  = f"{col}_norm"

        s = adata.obs[col]
        adata.obs[norm_col] = (s - s.min()) / (s.max() - s.min())
        print(f"Normalized '{col}' → '{norm_col}'")

        adata.obs[norm_col] = adata.obs[norm_col].where(
            adata.obs[norm_col] > norm_threshold
        )
        n_passing = (adata.obs[norm_col] >= norm_threshold).sum()
        n_total   = len(adata.obs[norm_col])
        print(f"Threshold applied → '{norm_col}' | {n_passing}/{n_total} spots pass ({100*n_passing/n_total:.1f}%)")
        score_cols_norm.append(norm_col)


    #plot normalized gfp score by itself...

    sq.pl.spatial_scatter(
        adata,
        img_res_key="hires",
        img_alpha = 0.4,
        color='tlscore_gfp_probe_norm',
        size=5.0,
        alpha=1.0,
        dpi=300,
        title='tlscore gfp probe norm',
    )
    plt.savefig("celltype_gfp_tlscore_" + args.label + "_norm.png",dpi=300,bbox_inches="tight",facecolor="white")

    sq.pl.spatial_scatter(
        adata,
        img_res_key="hires",
        img_alpha = 0.4,
        color='tltumscore_C11orf95_relafus1_norm',
        size=5.0,
        alpha=1.0,
        dpi=300,
        title='tlscore C11orf95_relafus1 norm',
    )

    plt.savefig("celltype_fusion_c11orf95_relafus1_tlscore_" + args.label + "_norm.png",
            dpi=300,bbox_inches="tight",facecolor="white")
       
    
    x_center = mean_x
    y_center = mean_y
    half_width = 1000
    half_height = 1000

    x_min = x_center - half_width
    y_min = y_center - half_height
    x_max = x_center + half_width
    y_max = y_center + half_height


    crop = (x_min, y_min, x_max, y_max)

    n = len(score_cols_norm)
    ncols = 2  # fewer columns → larger panels
    nrows = math.ceil(n / ncols)
    sq.pl.spatial_scatter(
        adata,
        img_res_key="hires",
        img_alpha = 0.4,
        color=score_cols_norm,
        size=5.0,
        alpha=1,
        crop_coord=crop,
        ncols=ncols,
        figsize=(3 * ncols, 3 * nrows),  # ~6 inches per panel
        wspace=0.01,  #0.15
        hspace=0.15,  #0.15
        dpi=300,
        vmin= 0.01,
        vmax= 1.0,
        title=score_cols_norm,
    )
    
    plt.savefig("celltype_grid_tlscore_gfp_zoom_" + args.label + "_" + str(x_center) + "_" + str(y_center) + "_norm.png",
            dpi=300,bbox_inches="tight",facecolor="white")


    #Define threshold and binary column

    score_col = "tlscore_gfp_probe_norm"
    base_col  = "tltumscore_C11orf95_relafus1_norm"  
    microglia = "tlscore_GFP_Immune_Microglia NN"
    threshold = np.nanpercentile(adata.obs[score_col], THRESHOLD)

    mask = adata.obs[score_col].values > threshold
    mask_microglia = gfp_immune_adata_tlscore_conf.obs["tlcelltype_subclass_assigned"] == microglia

    #Get scaled coordinates
    library_id = list(adata.uns["spatial"].keys())[0]
    scale = adata.uns["spatial"][library_id]["scalefactors"]["tissue_hires_scalef"]
    coords = adata.obsm["spatial"] * scale

    #Get scaled coordinates for microglia
    library_id_microglia = list(gfp_immune_adata_tlscore_conf.uns["spatial"].keys())[0]
    scale_microglia = gfp_immune_adata_tlscore_conf.uns["spatial"][library_id]["scalefactors"]["tissue_hires_scalef"]
    coords_immunity = gfp_immune_adata_tlscore_conf.obsm["spatial"] * scale

    #Base layer: continuous score of your marker
    fig, ax = plt.subplots(figsize=(10, 10))

    sq.pl.spatial_scatter(
        adata,
        img_res_key="hires",
        img_alpha = 0.4,
        color=base_col,
        ax=ax,
        title=f"{base_col} + GFP positive overlay",
    )

    #Overlay: solid red only where gfp_positive
    ax.scatter(
        coords[mask, 0],
        coords[mask, 1],
        c="red",
        s=8,
        alpha=1,
        zorder=5,
        label=f"GFP+ (n={mask.sum()})",
    )

    #Overlay: solid orange only where gfp_positive + immune
    ax.scatter(
        coords_immunity[:, 0],
        coords_immunity[:, 1],
        c="orange",
        s=8,
        alpha=1,
        zorder=5,
        label=f"GFP+ immune (n={len(gfp_immune_adata_tlscore_conf.obs)})",
    )


    #Overlay: solid yellow only where gfp_positive + microglia + immune
    ax.scatter(
        coords_immunity[mask_microglia, 0],
        coords_immunity[mask_microglia, 1],
        c="cyan",
        s=8,
        alpha=1,
        zorder=5,
        label=f"GFP+ immune + microglia (n={mask_microglia.sum()})",
    )
    ax.legend(fontsize=10, framealpha=0.7)

    plt.savefig("gfp_on_fusion_tlscore_" + args.label + ".png", dpi=300, bbox_inches="tight",facecolor="white")
    plt.clf()

    #zoom in

    
    x_center = mean_x  #6140
    y_center = mean_y #8309
    half_width = 1000
    half_height = 1000
    x_min = x_center - half_width
    y_min = y_center - half_height
    x_max = x_center + half_width
    y_max = y_center + half_height
    crop = (x_min, y_min, x_max, y_max)   # ← fixed

    score_col = "tlscore_gfp_probe_norm"
    base_col  = "tltumscore_C11orf95_relafus1_norm"
    microglia = "tlscore_GFP_Immune_Microglia NN"

    threshold = np.nanpercentile(adata.obs[score_col], THRESHOLD)
    mask = adata.obs[score_col].values > threshold
    mask_microglia = gfp_immune_adata_tlscore_conf.obs["tlcelltype_subclass_assigned"] == microglia

    library_id = list(adata.uns["spatial"].keys())[0]
    scale = adata.uns["spatial"][library_id]["scalefactors"]["tissue_hires_scalef"]
    # Try raw coords first (no scale factor)
    coords = adata.obsm["spatial"]

    library_id_microglia = list(gfp_immune_adata_tlscore_conf.uns["spatial"].keys())[0]
    scale_microglia = gfp_immune_adata_tlscore_conf.uns["spatial"][library_id]["scalefactors"]["tissue_hires_scalef"]
    coords_microglia = gfp_immune_adata_tlscore_conf.obsm["spatial"]


    fig, ax = plt.subplots(figsize=(10, 10))
    sq.pl.spatial_scatter(
        adata,
        img_res_key="hires",
        img_alpha = 0.4,
        color=base_col,
        crop_coord=crop,
        ax=ax,
        title=f"{base_col} + GFP positive overlay",
    )

    # Squidpy maps crop region → 0-1020 pixel space
    # So we need: (coord - crop_min) * scale → image pixels → squidpy axis
    coords_raw = adata.obsm["spatial"]   # raw spatial coords

    # Filter to spots within the crop window first
    in_crop = (
        (coords_raw[:, 0] >= x_min) & (coords_raw[:, 0] <= x_max) &
        (coords_raw[:, 1] >= y_min) & (coords_raw[:, 1] <= y_max)
    )

    combined_mask = mask & in_crop

    # Then shift and scale
    coords_crop = coords_raw.copy()
    coords_crop[:, 0] = (coords_raw[:, 0] - x_min) * scale
    coords_crop[:, 1] = (coords_raw[:, 1] - y_min) * scale

    ax.scatter(
        coords_crop[combined_mask, 0],
        coords_crop[combined_mask, 1],
        c="red",
        s=10,
        alpha=1,
        zorder=5,
        label=f"GFP+ (n={combined_mask.sum()})",
    )

    coords_immunity_raw = gfp_immune_adata_tlscore_conf.obsm["spatial"]

    in_crop_immunity = (
        (coords_immunity_raw[:, 0] >= x_min) & (coords_immunity_raw[:, 0] <= x_max) &
        (coords_immunity_raw[:, 1] >= y_min) & (coords_immunity_raw[:, 1] <= y_max)
    )

    combined_mask_microglia = mask_microglia & in_crop_immunity


    # Then shift and scale
    coords_crop_immunity = coords_immunity_raw.copy()
    coords_crop_immunity[:, 0] = (coords_immunity_raw[:, 0] - x_min) * scale
    coords_crop_immunity[:, 1] = (coords_immunity_raw[:, 1] - y_min) * scale

    ax.scatter(
        coords_crop_immunity[in_crop_immunity, 0],
        coords_crop_immunity[in_crop_immunity, 1],
        c="orange",
        s=10,
        alpha=1,
        zorder=5,
        label=f"GFP+ immune (n={in_crop_immunity.sum()})",
    )

    ax.scatter(
        coords_crop_immunity[combined_mask_microglia, 0],
        coords_crop_immunity[combined_mask_microglia, 1],
        c="cyan",
        s=10,
        alpha=1,
        zorder=5,
        label=f"GFP+ immune + microglia (n={combined_mask_microglia.sum()})",
    )

    ax.legend(fontsize=10)
    plt.savefig("gfp_on_fusion_tlscore_zoom_" + args.label + str(x_center) + "_" + str(y_center) + ".png", dpi=300, bbox_inches="tight",facecolor="white")
    plt.clf()
    

    #Isolate position for gfp spots in fusion scoring distribution

    score_col = "tlscore_gfp_probe_norm"
    base_col  = "tltumscore_C11orf95_relafus1_norm"


    threshold = np.nanpercentile(adata.obs[score_col], THRESHOLD)
    mask = adata.obs[score_col].values >= threshold
    

    #isolate fusion signature scores at location with GFP
    gfp_positive_base_scores = adata.obs[base_col].values[mask]
    print(gfp_positive_base_scores)
    all_base_scores = adata.obs[base_col].values
    print(all_base_scores)

    print(f"Number of GFP+ spots: {mask.sum()}")
    print(f"GFP+ {base_col} scores — mean: {gfp_positive_base_scores.mean():.3f}, "
        f"median: {np.median(gfp_positive_base_scores):.3f}, "
        f"min: {gfp_positive_base_scores.min():.3f}, "
        f"max: {gfp_positive_base_scores.max():.3f}")


    # Plot distribution
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Panel 1: histogram overlay
    axes[0].hist(all_base_scores, bins=50, color="steelblue", alpha=0.6, label="All fusion spots", density=True)
    axes[0].hist(gfp_positive_base_scores, bins=30, color="red", alpha=0.7, label="GFP+ spots", density=True)
    axes[0].axvline(np.median(gfp_positive_base_scores), color="red", linestyle="--", linewidth=1.5, label="GFP+ median")
    axes[0].axvline(np.median(all_base_scores), color="steelblue", linestyle="--", linewidth=1.5, label="All fusion median")
    axes[0].set_xlabel(base_col)
    axes[0].set_ylabel("Density")
    axes[0].set_title(f"{base_col} distribution\nGFP+ vs all spots")
    axes[0].legend(fontsize=9)

    # Panel 2: ECDF — shows exactly where GFP+ spots sit in the overall distribution
    def ecdf(data):
        x = np.sort(data)
        y = np.arange(1, len(x) + 1) / len(x)
        return x, y

    x_all, y_all = ecdf(all_base_scores)
    x_gfp, y_gfp = ecdf(gfp_positive_base_scores)

    axes[1].plot(x_all, y_all, color="steelblue", label="All fusion spots")
    axes[1].plot(x_gfp, y_gfp, color="red", label="GFP+ spots")
    axes[1].set_xlabel(base_col)
    axes[1].set_ylabel("Cumulative fraction")
    axes[1].set_title(f"ECDF — where GFP+ spots sit\nin {base_col} distribution")
    axes[1].legend(fontsize=9)

    plt.tight_layout()
    plt.savefig("gfp_in_fusion_distribution_" + args.label + ".png", dpi=300, bbox_inches="tight",facecolor="white")
    plt.clf()

    # Percentile rank of GFP+ spots within overall distribution
    percentile_rank = (all_base_scores < np.median(gfp_positive_base_scores)).mean() * 100
    print(f"GFP+ median {base_col} score is at the {percentile_rank:.1f}th percentile of all spots")


    #statistical significance gfp subset vs fusion
    
    # Get scores for non-GFP+ spots (complement)
    gfp_negative_base_scores = adata.obs[base_col].values[~mask]

    # Mann-Whitney U test (one-sided: testing if GFP+ scores are higher)
    stat, pval = mannwhitneyu(
        gfp_positive_base_scores,
        gfp_negative_base_scores,
        alternative="two-sided"
    )

    # Effect size: rank-biserial correlation
    n1, n2 = len(gfp_positive_base_scores), len(gfp_negative_base_scores)
    rank_biserial = 1 - (2 * stat) / (n1 * n2)

    print(f"\n── Mann-Whitney U Test ──────────────────────────")
    print(f"GFP+ spots      : n={n1}, median={np.median(gfp_positive_base_scores):.4f}")
    print(f"GFP- spots      : n={n2}, median={np.median(gfp_negative_base_scores):.4f}")
    print(f"U statistic     : {stat:.2f}")
    print(f"p-value         : {pval:.4e}")
    print(f"Effect size (r) : {rank_biserial:.4f}  (-1=GFP- higher, +1=GFP+ higher)")
    print(f"────────────────────────────────────────────────")

    # NNLS deconvolve assignment

    #| Aspect         |	NNLS deconvolution                      | tl.score_genes                              |
    #|----------------|-----------------------------------------|---------------------------------------------|
    #| Goal           | Estimate cell-type composition	        | Measure signature activity per cell         |
    #| Input          | Reference profiles + mixture expression |	Gene list + background genes          |
    #| Output         | Fractions or abundance-like weights     | Relative score per cell                     |
    #| Interpretation | “How much of each cell type is present?"| “How active is this gene set in this cell?” |
    #| Best use       | Bulk / spot deconvolution               | Single-cell or spatial signature scoring    |


    complete_markers = marker_allen.copy()
    complete_markers['C11orf95_relafus1'] = tumor['C11orf95_relafus1']

    prop_df = nnls_deconvolve(adata, complete_markers) #nnls_deconvolve(adata, marker_allen)
    adata.obs = adata.obs.join(prop_df)

    gfp_adata = adata[adata.obs['nnls_gfp_probe'] > 0]
    gfp_scores = gfp_adata.obs['nnls_gfp_probe'].values

    select_scores = ["nnls_" + k for k in marker_allen.keys()]

    scores_df = gfp_adata.obs[select_scores]
    scores_df = scores_df.drop(columns=["nnls_gfp_probe"])
    scores_df["predicted_celltype"] = scores_df.idxmax(axis=1)
    counts = scores_df["predicted_celltype"].value_counts()
    gfp_adata.obs["nnls_predicted_celltype"] = scores_df["predicted_celltype"]
    total = counts.sum()
    pcts = counts / total * 100

    fig, axes = plt.subplots(1, 2, figsize=(18, 6))

    counts.plot(
        ax=axes[0],
        kind="bar",
        color="steelblue"
        )

    axes[0].set_xlabel("predicted cell type highest score")
    axes[0].set_ylabel("Count")
    axes[0].set_title("Most frequent max-scoring cell type at GFP probe spot")

    cumulative = 0
    colors = plt.cm.tab20.colors

    for i, (celltype, pct) in enumerate(pcts.items()):
        axes[1].bar(0, pct, bottom=cumulative, color=colors[i % len(colors)], 
           width=0.5, label=f"{celltype} ({pct:.1f}%)")
        # Only label segments large enough to read
        if pct > 2:
            axes[1].text(
                0, cumulative + pct / 2,
                f"{celltype}\n{pct:.1f}%",
                ha="center", va="center",
                fontsize=7, color="white", fontweight="bold",
            )
        cumulative += pct

    axes[1].set_xlim(-0.5, 0.5)
    axes[1].set_ylim(0, 100)
    axes[1].set_xticks([])
    axes[1].set_ylabel("Percentage of GFP+ spots (%)")
    axes[1].set_title("GFP+ spots assigned to\nmax-scoring cell type")
    axes[1].legend(
        loc="upper right",
        bbox_to_anchor=(2.0, 1),
        fontsize=8,
        title="Cell type",
    )



    plt.tight_layout()
    plt.savefig("spatial_gfp_spots_top_score_tlscore_" + args.label + "_bar_nnls.png",dpi=300,bbox_inches="tight",facecolor="white")
    plt.clf()


    #plot NNLS gfp probes
    adata_filtered = adata[adata.obs["nnls_gfp_probe"] > THRESHOLD/100].copy()
    
    sq.pl.spatial_scatter(
        adata_filtered,
        img_res_key="hires",
        img_alpha = 0.4,
        color='nnls_gfp_probe',
        size=5.0,
        alpha=1.0,
        dpi=300,
        title='nnls gfp probe'
    )

    plt.savefig("celltype_gfp_nnls_" + args.label + ".png",dpi=300,bbox_inches="tight",facecolor="white")

    

    #plot NNLS tumor
    adata_filtered = adata[adata.obs["nnls_C11orf95_relafus1"] > THRESHOLD/100].copy()

    sq.pl.spatial_scatter(
        adata_filtered,
        img_res_key="hires",
        img_alpha = 0.4,
        color='nnls_gfp_probe',
        size=5.0,
        alpha=1.0,
        dpi=300,
        title='nnls C11orf95_relafus1'
    )

    plt.savefig("celltype_C11orf95_relafus1_nnls_" + args.label + ".png",dpi=300,bbox_inches="tight",facecolor="white")



    #subclassify immmune from nnls
    '''
    gfp_immune_adata = gfp_adata[gfp_adata.obs["nnls_predicted_celltype"] == "nnls_Immune"]
    #marker_allen_subclass["Immune"]
    for ct, genes in marker_allen_subclass["Immune"].items():
        #n = score_geneset(adata, genes, ct)
        sc.tl.score_genes(
            gfp_immune_adata,
            gene_list=genes,
            score_name= f"nnls_Immune_{ct}",
            use_raw=False
        )

    score_cols = [c for c in gfp_immune_adata.obs.columns if c.startswith("nnls_Immune_")]
    scores_df = gfp_immune_adata.obs[score_cols].copy()

    gfp_immune_adata.obs["nnls_Immune_subclass_assigned"] = scores_df.idxmax(axis=1)
    gfp_immune_adata.obs["nnls_Immune_subclass_score"] = scores_df.max(axis=1)
    counts = gfp_immune_adata.obs["nnls_Immune_subclass_assigned"].value_counts()

    ax = counts.plot(
        kind="bar",
        figsize=(12, 6))

    ax.set_xlabel("subclass type highest score")
    ax.set_ylabel("counts")
    ax.set_title("Most frequent max-scoring subclass at GFP probe spot")
    plt.tight_layout()

    plt.savefig("spatial_gfp_spots_top_score_tlscore_" + args.label + "_bar_nnls_immune_subclass.png",dpi=300,bbox_inches="tight")
    '''



if __name__ == '__main__':
    main()
