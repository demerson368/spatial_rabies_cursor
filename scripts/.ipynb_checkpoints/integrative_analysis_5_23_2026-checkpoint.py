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
                          spot_size=1.4, alpha=0.7, min_threshold=0.3,
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
    ax.imshow(img)
    ax.axis("off")
    
    legend_handles = []
    
    for i, ct in enumerate(celltypes):
        if ct not in adata.obs.columns:
            print(f"Warning: {ct} not found in adata.obs, skipping")
            continue
        
        scores = adata.obs[ct].values
        
        # Only plot spots above threshold percentile
        threshold = np.percentile(scores, min_threshold * 100)
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
    threshold = np.percentile(scores, min_percentile)
    mask = scores > threshold

    # Normalize scores for colormap
    vmin, vmax = scores[mask].min(), np.percentile(scores[mask], min_percentile)
    norm = plt.Normalize(vmin=vmin, vmax=vmax)
    colors = plt.get_cmap(cmap)(norm(scores[mask]))

    fig, ax = plt.subplots(figsize=figsize)
    ax.imshow(img)

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
        ax.imshow(img)
        ax.set_title(f"{ct}\n(no expression)", fontsize=8)
        ax.axis("off")
        return


    # Only show spots above threshold
    threshold = np.percentile(scores, min_percentile)
    mask = scores > threshold

    # Normalize scores for colormap
    vmin, vmax = scores[mask].min(), np.percentile(scores[mask], min_percentile)
    norm = plt.Normalize(vmin=vmin, vmax=vmax)
    colors = plt.get_cmap(cmap)(norm(scores[mask]))

    #fig, ax = plt.subplots(figsize=figsize)
    ax.imshow(img)

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
            thresh = np.percentile(s, score_pct)
            votes[:, j] += (s >= thresh).astype(int)

        # Method B — weighted score
        if f"wscore_{ct}" in adata.obs:
            s = adata.obs[f"wscore_{ct}"].values
            thresh = np.percentile(s, score_pct)
            votes[:, j] += (s >= thresh).astype(int)

        # Method C — NNLS proportion
        if f"nnls_{ct}" in adata.obs:
            s = adata.obs[f"nnls_{ct}"].values
            thresh = np.percentile(s, score_pct)
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
    parser.add_argument('min_pct',type=float,help="score must exceed minimum pct threshold to be visible")
    args = parser.parse_args()

    
    adata = sc.read_h5ad(args.adata_name)

    marker_table_allen,marker_allen = marker_unique_excel('/research/groups/mackgrp/home/demerson/AIB/yao_et_al_Nature_2023.xlsx','supertype_annotation')
    #remove ct with < 3 genes in adata
    marker_allen = minimum_threshold_gene_filter(adata,marker_allen)

    marker_allen['gfp_probe'] = ['Gfp', 'GFP', 'EGFP', 'eGFP', 'gfp', 'egfp','EGFP-L10a', 'Gfp-L10a', 'td-Tomato', 'mCherry']

    for ct, genes in marker_allen.items():
        #n = score_geneset(adata, genes, ct)
        sc.tl.score_genes(
            adata,
            gene_list        = genes,
            score_name       = f"tlscore_{ct}",
            use_raw          = True,
            n_bins           = 25,
            ctrl_size        = min(len(genes) * 5, 200),
        )
        #adata.obs[f"wscore_{ct}"] = weighted_score(adata, genes)
        #prop_df = nnls_deconvolve(adata, marker_allen)
        #adata.obs = adata.obs.join(prop_df)
        


    celltype_scores = list(marker_allen.keys())

    #Inspecting all GFP probed spots (score > 0) and analyzing which class has the highest score excluding GFP_probe

    gfp_scores = adata.obs['gfp_probe'].values
    gfp_adata = adata[gfp_scores > 0]

    scores_df = gfp_adata.obs[marker_allen.keys()]
    scores_df = scores_df.drop(columns=["gfp_probe"])
    scores_df["predicted_celltype"] = scores_df.idxmax(axis=1)
    counts = scores_df["predicted_celltype"].value_counts()


    ax = counts.plot(
        kind="bar",
        figsize=(12, 6)
    )

    ax.set_xlabel("cell type highest score")
    ax.set_ylabel("Frequency")
    ax.set_title("Most frequent max-scoring cell type at GFP probe spot")
    plt.tight_layout()
    plt.savefig("spatial_gfp_spots_top_score_" + args.label + "_bar.png",dpi=300,bbox_inches="tight")


    # Spatial feature plots — one panel per cell type 

    ncols = 4
    nrows = -(-len(celltype_scores) // ncols)   # ceiling division

    fig, axes = plt.subplots(nrows, ncols,figsize=(ncols * 4, nrows * 4))
    axes = axes.flatten()
    plot_type = 2

    for i, ct in enumerate(celltype_scores):
        if plot_type == 1:
            sc.pl.spatial(
                adata,
                color=ct,
                ax=axes[i],
                show=False,
                title=ct,
                cmap="magma",
                size=1.4,
                vmin=0,
                vmax=np.percentile(adata.obs[ct], args.min_pct),  # %5 instead clip top 1% for viz
                colorbar_loc=None,
            )
        else:
            fig = plot_spatial_with_border_v2(
                adata, ct,axes[i],fig,
                border_color="lime",      # neon green
                spot_size=10,
                min_percentile=args.min_pct)

        axes[i].set_xlabel("")
        axes[i].set_ylabel("")

    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)

    plt.suptitle("Cell type module scores allen — " + args.label, fontsize=14, y=1.01)
    plt.tight_layout()
    plt.savefig("spatial_celltype_scores_" + args.label  + "_" + str(args.min_pct) + "_v2.png", dpi=150, bbox_inches="tight")

    ## Enhance Brightness / Contrast of the Tissue Image
    img_key = "hires"  # or "lowres"
    library_id = list(adata.uns["spatial"].keys())[0]

    # Get the image
    img = adata.uns["spatial"][library_id]["images"][img_key].copy()

    img_bright = exposure.equalize_adapthist(img, clip_limit=0.02) #0.03
    adata.uns["spatial"][library_id]["images"][img_key] = img_bright
    #ct = "gfp_probe"
    #sc.pl.spatial(adata, color=ct, cmap="magma", img_key=img_key)
  
    celltypes_to_overlay = [
        "IT-ET Glut", "CTX-MGE GABA", "CTX-CGE GABA",
        "Astro-Epen", "OPC-Oligo", "Immune","gfp_probe"
    ]

    fig = plot_celltype_overlay(
        adata,
        celltypes=celltypes_to_overlay,
        img_key="hires",
        spot_size=4,          # increase if spots are hard to see originally 8
        alpha=0.5,
        min_threshold=args.min_pct/100,    # only show top 50% expressing spots per cell type
    )

    fig.savefig("celltype_overlay_" + args.label + "_" + str(args.min_pct) + "_v2.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    
    ## Side-by-side: Individual + Overlay Together

    fig, axes = plt.subplots(1, len(celltypes_to_overlay) + 1,
                          figsize=((len(celltypes_to_overlay) + 1) * 4, 4))

    # Individual plots
    for i, ct in enumerate(celltypes_to_overlay):
        sc.pl.spatial(adata, color=ct, ax=axes[i], show=False,
                  title=ct, cmap="magma", colorbar_loc=None,
                  vmax=np.percentile(adata.obs[ct], args.min_pct))

    # Overlay in last panel — draw manually into the existing axis
    library_id = list(adata.uns["spatial"].keys())[0]
    img = adata.uns["spatial"][library_id]["images"]["hires"]
    scale = adata.uns["spatial"][library_id]["scalefactors"]["tissue_hires_scalef"]
    coords = adata.obsm["spatial"] * scale

    axes[-1].imshow(img)
    for i, ct in enumerate(celltypes_to_overlay):
        scores = adata.obs[ct].values
        mask = scores > np.percentile(scores, args.min_pct)
        axes[-1].scatter(coords[mask, 0], coords[mask, 1],
                     c=[cm.tab20.colors[i % 20]], s=2,
                     alpha=0.5, label=ct, linewidths=0)
    axes[-1].legend(fontsize=6, loc="upper right", bbox_to_anchor=(1.4, 1))
    axes[-1].set_title("Overlay")
    axes[-1].axis("off")

    plt.tight_layout()
    fig.savefig("celltype_ind_plus_overlay_" + args.label + "_" + str(args.min_pct) + "_v2.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    # Dot plot — top marker genes per cell type across tissue domains 
    #Optional: bin spots by a spatial cluster or annotation if you have one 
    #e.g. adata.obs['spatial_domain'] from Leiden clustering on spatial graph

    if "leiden" not in adata.obs.columns:
        if "X_pca_harmony" in adata.obsm:
            use_rep = "X_pca_harmony"
        elif "X_scvi" in adata.obsm:
            use_rep = "X_scvi"
        elif "X_umap" in adata.obsm:
            use_rep = "X_umap"
        else:
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
            sc.pp.highly_variable_genes(adata, n_top_genes=3000)
            sc.pp.pca(adata)
            use_rep = "X_pca"
            sc.pp.neighbors(adata, use_rep="X_pca")   # or use_rep="X_pca_harmony"
            sc.tl.leiden(adata, resolution=0.5, key_added="leiden")

    # Pick top 3 marker genes per cell type that are actually in the data
    top_genes_per_ct = {}
    for ct, genes in marker_allen.items():
        present = [g for g in genes if g in adata.var_names][:3]
        if present:
            top_genes_per_ct[ct] = present

    all_top_genes = [g for genes in top_genes_per_ct.values() for g in genes]
    # deduplicate while preserving order
    seen = set()
    plot_genes = [g for g in all_top_genes if not (g in seen or seen.add(g))]

    sc.pl.dotplot(
        adata,
        var_names=top_genes_per_ct,   # dict groups genes by cell type in plot
        groupby="leiden",
        use_raw=False,
        cmap="Blues",
        dot_max=0.8,
        standard_scale="var",          # scale each gene 0–1 across clusters
        figsize=(min(len(plot_genes) * 0.45 + 2, 28), 5),
        title="Top marker genes per cell type × spatial cluster",
        save="spatial_celltype_scores_" + args.label + ".png",
    )

    # Heatmap of module scores across spatial clusters
    
    score_df = (
        adata.obs[celltype_scores + ["leiden"]]
        .groupby("leiden")[celltype_scores]
        .mean()
    )

    # Z-score across cell types for each cluster
    score_z = score_df.apply(lambda x: (x - x.mean()) / (x.std() + 1e-8), axis=1)

    fig, ax = plt.subplots(figsize=(max(len(celltype_scores) * 0.55, 12), 6))
    sns.heatmap(
        score_z.T,
        cmap="RdBu_r",
        center=0,
        linewidths=0.4,
        linecolor="white",
        ax=ax,
        cbar_kws={"label": "Z-score"},
        yticklabels=True,
        xticklabels=True,
    )
    ax.set_xlabel("Leiden spatial cluster")
    ax.set_ylabel("Cell type")
    ax.set_title("Cell type module score Z-scores per spatial cluster")
    plt.xticks(rotation=45, ha="right", fontsize=9)
    plt.yticks(fontsize=8)
    plt.tight_layout()
    plt.savefig("celltype_score_heatmap_" + args.label + ".png", dpi=150, bbox_inches="tight")

    adata.write_h5ad(args.adata_name)


if __name__ == '__main__':
    main()
