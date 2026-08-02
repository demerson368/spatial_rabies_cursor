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
import squidpy as sq
import math
from collections import defaultdict
from PIL import Image
from skimage import transform
import cv2
from sklearn.preprocessing import MinMaxScaler
import cell2location
from cell2location.utils.filtering import filter_genes
from cell2location.models import RegressionModel
from cell2location.models import Cell2location

def IF_registrar_coord(adata,scale_factor = 0.5,select_channel = "red"):
    # Get the H&E reference image

    he_image = adata.uns["spatial"][library_id]["images"]["hires"]

    # Option A: Manual alignment with known scale factor

    if_resized = transform.resize(
        if_image,
        (he_image.shape[0], he_image.shape[1])
    )


    if select_channel == "red":
        if_resized = if_resized[:, :, 0]

    # Extract IF Intensity per Spot
    coords = adata.obsm["spatial"]   # (n_spots, 2) array of x,y coordinates
    scale  = adata.uns["spatial"][library_id]["scalefactors"]["tissue_hires_scalef"]

    # Extract mean IF intensity in a window around each spot
    spot_radius = 8   # pixels, adjust to spot size

    if_intensities = []

    for x, y in coords:
        # Scale coordinates to image space
        xi = int(x * scale)
        yi = int(y * scale)

        # Extract patch around spot
        patch = if_resized[
            max(0, yi - spot_radius):yi + spot_radius,
            max(0, xi - spot_radius):xi + spot_radius
        ]

        if_intensities.append(patch.mean())

    # Add to adata.obs
    adata.obs["IF_intensity"] = if_intensities
    #print(adata.obs["IF_intensity"].describe())

    return if_intensities


def marker_unique_excel(exl_file,sheet_name,cluster = "supertype.markers.combo",class_label = "class_label"):

    """
        split marker unique genes into class. Based on Allen Mouse Brain.

        /research/groups/mackgrp/home/demerson/CellTypeEstimate/db/mm.brain.v2.xlsx

        class_label = Class for mm.brain.v2.xlsx
        subclass_label = subClass for mm.brain.v2.xlsx
        sheet = Sheet1 for  mm.brain.v2.xlsx
        cluster = geneSymbolmore1 for  mm.brain.v2.xlsx


    """
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
        if f"tlscore_{ct}" in adata.obs:
            s = adata.obs[f"tlscore_{ct}"].values
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


def plot_celltype_overlay(adata, celltypes, img_key="hires", 
                          spot_size=1.4, alpha=0.7, min_threshold=0.3,
                          figsize=(10, 10),border_color="white",colors = cm.tab20.colors,crop_coord=None):
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
    
    # ── Resolve zoom from crop_coord if provided ──────────────────────────────
    if crop_coord is not None:
        x_min, y_min, x_max, y_max = crop_coord
        xlim = (x_min, x_max)
        ylim = (y_max, y_min)   # inverted: imshow origin is top-left

        # Crop the image array to avoid rendering the full tissue
        img = img[y_min:y_max, x_min:x_max]

        # Shift coordinates so they align with the cropped image
        coords = coords.copy()
        coords[:, 0] -= x_min
        coords[:, 1] -= y_min

        # Recompute limits in cropped space
        xlim = (0, x_max - x_min)
        ylim = (y_max - y_min, 0)   # still inverted

    
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
        threshold = np.percentile(scores, min_threshold * 100)
        mask = scores > threshold
        
        # Normalize scores for alpha/size scaling
        s_norm = (scores[mask] - scores[mask].min()) / (scores[mask].max() - scores[mask].min() + 1e-8)
        
        color = colors[i % len(colors)]
        

        # ── Draw border layer first (slightly larger, neon color) ─────────────────
        scatter = ax.scatter(
            coords[mask, 0],
            coords[mask, 1],
            s=spot_size * (1 + s_norm * 3), #spot_size * 1.6, #spot_size * (1 + s_norm * 3),        # slightly bigger than fill
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

def main():

    """
        create or import celltype on spots either direct calling or transfer/modeling:

        tlscore - score ppn spot
        tlscore_norm - score on spot. Normalized
        NNLS - deconvolution on spot.  set to max fraction per spot
        cellestimate - import cellestimate calls per spot
        predicted.id - import celltype transfer from scRNA to spot
        cell2location - import cell2location calls

        add features to adata and export loupe ../for_loupe/*.csv 


    """


    parser = argparse.ArgumentParser()
    parser.add_argument('--h5ad', type=str, help="h5ad file spatial output")
    parser.add_argument('--label', type=str, help="spatial sample label")
    parser.add_argument('--ref_tag', type=str, default="v2", help="ref tag for tagged outputs")
    parser.add_argument('--cohort', type=str, default="Tuyu_Apr2026", help="processed data cohort subdir")
    parser.add_argument('--cellestimate_cf', type=str, default="Yao_et_al_mouse_format", help="cellestimate filename suffix")
    parser.add_argument('--h5ad_ext', type=str, default=None, help="extended h5ad output path")
    parser.add_argument('--loupe_out', type=str, default=None, help="Loupe CSV output path")
    parser.add_argument('--c2l_h5ad', type=str, default=None, help="cell2location sp.h5ad input path")
    parser.add_argument('--transfer_txt', type=str, default=None, help="label transfer txt input path")
    parser.add_argument('--cellestimate_txt', type=str, default=None, help="cellestimate txt input path")
    parser.add_argument('--marker_excel', type=str, required=True, help="Allen/mouse brain marker Excel path")
    parser.add_argument('--marker_sheet', type=str, default="supertype_annotation", help="Excel sheet name")
    parser.add_argument('--marker_cluster_col', type=str, default="supertype.markers.combo", help="gene list column")
    parser.add_argument('--marker_class_col', type=str, default="class_label", help="cell class column")
    args = parser.parse_args()

    adata = sc.read_h5ad(args.h5ad)
    label = args.label
    ref_tag = args.ref_tag
    cohort = args.cohort
    processed_dir = f"../data/processed/{cohort}"
    loupe_out = args.loupe_out or f"../for_loupe/spatial_{label}{'' if ref_tag == 'v2' else f'_{ref_tag}'}_for_loupe.csv"
    h5ad_ext = args.h5ad_ext or (
        f"{processed_dir}/spatial_{label}{'' if ref_tag == 'v2' else f'_{ref_tag}'}_ext.h5ad"
    )
    c2l_h5ad = args.c2l_h5ad or f"../results_cell2location/c2l_spatial_{label}_{ref_tag}/sp.h5ad"
    transfer_txt = args.transfer_txt or (
        f"{processed_dir}/spatial_{label}{'' if ref_tag == 'v2' else f'_{ref_tag}'}_predicted.id_transfer.txt"
    )
    cellestimate_txt = args.cellestimate_txt or (
        f"{processed_dir}/spatial_{label}{args.cellestimate_cf}_cellType.txt"
    )

    print(
        "Loading marker atlas:",
        f"file={args.marker_excel}",
        f"sheet={args.marker_sheet}",
        f"cluster={args.marker_cluster_col}",
        f"class={args.marker_class_col}",
        sep=" ",
    )
    marker_table_allen, marker_allen = marker_unique_excel(
        args.marker_excel,
        args.marker_sheet,
        cluster=args.marker_cluster_col,
        class_label=args.marker_class_col,
    )

    #remove ct with < 3 genes in adata
    marker_allen = minimum_threshold_gene_filter(adata,marker_allen)

    marker_allen['gfp_probe'] = ['Gfp', 'GFP', 'EGFP', 'eGFP', 'gfp', 'egfp','EGFP-L10a', 'Gfp-L10a', 'td-Tomato', 'mCherry']
    tumor = {}
    instance = pd.read_excel("../reference/93_core_gene_set_commonly_upregulated_in_human_and_mouse_C11orf95_RELAFUS1_driven_ependymoma.xlsx",sheet_name="Sheet1",header=None)
    tumor['C11orf95_relafus1'] = instance[0].str.capitalize().tolist()


    for ct, genes in marker_allen.items():
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
   
    # create normalized scores
    score_cols = []
    for ct, genes in marker_allen.items():
        score_cols.append('tlscore_' + ct)
    score_cols.append('tltumscore_C11orf95_relafus1')

    norm_threshold = 0.07
    score_cols_norm = []
    for col in score_cols:
        if col not in adata.obs.columns:
            print(f"Skipping '{col}': not found in adata.obs")
            continue#

        norm_col  = f"norm_{col}"#

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

    # Assign prediction based on score max (unnormalized)
    score_cols = [c for c in adata.obs.columns if c.startswith("tlscore_")]
    scores = adata.obs[score_cols].copy()

    top2 = np.sort(scores.values, axis=1)[:, -2:]
    margin = top2[:, 1] - top2[:, 0]   # best minus second-best
    adata.obs["tlcelltype_assigned"] = scores.idxmax(axis=1)
    adata.obs["tlcelltype_score"] = scores.max(axis=1)
    adata.obs["tlcelltype_margin"] = margin

    # Require both absolute score AND separation from runner-up
    low_conf = (adata.obs["tlcelltype_score"] < 0.1) | (adata.obs["tlcelltype_margin"] < 0.05)
    adata.obs.loc[low_conf, "tlcelltype_assigned"] = "Unknown"

    
    scaler = MinMaxScaler()

    score_cols = []
    for ct, genes in marker_allen.items():
        if ct != 'gfp_probe':
            score_cols.append('tlscore_' + ct)
    scores_df = adata.obs[score_cols].copy()

    scores_norm = pd.DataFrame(
        scaler.fit_transform(scores_df),
        index   = scores_df.index,
        columns = score_cols #[c.replace("score_", "") for c in score_cols]
    )

    # Assign based on normalized scores
    min_threshold = 0.3

    max_norm_scores = scores_norm.max(axis=1)

    adata.obs["celltype_assigned_norm"] = np.where(
        max_norm_scores >= min_threshold,
        scores_norm.idxmax(axis=1),
        "Unassigned"
    )

    loupe_export = adata.obs[[
         "tlcelltype_assigned",
         "celltype_assigned_norm"
    ]].reset_index()
    loupe_export.columns = ["Barcode", "tlcelltype_assigned", "celltype_assigned_norm"]

    filename = loupe_out
    loupe_export.to_csv(filename, index=False)
    print(f"Exported {len(loupe_export)} spots")


    #NNLS
    prop_df = nnls_deconvolve(adata, marker_allen)

    new_cols = prop_df.columns.difference(adata.obs.columns)

    if len(new_cols) > 0:
        adata.obs = adata.obs.join(prop_df[new_cols])
    else:
        print("All columns already exist in adata.obs, skipping join.")

    #total cells
    select_scores = ["nnls_" + k for k in marker_allen.keys()]

    scores_df = adata.obs[select_scores]
    scores_df["predicted_celltype"] = scores_df.idxmax(axis=1)
    adata.obs["nnls_celltype_assigned"] = scores_df["predicted_celltype"]

    loupe_export = adata.obs[[
        "tlcelltype_assigned",
        "celltype_assigned_norm",
        "nnls_celltype_assigned",
    ]].reset_index()
    loupe_export.columns = ["Barcode", "tlcelltype_assigned", "celltype_assigned_norm",
                        "nnls_celltype_assigned"] #"norm_tltumscore_C11orf95_relafus1","norm_tlscore_gfp_probe"]

    col = "norm_tltumscore_C11orf95_relafus1"

    loupe_export[col] = np.where(
        adata.obs[col].values > 0.07,
        col,
        "not fusion"
    )

    col2 = "norm_tlscore_gfp_probe"
    loupe_export[col2] = np.where(
        adata.obs[col2].values > 0.07,   # adjust threshold as needed
        col2,
        "not gfp"
    )


    loupe_export.to_csv(loupe_out, index=False)

    # Add celltypes from spot identification from celltypeestimate
    celltype_est = pd.read_csv(cellestimate_txt, sep="\t", engine="python")

    celltype_est = celltype_est.rename(columns={
        "Unnamed: 0": "Barcode",
        "cell_type": "cellestimate"
    })

    loupe_export = loupe_export.merge(
        celltype_est,
        on="Barcode",
        how="left",
        suffixes=("", "_right")
    )

    ## test adding labels from label transfer from ref (merged_dJ33_J51_D87_J58_J59_J61_D85_J52_J53)

    predicted_id_df = pd.read_csv(transfer_txt, sep="\t", engine="python")

    predicted_id_df = predicted_id_df.rename(columns={
        "Unnamed: 0": "Barcode",
        "predicted.id": "predicted.id"
    })

    loupe_export = loupe_export.merge(
        predicted_id_df,
        on="Barcode",
        how="left",
        suffixes=("", "_right")
    )


    ## Add cell2location
    adata_cell_2_location = sc.read_h5ad(c2l_h5ad)

    cell2loc_series = adata_cell_2_location.obs["dominant_celltype"].rename("cell2location")

    # Map onto loupe_export by barcode — unmatched barcodes get NaN
    loupe_export["cell2location"] = loupe_export["Barcode"].map(cell2loc_series)

    adata.obs["cell2location"] = loupe_export["cell2location"]
    adata.obs["predicted.id"] = loupe_export["predicted.id"].astype(str)
    adata.obs["predicted.id"] = adata.obs["predicted.id"].fillna("Unknown").astype(str)
    adata.obs["cellestimate"] = loupe_export["cellestimate"]


    # Fix all object dtype columns before writing
    for col in adata.obs.columns:
        if adata.obs[col].dtype == object:
            adata.obs[col] = adata.obs[col].fillna("Unknown").astype(str)

    # Also fix categorical columns with mixed types
    for col in adata.obs.columns:
        if hasattr(adata.obs[col], 'cat'):
            adata.obs[col] = adata.obs[col].astype(str)


    #rewrite adata with 6 added obs fields:
    #
    #cell2location
    #predicted.id
    #cellestimate
    #nnls_celltype_assigned
    #celltype_assigned_norm
    #tlcelltype_assigned

    adata.write(h5ad_ext)

    loupe_export.to_csv(loupe_out, index=False, na_rep="Unknown")


if __name__ == '__main__':
    main()

