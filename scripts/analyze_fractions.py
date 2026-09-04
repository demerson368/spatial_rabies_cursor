import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq

def main():

    combined_label = "Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um_merged_nontumor_v1_broad_nontumor"
    path = f"../results_cell2location/c2l_spatial_{combined_label}/sp.h5ad"

    adata = sc.read_h5ad(path)

    factor_names = list(adata.uns["mod"]["factor_names"])
    m = adata.obsm["means_cell_abundance_w_sf"]
    abund = m.copy()

    #if isinstance(m, pd.DataFrame):
    #    abund = m.copy()
    #if not abund.index.equals(adata.obs_names):
    #    abund.index = adata.obs_names
    #else:
    #    abund = pd.DataFrame(np.asarray(m), index=adata.obs_names, columns=factor_names)

    
    # Rank 1 and rank 2 per spot
    rank1 = abund.idxmax(axis=1).rename("dominant_celltype")
    rank1_abund = abund.max(axis=1).rename("dominant_abundance")

    # Second-largest: mask top value, then argmax again
    abund_masked = abund.copy()
    for ct in abund.columns:
        is_top = rank1.eq(ct)
        abund_masked.loc[is_top, ct] = -np.inf

    rank2 = abund_masked.idxmax(axis=1).rename("secondary_celltype")
    rank2_abund = abund_masked.max(axis=1).rename("secondary_abundance")
    summary = pd.concat([rank1, rank1_abund, rank2, rank2_abund], axis=1)

    # All GABA-dominant locations
    gaba_spots = summary[summary["dominant_celltype"] == "meanscell_abundance_w_sf_GABA"].copy()
    gaba_spots["margin"] = gaba_spots["dominant_abundance"] - gaba_spots["secondary_abundance"]
    #gaba_spots.head(10)

    print(gaba_spots["secondary_celltype"].value_counts())

    print("mean primary gaba: ",np.nanmean(gaba_spots["dominant_abundance"]))
    print("mean secondary for gaba: ",np.nanmean(gaba_spots["secondary_abundance"]))
    
    

    #Spatial view (optional)
    #adata.obs["secondary_celltype"] = summary["secondary_celltype"]
    #adata.obs["c2l_margin"] = summary["dominant_abundance"] - summary["secondary_abundance"]
    #gaba_mask = adata.obs["dominant_celltype"] == "GABA"  # or use adata.obs if already present
    #adata_gaba = adata[gaba_mask].copy()
    #sq.pl.spatial_scatter(adata_gaba, color="secondary_celltype", size=1.5)


if __name__ == '__main__':
    main()
