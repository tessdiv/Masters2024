r# Script for going from raw data -> first round of all cells leiden clustering

from glob import glob
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

import anndata
import scanpy as sc

adata = anndata.read_h5ad("/home/tessdiv/limitedAccessBladder/NOBACKUP/IMC_BCG/Trine_final/Trine_no_miss_masks/results/BCG_all.h5ad")

sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True,var_type='protein')
adata.layers['raw'] = adata.X.copy()
adata = adata[(adata.obs['area'] > 1)]
#NEW FILTER STEP ADDED FOR REMOVING ROIs WITH VERY FEW CELLS 
# Count cells for each ROI
cells = adata.obs.groupby(["roi"]).size().rename("cell_count")
# Join counts back to adata.obs
adata.obs = adata.obs.join(cells, on='roi')
adata = adata[adata.obs['cell_count'] >= 400]

adata.var['markers'] = [marker.split('(')[0] if (len(marker.split('(')) == 2) else marker  for marker in adata.var.index.to_list()]
adata.var['channel'] = adata.var.index.copy()

adata.obs['total_counts'] = adata.X.sum(1)
adata.var['total_signal'] = adata.X.sum(0)


# Load the CSV file with mob_id and sample mapping
csv_file_path = '/home/tessdiv/limitedAccessBladder/NOBACKUP/IMC_BCG/clin_info.csv'
mob_id_mapping = pd.read_csv(csv_file_path,sep=';')

# Merge andata with mob_id_mapping based on the 'sample' column
adata.obs['mob_id'] = adata.obs['roi'].map(mob_id_mapping.set_index('roi')['mob_id'])

# Now, andata_with_mob_id contains both the 'mob_id' and 'sample' columns
adata.obs['pre_post'] = adata.obs['mob_id'].str.split('_').str.get(1)
adata.obs['mob'] = adata.obs['mob_id'].str.split('_').str.get(0)
adata.obs['mob']

adata[adata.obs.mob_id.isna()]

# Generate unique cell IDs
num_cells = adata.n_obs  # Number of cells in the AnnData object
cell_ids = ["cell_" + str(i) for i in range(num_cells)]

# Add the cell IDs as a new column to the obs DataFrame
adata.obs['cell_id'] = pd.Series(cell_ids, index=adata.obs.index)


#log normalize 
sc.pp.log1p(adata)
adata.layers['log'] = adata.X.copy()

#z scale the data
sc.pp.scale(adata)
adata.layers['z_scale'] = adata.X.copy()


sc.pp.highly_variable_genes(adata)

print("highly variable genes")

use_highly_variable=['CD14(Nd143)', 'CD16(Sm147)', 'panCytokeratin(Nd148)', 'XCR1(Sm149)',
       'FOXP3(Gd155)', 'CD4(Gd156)', 'CD68(Tb159)', 'NKp46(Dy161)',
       'CD8a(Dy162)', 'CD3(Yb170)', 'CD66b(Yb174)', 'GATA3(Yb171)', 'EpCAM(Nd144)']

print(use_highly_variable)

adata.var["highly_variable"] = [True if ele in use_highly_variable else False for ele in adata.var_names]

# Set the random seed for reproducibility
np.random.seed(1)  

#calculate  principle components
sc.tl.pca(adata, use_highly_variable=True)
sc.external.pp.harmony_integrate(adata, 'mob_id')

sc.pp.neighbors(adata,use_rep='X_pca_harmony',n_neighbors = 60)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=1)



# Rank most differentially expressed genes between clusters
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon') #Find out if you can change the layers here to for instance the log. Might be more informative
#Try different tests as well. T-test over dispersion(‘t-test_overestim_var’); (wilcoxon), log regression. Try as well on the raw signal. WIlcoxon better than t-test on raw
sc.pl.rank_genes_groups(adata, n_genes=12, sharey=False, save='filter.pdf')


sc.pl.umap(adata, color=['leiden'],legend_loc='on data',save='filter.pdf')

adata.write("/home/tessdiv/limitedAccessBladder/NOBACKUP/IMC_BCG/Tessa_test/final_analysis/May3/filter.h5ad") 