## Performing StabMAP and reducedMNN to Align seqFISH and scRNA Datasets
To align all seqFISH cells that passed quality control with the downsampled scRNA cells from the extended gastrulation atlas, we applied StabMAP15 followed by reducedMNN16. For StabMAP, we used the top 30 principal components calculated for both the scRNA and seqFISH datasets, as described above, with projectALL = TRUE. The resulting StabMAP embedding was then reweighted to ensure equal contribution from both datasets using reWeightEmbedding(). To correct for any remaining technical differences between the seqFISH and scRNA datasets, we applied reducedMNN_batchFactor() with k = 10, setting the batch factor to either seqFISH or scRNA as appropriate.

## UMAP Generation and Joint Clustering for Spatiotemporal Atlas
A low-dimensional UMAP was generated to visualize both scRNA and seqFISH cells within the same space, using runUMAP applied to the batch-corrected StabMAP embedding. Joint clustering and subsequent sub-clustering were performed on this embedding with FindNeighbors (dims = 1:60) and FindClusters (resolution = 1 for clustering; resolution = 3 for sub-clustering). As a quality control step, clusters containing more than 98% or fewer than 2% of either seqFISH or scRNA cells were labeled as ‘poor joint clusters.’

## Label Transfer across scRNA and seqFISH Datasets
We applied K-nearest neighbors (K = 5) to classify the seqFISH cells based on cell type annotations from the extended mouse atlas, incorporating both extended and original atlas labels, as well as embryonic stage and subdissection labels. Additionally, we used the same K-nearest neighbors approach (K = 5) to classify the scRNA cells according to the cell types assigned to the seqFISH cells.

## Gene Expression Imputation for seqFISH Cells
We leveraged the scRNA atlas data, which provides whole transcriptome counts, to impute full gene expression profiles for the seqFISH-resolved cells. For each seqFISH-resolved cell, we identified the K-nearest neighbors (K = 5) and calculated the mean expression vector across all genes from the extended mouse atlas dataset.


### StabMAP output can be downloaded from this link: http://shiny.maths.usyd.edu.au/SpatiotemporalMouseAtlas/

The output from 1_StabMAP_scRNA_seqFISH.ipynb (StabMAP_allstages_postQC_reweighted_MNNcorrected_original_counts_without_neigh.Rds) which contains the low-dimensional embedding 
for all seqFISH and scRNA cells in the spatiotemporal atlas can be downloaded from the front page of the interactive exploration Shiny app above.
