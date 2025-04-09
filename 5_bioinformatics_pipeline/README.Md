## Bioinformatics pipeline to project in vitro models into a spatiotemporal framework

To project additional scRNA-seq datasets into the spatiotemporal atlas, we used StabMAP and reducedMNN. Principal components were calculated for the gastruloid cells after performing rPCA integration using Seurat v4 as described above in the ‘Calculating Principal Components’ sections. Briefly, multiple rounds of rPCA integration were performed, first across samples within each gastruloid timepoint and then across timepoints. Next, to align all gastruloid cells with the spatiotemporal atlas cells, we applied StabMAP followed by reducedMNN. For StabMAP, we used the top 30 principal components calculated for the scRNA, seqFISH and query gastruloid datasets, as described above, with projectALL = TRUE. The resulting StabMAP embedding was then reweighted to ensure equal contribution from all three datasets using reWeightEmbedding(). To correct for any remaining technical differences between the seqFISH, scRNA and gastruloid datasets, we applied reducedMNN_batchFactor() with k = 10, setting the batch factor to either seqFISH, scRNA or gastruloid as appropriate. Cell type labels, and imputed AP and DV values for the gastruloid cells were assigned as described in the ‘Label Transfer across scRNA and seqFISH Datasets” and ‘Gene Expression and AP/DV Coordinate Imputation for seqFISH and scRNA Cells’ sections above. 

### Statistical Analysis of Cell Type Distribution along the AP Axis

To quantitatively assess differences in cell type distribution between query cell and scRNA cells along the anteroposterior (AP) axis, we implemented two complementary statistical tests into the bioinformatics pipeline. 
1.	Kolmogorov-Smirnov (KS) Test: To determine whether the distribution of cells along the AP axis differs significantly between query and scRNA-seq cells, we applied a KS test. This test provides a statistical measure of divergence in AP positioning between the two datasets. 
2.	Jaccard Index Analysis: To further assess similarity in cell type composition, we calculated Jaccard indices across discrete bins spanning the AP axis. This analysis identifies regions long the AP axis where cell type composition is similar or divergent in query versus scRNA datasets.



### Spatial QC in the Bioinformatics Pipeline

To construct the spatiotemporal atlas, we integrated seqFISH and scRNA-seq datasets into a shared low-dimensional embedding and performed high-resolution clustering (see UMAP Generation and Joint Clustering for Spatiotemporal Atlas sub-section in scRNA-seq Analyses section). Poor joint clusters were identified as clusters where >98% of cells originated from only one dataset (see Supplementary Figure 3A). scRNA-seq cells in these clusters lacked sufficient neighboring seqFISH cells, making spatial imputation (e.g., AP or DV values) unreliable. This discrepancy arises because the scRNA dataset includes tissues absent from seqFISH slices (e.g., mature yolk sac endothelium and mesothelium). When projecting additional query datasets, such as gastruloids, onto the spatiotemporal atlas, poor joint cluster label information is transferred to the query cells. Query cells receiving the poor joint cluster label = TRUE are identified as those cells failing spatial QC (Figure 4C), as their spatial information is unreliable due to alignment with scRNA-seq cells that lack spatial reference points in seqFISH. For instance, in Figure 4, some gastruloid cardiomyocytes and yolk sac endothelium fail spatial QC because their scRNA-seq counterparts lack seqFISH equivalents. 

### Data required can be downloaded from this link: http://shiny.maths.usyd.edu.au/SpatiotemporalMouseAtlas/

To run StabMAP, you need to provide principal components and logcounts for seqFISH, scRNA and query cells. 
Please find  links to the PC and logcounts files for seqFISH and scRNA datasets at the link above.

Principal Components for seqFISH and scRNA cells:

seqFISH_pca.rds
atlas_pca.rds

Logcounts matrices for seqFISH and scRNA cells:

atlas_logcounts_ss.rds
seqFISH_logcounts_ss.rds
