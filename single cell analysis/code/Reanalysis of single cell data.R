# This code is for re-analysis of the single-cell data from Hauser et al., iScience 2020.
# Written by Seungmin Han (23rd Nov 2022)

library(scater)
library(scran)
library(Seurat)

# Preparing data object
# The data from Hauser et al., iScience 2020. (https://doi.org/10.6084/m9.figshare.13157687)
Embryo_Seu <- readRDS("Embryonic SMG Integrated.rds")
Embryo_Seu@meta.data$CellType <- as.character(Embryo_Seu@meta.data$CellType)
# The data from Hauser et al., iScience 2020. (https://doi.org/10.6084/m9.figshare.13157726)
Postnatal_Seu <- readRDS("Postnatal SMG Integrated.rds")
Postnatal_Seu@meta.data$CellType <- as.character(Postnatal_Seu@meta.data$CellType)
tmp_idx <- Embryo_Seu@meta.data$CellType %in% c('Myoepithelial', 'End bud', 'Basal duct', 'Krt19+ duct')
epi_embryo_Seu <- Embryo_Seu[,tmp_idx]

tmp_idx <- Postnatal_Seu@meta.data$CellType %in% c('Acinar', 'Ascl3+ duct', 'Basal duct', 'Bpifa2+', 'Bpifa2+ Proacinar',
                                                   'GCT', 'Intercalated duct', 'Krt19+ duct', 'Mitotic cells', 'Myoepithelial',
                                                   'Smgc+', 'Smgc+ Proacinar', 'Striated duct')
epi_Postnatal_Seu <- Postnatal_Seu[,tmp_idx]
tmp_idx <- epi_Postnatal_Seu@meta.data$stage %in% c('P1')
epi_P1_Seu <- epi_Postnatal_Seu[,tmp_idx]
DefaultAssay(object = epi_embryo_Seu) <- "RNA"
DefaultAssay(object = epi_P1_Seu) <- "RNA"
epi_embryo_sce <- as.SingleCellExperiment(epi_embryo_Seu)
epi_P1_sce <- as.SingleCellExperiment(epi_P1_Seu)
common_genes <- intersect( rownames(epi_embryo_Seu), rownames(epi_P1_Seu) )
tmp_idx <- rownames(epi_embryo_sce) %in% common_genes
epi_embryo_sce <- epi_embryo_sce[tmp_idx,]
tmp_idx <- rownames(epi_P1_sce) %in% common_genes
epi_P1_sce <- epi_P1_sce[tmp_idx,]
epi_P1_sce <- epi_P1_sce[rownames(epi_embryo_sce),]

colData(epi_embryo_sce)$RNA_snn_res.0.6 <- NULL
colData(epi_embryo_sce)$RNA_snn_res.0.9 <- NULL
colData(epi_embryo_sce)$RNA_snn_res.0.5 <- NULL
colData(epi_embryo_sce)$integrated_snn_res.0.6 <- NULL
colData(epi_P1_sce)$RNA_snn_res.0.9 <- NULL
colData(epi_P1_sce)$integrated_snn_res.0.9 <- NULL
colData(epi_P1_sce)$celltype.stage.sex <- NULL
colData(epi_P1_sce)$integrated_snn_res.0.6 <- NULL
colData(epi_P1_sce)$sex <- NULL
epi_EmbryoP1_sce <- cbind(epi_embryo_sce, epi_P1_sce)
assays(epi_EmbryoP1_sce)$logcounts <- NULL
reducedDims(epi_EmbryoP1_sce)$PCA <- NULL
reducedDims(epi_EmbryoP1_sce)$UMAP <- NULL

# Subsetting the data into E14 & E16
tmp_idx <- epi_EmbryoP1_sce$stage %in% c('E14', 'E16')
epi_E14E16_sce <- epi_EmbryoP1_sce[,tmp_idx]

# Normalisation by the deconvolution method
set.seed(1000)
clusters <- quickCluster(epi_E14E16_sce, min.mean=0.1, method="igraph")
epi_E14E16_sce <- computeSumFactors(epi_E14E16_sce, cluster=clusters, min.mean=0.1)
epi_E14E16_sce <- normalize(epi_E14E16_sce)

# Cell cycle analysis
library(org.Mm.eg.db)
ensembl <- mapIds(org.Mm.eg.db, keys=rownames(epi_E14E16_sce), keytype="SYMBOL", column="ENSEMBL")
set.seed(100)
mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", 
                                package="scran"))
assignments <- cyclone(epi_E14E16_sce, mm.pairs, gene.names=ensembl, assay.type="logcounts")
options(na.action="na.omit")
design <- model.matrix(~ G1 + G2M, assignments$score)
tmp_NAidx <- union(which(is.na(assignments$score$G1)), which(is.na(assignments$score$G2M)))
tmp_idx <- setdiff(c(1:dim(epi_E14E16_sce)[2]),tmp_NAidx)
tmp_sce <- epi_E14E16_sce[,tmp_idx]
fit.block <- trendVar(epi_E14E16_sce, design=design, parametric=TRUE, use.spikes=NA)
dec.block <- decomposeVar(epi_E14E16_sce, fit.block)
library(limma)
epi_E14E16_sce.block <- epi_E14E16_sce
assay(epi_E14E16_sce.block, "corrected") <- removeBatchEffect(
  logcounts(epi_E14E16_sce), covariates=design[,-1])
epi_E14E16_sce.block <- denoisePCA(epi_E14E16_sce.block, technical=dec.block, 
                                   assay.type="corrected")
epi_E14E16_sce$G1score <- epi_E14E16_sce.block$G1score <- assignments$score$G1[tmp_idx]
epi_E14E16_sce$Sscore <- epi_E14E16_sce.block$Sscore <- assignments$score$S[tmp_idx]
epi_E14E16_sce$G2Mscore <- epi_E14E16_sce.block$G2Mscore <- assignments$score$G2M[tmp_idx]
epi_E14E16_sce$phases <- epi_E14E16_sce.block$phases <- assignments$phases[tmp_idx]
fit <- trendVar(epi_E14E16_sce, parametric=TRUE, use.spikes=NA) 
epi_E14E16_sce <- denoisePCA(epi_E14E16_sce, technical=fit$trend)
fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))
out <- plotReducedDim(epi_E14E16_sce, use_dimred="PCA", ncomponents=2, colour_by="G1score", 
                      size_by="G2Mscore") + fontsize + ggtitle("Before removal")
out2 <- plotReducedDim(epi_E14E16_sce.block, use_dimred="PCA", ncomponents=2, 
                       colour_by="G1score", size_by="G2Mscore") + fontsize + 
  ggtitle("After removal")
multiplot(out, out2, cols=2)

# Dimension reduction
set.seed(1000)
sepi_E14E16_sce.block <- runTSNE(epi_E14E16_sce.block, use_dimred="PCA", perplexity=50)
epi_E14E16_sce.block <- runUMAP(epi_E14E16_sce.block, use_dimred="PCA")

# Removing batch effect using Harmony
library(harmony)
epi_E14E16_sce.block_v2 <- epi_E14E16_sce.block %>% RunHarmony(c("stage"), plot_convergence = TRUE)
set.seed(1000)
epi_E14E16_sce.block_v2 <- runTSNE(epi_E14E16_sce.block_v2, use_dimred="HARMONY", perplexity=50)
epi_E14E16_sce.block_v2 <- runUMAP(epi_E14E16_sce.block_v2, use_dimred="HARMONY")

# Clustering
k_values <- c(5:30)
all_modularity <- vector(mode="double", length=0)
all_cluster_no <- vector(mode="integer", length=0)
for (k_tmp in k_values) {
  snn.gr <- buildSNNGraph(epi_E14E16_sce.block_v2, k = k_tmp, use.dimred="HARMONY")
  cluster.out <- igraph::cluster_walktrap(snn.gr)
  my.clusters <- cluster.out$membership
  key_name1 <- paste0("cluster_E14E16_CC_k",as.character(k_tmp))
  colData(epi_E14E16_sce.block_v2)[key_name1] <- factor(my.clusters)
  all_modularity <- append( all_modularity, igraph::modularity(cluster.out) )
  all_cluster_no <- append( all_cluster_no, dim(table(my.clusters)) )
}
louvain_E14E16_CC_v2_df <- data.frame(k_value = k_values, cluster_no = all_cluster_no, modularity = all_modularity, stringsAsFactors = FALSE)

# Cell type annotations
colData(epi_E14E16_sce.block_v2)$celltype_SH <- vector(mode='character', length=dim(epi_E14E16_sce.block_v2)[2])
tmp_idx <- epi_E14E16_sce.block_v2$cluster_E14E16_CC_k7 %in% c('6', '11')
epi_E14E16_sce.block_v2$celltype_SH[tmp_idx] <- c('A. End_bud')
tmp_idx <- epi_E14E16_sce.block_v2$cluster_E14E16_CC_k7 %in% c('2')
epi_E14E16_sce.block_v2$celltype_SH[tmp_idx] <- c('B. Myoepithelial')
tmp_idx <- epi_E14E16_sce.block_v2$cluster_E14E16_CC_k7 %in% c('1')
epi_E14E16_sce.block_v2$celltype_SH[tmp_idx] <- c('C. Basal_duct')
tmp_idx <- epi_E14E16_sce.block_v2$cluster_E14E16_CC_k7 %in% c('5', '7', '10')
epi_E14E16_sce.block_v2$celltype_SH[tmp_idx] <- c('D. Proacinal')
tmp_idx <- epi_E14E16_sce.block_v2$cluster_E14E16_CC_k7 %in% c('4')
epi_E14E16_sce.block_v2$celltype_SH[tmp_idx] <- c('E. Luminal_duct')
tmp_idx <- epi_E14E16_sce.block_v2$cluster_E14E16_CC_k7 %in% c('3', '8', '9', '12')
epi_E14E16_sce.block_v2$celltype_SH[tmp_idx] <- c('F. EMT')

# Filtering out EMT cells
tmp_idx <- epi_E14E16_sce.block_v2$celltype_SH %in% c('F. EMT')
epi_E14E16_sce.block_v3 <- epi_E14E16_sce.block_v2[,!tmp_idx]

# Defining marker genes to run monocle
tmp_idx <- epi_E14E16_sce.block_v3$stage %in% c('E16')
epi_E16_sce.block_tmp <- epi_E14E16_sce.block_v3[,tmp_idx]
my.clusters <- colData(epi_E16_sce.block_tmp)$celltype_SH
markers_sub <- findMarkers(epi_E16_sce.block_tmp, my.clusters, direction="up")
allmarkers_E16 <- c()
allcelltypenames <- names( table(epi_E16_sce.block_tmp$celltype_SH) )
allcelltypenames <- setdiff(allcelltypenames, c('A. End_bud'))
for ( typename in allcelltypenames ) {
  marker_sub.set <- markers_sub[[typename]] 
  tmp_idx <- ( (marker_sub.set$FDR < 0.01) & (marker_sub.set$summary.logFC > 0.445) )
  tmp_genes <- rownames(marker_sub.set)[tmp_idx] 
  allmarkers_E16 <- union(allmarkers_E16, tmp_genes)   
}

# Running monocle
rowData(epi_E14E16_sce.block_v3)$gene_short_name <- rownames(epi_E14E16_sce.block_v3)
my_colData <- colnames( colData(epi_E14E16_sce.block_v3) )
epi_E14E16_cds.block <- convertTo( epi_E14E16_sce.block_v3, type=("monocle") )
epi_E14E16_cds.block <- estimateSizeFactors(epi_E14E16_cds.block)
epi_E14E16_cds.block <- estimateDispersions(epi_E14E16_cds.block)
epi_E14E16_cds.block <- setOrderingFilter(epi_E14E16_cds.block, allmarkers_E16)
epi_E14E16_cds.block <- reduceDimension(epi_E14E16_cds.block, max_components = 2, reduction_method = 'DDRTree', residualModelFormulaStr = "~stage", norm_method = 'log', verbose = F)
epi_E14E16_cds.block <- orderCells(epi_E14E16_cds.block, num_paths = 6)
pData(epi_E14E16_cds.block)$sample <- pData(epi_E14E16_cds.block)$celltype_SH
pData(epi_E14E16_cds.block)$sample_name <- NULL
plot_cell_trajectory(epi_E14E16_cds.block, color_by = "State") + 
  theme(aspect.ratio=1) + 
  scale_color_manual( values=c('#E41A1C', '#377EB8', '#4DAF4A', '#984EA3', '#FF7F00', '#FFDD33', '#A65628') )
plot_cell_trajectory(epi_E14E16_cds.block, color_by = "stage") +
  theme(aspect.ratio=1)
epi_E14E16_cds.block <- orderCells(epi_E14E16_cds.block, num_paths = 6, root_state = 6)
plot_cell_trajectory(epi_E14E16_cds.block, color_by = "Pseudotime") +
  theme(aspect.ratio=1)
