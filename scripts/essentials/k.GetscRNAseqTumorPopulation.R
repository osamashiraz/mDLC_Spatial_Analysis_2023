# load("../inputs/MDLC-1_SeuratObj.Rdata")

# # EPI CELLS
# MDLC1_scRNAseq_epi = subset(AL14, CellTypes_new %in% c("Ductal", "Lobular"))
# MDLC1_scRNAseq_epi = NormalizeData(MDLC1_scRNAseq_epi,verbose = T)
# MDLC1_scRNAseq_epi = FindVariableFeatures(MDLC1_scRNAseq_epi,verbose = T)
# MDLC1_scRNAseq_epi = ScaleData(MDLC1_scRNAseq_epi,features = rownames(MDLC1_scRNAseq_epi@assays$RNA@counts), verbose = T)
# dims = 1:20
# set.seed(123); MDLC1_scRNAseq_epi  = RunPCA(MDLC1_scRNAseq_epi,verbose = T) %>%
#   RunTSNE(dims=dims,verbose=T) %>%
#   FindNeighbors(dims=dims,verbose=T)
# MDLC1_scRNAseq_epi = FindClusters(MDLC1_scRNAseq_epi, resolution = c(0.1,0.2 ,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), verbose=T)
# set.seed(123); MDLC1_scRNAseq_epi = RunUMAP(MDLC1_scRNAseq_epi, dims = dims)

# save(MDLC1_scRNAseq_epi, file = "../inputs/MDLC-1_tumorOnly_SeuratObj.Rdata")


library(Seurat)

load("../inputs/MDLC-1_tumorOnly_SeuratObj.Rdata")

MDLC1_scRNAseq_epi = ScaleData(MDLC1_scRNAseq_epi, features = rownames(MDLC1_scRNAseq_epi))

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

MDLC1_scRNAseq_epi <- CellCycleScoring(MDLC1_scRNAseq_epi, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)

cell_cols = c("Lobular"="red","Ductal"="blue", "Endo"="#6D4C41","Fibro"="#A1887F",
              "Smooth Muscle"="#D7CCC8","Macro"="#FF6F00","T cells"="#FDD835",
              "Pre-invasive"="#7E57C2","Non-invasive"="#B39DDB")

cols = c("lightgray", "#ca9446", "#543005")

