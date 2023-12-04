load("../inputs/TCGA_ILCvNST.RData")

library(DESeq2)
ILC_IDC_res = as.data.frame(ILC_IDC_res) # ILC_IDC_res

table(ILC_IDC_res$log2FoldChange > 0)

TCGA_ILC = intersect(rownames(subset(ILC_IDC_res, log2FoldChange > 1.25 & padj <= 0.05)),
                     rownames(mDLC_DSP_log2))
TCGA_NST = intersect(rownames(subset(ILC_IDC_res, log2FoldChange < -1.25  & padj <= 0.05)),
                     rownames(mDLC_DSP_log2))


lobular_all = subset(LvsD_ALL_er, logFC > 0.25 & adj.P.Val <= 0.3)$label
ductal_all = subset(LvsD_ALL_er, logFC < -0.25 & adj.P.Val <= 0.3)$label