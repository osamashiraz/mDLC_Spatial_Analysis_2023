#  GENESET SELECTION & FILTERING BY DEGS
gset_ = c(Hallmark$genesets, KEGG$genesets[c("Adherens Junction")])
names(gset_) = c(paste0("Hallmark - ", names(Hallmark$genesets)), "KEGG - Adherens Junction")

# db_group = c(paste0("Hallmark - ", names(Hallmark$genesets)), "KEGG - Adherens Junction")

# # remove CDH1 - gives sig ER signaling with LvsD_ALL_er_sig$genes
gset_ = sapply(gset_, function(x){
  x = setdiff(x, "CDH1")
  return(intersect(x, rownames(mDLC_DSP_log2)))
})

degs = subset(LvsD_ALL_er, absLogFC > 0.25 & P.Value < 0.1)$label 

gset_sub = sapply(gset_, function(x){
  return(intersect(x, degs))
})

names(sapply(gset_sub, length))[sapply(gset_sub, length) >= 3] -> selectPathways

gset_sub = gset_sub[selectPathways]
length(gset_sub)


# SAMPLE SELECTION
toRemove = subset(ROI_annot_clean, SampleID %in% "2460" & SegmentLabel %in% "Ductal") %>% rownames()
sampleIDs = colnames(mDLC_DSP_log2)  # ALL SAMPLES
sampleIDs = setdiff(sampleIDs, toRemove)

# RUN GSVA
library(GSVA)
sig_score = gsva(expr = as.matrix(mDLC_DSP_log2[, sampleIDs]), 
                 gset.idx.list = gset_sub,
                 method="gsva", kcdf = "Poisson", parallel.sz = 1)

EA_DSP_OUT = enrichmentAnalysis(sig_score, 0.15, subset(ROI_annot[colnames(sig_score), ], SegmentLabel %in% "Ductal") %>% rownames(),
                                subset(ROI_annot[colnames(sig_score), ], SegmentLabel %in% "Lobular") %>% rownames())

EA_DSP_df = as.data.frame(EA_DSP_OUT)
EA_DSP_df$pathway = names(EA_DSP_OUT$group)
EA_DSP_df$db = sapply(EA_DSP_df$pathway, FUN = function(x){
  strsplit(x, split = " - ")[[1]][1]
})

EA_DSP_df = EA_DSP_df[,c("pathway","db","lobular", "ductal", "estimate","group","pvalues", "pvalues_simple", "fdr", "fdr_simple")]

estimate = EA_DSP_OUT$estimate
group = EA_DSP_OUT$group
pvalues_cap = EA_DSP_OUT$pvalues_cap
pvalues = EA_DSP_OUT$pvalues
pvalues_simple = EA_DSP_OUT$pvalues_simple
fdr = EA_DSP_OUT$fdr
fdr_simple = EA_DSP_OUT$fdr_simple
fdr_cap = EA_DSP_OUT$fdr_cap

names(pvalues)[pvalues <= 0.05 & group != "None"] -> SigPathways



bot_annot2 = function(sampleIDs){HeatmapAnnotation(annotation_name_side = "right",
                                                   Patient = patient_map[ROI_annot[sampleIDs, "SampleID"]],
                                                   Histology = ROI_annot[sampleIDs, "SegmentLabel"],
                                                   # PAM50 = PAM50.subtype$subtype[sampleIDs],
                                                   #ROI_ID = as.factor(ROI_annot[sampleIDs, "ROILabel"]),
                                                   annotation_name_gp = gpar(fontsize = 10, fontface="bold", col="#616161"),
                                                   col = annot_col)}

# indx <- sapply(gset_sub, length)
# res <- as.data.frame(do.call(rbind,lapply(gset_sub, `length<-`,
#                                           max(indx))))
# res[is.na(res)] = " "
# 
# res_df = as.data.frame(res)
# res_df = data.frame(pathway = names(gset_sub), res_df)
# res_df[1:5, ]
# rownames(res_df) = names(gset_sub)
# 
# EA_DSP_df$hitCounts = sapply(gset_sub[EA_DSP_df$pathway], length)
# EA_DSP_df$hits = res_df[match(EA_DSP_df$pathway, res_df$pathway), -1]
# write.table(EA_DSP_df, file = "../outputs/hallmarksignatures DSP.tsv", quote = F, row.names = F, sep = "\t")
