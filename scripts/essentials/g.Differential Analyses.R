


# LIMMA DIFFERENTIAL ANALYSIS - PAN-PATIENT ALL ER+ -------------------------

toRemove = subset(ROI_annot, SampleID %in% "MDLC3" & SegmentLabel %in% "Ductal") %>% rownames()
sampleIDs = colnames(mDLC_DSP_log2)  # ALL SAMPLES
sampleIDs = setdiff(sampleIDs, toRemove)

## REMOVE COMMENTS TO RERUN

# mat_limma = as.matrix(mDLC_DSP_log2[, sampleIDs])

# # DESIGN
# Histology = factor(ROI_annot[sampleIDs, "SegmentLabel"])
# PatientID = factor(ROI_annot[sampleIDs, "SampleID"])
# design <- model.matrix(~0 + Histology + PatientID)
# dim(design)
# 
# colnames(design)
# colnames(design)[c(1,2)] <- c("Ductal","Lobular")
# 

# # RUN LIMMA
# library(limma)
# 
# v = voom(mat_limma, design)
# 
# LvsD_ALL_er = run_Limma(fit = lmFit(v), mat = mat_limma,
#                          contrast = makeContrasts(Lobular-Ductal, levels=design))
# LvsD_ALL_er$gene = rownames(LvsD_ALL_er)
# LvsD_ALL_er[1:5, ]
# 
# # SAVE RESULTS
# write.csv(LvsD_ALL_er, file = "../outputs/ERpos_LvsD_Limma.csv", row.names = F, quote = F)



LvsD_ALL_er = read.csv("../outputs/ERpos_LvsD_Limma.csv")
rownames(LvsD_ALL_er) = LvsD_ALL_er$label
LvsD_ALL_er_sig = subset(LvsD_ALL_er, adj.P.Val <= 0.3 & abs(logFC) > 0.25)


# ----------------------------------------------------------------------------------------------------









# LIMMA DIFFERENTIAL ANALYSIS - PATIENT WISE ---------------------------------------------------------------------------



## MDLC-1 DEG ANALYSIS ---------------------------------------------------------------------------

# sampleIDs = subset(ROI_annot[colnames(mDLC_DSP_log2),], SampleID %in% "MDLC1") %>% rownames()  # MDLC1
# mat_limma = mDLC_DSP_log2[, sampleIDs]
# 
# # DESIGN
# Histology = factor(ROI_annot[sampleIDs, "SegmentLabel"])
# design <- model.matrix(~0 + Histology)
# dim(design)
# 
# colnames(design)
# colnames(design)[c(1,2)] <- c("Ductal","Lobular")
# 
# v = voom(mat_limma, design)
# LvsD_MDLC1 = run_Limma(fit = lmFit(v), mat = mat_limma, 
#                      contrast = makeContrasts(Lobular-Ductal, levels=design))
# 
# LvsD_MDLC1[LvsD_MDLC1$label %in% "CDH1",]

# write.csv(LvsD_MDLC1, file = "../outputs/MDLC1_LvsD_Limma.csv", row.names = F, quote = F)


LvsD_MDLC1 = read.csv("../outputs/MDLC1_LvsD_Limma.csv")
rownames(LvsD_MDLC1) = LvsD_MDLC1$label
LvsD_MDLC1_sig = subset(LvsD_MDLC1, adj.P.Val <= 0.3 & abs(logFC) > 0.25)







## MDLC-2 DEG ANALYSIS ---------------------------------------------------------------------------


# sampleIDs = subset(ROI_annot[colnames(mDLC_DSP_log2),], SampleID %in% "MDLC2") %>% rownames()  # MDLC2
# mat_limma = mDLC_DSP_log2[, sampleIDs]
# 
# # DESIGN
# Histology = factor(ROI_annot[sampleIDs, "SegmentLabel"])
# design <- model.matrix(~0 + Histology)
# dim(design)
# 
# colnames(design)
# colnames(design)[c(1,2)] <- c("Ductal","Lobular")
# 
# v = voom(mat_limma, design)
# LvsD_MDLC2 = run_Limma(fit = lmFit(v), mat = mat_limma, 
#                      contrast = makeContrasts(Lobular-Ductal, levels=design))
# 
# LvsD_MDLC2[LvsD_MDLC2$label %in% "CDH1",]
# 
# write.csv(LvsD_MDLC2, file = "../out/Limma/MDLC2_LvsD_Limma.csv", row.names = F, quote = F)

LvsD_MDLC2 = read.csv("../outputs/MDLC2_LvsD_Limma.csv")
rownames(LvsD_MDLC2) = LvsD_MDLC2$label
LvsD_MDLC2_sig = subset(LvsD_MDLC2, adj.P.Val <= 0.3 & abs(logFC) > 0.25)






## MDLC-3 DEG ANALYSIS ---------------------------------------------------------------------------

# sampleIDs = subset(ROI_annot[colnames(mDLC_DSP_log2),], SampleID %in% "MDLC3") %>% rownames()  # MDLC3
# mat_limma = mDLC_DSP_log2[, sampleIDs]
# 
# # DESIGN
# Histology = factor(ROI_annot[sampleIDs, "SegmentLabel"])
# design <- model.matrix(~0 + Histology)
# dim(design)
# 
# colnames(design)
# colnames(design)[c(1,2)] <- c("Ductal","Lobular")
# 
# v = voom(mat_limma, design)
# LvsD_MDLC3 = run_Limma(fit = lmFit(v), mat = mat_limma, 
#                      contrast = makeContrasts(Lobular-Ductal, levels=design))
# 
# LvsD_MDLC3[LvsD_MDLC3$label %in% "CDH1",]
# 
# write.csv(LvsD_MDLC3, file = "../out/Limma/MDLC3_LvsD_Limma.csv", row.names = F, quote = F)

LvsD_MDLC3 = read.csv("../outputs/MDLC3_LvsD_Limma.csv")
rownames(LvsD_MDLC3) = LvsD_MDLC3$label
LvsD_MDLC3_sig = subset(LvsD_MDLC3, adj.P.Val <= 0.3 & abs(logFC) > 0.25)




