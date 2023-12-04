
# GENESET SELECTION & FILTERING BY DEGS
gset_ = c(GO_BP$genesets, GO_CC$genesets, GO_MF$genesets, Hallmark$genesets, KEGG$genesets, REACTOME$genesets)


# FILTER FOR SPECIFIC BREAST CANCER RELATED THEMES
gset_ = c(GO_BP$genesets, GO_CC$genesets, GO_MF$genesets, Hallmark$genesets, KEGG$genesets, REACTOME$genesets)
db_group = c(paste0("GP_BP - ", names(GO_BP$genesets)), paste0("GP_CC - ", names(GO_CC$genesets)),
             paste0("GP_MF - ", names(GO_MF$genesets)), paste0("Hallmark - ", names(Hallmark$genesets)),
             paste0("KEGG - ", names(KEGG$genesets)), paste0("REACTOME - ", names(REACTOME$genesets)))
gset_ = gset_[grep(x = names(gset_), pattern = "Collagen|ECM|Extracellular Matrix|Ncam1|Immunoglobulin|Mediated Immunity|Leptin|Glucose Metaboli|Gluconeogenesis|Mitotic|Cell Cycle|Senescence|Proteoglycan|Tp53|Proteosome|Mhc Class I|Antigen Processing|Oxidative Phosphorylation|Electron Transport Chain|Respiration|Sumoylation|Acetyltransferase|Oxidative Stress|Double Strand|DNA Repair|Mtor|Phosphatidylinositol|Protein Maturation|Unfolded Protein Response|Mapk Targets")]

db_group = db_group[grep(x = db_group, pattern = "Collagen|ECM|Extracellular Matrix|Ncam1|Immunoglobulin|Mediated Immunity|Leptin|Glucose Metaboli|Gluconeogenesis|Mitotic|Cell Cycle|Senescence|Proteoglycan|Tp53|Proteosome|Mhc Class I|Antigen Processing|Oxidative Phosphorylation|Electron Transport Chain|Respiration|Sumoylation|Acetyltransferase|Oxidative Stress|Double Strand|DNA Repair|Mtor|Phosphatidylinositol|Protein Maturation|Unfolded Protein Response|Mapk Targets", value = F)]
names(db_group) = names(gset_)

# # remove CDH1
gset_ = sapply(gset_, function(x){
  x = setdiff(x, "CDH1")
  return(intersect(x, rownames(mDLC_DSP_log2)))
})

degs = subset(LvsD_ALL_er, absLogFC > 0.25 & P.Value < 0.1)$label 


# FILTER GENESETS BY DEGS
gset_sub = sapply(gset_, function(x){
  return(intersect(x, degs))
})

names(sapply(gset_sub, length))[sapply(gset_sub, length) >= 3] -> selectPathways

gset_sub = gset_sub[selectPathways]
length(gset_sub)


# PERFORM GSVA
toRemove = subset(ROI_annot_clean, SampleID %in% "2460" & SegmentLabel %in% "Ductal") %>% rownames()
sampleIDs = colnames(mDLC_DSP_log2)  # ALL SAMPLES
sampleIDs = setdiff(sampleIDs, toRemove)

library(GSVA)
sig_score = gsva(expr = as.matrix(mDLC_DSP_log2[, sampleIDs]), 
                 gset.idx.list = gset_sub,
                 method="gsva", kcdf = "Poisson", parallel.sz = 1)





#  SIGNIFICANCE TESTING 

EA_DSP_OUT = enrichmentAnalysis(sig_score, 0.15, subset(ROI_annot[colnames(sig_score), ], SegmentLabel %in% "Ductal") %>% rownames(),
                                subset(ROI_annot[colnames(sig_score), ], SegmentLabel %in% "Lobular") %>% rownames())

EA_DSP_df = as.data.frame(EA_DSP_OUT)
EA_DSP_df$pathway = names(EA_DSP_OUT$group)
EA_DSP_df$db = sapply(db_group[selectPathways], FUN = function(x){
  strsplit(x, split = " - ")[[1]][1]
})

table(EA_DSP_df$db)
EA_DSP_df = EA_DSP_df[,c("pathway","db","lobular", "ductal", "estimate","group","pvalues", "pvalues_simple", "fdr", "fdr_simple")]


estimate = EA_DSP_OUT$estimate
group = EA_DSP_OUT$group
pvalues_cap = EA_DSP_OUT$pvalues_cap
pvalues = EA_DSP_OUT$pvalues
pvalues_simple = EA_DSP_OUT$pvalues_simple
fdr = EA_DSP_OUT$fdr
fdr_simple = EA_DSP_OUT$fdr_simple
fdr_cap = EA_DSP_OUT$fdr_cap

table(pvalues <= 0.05)
names(pvalues)[pvalues <= 0.05 & group != "None"] -> SigPathways
length(SigPathways)

group_sub = group[group %notin% "None"]
pvalues_sub = pvalues[group %notin% "None"]

all_pathays = unique(names(group_sub[names(pvalues_sub)[pvalues_sub < 0.05]]))
length(all_pathays)
pvalues_sub[grep(x = names(pvalues_sub), pattern = "Adherens Junction", value = T)]

myList = gset_sub[unique(all_pathays)] # p < 0.05 


# BACKGROUND PERMUTATIONS - FILTERING

load("../outputs/5000_bg_permutations_select_mSIG.Rdata")

bg_sig = rowSums(pvalues_bg_repeats_5000 <= 0.05, na.rm = T)
bg_pct = 100*bg_sig/5000
bg_pct_simple = ifelse(bg_pct > 5, 1, 0)

table(bg_pct_simple)

bg_pct_annots = function(ids){
  HeatmapAnnotation(which = "row", `BG > 5%` = bg_pct_simple[ids], 
                    col = list(`BG > 5%` = c("1" = "red", "0" = "gray")))
}

EA_DSP_df$bg_pct = bg_pct[EA_DSP_df$pathway]

bg_noise = names(bg_pct_simple)[bg_pct_simple == 1]
myList = myList[names(myList) %notin% bg_noise]


# JACCARD SIMILARITY

jacMat = sapply(myList, function(x){
  return(sapply(myList, function(y) length(intersect(x, y))/(length(c(x,y))-length(intersect(x, y))), simplify = F))
})
mode(jacMat)="numeric"

colnames(jacMat) = names(myList)
rownames(jacMat) = names(myList)



# DEFINE SIGNATURE CLUSTERS BASED ON STRONG JI SIMILARITY

similarityCutoff = 0.5
memberShipCutoff = 2
jacMat_sub = jacMat[rowSums(jacMat>=similarityCutoff) >= memberShipCutoff, rowSums(jacMat >=similarityCutoff) >= memberShipCutoff]

set.seed(123); clusterObj_sub = hclust(dist(jacMat_sub), method = "ward.D2")

set.seed(123); memberShip_sub = cutree(clusterObj_sub, k = 20)

memberShip_sub[memberShip_sub %in% c(1,8)] = 1 # immune
memberShip_sub[memberShip_sub %in% c(2,7,9,11)] = 2 # cell cyle
memberShip_sub[memberShip_sub %in% c(3,10)] = 3 # oxidative stress
memberShip_sub[memberShip_sub %in% c(4,6)] = 4 # Leptin & glucose
memberShip_sub[memberShip_sub %in% c(14,17)] = 6 # Collagen
memberShip_sub[grep(x = names(memberShip_sub), pattern = "Collagen Metabolic|Collagen Biosynthetic", value = T)] = 7 # Collagen Metabolism
memberShip_sub[grep(x = names(memberShip_sub), pattern = "B Cell Mediated", value = T)] = 8 # Immune
memberShip_sub[grep(x = names(memberShip_sub), pattern = "Protein Maturation", value = T)] = 9 # Protein
memberShip_sub[grep(x = names(memberShip_sub), pattern = "Proteoglycan", value = T)] = 10 # Metabolism
memberShip_sub[memberShip_sub %in% c(5,13)] = 5 # cell cycle senescence
memberShip_sub[memberShip_sub %in% c(5,13)] = 5 # cell cycle senescence

memberShip_sub_new = memberShip_sub*0
memberShip_sub_new[memberShip_sub %in% 2] = 1 # cell cycle
memberShip_sub_new[memberShip_sub %in% 5] = 2 # senesence
memberShip_sub_new[memberShip_sub %in% 12] = 3 # DBS repair
memberShip_sub_new[memberShip_sub %in% 1] = 4 # Immune T cells
memberShip_sub_new[memberShip_sub %in% 15] = 5 # Immune Immunoglobublin
memberShip_sub_new[memberShip_sub %in% 8] = 6 # Immune B cells
memberShip_sub_new[memberShip_sub %in% 6] = 7 # Collagen
memberShip_sub_new[memberShip_sub %in% 3] = 8 # Metabolism -> Oxidative stress
memberShip_sub_new[memberShip_sub %in% 4] = 9 # Metabolism -> Leptin/glucose
memberShip_sub_new[memberShip_sub %in% 7] = 10 # Metabolism -> Collagen metabolism
memberShip_sub_new[memberShip_sub %in% 10] = 11 # Metabolism -> Proteoglycan metabolism
memberShip_sub_new[memberShip_sub %in% 20] = 12 # Oncogeneic Signaling - MAPK/ERk
memberShip_sub_new[memberShip_sub %in% 19] = 13 # Oncogeneic Signaling - mTOR
memberShip_sub_new[memberShip_sub %in% 9] = 14 # Protein maturation & modification - Protein Maturation
memberShip_sub_new[memberShip_sub %in% 16] = 15 # Protein maturation & modification - Acetyletransferase
memberShip_sub_new[memberShip_sub %in% 18] = 16 # Phosphotidylinositol


# set.seed(123); sigClust_HT = Heatmap(as.matrix(jacMat_sub), 
#                                       row_split = memberShip_sub, column_split = memberShip_sub,
#                                       row_order = clusterObj_sub$order, column_order = clusterObj_sub$order,
#                                       cluster_column_slices = F, cluster_columns = F, 
#                                       cluster_rows = F, show_heatmap_legend = F,
#                                       row_names_gp = gpar(fontsize = 8, fontface = "bold"),
#                                       col = colFun3, border = T, row_gap = unit(0,"cm"), column_gap = unit(0,'cm'),
#                                       top_annotation = enrichedInAnnot(colnames(jacMat_sub), label_pval = F),
#                                       left_annotation = enrichedInAnnot(colnames(jacMat_sub), 
#                                                                         which = "row", label_pval = F),
#                                       cluster_row_slices = F, show_row_names = T, show_column_names = F)
# sigClust_HT



# DEFINE RECURRING BIOLOGICAL THEMES/CLUSTERS
biological_clusters = ifelse(memberShip_sub_new %in% c(1,2,3), "Cell Cycle",
                             ifelse(memberShip_sub_new %in% c(4,5,6), "Immune",
                                    ifelse(memberShip_sub_new %in% c(7), "Collagen & ECM",
                                           ifelse(memberShip_sub_new %in% c(8,9,10,11), "Metabolism",
                                                  ifelse(memberShip_sub_new %in% c(12,13), "Oncogeneic Signaling",
                                                         ifelse(memberShip_sub_new %in% c(14,15), "Protein Maturation",
                                                                ifelse(memberShip_sub_new == 16, "Phosphatidylinositol",
                                                                       memberShip_sub_new)))))))

names(biological_clusters) = NULL
tbl = table(biological_clusters)

biological_clusters = factor(biological_clusters, levels = names(tbl)[order(tbl, decreasing = T)])
biological_clusters = factor(biological_clusters, 
                             levels = c("Cell Cycle","Immune", # Immune
                                        "Collagen & ECM", "Metabolism", # METABOLISM
                                        "Oncogeneic Signaling", # Oncogenic Pathways
                                        "Protein Maturation", "Phosphatidylinositol"))



## mean signature score


biological_themes_df = data.frame(cluster = memberShip_sub_new)
biological_themes_df$pathway = rownames(biological_themes_df)
rownames(biological_themes_df) = NULL
biological_themes_df$signature_clusters = 
  ifelse(biological_themes_df$cluster == 1, "L1: Cell Cycle",
         ifelse(biological_themes_df$cluster == 2, "L2: TP53 & Senescence",
                ifelse(biological_themes_df$cluster == 3, "D1: DSB/HR",
                       ifelse(biological_themes_df$cluster == 4, "D3: Lymphocyte Immunity",
                              ifelse(biological_themes_df$cluster == 5, "L4: Immunoglobulin",
                                     ifelse(biological_themes_df$cluster == 6, "L3: B-Cell Immunity",
                                            ifelse(biological_themes_df$cluster == 7, "L5: Collagen & ECM",
                                                   ifelse(biological_themes_df$cluster == 8, "D4: Oxidative Stress",
                                                          ifelse(biological_themes_df$cluster == 9, "L6: Leptin Signaling/Glucose Metabolism",
                                                                 ifelse(biological_themes_df$cluster == 10, "L7: Collagen Metabolism",
                                                                        ifelse(biological_themes_df$cluster == 11, "L8: Proteoglycan Metabolism",
                                                                               ifelse(biological_themes_df$cluster == 12, "L9: MAPK/ERK",
                                                                                      ifelse(biological_themes_df$cluster == 13, "D5: mTOR",
                                                                                             ifelse(biological_themes_df$cluster == 14, "D7: Protein Regulation",
                                                                                                    ifelse(biological_themes_df$cluster == 15, "D6: Acetyltransferase",
                                                                                                           ifelse(biological_themes_df$cluster == 16, "D8: Phosphatidylinositol",NA))))))))))))))))
biological_themes_df[grep(x = biological_themes_df$pathway, pattern = "Antigen"), "signature_clusters"] = "D2: Antigen Presentation" 

clusters = unique(biological_themes_df$signature_clusters)

clusters_scores = matrix(NA, nrow = length(clusters), ncol = ncol(sig_score))
rownames(clusters_scores) = clusters
colnames(clusters_scores) = colnames(sig_score)

#i = themes[1]
for (i in clusters){
  pathways = subset(biological_themes_df, signature_clusters == i)$pathway
  clusters_scores[i,] = colMeans(sig_score[pathways,])
}




# indx <- sapply(gset_sub, length)
# res <- as.data.frame(do.call(rbind,lapply(gset_sub, `length<-`,
#                                           max(indx))))
# res[is.na(res)] = " "
# 
# res_df = as.data.frame(res)
# res_df = data.frame(pathway = names(gset_sub), res_df)
# res_df[1:5, ]
# 
# 
# rownames(biological_themes_df) = biological_themes_df$pathway
# EA_DSP_df_sub = subset(EA_DSP_df, pathway %in% biological_themes_df$pathway)
# EA_DSP_df_sub$Clusters = biological_themes_df[EA_DSP_df_sub$pathway, "signature_clusters"]
# EA_DSP_df_sub$Clusters[is.na(EA_DSP_df_sub$Clusters)] = ""
# 
# EA_DSP_df_sub = EA_DSP_df_sub[,c("db","pathway","Clusters","lobular", "ductal","bg_pct", "estimate","group","pvalues", "pvalues_simple", "fdr", "fdr_simple", "bg_pct")]
# EA_DSP_df_sub$hitCounts = sapply(gset_sub[EA_DSP_df_sub$pathway], length)
# EA_DSP_df_sub$hits = res_df[match(EA_DSP_df_sub$pathway, res_df$pathway), -1]
# 
# write.table(EA_DSP_df_sub, file = "../outputs/biologicalThemes_DSP.tsv", quote = F, row.names = F, sep = "\t")



