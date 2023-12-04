<style type="text/css">
.main-container {
  max-width: 1800px;
  margin-left: auto;
  margin-right: auto;
}
</style>

# DISCLAIMER

This rmarkdown is produced to reproduce the bioinformatics analysis in
“Spatial molecular profiling of mixed invasive ductal-lobular breast
cancers reveals heterogeneity in intrinsic molecular subtypes, oncogenic
signatures, and mutations”. Run the script in sequential order. Make
sure packages listed in the session info text file are installed.
Changes in R and/or its package versions may result in changes in the
expected output.

# INSTALL PACKAGES

# PREPARATION

``` r
# PREPARE ENVIRONMENT
source("./essentials/a.Load_packages_and_functions.R")

# LOAD INPUT DSP DATA AND ANNOTATIONS
source("./essentials/b.Load_dataset.R")
```

#………………………………………………………………………………………………….

# A. PREPARE QC PLOTS

``` r
source("./essentials/c.Prepare_QC_plots.R")
```

## Sup Figure 1A - Sequencing Saturation

``` r
supFig1A
```

![](2023-09-28---Main_files/figure-markdown_github/unnamed-chunk-5-1.png)

``` r
# ggsave(plot = supFig1A, file = "../outputs/Sup Fig.1A - Sequencing Saturation.png", dpi = 800, width=5, height=3.5)
```

## Sup Figure 1B - Signal to Noise Ratio

``` r
supFig1B
```

![](2023-09-28---Main_files/figure-markdown_github/unnamed-chunk-6-1.png)

``` r
# ggsave(plot = supFig1B, file = "../outputs/Sup Fig.1B - Signal to Noise Ratio.png", dpi = 800, width=5, height=3.5)
```

## Sup Figure 1C - Count Artifacts Heatmap

``` r
supFig1C
```

![](2023-09-28---Main_files/figure-markdown_github/unnamed-chunk-7-1.png)

``` r
png(filename = "../outputs/Sup Fig.1C - Count Artifacts Heatmap.png", res = 200, width=1500, height=1000)
supFig1C
dev.off()
```

    ## png 
    ##   2

## Figure 1E - E-cadherin Intensity in Ductal vs Lobular Regions

``` r
ecad_imagej$group = paste0(ecad_imagej$Patient, " - ", ecad_imagej$Label)

ecad_svg = ggplot(ecad_imagej, aes(Label, Ecadherin.normalized)) + geom_bar(stat = "identity", aes(fill = Label)) + 
  myTheme + scale_fill_manual("Histology", values = annot_col$Histology)+ facet_grid(.~Patient) + xlab("")+
  ylab("Normalized Intensity") + ylim(c(0,3)) + 
  theme(axis.text.y = element_text(size = 15, color = "#7F7F7F"), 
        axis.text.x = element_text(size = 10, color = c("blue", "red"), hjust = c(0.65,0.45)))

ecad_svg
```

![](2023-09-28---Main_files/figure-markdown_github/unnamed-chunk-8-1.png)

``` r
# ggsave(plot = ecad_svg, file = "../outputs/Fig.1e - Ecadherin Intensites.png", dpi = 1000,  width=5, height=2.5)
```

## - - - - - - - - - - - - - - - - - - - - - - - - - - -

# B. EXPLORATORY TRANSCRIPTOMIC ANALYSIS

## I. SAMPLE-SAMPLE CLUSTERING - TSNE

``` r
source("./essentials/d.Generate_tsne.R")

plt_df_top15 = mkplotDF(tsne_res_top15$Y, sampleIDs)

supFig1D = ggplot(plt_df_top15,  aes(x = DIM_1, y = DIM_2, color = Histology, shape = Patient, 
                                     label = labels)) +
  myTheme + geom_point(size = 3, alpha = 0.8) + 
  scale_color_manual(values = c("Lobular"="red","Ductal"="blue"))
```

### Sup Figure 1D - TSNE

``` r
supFig1D
```

![](2023-09-28---Main_files/figure-markdown_github/unnamed-chunk-10-1.png)

``` r
# ggsave(plot = supFig1D, file = "../outputs/Sup Fig.1D - TSNE Plot.png", dpi = 800, width=5, height=3.5)
```

## II. CLINICAL MARKER EXPRESSION

``` r
goi = c("CDH1","ESR1","ERBB2","PGR")

df = data.frame(scale(t(mDLC_DSP_log2[goi, sampleIDs])),
                hist = ROI_annot[mDLC_DSP_log2[,sampleIDs] %>% colnames(), "SegmentLabel"],
                #ID = ROI_annot[mDLC_DSP_log2[, sampleIDs] %>% colnames(), "ROILabel"],
                Patient = patient_map[ROI_annot[mDLC_DSP_log2[,sampleIDs] %>% colnames(), "SampleID"]])
df_melt = melt(df)

SupFig1E = ggplot(df_melt, aes(hist, value, color = hist)) + geom_boxplot(alpha = 0.5, outlier.color = NA) + 
  ggbeeswarm::geom_quasirandom(alpha = 0.5, size = 4, shape = 20) + 
  myTheme + scale_color_manual("Histology", values = annot_col$Histology) + 
  facet_wrap(.~variable + Patient, ncol = 3) + 
  stat_compare_means(comparisons = list(c("Lobular","Ductal")), label = "p.signif", label.y = 1.5,
                     color = 'red', size = 6) + ylim(c(-2,2))
```

### Sup Figure 1E - Clinical Marker Expression

``` r
SupFig1E
```

![](2023-09-28---Main_files/figure-markdown_github/unnamed-chunk-12-1.png)

``` r
# ggsave(plot = SupFig1E, file = "../outputs/Sup Fig.1E - Clinical Marker Expression.png", dpi = 800, width=6, height=10)
```

## III. EMT SIGNATURE ANALYSIS

``` r
EMT_markers = c("FOXC1", "VIM", "CAV1", "CDH2", "CDH1", "SOX4", "SNAI2", "ZEB1","ZEB2", "TWIST1", "SNAI1") # adapted from pmid: 26510554

EMT_markers_df = t(scale(t(mDLC_DSP_log2)))[EMT_markers, ]
EMT_markers_df = melt(EMT_markers_df, factorsAsStrings = T)

EMT_markers_df$value[EMT_markers_df$value > 2] = 2
EMT_markers_df$value[EMT_markers_df$value < -2] = -2

EMT_markers_df$ID = as.character(EMT_markers_df$Var2)
EMT_markers_df$hist = ROI_annot[EMT_markers_df$ID, "SegmentLabel"]
EMT_markers_df$Patient = ROI_annot[EMT_markers_df$ID, "BIOS_ID"]

SupFig1F = ggplot(EMT_markers_df, aes(Var1, value, fill = hist)) + 
  geom_boxplot(alpha = 0.2, outlier.color = NA) + 
  geom_jitter(alpha = 1, size = 1, shape = 20, color = "black") + 
  myTheme + scale_fill_manual("Histology", values = annot_col$Histology) + 
  facet_wrap(.~Patient) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5)) + 
  xlab("") + ylab("Expression")
```

### Sup Figure 1F - EMT Markers Expression

``` r
SupFig1F
```

![](2023-09-28---Main_files/figure-markdown_github/unnamed-chunk-14-1.png)

``` r
# ggsave(file = "../outputs/Sup Fig.1F - EMT Markers.png", bg = "white", dpi = 1000,  plot = SupFig1F, width=10, height=3)
```

## IV. TRANSCRIPTOMIC LANDSCAPE HEATMAP

### a. PERFORM CONSENSUS CLUSTERING

### b. PERFORM INTRINSIC MOLECULAR SUBTYPES

``` r
source("./essentials/f.Perform_Molecular_Subtyping.R")
```

``` r
# PAM50 PANEL
bot_annot_side = function(sampleIDs, side){
  HeatmapAnnotation(which = side, #name = "Annotations",
                    Patient = ROI_annot[sampleIDs, "BIOS_ID"],
                    Histology = ROI_annot[sampleIDs, "SegmentLabel"], 
                    Clusters = as.numeric(order_df[sampleIDs,]$Cluster), #direction = "horizontal",
                    gap = unit(0,"mm"), border = F, #gp = gpar(col = "black"), 
                    simple_anno_size = unit(5,"mm"),
                    col = c(annot_col, list(Clusters = cluster_cols)))}


pam50score_mat = t(PAM50.subtype$subtype.proba)

set.seed(123);  pam50_heatmap = Heatmap(pam50score_mat[-5,], col = score_colors, name = "PAM50 score", 
                                        #width = unit(8, "cm"), 
                                        height = unit(2, "cm"),
                                        show_column_names = F, show_row_names = T, 
                                        border = T, cluster_row_slices = F, 
                                        cluster_column_slices = F, row_title = " ", 
                                        column_labels = ROI_annot[sampleIDs, "ROILabel"],
                                        border_gp = gpar(lwd = 2),column_gap = unit(0.1,"cm"),
                                        row_names_gp = gpar(fontsize = 10, fontface = "bold"),
                                        heatmap_legend_param = list(legend_direction = "horizontal"),
                                        show_row_dend = F, show_column_dend = F)



# TOP VARIABLE GENE PANLE
topVarMat = t(scale(t(mDLC_DSP_log2[select_features[1:100], rownames(order_df)])))
dim(topVarMat)
```

    ## [1] 100  26

``` r
ht_opt$TITLE_PADDING = unit(c(4, 4), "points")

set.seed(123);  topVar_heatmap = Heatmap(topVarMat, col = colFun, name = "mRNA Exp", 
                                         height = unit(8, "cm"),
                                         show_column_names = F, show_row_names = F,  
                                         cluster_column_slices = T, cluster_row_slices = T, 
                                         row_labels = ROI_annot[sampleIDs, "ROILabel"],
                                         border = T, column_gap = unit(0.1,"cm"), 
                                         column_dend_height = unit(0.8, "cm"), border_gp = gpar(lwd = 2),
                                         column_km = 7, column_km_repeats = 100, #
                                         column_title = " ", column_split = order_df$IDs,
                                         column_names_gp = gpar(fontsize = 10, fontface = "bold"),
                                         show_row_dend = F, show_column_dend = T, 
                                         top_annotation = bot_annot_side(colnames(topVarMat), 
                                                                         side = "column"))
```

### Figure 1A - Transcriptomic Clases and Molecular Subtypes

``` r
draw(topVar_heatmap %v% pam50_heatmap, ht_gap = unit(0.1, "cm"))
```

![](2023-09-28---Main_files/figure-markdown_github/unnamed-chunk-18-1.png)

``` r
png("../outputs/Fig.2A - Transcriptomic Clases and Molecular Subtypes.png", res = 1000, 
    width=5, height=10, units = "in")
draw(topVar_heatmap %v% pam50_heatmap, ht_gap = unit(0.1, "cm"))
dev.off()
```

    ## png 
    ##   2

## - - - - - - - - - - - - - - - - - - - - - - - - - - -

# C. DEG ANALYSIS

``` r
source("./essentials/g.Differential Analyses.R")
```

## I. PAN-PATIENT DEG PLOTS

``` r
# PAN-PATIENT DEG VOLCANO PLOTS
Fig2B_nolabs = volcanoPlt(res = LvsD_ALL_er, sig_res = LvsD_ALL_er_sig, pcut = 0.3, fc_cut = 0.25, topN = 0, annotsize=0) + xlim(c(-1.25, 1.25))

Fig2B = volcanoPlt(res = LvsD_ALL_er, sig_res = LvsD_ALL_er_sig, pcut = 0.3, fc_cut = 0.25, topN = 7, annotsize=3) + xlim(c(-1.25, 1.25))

# PAN-PATIENT DEG HEATMAP
labs = paste0(c("Ductal","Lobular"), " (N=", table(LvsD_ALL_er_sig$logFC > 0), ")")
goi_group = ifelse(LvsD_ALL_er_sig$logFC > 0, labs[2], labs[1])

sampleIDs = colnames(mDLC_DSP_log2)
sampleIDs = setdiff(colnames(mDLC_DSP_log2), toRemove)

LvD_All_er_ht = makeHeatmap(LvsD_ALL_er_sig$label, LvsD_ALL_er_sig, goi_group, 
                            sampleIDs, show_column_names = F)
```

### Figure 2B - Lobular vs Ductal DEG Heatmap

``` r
draw(LvD_All_er_ht[[1]] + LvD_All_er_ht[[2]] + LvD_All_er_ht[[3]], merge_legend = TRUE, gap = unit(0, "mm"), 
     heatmap_legend_side = "left")
```

![](2023-09-28---Main_files/figure-markdown_github/unnamed-chunk-21-1.png)

``` r
png("../outputs/Fig.2B - DEG Heatmap.png", res = 1000, width=3, height=6, units = "in")
draw(LvD_All_er_ht[[1]] + LvD_All_er_ht[[2]] + LvD_All_er_ht[[3]], merge_legend = TRUE, gap = unit(0, "mm"), heatmap_legend_side = "left")
dev.off()
```

    ## png 
    ##   2

### Figure 2B - Lobular vs Ductal DEG Volcano Plot

``` r
Fig2B
```

![](2023-09-28---Main_files/figure-markdown_github/unnamed-chunk-22-1.png)

``` r
# ggsave(plot = Fig2B, file = "../outputs/Fig.2B - DEG Volcano Plot.png", bg = "white",  dpi = 1000, width=6, height=3)

# ggsave(plot = Fig2B_nolabs, file = "../outputs/Fig.2B - DEG Volcano Plot - no labels.png", bg = "white", dpi = 1000, width=6, height=3)
```

## II. PATIENT WISE DEG PLOTS

``` r
# PATIENT-WISE ANALYSIS
## MDLC-1 ANALYSIS
SupFig2A = volcanoPlt(res = LvsD_MDLC1, sig_res = LvsD_MDLC1_sig, pcut = 0.3, fc_cut = 0.25, topN = 10, annotsize=3)

## MDLC-2 ANALYSIS
SupFig2B = volcanoPlt(res = LvsD_MDLC2, sig_res = LvsD_MDLC2_sig, pcut = 0.3, fc_cut = 0.25, topN = 10, annotsize=3)

## MDLC-3 ANALYSIS
SupFig2C = volcanoPlt(res = LvsD_MDLC3, sig_res = LvsD_MDLC3_sig, pcut = 0.3, fc_cut = 0.25, topN = 10, annotsize=3)



# Overlap of Patient DEGs with Pan-Patient DEGs
lobular_all = subset(LvsD_ALL_er, logFC > 0.25 & adj.P.Val <= 0.3)$label
ductal_all = subset(LvsD_ALL_er, logFC < -0.25 & adj.P.Val <= 0.3)$label

lobular_MDLC1 = subset(LvsD_MDLC1, logFC > 0.25 & adj.P.Val <= 0.3)$label
ductal_MDLC1 = subset(LvsD_MDLC1, logFC < -0.25 & adj.P.Val <= 0.3)$label

lobular_MDLC2 = subset(LvsD_MDLC2, logFC > 0.25 & adj.P.Val <= 0.3)$label
ductal_MDLC2 = subset(LvsD_MDLC2, logFC < -0.25 & adj.P.Val <= 0.3)$label

lobular_MDLC3 = subset(LvsD_MDLC3, logFC > 0.25 & adj.P.Val <= 0.3)$label
ductal_MDLC3 = subset(LvsD_MDLC3, logFC < -0.25 & adj.P.Val <= 0.3)$label

library("ggvenn")

# use list as input 
HL <-list(lobular_all = lobular_all, lobular_MDLC1 = lobular_MDLC1,
         lobular_MDLC2 = lobular_MDLC2, lobular_MDLC3 = lobular_MDLC3)

# create customised venn diagram
SupFig2A_venn_l = ggvenn::ggvenn(HL[c(1,2)],show_elements=F, stroke_color="Red", stroke_linetype="solid")
SupFig2B_venn_l = ggvenn::ggvenn(HL[c(1,3)],show_elements=F,stroke_color="Red", stroke_linetype="solid")
SupFig2C_venn_l = ggvenn::ggvenn(HL[c(1,4)],show_elements=F,stroke_color="Red", stroke_linetype="solid")


HD <-list(ductal_all = ductal_all, ductal_MDLC1 = ductal_MDLC1,
         ductal_MDLC2 = ductal_MDLC2, ductal_MDLC3 = ductal_MDLC3)

SupFig2A_venn_d = ggvenn::ggvenn(HD[c(1,2)],show_elements=F, stroke_color="Red", stroke_linetype="solid")
SupFig2B_venn_d = ggvenn::ggvenn(HD[c(1,3)],show_elements=F,stroke_color="Red", stroke_linetype="solid")
SupFig2C_venn_d = ggvenn::ggvenn(HD[c(1,4)],show_elements=F,stroke_color="Red", stroke_linetype="solid")
```

### Sup Figure 2 - Patient Wise DEG Volcano Plots

``` r
SupFig2A
```

![](2023-09-28---Main_files/figure-markdown_github/unnamed-chunk-24-1.png)

``` r
# ggsave(plot = SupFig2A, file = "../outputs/Sup Fig.2A - MDLC1 DEG Volcano Plot.png", bg = "white", dpi = 1000, width=6, height=7)

SupFig2B
```

![](2023-09-28---Main_files/figure-markdown_github/unnamed-chunk-24-2.png)

``` r
# ggsave(plot = SupFig2B, file = "../outputs/Sup Fig.2B - MDLC2 DEG Volcano Plot.png",  bg = "white", dpi = 1000, width=6, height=7)

SupFig2C
```

![](2023-09-28---Main_files/figure-markdown_github/unnamed-chunk-24-3.png)

``` r
# ggsave(plot = SupFig2C, file = "../outputs/Sup Fig.2C - MDLC3 DEG Volcano Plot.png",  bg = "white", dpi = 1000, width=6, height=7)
```

### Sup Figure 2 - Overlap of Patient DEGs with Pan-Patient DEGs - lobular

``` r
SupFig2A_venn_l
```

![](2023-09-28---Main_files/figure-markdown_github/unnamed-chunk-25-1.png)

``` r
# ggsave(plot = SupFig2A_venn_l, file = "../outputs/Sup Fig.2A - MDLC1 Lobular DEG Overlap.png", bg = "white", dpi = 1000, width=6, height=7)

SupFig2B_venn_l
```

![](2023-09-28---Main_files/figure-markdown_github/unnamed-chunk-25-2.png)

``` r
# ggsave(plot = SupFig2B_venn_l, file = "../outputs/Sup Fig.2B - MDLC2 Lobular DEG Overlap.png", bg = "white", dpi = 1000, width=6, height=7)

SupFig2C_venn_l
```

![](2023-09-28---Main_files/figure-markdown_github/unnamed-chunk-25-3.png)

``` r
# ggsave(plot = SupFig2C_venn_l, file = "../outputs/Sup Fig.2C - MDLC3 Lobular DEG Overlap.png", bg = "white", dpi = 1000, width=6, height=7)
```

### Sup Figure 2 - Overlap of Patient DEGs with Pan-Patient DEGs - ductal

``` r
SupFig2A_venn_d
```

![](2023-09-28---Main_files/figure-markdown_github/unnamed-chunk-26-1.png)

``` r
# ggsave(plot = SupFig2A_venn_d, file = "../outputs/Sup Fig.2A - MDLC1 Ductal DEG Overlap.png",  bg = "white", dpi = 1000, width=6, height=7)

SupFig2B_venn_d
```

![](2023-09-28---Main_files/figure-markdown_github/unnamed-chunk-26-2.png)

``` r
# ggsave(plot = SupFig2B_venn_d, file = "../outputs/Sup Fig.2B - MDLC2 Ductal DEG Overlap.png", bg = "white", dpi = 1000, width=6, height=7)

SupFig2C_venn_d
```

![](2023-09-28---Main_files/figure-markdown_github/unnamed-chunk-26-3.png)

``` r
# ggsave(plot = SupFig2C_venn_d, file = "../outputs/Sup Fig.2C - MDLC3 Ductal DEG Overlap.png",  bg = "white", dpi = 1000, width=6, height=7)
```

## III. DEG OVERLAP WITH PURE ILC VS IDC

``` r
source("./essentials/h.DEG_overlap_with_TCGA.R")
```

### Figure 2C - DEG OVERLAP WITH PURE ILC VS IDC

``` r
intersect(lobular_all, TCGA_ILC)
```

    ## [1] "KLK11"   "BTG2"    "SHROOM1" "PDK4"    "KLK10"   "NRAP"

``` r
lob_int = seqsetvis::ssvFeatureVenn(list(`DSP Lobular` = lobular_all, `TCGA ILC` = TCGA_ILC),
                                    circle_colors = c("red", "pink"), counts_txt_size = 0) + 
  theme(legend.text = element_text(size = 12, face = "bold"))

lob_int
```

![](2023-09-28---Main_files/figure-markdown_github/unnamed-chunk-28-1.png)

``` r
# ggsave(file = "../outputs/Fig.2C - mDLC lobular vs TCGA ILC DEG Overlap.png", dpi = 800,  bg = "white", plot = lob_int, width=4, height=4)



intersect(ductal_all, TCGA_NST)
```

    ## [1] "DCD"      "CDH1"     "SERPINA1" "SCD"

``` r
duc_int = seqsetvis::ssvFeatureVenn(list(`DSP Ductal` = ductal_all, `TCGA IDC` = TCGA_NST),
                                    circle_colors = c("blue","cyan"),counts_txt_size = 0) + 
  theme(legend.text = element_text(size = 12, face = "bold"))

duc_int
```

![](2023-09-28---Main_files/figure-markdown_github/unnamed-chunk-28-2.png)

``` r
# ggsave(file = "../outputs/Fig.2C - mDLC ductal vs TCGA NST DEG Overlap.png", dpi = 800, bg = "white", plot = duc_int, width=4, height=4)
```

## - - - - - - - - - - - - - - - - - - - - - - - - - - -

# D. PATHWAY ANALYSIS

##I. LOAD MSIG GSETS

``` r
# library(hypeR)
# KEGG     <- msigdb_gsets(species="Homo sapiens", category="C2", subcategory="CP:KEGG", clean = T)
# REACTOME <- msigdb_gsets(species="Homo sapiens", category="C2", subcategory="CP:REACTOME", clean = T)
# GO_BP <- hypeR::msigdb_gsets(species="Homo sapiens", category="C5", subcategory="BP",clean = T)
# GO_MF <- hypeR::msigdb_gsets(species="Homo sapiens", category="C5", subcategory="MF",clean = T)
# GO_CC <- hypeR::msigdb_gsets(species="Homo sapiens", category="C5", subcategory="CC",clean = T)
# Hallmark <- hypeR::msigdb_gsets(species="Homo sapiens", category="H",clean = T)
# 
# save(list = c("KEGG", "REACTOME", "GO_BP", "GO_MF", "GO_CC", "Hallmark"), file = "../inputs/hypeR_msigDB.Rdata")

load("../inputs/hypeR_msigDB.Rdata")
```

## II. HALLMARK SIGNATURE ANALYSIS

### Figure 2D - Differential Hallmark Pathways

``` r
library(ComplexHeatmap)
cols = c("#003d30", "#7dc8bf","white", "#ca9446", "#543005")
colFun4 = circlize::colorRamp2(c(-1, -0.25, 0, 0.25, 1), cols, space = 'RGB')
column_names = colnames(sig_score)

names(pvalues)[pvalues <= 0.05 & group != "None"] -> SigPathways

set.seed(123);  gsva_heatmap = Heatmap(t((t(sig_score[SigPathways,]))), 
                                       col = colFun4, show_column_names = F, show_row_names = T,  
                                       name = "GSVA", cluster_column_slices = F,cluster_row_slices = F,
                                       border = T, row_gap = unit(0,"mm"), column_gap = unit(0,'mm'), 
                                       column_labels = ROI_annot[colnames(sig_score), "SegmentLabel"],
                                       row_names_gp = gpar(fontsize = 8, fontface = "bold"),
                                       column_split = ROI_annot[colnames(sig_score), "SegmentLabel"],
                                       top_annotation = bot_annot2(colnames(sig_score)))
                                       # top_annotation = bot_annot2(colnames(sig_score)))
gsva_heatmap
```

![](2023-09-28---Main_files/figure-markdown_github/unnamed-chunk-31-1.png)

``` r
png("../outputs/Fig.2D - Differential Hallmark Pathways.png", res = 1000, 
    width=8, height=6, units = "in")
draw(gsva_heatmap, ht_gap = unit(0.1, "cm"))
dev.off()
```

    ## png 
    ##   2

## II. RECURRING BIOLOGICAL THEMES

##### Figure 2E - Biological Signatures (left panel - jaccard similarity heatmap)

``` r
set.seed(123); bioThemes_HT = Heatmap(as.matrix(jacMat_sub), name = "Jaccard Similarity",
                                      row_split = biological_clusters, column_split = biological_clusters, 
                                      row_title_rot = 0, column_title = " ",
                                      cluster_columns = T, cluster_rows = T, cluster_row_slices = F,
                                      cluster_column_slices = F, 
                                      col = circlize::colorRamp2(c(0,0.5,1), 
                                                                 colors = c("white","#ca9446", "#543005"),
                                                                 space = 'RGB'),
                                      border = T, row_gap = unit(0,"mm"), column_gap = unit(0,'mm'), 
                                      # right_annotation = bg_pct_annots(colnames(jacMat_sub)),
                                      border_gp = gpar(col = "gray", lty = "dotted",alpha = 0.5), use_raster = F,
                                      top_annotation = enrichedInAnnot(colnames(jacMat_sub), label_pval = F),
                                      row_names_gp = gpar(fontsize = 8, fontface = "bold"),
                                      left_annotation = enrichedInAnnot(colnames(jacMat_sub), 
                                                                        label_pval = F, which = "row"),
                                      show_row_names = F, show_column_names = F)

bioThemes_HT
```

![](2023-09-28---Main_files/figure-markdown_github/unnamed-chunk-33-1.png)

``` r
png("../outputs/Fig.2E - Recurring Biological Themes (left panel).png", res = 1000, 
    width=8, height=6, units = "in")
draw(bioThemes_HT, ht_gap = unit(0.1, "cm"))
dev.off()
```

    ## png 
    ##   2

##### Figure 2E - Biological Signatures (right panel - mean score)

``` r
cols = c("#003d30", "#7dc8bf","white", "#ca9446", "#543005")
colFun4 = circlize::colorRamp2(c(-1, -0.25, 0, 0.25, 1), cols, space = 'RGB')

bot_annot2 = function(sampleIDs){
  HeatmapAnnotation(annotation_name_side = "right",
                    Patient = patient_map[ROI_annot[sampleIDs, "SampleID"]],
                    Histology = ROI_annot[sampleIDs, "SegmentLabel"],
                    annotation_name_gp = gpar(fontsize = 10, fontface="bold", col="#616161"),
                    col = annot_col)}


order_ = c("Immune", "Cell Cycle", "Metabolism", "Metabolism", "Metabolism", "Metabolism", "Immune", "Immune",
           "Cell Cycle", "Protein Maturation", "Cell Cycle", "Collagen & ECM", "Immune", "Protein Maturation",
           "Phosphatidylinositol", "Oncogeneic Signaling","Oncogeneic Signaling")
order_ = factor(order_, levels = levels(biological_clusters))

set.seed(123); mean_score_HT = Heatmap(clusters_scores, col = colFun4, 
                                row_split = order_, 
                                row_title_rot = 0, cluster_row_slices = F, cluster_column_slices = F,
                                column_split = ROI_annot[colnames(clusters_scores), "SegmentLabel"],
                                #border_gp = gpar(col = "gray", lty = "dotted",alpha = 0.5), use_raster = F,
                                name = 'Mean ES', border = T, row_gap = unit(1,"mm"), column_gap = unit(0,'mm'),
                                top_annotation = bot_annot2(colnames(clusters_scores)), show_column_names = F)



mean_score_HT
```

![](2023-09-28---Main_files/figure-markdown_github/unnamed-chunk-34-1.png)

``` r
png("../outputs/Fig.2E - Recurring Biological Themes (right panel).png", res = 1000, 
    width=8, height=6, units = "in")
draw(mean_score_HT, ht_gap = unit(0.1, "cm"))
dev.off()
```

    ## png 
    ##   2

### Sup Figure 4 - Genes Associated with Differential Biological Themes

``` r
goi_d4 = unique(unlist(gset_sub[subset(biological_themes_df, 
                                       signature_clusters %in% "D4: Oxidative Stress")$pathway]))
goi_d5 = unique(unlist(gset_sub[subset(biological_themes_df, 
                                       signature_clusters %in% "D5: mTOR")$pathway]))
goi_d6 = unique(unlist(gset_sub[subset(biological_themes_df, 
                                       signature_clusters %in% "D6: Acetyltransferase")$pathway]))
goi_d7 = unique(unlist(gset_sub[subset(biological_themes_df, 
                                       signature_clusters %in% "D7: Protein Regulation")$pathway]))
goi_d8 = unique(unlist(gset_sub[subset(biological_themes_df, 
                                       signature_clusters %in% "D8: Phosphatidylinositol")$pathway]))
goi_l2 = unique(unlist(gset_sub[subset(biological_themes_df, 
                                       signature_clusters %in% "L2: TP53 & Senescence")$pathway]))
goi_l5 = unique(unlist(gset_sub[subset(biological_themes_df, 
                                       signature_clusters %in% "L5: Collagen & ECM")$pathway]))
goi_l7 = unique(unlist(gset_sub[subset(biological_themes_df, 
                                       signature_clusters %in% "L7: Collagen Metabolism")$pathway]))
goi_l9 = unique(unlist(gset_sub[subset(biological_themes_df, 
                                       signature_clusters %in% "L9: MAPK/ERK")$pathway]))

# repeat the below with all above gois to make supplementary figure 4 panels B-I
df = data.frame(scale(t(mDLC_DSP_log2[goi_l9, sampleIDs])),
                hist = ROI_annot[mDLC_DSP_log2[,sampleIDs] %>% colnames(), "SegmentLabel"],
                Patient = patient_map[ROI_annot[mDLC_DSP_log2[,sampleIDs] %>% colnames(), "SampleID"]])
df_melt = melt(df)

ggplot(df_melt, aes(hist, value, color = hist)) + geom_boxplot(alpha = 0.5, outlier.color = NA) + 
  ggbeeswarm::geom_quasirandom(alpha = 0.5, size = 4, shape = 20) + 
  myTheme + scale_color_manual("Histology", values = annot_col$Histology) + 
  facet_wrap(.~variable, ncol = 3) + 
  stat_compare_means(comparisons = list(c("Lobular","Ductal")), label = "p.signif", label.y = 1.5,
                     color = 'red', size = 6) + ylim(c(-2,2)) + ylab("Gene Expression") + xlab(" ")
```

![](2023-09-28---Main_files/figure-markdown_github/unnamed-chunk-35-1.png)

## - - - - - - - - - - - - - - - - - - - - - - - - - - -

# E. SINGLE CELL RNASEQ ANALYSIS

``` r
source("./essentials/k.GetscRNAseqTumorPopulation.R")
```

## I. IDENTIFICATION OF DUCTAL AND LOBULAR CLUSTERS

``` r
# EXPLORING TRANSCRIPTOMIC CLUSTERS, CDH1 EXP DIFFERENCES AND CELL PHASE DIFFERENCES

tsnePlot = DimPlot(MDLC1_scRNAseq_epi, reduction = 'tsne', group.by = "RNA_snn_res.0.1", 
                   cols = c("#AED581", "#8BC34A"), label = F, pt.size = 2) + myTheme_umap

cellPhase_plt = DimPlot(MDLC1_scRNAseq_epi, reduction = 'tsne', group.by = "Phase", label = F, pt.size = 2, 
        cols = c("G1"="#FFE082", "G2M"="#EF9A9A", "S"="#C5E1A5")) + myTheme_umap

CDH1_exp_plt = FeaturePlot(MDLC1_scRNAseq_epi, features = c("CDH1"), slot = "scale.data",
                           pt.size = 2, cols = c("gray", "#ca9446", "#543005"), min.cutoff = -1, max.cutoff = 2,
                           reduction = "tsne") & myTheme_umap


CDH1_df = data.frame(exp = MDLC1_scRNAseq_epi@assays$RNA@scale.data["CDH1",] %>% as.numeric(), ID = colnames(MDLC1_scRNAseq_epi))
CDH1_df$Cluster = MDLC1_scRNAseq_epi$RNA_snn_res.0.1

CDH1_exp_barplt = ggplot(CDH1_df, aes(Cluster, exp, color = Cluster)) + 
  geom_boxplot(alpha = 1, outlier.color = NA, show.legend = F) + 
  ggbeeswarm::geom_quasirandom(alpha = 0.5, size = 4, shape = 20) + 
  myTheme + scale_color_manual("Cluster", values = c("#AED581", "#8BC34A")) + xlab(" ") + ylab("CDH1 Expression") + 
  stat_compare_means(comparisons = list(c("0","1")), label = "p.signif", color = 'red', size = 6) + ylim(c(-1,4.4))
```

### Figure 4A - Transcriptomic Clusters

``` r
fig4A = ggarrange(plotlist = list(tsnePlot,cellPhase_plt,CDH1_exp_plt, CDH1_exp_barplt), nrow = 1)
fig4A
```

![](2023-09-28---Main_files/figure-markdown_github/unnamed-chunk-38-1.png)

``` r
# ggsave(plot = fig4A, file = "../outputs/Fig.4A - scRNA Transcriptomic Clusters.png",bg = "white",  dpi = 800, width=12, height=3)
```

## II. DSP SIGNATURE BASED DECONVOLUTION

``` r
# IDENTIFY DUCTAL AND LOBULAR CELLS USING DSP MDLC-1 DUCTAL VS LOBULAR DEGs
LvsD_MDLC1 = read.csv("../outputs/MDLC1_LvsD_Limma.csv")
rownames(LvsD_MDLC1) = LvsD_MDLC1$label
LvsD_MDLC1 = LvsD_MDLC1[order(LvsD_MDLC1$absLogFC, decreasing = T),]
lobular_MDLC1_res = subset(LvsD_MDLC1, logFC > 0)
ductal_MDLC1_res = subset(LvsD_MDLC1, logFC < 0)


Ductal_DSP_signature = FeatureModulePlot(MDLC1_scRNAseq_epi, head(rownames(ductal_MDLC1_res), 100), 
                                         "Ductal_Cells", cols = c("gray", "#ca9446", "#543005"), label = F, 
                                         reduction = "tsne") & myTheme_umap
  
Lobular_DSP_signature = FeatureModulePlot(MDLC1_scRNAseq_epi, head(rownames(lobular_MDLC1_res), 100), 
                                          "Lobular_Cells", cols = c("gray", "#ca9446", "#543005"), 
                                          label = F, reduction = "tsne")
```

### Figure 4B - Deconvuluted Ductal and Lobular Cells Using DSP

``` r
fig4B = ggarrange(plotlist = list(Ductal_DSP_signature, Lobular_DSP_signature), nrow = 1)
fig4B
```

![](2023-09-28---Main_files/figure-markdown_github/unnamed-chunk-40-1.png)

``` r
# ggsave(plot = fig4B, file = "../outputs/Fig.4B - scRNA Ductal and Lobular Populations - DSP Deconvulution.png", bg = "white", dpi = 800, width=8, height=5)
```

## III. LOBULAR VS DUCTAL MARKER GENES

``` r
# # FIND DIFFERENTIALLY EXPRESSED GENES - MARKER GENES
# MDLC1_scRNAseq_epi_markers <- FindMarkers(MDLC1_scRNAseq_epi, ident.1 = "Lobular", ident.2 = "Ductal", 
#                                     group.by = "CellTypes_new", only.pos = F, logfc.threshold = 0.25)
# 
# MDLC1_scRNAseq_epi_markers$log10Padj = -log10(MDLC1_scRNAseq_epi_markers$p_val_adj)
# MDLC1_scRNAseq_epi_markers$absLogFC = abs(MDLC1_scRNAseq_epi_markers$avg_log2FC)
# MDLC1_scRNAseq_epi_markers$label = rownames(MDLC1_scRNAseq_epi_markers)
# max(MDLC1_scRNAseq_epi_markers$p_val_adj)
# 
# write.table(MDLC1_scRNAseq_epi_markers, file = "../out/outputs/MDLC1_scRNAseq_DEG_LvsD.tsv",
#             quote = F, row.names = F, sep = "\t")
```

``` r
MDLC1_scRNAseq_epi_markers = read.delim("../outputs/MDLC1_scRNAseq_DEG_LvsD.tsv")
rownames(MDLC1_scRNAseq_epi_markers) = MDLC1_scRNAseq_epi_markers$label

fc_cut = 0.5
pcut = 0.05
topN = 10

MDLC1_scRNAseq_markers_sig = subset(MDLC1_scRNAseq_epi_markers, p_val_adj <= pcut)
table(MDLC1_scRNAseq_markers_sig$avg_log2FC > 0)
```

    ## 
    ## FALSE  TRUE 
    ##   445   275

``` r
MDLC1_scRNAseq_markers_sig$label = rownames(MDLC1_scRNAseq_markers_sig)

sig_res = subset(MDLC1_scRNAseq_markers_sig, abs(avg_log2FC) > fc_cut & p_val_adj <= pcut)

topUP = subset(MDLC1_scRNAseq_markers_sig, avg_log2FC > fc_cut & p_val_adj <= pcut)
dim(topUP)
```

    ## [1] 145   8

``` r
topUP = topUP[order(topUP$absLogFC, topUP$log10Padj, decreasing = T),] %>% head(topN)

topDn = subset(MDLC1_scRNAseq_markers_sig, avg_log2FC < -fc_cut & p_val_adj <= pcut)
dim(topDn)
```

    ## [1] 120   8

``` r
topDn = topDn[order(topDn$absLogFC, topDn$log10Padj, decreasing = T),] %>% head(topN)


scRNAseq_DEG_plot = ggplot(MDLC1_scRNAseq_epi_markers, aes(x = avg_log2FC, y = log10Padj, 
                                                           color = factor(avg_log2FC>0))) + 
  myTheme + geom_point(alpha=0.8, aes(size = log10Padj)) + 
  scale_color_manual("Log2 FC",values = c("blue", "red"), labels = c("Ductal", "Lobular")) + 
  guides(color = guide_legend(override.aes = list(size=6)), size="none") + xlab("Log2 FC") + 
  geom_hline(yintercept = -log10(0.05), linetype="dashed") +
  geom_vline(xintercept = 0.2, linetype="dashed") + 
  geom_vline(xintercept = -0.2, linetype="dashed") + 
  ggrepel::geom_label_repel(data = rbind(topUP, topDn), max.overlaps = 100, fontface = "bold", color = "black",
                            aes(label=label), show.legend = F, box.padding = 0.5) + 
  ggplot2::annotate(geom = "text", x = 2, y = 0, size=4, fontface="bold",color="darkgray", 
                    label = paste0("Lobular (",as.numeric(table(sig_res$avg_log2FC>0)[2]),")")) + 
  ggplot2::annotate(geom = "text", x = -2.5, y = 0, size=4, fontface="bold",color="darkgray",
                    label= paste0("Ductal (",as.numeric(table(sig_res$avg_log2FC>0)[1]),")")) + 
  theme(legend.position = "Dn")
```

### Sup Figure 6B - scRNAseq Ductal vs Lobular DEGs Volcano Plot

``` r
scRNAseq_DEG_plot
```

![](2023-09-28---Main_files/figure-markdown_github/unnamed-chunk-43-1.png)

``` r
# ggsave(plot = scRNAseq_DEG_plot,  file = "../outputs/Sup Fig.6B - scRNAseq Ductal vs Lobular DEG Volcano Plot.png", bg = "white", dpi = 800, width=8, height=5)
```

## IV. DEG OVERLAP: scRNAseq vs DSP

``` r
MDLC1_scRNAseq_epi_markers$group = ifelse(MDLC1_scRNAseq_epi_markers$avg_log2FC > 0, "Lobular", "Ductal")
LvsD_MDLC1_sig$group = ifelse(LvsD_MDLC1_sig$logFC > 0, "Lobular","Ductal")
table(LvsD_MDLC1_sig$group)
```

    ## 
    ##  Ductal Lobular 
    ##     522     600

``` r
library(ComplexHeatmap)
features = intersect(LvsD_MDLC1_sig$label, MDLC1_scRNAseq_epi_markers$label)
mat = MDLC1_scRNAseq_epi@assays$RNA@scale.data[features, ]

# set.seed(123); Heatmap(t(scale(t(mat))), show_column_names = F,show_row_names = F, cluster_row_slices = F, 
#                        name = "mRNA", column_split = MDLC1_scRNAseq_epi$CellTypes_new, show_row_dend = F, 
#                        show_column_dend = F, col = colFun, 
#                        top_annotation = HeatmapAnnotation(Histology = MDLC1_scRNAseq_epi$CellTypes_new,
#                                                           col = list(Histology = annot_col$Histology)))


scRNA_lobular = subset(MDLC1_scRNAseq_markers_sig, avg_log2FC > 0)$label
scRNA_ductal = subset(MDLC1_scRNAseq_markers_sig, avg_log2FC < 0)$label

DSP_lobular = subset(LvsD_MDLC1, logFC > 0.25 & adj.P.Val <= 0.3)$label
DSP_ductal = subset(LvsD_MDLC1, logFC < -0.25 & adj.P.Val <= 0.3)$label


intersect(DSP_lobular, scRNA_lobular)
```

    ##  [1] "AGR2"     "MBOAT2"   "GFRA1"    "TMEM150C" "CBX6"     "C3orf52" 
    ##  [7] "MAST4"    "BMP2K"    "KLHDC7A"  "BEX3"     "ISOC1"    "ANKH"    
    ## [13] "CLSTN2"   "SEC14L2"  "PRKAA2"   "SASH1"    "MSX2"     "IL1R1"   
    ## [19] "GALNT6"   "MAP3K1"   "ZSCAN18"  "ANXA5"    "DMKN"     "AR"      
    ## [25] "PDK4"     "HPX"      "TM4SF1"   "QDPR"     "HK2"      "DGCR6L"  
    ## [31] "ALCAM"    "SLC40A1"  "KIF13B"   "THRB"     "SH3BGRL"  "GSTP1"   
    ## [37] "LRRFIP1"  "HMGCS2"   "GLCE"     "ALDH3B2"  "SOCS2"    "HINT1"   
    ## [43] "COMT"     "IGSF21"   "TNFSF10"  "RTN4"     "ZFP36L2"  "TSPAN13" 
    ## [49] "FOXP1"    "PRR15"    "PERP"     "GOLM1"    "CD9"      "CYBRD1"  
    ## [55] "LDLRAD4"  "IL6ST"    "PRKAA1"   "SHROOM1"  "CYP4X1"   "SNRPN"   
    ## [61] "SAR1A"    "TPBG"     "CSNK1A1"  "ACSL3"    "SLC1A1"   "C15orf48"
    ## [67] "XBP1"     "SNRPD3"   "DDX17"    "NFIA"     "PNRC2"    "PDIA3"   
    ## [73] "SRSF5"    "SCARB2"   "MACF1"    "SCUBE2"   "PTMS"     "DHCR24"  
    ## [79] "LIMA1"    "MAGED2"   "CMTM6"    "C5orf15"  "CD63"     "CAMK2N1" 
    ## [85] "BTG2"     "SUB1"

``` r
DSP_scRNASeq_overlap_lobular = seqsetvis::ssvFeatureVenn(list(`DSP` = DSP_lobular, 
                                                              `scRNAseq` = scRNA_lobular), circle_colors = c("red", "pink")) + theme(legend.text = element_text(size = 12, face = "bold"))



intersect(DSP_ductal, scRNA_ductal)
```

    ##  [1] "AMIGO2"    "WNT11"     "HPGD"      "SLC28A3"   "S100G"     "PYDC1"    
    ##  [7] "PCDH7"     "HOXA9"     "NTRK2"     "SERPING1"  "WFDC2"     "SCRN1"    
    ## [13] "CNFN"      "SYT7"      "CHI3L2"    "CD109"     "RGL3"      "SP5"      
    ## [19] "SCARA3"    "CITED4"    "CD44"      "SLC7A2"    "CEACAM1"   "IRX1"     
    ## [25] "SKAP2"     "VPS50"     "ITGB8"     "HEPACAM2"  "S100A6"    "SLC6A4"   
    ## [31] "SULF2"     "TINAGL1"   "AK5"       "BCKDK"     "DCLK1"     "GSC"      
    ## [37] "EPHA3"     "RAB31"     "ADGRF1"    "TMC5"      "S100A9"    "COX7A1"   
    ## [43] "PIGR"      "PROM1"     "ITPRIPL2"  "CXCL16"    "ANKRD36C"  "IL27RA"   
    ## [49] "CBR3"      "MAN1A1"    "SLC16A1"   "ARPC1B"    "RAB11FIP1" "PCED1B"   
    ## [55] "METRN"     "ST6GAL1"   "SULF1"     "TMEM37"    "PLP2"      "GTF2H2C"  
    ## [61] "ERBB4"     "PGGHG"     "RERG"      "YPEL3"     "PYCARD"    "TRPS1"    
    ## [67] "NAT1"      "CLU"       "PYCR2"     "CRABP2"    "ZG16B"     "SELENOH"  
    ## [73] "S100A13"   "COX6C"     "RAC1"      "ZYX"       "CNDP2"

``` r
DSP_scRNASeq_overlap_ductal = seqsetvis::ssvFeatureVenn(list(`DSP` = DSP_ductal, 
                                                             `scRNAseq` = scRNA_ductal), circle_colors = c("blue","cyan")) +  theme(legend.text = element_text(size = 12, face = "bold"))
```

### Sup Figure 6C - DSP vs scRNAseq DEG Overlap

``` r
DSP_scRNASeq_overlap_lobular
```

![](2023-09-28---Main_files/figure-markdown_github/unnamed-chunk-45-1.png)

``` r
# ggsave(file = "../outputs/Sup Fig.6C - DSP lobular vs scRNAseq Lobular DEG Overlap.png", dpi = 800,  bg = "white", plot = DSP_scRNASeq_overlap_lobular, width=4, height=4)

DSP_scRNASeq_overlap_ductal
```

![](2023-09-28---Main_files/figure-markdown_github/unnamed-chunk-45-2.png)

``` r
# ggsave(file = "../outputs/Sup Fig.6C - DSP ductal vs scRNAseq ductal DEG Overlap.png", dpi = 800,  bg = "white", plot = DSP_scRNASeq_overlap_ductal, width=4, height=4)
```

## V. SCRNA TUMOR HETEROGENEITY

### a. IDENTIFICATION OF ROBUST CLUSTERS

``` r
library(clustree)

clustree(MDLC1_scRNAseq_epi, prefix = "RNA_snn_res.")
```

![](2023-09-28---Main_files/figure-markdown_github/unnamed-chunk-46-1.png)

### b. DEFINING SUBCLUSTERS

``` r
MDLC1_scRNAseq_epi$seurat_clusters = as.character(MDLC1_scRNAseq_epi$RNA_snn_res.1)
MDLC1_scRNAseq_epi$RNA_snn_res.1_new = ifelse(MDLC1_scRNAseq_epi$seurat_clusters == "5", "D1", 
                                    ifelse(MDLC1_scRNAseq_epi$seurat_clusters == "2", "D2", 
                                           ifelse(MDLC1_scRNAseq_epi$seurat_clusters == "0", "D3", 
                                                  ifelse(MDLC1_scRNAseq_epi$seurat_clusters == "4", "L1", 
                                                         ifelse(MDLC1_scRNAseq_epi$seurat_clusters == "3", "L2", 
                                                                ifelse(MDLC1_scRNAseq_epi$seurat_clusters == "1", "L3", MDLC1_scRNAseq_epi$seurat_clusters
                                                                ))))))

cluster_cols = c("D1" = "#90CAF9", "D2" = "#2196F3", "D3" = "#0D47A1", "L1" = "#EF9A9A", "L2" = "#F44336", "L3" = "#B71C1C")


MDLC1_scRNAseq_epi$RNA_snn_res.1_new = factor(MDLC1_scRNAseq_epi$RNA_snn_res.1_new, levels = c("D1","D2","D3","L1","L2","L3"), labels = c("D1","D2","D3","L1","L2","L3"))

subCluster_Tnse = DimPlot(MDLC1_scRNAseq_epi, reduction = 'tsne', label = F, 
        group.by = "RNA_snn_res.1_new", cols = cluster_cols, pt.size = 2) + myTheme_umap +
  theme(legend.direction="horizontal")
```

### c. COMPUTE MARKER GENES

``` r
source("./essentials/l.subcluster_markerGenes.R")
```

### Figure 4C - All Marker Genes

``` r
subCluster_Tnse
```

![](2023-09-28---Main_files/figure-markdown_github/unnamed-chunk-49-1.png)

``` r
# ggsave(plot = subCluster_Tnse,  file = "../outputs/Fig.4C - scRNA subcluster tsne.png", dpi = 800, width=4, height=2)


draw(Markers_Ht, ht_gap = unit(0.1, "cm"))
```

![](2023-09-28---Main_files/figure-markdown_github/unnamed-chunk-49-2.png)

``` r
png("../outputs/Fig.4C - scRNA subcluster all marker genes.png", res = 200, 
    width=3, height=2, units = "in")
draw(Markers_Ht, ht_gap = unit(0.1, "cm"))
dev.off()
```

    ## png 
    ##   2

### Sup Figure 6D - Top 5 Marker Genes

``` r
draw(topMarkers_Ht, ht_gap = unit(0.1, "cm"))
```

![](2023-09-28---Main_files/figure-markdown_github/unnamed-chunk-50-1.png)

``` r
png("../outputs/Sup Fig.6D - scRNA subcluster top marker genes.png", res = 1000, 
    width=5, height=4, units = "in")
draw(topMarkers_Ht, ht_gap = unit(0.1, "cm"))
dev.off()
```

    ## png 
    ##   2

## VI. MAPPING SUBCLUSTERS TO SPATIAL ROIs

``` r
sampleIDs = subset(ROI_annot[colnames(mDLC_DSP_log2),], SampleID %in% "TP19") %>% rownames()  # MDLC1/TP19


topN = subset(subcluster_markers_sig, gene %in% rownames(mDLC_DSP_log2))
subcluster_markers = split(topN$gene, topN$cluster)


library(GSVA)
subcluster_ROI_score = gsva(expr = as.matrix(mDLC_DSP_log2[, sampleIDs]), 
                            gset.idx.list = subcluster_markers,
                            method = "ssgsea", kcdf="Gaussian")
```

    ## Estimating ssGSEA scores for 6 gene sets.
    ## [1] "Calculating ranks..."
    ## [1] "Calculating absolute values from ranks..."
    ##   |                                                                              |                                                                      |   0%  |                                                                              |============                                                          |  17%  |                                                                              |=======================                                               |  33%  |                                                                              |===================================                                   |  50%  |                                                                              |===============================================                       |  67%  |                                                                              |==========================================================            |  83%  |                                                                              |======================================================================| 100%
    ## 
    ## [1] "Normalizing..."

``` r
mat = t(scale(t(subcluster_ROI_score)))
colnames(mat) = gsub(x = colnames(mat), pattern = "2460.2459.TP19.[|]|[.]", replacement = "")

idx = apply(mat, MARGIN = 2, which.max)
a=rownames(mat)[idx]
names(a) = colnames(mat)
a
```

    ## 009|Lobular  005|Ductal  006|Ductal 008|Lobular  007|Ductal 010|Lobular 
    ##        "L2"        "D3"        "D2"        "L3"        "D1"        "L1"

``` r
set.seed(123); Fig4D = Heatmap(t(mat), name = "GSVA", col = colFun,border = T, show_row_dend = F, 
                               show_column_dend = F)
```

### Figure 4D - Subcluster ROI mapping

``` r
draw(Fig4D, ht_gap = unit(0.1, "cm"))
```

![](2023-09-28---Main_files/figure-markdown_github/unnamed-chunk-52-1.png)

``` r
png("../outputs/Fig.4D - Subcluster_ROI Mapping.png", res = 1000, 
    width=5, height=3.5, units = "in")
draw(Fig4D, ht_gap = unit(0.1, "cm"))
dev.off()
```

    ## png 
    ##   2

## VII. PATHWAY ANALYSIS

### a. VALIDATING DSP HALLMARK PATHWAYS

### Sup Figure 7A - Hallmark Pathways scRNAseq

``` r
draw(hallmark_gsva_heatmap_scRNA, heatmap_legend_side = "left")
```

![](2023-09-28---Main_files/figure-markdown_github/unnamed-chunk-54-1.png)

``` r
png("../outputs/Sup Fig.7A - Hallmark Pathways scRNAseq.png", res = 1000, 
    width=7, height=4, units = "in")
draw(hallmark_gsva_heatmap_scRNA, heatmap_legend_side = "left")
dev.off()
```

    ## png 
    ##   2

### Sup Figure 7B - Subcluster Hallmark Pathway Scores

``` r
library(dplyr)

sig_score_df$Cluster = as.character(MDLC1_scRNAseq_epi@active.ident[sig_score_df$Cell])

sig_score_df_mean = subset(sig_score_df, Pathway %in% myDSP_sig)[,c("Pathway","Cluster","Score")] %>% 
  group_by(Cluster, Pathway) %>% 
  summarise_at(.vars = "Score",.funs = median) %>% as.data.frame()

sig_score_df_mean = subset(sig_score_df_mean, Pathway %in% myDSP_sig)
sig_score_df_mean$Pathway = gsub(x = sig_score_df_mean$Pathway, pattern = "Hallmark - | KEGG -", replacement = "")

hallmark_subcluster_barplot_scRNA = ggplot(sig_score_df_mean, aes(x = Cluster, y = Score, fill = Cluster)) + 
  geom_bar(stat = "identity") + scale_fill_manual(values = cluster_cols)+
  myTheme+facet_wrap(.~reorder(Pathway, Cluster), ncol = 2)


hallmark_subcluster_barplot_scRNA
```

![](2023-09-28---Main_files/figure-markdown_github/unnamed-chunk-55-1.png)

``` r
# ggsave(plot = hallmark_subcluster_barplot_scRNA, file = "../outputs/Sup Fig.7B - Subcluster Hallmark Pathway Scores.png", dpi = 800, width=6, height=8)
```

### b. VALIDATING DSP BIOLOGICAL THEMES

### Sup Figure 7C - Biological Theme Signatures

``` r
draw(biologicalTheme_heatmap_scRNA, heatmap_legend_side = "left")
```

![](2023-09-28---Main_files/figure-markdown_github/unnamed-chunk-57-1.png)

``` r
png("../outputs/Sup Fig.7C - Biological Themes scRNAseq.png", res = 1000, 
    width=9, height=6, units = "in")
draw(biologicalTheme_heatmap_scRNA, heatmap_legend_side = "left")
dev.off()
```

    ## png 
    ##   2

### Sup Figure 7D - Biological Theme - Mean Signature Score

``` r
draw(biologicalTheme_heatmap_scRNA_mean, heatmap_legend_side = "left")
```

![](2023-09-28---Main_files/figure-markdown_github/unnamed-chunk-58-1.png)

``` r
png("../outputs/Sup Fig.7D - Biological Themes scRNAseq Mean Scores.png", res = 1000, 
    width=6, height=4, units = "in")
draw(biologicalTheme_heatmap_scRNA_mean, heatmap_legend_side = "left")
dev.off()
```

    ## png 
    ##   2

## - - - - - - - - - - - - - - - - - - - - - - - - - - -

# SESSION INFO

``` r
sessionInfo()
```

    ## R version 4.3.1 (2023-06-16 ucrt)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 11 x64 (build 22621)
    ## 
    ## Matrix products: default
    ## 
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_United States.utf8 
    ## [2] LC_CTYPE=English_United States.utf8   
    ## [3] LC_MONETARY=English_United States.utf8
    ## [4] LC_NUMERIC=C                          
    ## [5] LC_TIME=English_United States.utf8    
    ## 
    ## time zone: America/New_York
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats4    grid      stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] clustree_0.5.0              ggraph_2.1.0               
    ##  [3] SeuratObject_4.1.4          Seurat_4.4.0               
    ##  [5] GSVA_1.48.3                 DESeq2_1.40.2              
    ##  [7] SummarizedExperiment_1.30.2 MatrixGenerics_1.12.3      
    ##  [9] matrixStats_1.0.0           GenomicRanges_1.52.0       
    ## [11] GenomeInfoDb_1.36.2         IRanges_2.34.1             
    ## [13] S4Vectors_0.38.1            ggvenn_0.1.10              
    ## [15] genefu_2.32.0               AIMS_1.32.0                
    ## [17] Biobase_2.60.0              BiocGenerics_0.46.0        
    ## [19] e1071_1.7-13                iC10_1.5                   
    ## [21] iC10TrainingData_1.3.1      impute_1.74.1              
    ## [23] pamr_1.56.1                 cluster_2.1.4              
    ## [25] biomaRt_2.56.1              survcomp_1.50.0            
    ## [27] prodlim_2023.08.28          survival_3.5-5             
    ## [29] ConsensusClusterPlus_1.64.0 Rtsne_0.16                 
    ## [31] limma_3.56.2                ComplexHeatmap_2.16.0      
    ## [33] reshape2_1.4.4              ggthemes_4.2.4             
    ## [35] gridExtra_2.3               ggpubr_0.6.0               
    ## [37] EnvStats_2.8.1              dplyr_1.1.2                
    ## [39] cowplot_1.1.1               ggplot2_3.4.3              
    ## [41] magrittr_2.0.3             
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] spatstat.sparse_3.0-2       bitops_1.0-7               
    ##   [3] httr_1.4.7                  webshot_0.5.5              
    ##   [5] RColorBrewer_1.1-3          doParallel_1.0.17          
    ##   [7] sctransform_0.4.0           tools_4.3.1                
    ##   [9] backports_1.4.1             utf8_1.2.3                 
    ##  [11] R6_2.5.1                    HDF5Array_1.28.1           
    ##  [13] uwot_0.1.16                 lazyeval_0.2.2             
    ##  [15] rhdf5filters_1.12.1         GetoptLong_1.0.5           
    ##  [17] sp_2.0-0                    withr_2.5.1                
    ##  [19] prettyunits_1.1.1           progressr_0.14.0           
    ##  [21] cli_3.6.1                   spatstat.explore_3.2-3     
    ##  [23] labeling_0.4.3              spatstat.data_3.0-1        
    ##  [25] ggridges_0.5.4              proxy_0.4-27               
    ##  [27] pbapply_1.7-2               Rsamtools_2.16.0           
    ##  [29] systemfonts_1.0.4           yulab.utils_0.0.8          
    ##  [31] svglite_2.1.2               parallelly_1.36.0          
    ##  [33] readxl_1.4.3                rstudioapi_0.15.0          
    ##  [35] RSQLite_2.3.1               generics_0.1.3             
    ##  [37] gridGraphics_0.5-1          shape_1.4.6                
    ##  [39] BiocIO_1.10.0               spatstat.random_3.1-6      
    ##  [41] ica_1.0-3                   car_3.1-2                  
    ##  [43] Matrix_1.6-1                ggbeeswarm_0.7.2           
    ##  [45] fansi_1.0.4                 abind_1.4-5                
    ##  [47] lifecycle_1.0.3             yaml_2.3.7                 
    ##  [49] carData_3.0-5               rhdf5_2.44.0               
    ##  [51] BiocFileCache_2.8.0         blob_1.2.4                 
    ##  [53] promises_1.2.1              seqsetvis_1.20.0           
    ##  [55] crayon_1.5.2                miniUI_0.1.1.1             
    ##  [57] lattice_0.21-8              beachmat_2.16.0            
    ##  [59] annotate_1.78.0             KEGGREST_1.40.0            
    ##  [61] magick_2.8.1                pillar_1.9.0               
    ##  [63] knitr_1.45                  rjson_0.2.21               
    ##  [65] future.apply_1.11.0         codetools_0.2-19           
    ##  [67] leiden_0.4.3                glue_1.6.2                 
    ##  [69] data.table_1.14.8           vctrs_0.6.3                
    ##  [71] png_0.1-8                   bootstrap_2019.6           
    ##  [73] cellranger_1.1.0            gtable_0.3.4               
    ##  [75] cachem_1.0.8                xfun_0.40                  
    ##  [77] mime_0.12                   S4Arrays_1.0.6             
    ##  [79] tidygraph_1.2.3             SingleCellExperiment_1.22.0
    ##  [81] iterators_1.0.14            pbmcapply_1.5.1            
    ##  [83] lava_1.7.2.1                ellipsis_0.3.2             
    ##  [85] fitdistrplus_1.1-11         ROCR_1.0-11                
    ##  [87] nlme_3.1-162                survivalROC_1.0.3.1        
    ##  [89] bit64_4.0.5                 RcppAnnoy_0.0.21           
    ##  [91] progress_1.2.2              filelock_1.0.2             
    ##  [93] UpSetR_1.4.0                irlba_2.3.5.1              
    ##  [95] vipor_0.4.5                 KernSmooth_2.23-21         
    ##  [97] colorspace_2.1-0            rmeta_3.0                  
    ##  [99] DBI_1.1.3                   tidyselect_1.2.0           
    ## [101] bit_4.0.5                   compiler_4.3.1             
    ## [103] curl_5.0.2                  rvest_1.0.3                
    ## [105] graph_1.78.0                xml2_1.3.5                 
    ## [107] plotly_4.10.2               DelayedArray_0.26.7        
    ## [109] rtracklayer_1.60.1          checkmate_2.2.0            
    ## [111] scales_1.2.1                lmtest_0.9-40              
    ## [113] rappdirs_0.3.3              goftest_1.2-3              
    ## [115] stringr_1.5.0               digest_0.6.33              
    ## [117] spatstat.utils_3.0-3        rmarkdown_2.24             
    ## [119] XVector_0.40.0              htmltools_0.5.6            
    ## [121] pkgconfig_2.0.3             sparseMatrixStats_1.12.2   
    ## [123] highr_0.10                  dbplyr_2.3.3               
    ## [125] fastmap_1.1.1               htmlwidgets_1.6.2          
    ## [127] rlang_1.1.1                 GlobalOptions_0.1.2        
    ## [129] shiny_1.7.5                 SuppDists_1.1-9.7          
    ## [131] DelayedMatrixStats_1.22.6   farver_2.1.1               
    ## [133] zoo_1.8-12                  jsonlite_1.8.7             
    ## [135] BiocParallel_1.34.2         mclust_6.0.0               
    ## [137] BiocSingular_1.16.0         RCurl_1.98-1.12            
    ## [139] GenomeInfoDbData_1.2.10     ggplotify_0.1.2            
    ## [141] patchwork_1.1.3             Rhdf5lib_1.22.1            
    ## [143] munsell_0.5.0               Rcpp_1.0.11                
    ## [145] viridis_0.6.4               reticulate_1.31            
    ## [147] stringi_1.7.12              MASS_7.3-60                
    ## [149] zlibbioc_1.46.0             plyr_1.8.8                 
    ## [151] parallel_4.3.1              listenv_0.9.0              
    ## [153] ggrepel_0.9.3               deldir_1.0-9               
    ## [155] graphlayouts_1.0.0          Biostrings_2.68.1          
    ## [157] splines_4.3.1               tensor_1.5                 
    ## [159] hms_1.1.3                   circlize_0.4.15            
    ## [161] locfit_1.5-9.8              igraph_1.5.1               
    ## [163] spatstat.geom_3.2-5         ggsignif_0.6.4             
    ## [165] ScaledMatrix_1.8.1          XML_3.99-0.14              
    ## [167] evaluate_0.21               tweenr_2.0.2               
    ## [169] httpuv_1.6.11               foreach_1.5.2              
    ## [171] polyclip_1.10-4             RANN_2.6.1                 
    ## [173] tidyr_1.3.0                 purrr_1.0.2                
    ## [175] scattermore_1.2             future_1.33.0              
    ## [177] clue_0.3-64                 ggforce_0.4.1              
    ## [179] rsvd_1.0.5                  eulerr_7.0.0               
    ## [181] broom_1.0.5                 xtable_1.8-4               
    ## [183] restfulr_0.0.15             later_1.3.1                
    ## [185] rstatix_0.7.2               viridisLite_0.4.2          
    ## [187] class_7.3-22                tibble_3.2.1               
    ## [189] memoise_2.0.1               beeswarm_0.4.0             
    ## [191] AnnotationDbi_1.62.2        GenomicAlignments_1.36.0   
    ## [193] globals_0.16.2              GSEABase_1.62.0

``` r
writeLines(capture.output(sessionInfo()), paste0("./",Sys.Date(),"-", "_Session Info.txt"))

write.csv(data.frame(x=(.packages())), file = "../inputs/required_packages.csv", row.names = F)
```
