library(reshape2)  # for melt
library(cowplot)   # for plot_grid

# Graph Q3 value vs negGeoMean of Negatives
ann_of_interest <- "region"
Stat_data <- 
  data.frame(row.names = colnames(mDLC_DSP),
             Annotation = ROI_annot[colnames(mDLC_DSP), "SegmentLabel"],
             Q3 = unlist(apply(mDLC_DSP_raw[, rownames(ROI_annot)], 2,
                               quantile, 0.75, na.rm = TRUE)),
             label = ROI_annot[colnames(mDLC_DSP), "ROILabel"],
             flag = ROI_annot[colnames(mDLC_DSP), "flagRemove"],
             `Sequencing Saturation` = ROI_annot[colnames(mDLC_DSP), "SequencingSaturation"],
             `AOI Surface Area` = ROI_annot[colnames(mDLC_DSP), "AOISurfaceArea"],
             NegProbe = as.numeric(mDLC_DSP_raw["NegProbe-WTX", rownames(ROI_annot)]))
Stat_data_m <- melt(Stat_data, measure.vars = c("Q3", "NegProbe"),
                    variable.name = "Statistic", value.name = "Value")

supFig1A <- ggplot(Stat_data,
                   aes(y = `Sequencing.Saturation`, x = `AOI.Surface.Area`, color = Annotation)) +
  geom_abline(intercept = 50, slope = 0, lty = "dashed", color = "red") +
  geom_vline(xintercept = 5000, lty = "dashed", color = "red") +
  geom_point(size=3) + guides(color = "none") + theme_bw() + 
  #ggrepel::geom_label_repel(aes(label = label), max.overlaps = 100) +
  ylim(c(0, 100)) + 
  scale_x_continuous(trans = "log2") + 
  theme(aspect.ratio = 1) +scale_color_manual(values = c("Lobular"="red","Ductal"="blue")) +
  labs(y = "Sequencing Saturation", x = "AOI Surface Area") + myTheme


supFig1B <- ggplot(Stat_data,
                   aes(x = NegProbe, y = Q3, color = Annotation)) +
  geom_abline(intercept = 0, slope = 1, lty = "dashed", color = "red") +
  geom_point(size=3) + guides(color = "none") + theme_bw() +
  scale_x_continuous(trans = "log2") + 
  scale_y_continuous(trans = "log2") + #ggrepel::geom_label_repel(aes(label = label)) +
  theme(aspect.ratio = 1) + scale_color_manual(values = c("Lobular"="red","Ductal"="blue")) +
  labs(x = "NegProbe Counts", y = "Q3 Counts") + myTheme


colfun = circlize::colorRamp2(breaks = c(0, 10, 15, 20, 35, 50), colors = c("white", "#84FFFF","#00E5FF","#00B8D4", 
                                                                            "#0277BD", "#0D47A1"))
set.seed(123); supFig1C = Heatmap(t(as.matrix(mDLC_DSP[sample(rownames(mDLC_DSP), size = 1000),])), 
                                  use_raster = T, col = colfun, show_column_names = F, name = "Q3 Counts",
                                  row_labels = ROI_annot[colnames(mDLC_DSP), "ROILabel"], 
                                  show_column_dend = F, show_row_dend = F,
                                  row_split = ROI_annot[colnames(mDLC_DSP), "SegmentLabel"],
                                  left_annotation = HeatmapAnnotation(which = "row",show_annotation_name = F,
                                                                      Annotation = ROI_annot[colnames(mDLC_DSP), "SegmentLabel"],
                                                                      col = list(Annotation = c("Lobular"="red", "Ductal" = "blue"))))

