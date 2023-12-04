library(ConsensusClusterPlus)

# select features
feature.vars = apply(dim_mat, MARGIN = 1, var)
feature.vars = sort(feature.vars,decreasing = T)
names(feature.vars[feature.vars>0.15]) %>% length()
topPCT = 15
select_features = names(sort(feature.vars, dec = T))[1:(length(feature.vars)*topPCT/100)]
length(select_features)


d = as.matrix(mDLC_DSP_log2[,])
d = d[select_features,]
d = sweep(d,1, apply(d,1,median,na.rm=T))
d[1:5,1:5]
class(d)

results = ConsensusClusterPlus(d, maxK = 8, reps = 50, pItem = 0.8, pFeature = 1, 
                               clusterAlg = "hc", distance = "pearson", seed = 123)
icl = calcICL(results)


sampleIDs = colnames(mDLC_DSP_log2)  # ALL SAMPLES

consensus_mat = results[[6]]
mat = consensus_mat$consensusMatrix

colnames(mat) = names(consensus_mat$consensusClass)
rownames(mat) = colnames(mat)

cluster_cols = consensus_mat$clrs[[1]]
names(cluster_cols) =  consensus_mat$consensusClass


# switch TP19 ductal and lobular clusters & 2460 at end
order_df = data.frame(Cluster = consensus_mat$consensusClass, Histology = ROI_annot[colnames(mat), "SegmentLabel"])
table(order_df$Cluster, order_df$Histology)
order_df_Cluster = order_df$Cluster
order_df$Cluster[order_df_Cluster %in% 1] = 4
order_df$Cluster[order_df_Cluster %in% 2] = 3
order_df$Cluster[order_df_Cluster %in% 3] = 2
order_df$Cluster[order_df_Cluster %in% 4] = 1
order_df$Cluster[order_df_Cluster %in% 5] = 6
order_df$Cluster[order_df_Cluster %in% 6] = 5
order_df$IDs = paste0(order_df$Cluster, "-", order_df$Histology)
order_df$IDs = factor(order_df$IDs, labels = sort(order_df$IDs), levels = sort(order_df$IDs))


ht_opt$ROW_ANNO_PADDING = unit(0.5, "mm") # remove padding between heatmap and annotation
score_colors = circlize::colorRamp2(c(0,0.5,1), colors = c("white","#ca9446", "#543005"), space = 'RGB')

