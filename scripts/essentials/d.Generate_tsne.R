library(Rtsne)

sampleIDs = colnames(mDLC_DSP_log2)  # ALL SAMPLES

# input matrix
dim_mat = mDLC_DSP_log2[, sampleIDs]
sort(colnames(dim_mat))
dim(dim_mat)
dim_mat[1:5,1:5]

# select features
feature.vars = apply(dim_mat, MARGIN = 1, var)
feature.vars = sort(feature.vars,decreasing = T)
names(feature.vars[feature.vars>0.15]) %>% length()
topPCT = 15
select_features = names(sort(feature.vars, dec = T))[1:(length(feature.vars)*topPCT/100)]
length(select_features)

## tsne
set.seed(123); tsne_res_top15 <- Rtsne(t(dim_mat[select_features,]), perplexity = 5, 
                                              check_duplicates = FALSE) 

dim_res = tsne_res_top15$Y # pca_res$x