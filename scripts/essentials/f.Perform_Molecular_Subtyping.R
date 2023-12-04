## 2. GeneFu
library(genefu)
data("pam50")
data("pam50.robust")
data("nkis")

# pam50 genes 
pam50.genes = rownames(pam50$centroids)
pam50.genes[pam50.genes %notin% rownames(mDLC_DSP_log2)]
pam50.genes[pam50.genes %in% "CDCA1"]="NUF2"
pam50.genes[pam50.genes %in% "KNTC2"]="NDC80"
pam50.genes[pam50.genes %in% "ORC6L"]="ORC6"

sampleIDs = colnames(mDLC_DSP_log2)
pam50_tpm = t(mDLC_DSP_log2[pam50.genes, sampleIDs])

colnames(pam50_tpm) = rownames(pam50$centroids)

PAM50.subtype <- molecular.subtyping(sbt.model="pam50", data = pam50_tpm, annot = annot.nkis, do.mapping = FALSE)
sort(PAM50.subtype$subtype)

save(PAM50.subtype, file = "../outputs/PAM50.subtype.Rdata")