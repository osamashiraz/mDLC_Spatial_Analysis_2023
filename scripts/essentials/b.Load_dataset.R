# ROI Annotations
ROI_annot = readxl::read_excel("../inputs/DPS Data.xlsx", sheet = "SegmentProperties") %>% as.data.frame()
dim(ROI_annot)
head(ROI_annot)
table(ROI_annot$SegmentLabel)
table(ROI_annot$SampleID)
rownames(ROI_annot) = ROI_annot$Sample_ID



# AOIs flagged for removal: AOI 004: Lobular, AOI 016 Lobular, AOI 016 Ductal - SEE GEOMX PPT
ROI_annot$flagRemove = ifelse(grepl(x = ROI_annot$Sample_ID, pattern = "004.[|].Lobular|016"), "Yes", "No")
ROI_annot_clean = subset(ROI_annot, flagRemove %in% "No")

# table(ROI_annot$SegmentLabel, ROI_annot$SampleID)
# gridExtra::grid.table(addmargins(table(ROI_annot_clean$SegmentLabel, ROI_annot_clean$SampleID)))
# gridExtra::grid.table(ROI_annot_clean[, c("ROILabel", "SegmentLabel", "SampleID")])

# RENAME PATIENT IDS
# patient_map = c("BIOS76976", "BIOS19209", "BIOS60674")
patient_map = c("MDLC-3", "MDLC-2", "MDLC-1")
names(patient_map) = c("2460", "2459", "TP19")

ROI_annot$BIOS_ID = ifelse(ROI_annot$SampleID %in% "2459", "MDLC-2",
                           ifelse(ROI_annot$SampleID %in% "TP19", "MDLC-1",
                                  ifelse(ROI_annot$SampleID %in% "2460",
                                         "MDLC-3",ROI_annot$SampleID)))



# LOAD COUNTS
mDLC_DSP_raw = readxl::read_excel("../inputs/DPS Data.xlsx", sheet = "TargetCountMatrix") %>% as.data.frame()
rownames(mDLC_DSP_raw) = mDLC_DSP_raw$TargetName
mDLC_DSP_raw = mDLC_DSP_raw[,-1]

mDLC_DSP = readxl::read_excel("../inputs/DPS Data.xlsx", sheet = "Q3 TargetCountMatrix") %>% as.data.frame()
mDLC_DSP[1:5, 1:5]
rownames(mDLC_DSP) = mDLC_DSP$Gene
mDLC_DSP = mDLC_DSP[,-1]
mDLC_DSP_log2 = log2(mDLC_DSP+1)

selectSamples = subset(ROI_annot, flagRemove %in% "No") %>% rownames() # remove flagged
mDLC_DSP_log2 = mDLC_DSP_log2[, selectSamples]





# E-CADHERIN IF INTENSITIES
ecad_imagej = read.delim("../inputs/E-cadherin IF intensities.tsv")
