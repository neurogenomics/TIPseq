## First batch TIP-seq analysis

## Peak files
# batch 05/01/22
tip_1_05_01_R1_ac <- ChIPseeker::readPeakFile("/Users/serachoi/Documents/EpiCompare_extra/peakfiles/TIP/S_1_05_01_22_R1.peaks.bed.stringent.bed", as="GRanges")
tip_1_05_01_R2_ac <- ChIPseeker::readPeakFile("/Users/serachoi/Documents/EpiCompare_extra/peakfiles/TIP/S_1_05_01_22_R2.peaks.bed.stringent.bed", as="GRanges")

# batch 28/01/22
tip_2_28_01_ac <- ChIPseeker::readPeakFile("/Users/serachoi/Documents/EpiCompare_extra/peakfiles/TIP/S_2_28_01_22_R1.peaks.bed.stringent.bed", as="GRanges")
tip_3_28_01_me <- ChIPseeker::readPeakFile("/Users/serachoi/Documents/EpiCompare_extra/peakfiles/TIP/S_3_28_01_22_R1.peaks.bed.stringent.bed", as="GRanges")
tip_4_28_01_ac <- ChIPseeker::readPeakFile("/Users/serachoi/Documents/EpiCompare_extra/peakfiles/TIP/S_4_28_01_22_R1.peaks.bed.stringent.bed", as="GRanges")
tip_5_28_01_me <- ChIPseeker::readPeakFile("/Users/serachoi/Documents/EpiCompare_extra/peakfiles/TIP/S_5_28_01_22_R1.peaks.bed.stringent.bed", as="GRanges")
tip_6_28_01_ac <- ChIPseeker::readPeakFile("/Users/serachoi/Documents/EpiCompare_extra/peakfiles/TIP/S_6_28_01_22_R1.peaks.bed.stringent.bed", as="GRanges")

# batch 03/02/22
tip_4_03_02_ac <- ChIPseeker::readPeakFile("/Users/serachoi/Documents/EpiCompare_extra/peakfiles/TIP/S_4_03_02_22_R1.peaks.bed.stringent.bed", as="GRanges")

# batch 18/02/22
tip_8_18_02_me <- ChIPseeker::readPeakFile("/Users/serachoi/Documents/EpiCompare_extra/peakfiles/TIP/S_8_18_02_22_R1.peaks.bed.stringent.bed", as="GRanges")
tip_9_18_02_ac <- ChIPseeker::readPeakFile("/Users/serachoi/Documents/EpiCompare_extra/peakfiles/TIP/S_9_18_02_22_R1.peaks.bed.stringent.bed", as="GRanges")
tip_10_18_02_me <- ChIPseeker::readPeakFile("/Users/serachoi/Documents/EpiCompare_extra/peakfiles/TIP/S_10_18_02_22_R1.peaks.bed.stringent.bed", as="GRanges")
tip_11_18_02_ac <- ChIPseeker::readPeakFile("/Users/serachoi/Documents/EpiCompare_extra/peakfiles/TIP/S_11_18_02_22_R1.peaks.bed.stringent.bed", as="GRanges")

## create peaklist
peaklist <- list(tip_1_05_01_R1_ac, tip_1_05_01_R2_ac,
                 tip_2_28_01_ac, tip_4_28_01_ac, tip_6_28_01_ac,
                 tip_4_03_02_ac,
                 tip_9_18_02_ac, tip_11_18_02_ac,
                 tip_3_28_01_me, tip_5_28_01_me,
                 tip_8_18_02_me, tip_10_18_02_me)

## names
my_label <- c("H3K27ac_Abcam.phase_1_05_jan_2022.S1_R1",
              "H3K27ac_Abcam.phase_1_05_jan_2022.S1_R2",
              "H3K27ac_Diagenode.phase_2_28_jan_2022.S2_R1",
              "H3K27ac_Abcam.phase_2_28_jan_2022.S4_R1",
              "H3K27ac_Abcam.phase_2_28_jan_2022.S6_R1",
              "H3K27ac_Abcam.phase_2_03_feb_2022.S4_R1",
              "H3K27ac_Abcam.phase_3_18_feb_2022.S9_R1",
              "H3K27ac_Abcam.phase_3_18_feb_2022.S11_R1",
              "H3K27me3_CellSignalling.phase_2_28_jan_2022.S3_R1",
              "H3K27me3_CellSignalling.phase_2_28_jan_2022.S5_R1",
              "H3K27me3_CellSignalling.phase_3_18_feb_2022.S8_R1",
              "H3K27me3_CellSignalling.phase_3_18_feb_2022.S10_R1")
names(peaklist) <- my_label

## reference files
encode_hg38 <- ChIPseeker::readPeakFile("/Users/serachoi/Documents/EpiCompare_extra/peakfiles/ENCODE/ENCODE_H3K27ac_hg38_ENCFF038DDS.bed", as="GRanges")
reference <- list("H3K27ac_ENCODE_reference"=encode_hg38)

## blacklist
data("hg38_blacklist")

## Run EpiCompare
EpiCompare(peakfiles = peaklist,
           genome_build = "hg38",
           blacklist = hg38_blacklist,
           reference = reference,
           picard_files = NULL,
           upset_plot = TRUE,
           stat_plot = TRUE,
           save_output = TRUE,
           chrmHMM_plot = TRUE,
           chrmHMM_annotation = "K562",
           chipseeker_plot = TRUE,
           enrichment_plot = TRUE,
           tss_plot = FALSE,
           interact = TRUE,
           output_dir = "/Users/serachoi/Desktop/TIP-seq")





