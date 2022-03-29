## First batch TIP-seq analysis

## Peak and picard files
# batch 05/01/22
tip_1_05_01_R1_ac <- ChIPseeker::readPeakFile("/Users/serachoi/Documents/EpiCompare_extra/peakfiles/TIP/S_1_05_01_22_R1.peaks.bed.stringent.bed", as="GRanges")
tip_1_05_01_R2_ac <- ChIPseeker::readPeakFile("/Users/serachoi/Documents/EpiCompare_extra/peakfiles/TIP/S_1_05_01_22_R2.peaks.bed.stringent.bed", as="GRanges")
tip_1_05_01_R1_ac_picard <- read.table("/Users/serachoi/Documents/EpiCompare_extra/picard/TIP/S_1_R1_05_01_2022.target.markdup.MarkDuplicates.metrics.txt", header = TRUE, fill = TRUE)
tip_1_05_01_R2_ac_picard <- read.table("/Users/serachoi/Documents/EpiCompare_extra/picard/TIP/S_1_R2_05_01_2022.target.markdup.MarkDuplicates.metrics.txt", header = TRUE, fill = TRUE)

# batch 28/01/22
tip_2_28_01_ac <- ChIPseeker::readPeakFile("/Users/serachoi/Documents/EpiCompare_extra/peakfiles/TIP/S_2_28_01_22_R1.peaks.bed.stringent.bed", as="GRanges")
tip_3_28_01_me <- ChIPseeker::readPeakFile("/Users/serachoi/Documents/EpiCompare_extra/peakfiles/TIP/S_3_28_01_22_R1.peaks.bed.stringent.bed", as="GRanges")
tip_4_28_01_ac <- ChIPseeker::readPeakFile("/Users/serachoi/Documents/EpiCompare_extra/peakfiles/TIP/S_4_28_01_22_R1.peaks.bed.stringent.bed", as="GRanges")
tip_5_28_01_me <- ChIPseeker::readPeakFile("/Users/serachoi/Documents/EpiCompare_extra/peakfiles/TIP/S_5_28_01_22_R1.peaks.bed.stringent.bed", as="GRanges")
tip_6_28_01_ac <- ChIPseeker::readPeakFile("/Users/serachoi/Documents/EpiCompare_extra/peakfiles/TIP/S_6_28_01_22_R1.peaks.bed.stringent.bed", as="GRanges")

tip_2_28_01_ac_picard <- read.table("/Users/serachoi/Documents/EpiCompare_extra/picard/TIP/S_2_R1.target.markdup.MarkDuplicates.metrics.txt", header = TRUE, fill = TRUE)
tip_3_28_01_me_picard <- read.table("/Users/serachoi/Documents/EpiCompare_extra/picard/TIP/S_3_R1.target.markdup.MarkDuplicates.metrics.txt", header = TRUE, fill = TRUE)
tip_4_28_01_ac_picard <- read.table("/Users/serachoi/Documents/EpiCompare_extra/picard/TIP/S_4_R1.target.markdup.MarkDuplicates.metrics.txt", header = TRUE, fill = TRUE)
tip_5_28_01_me_picard <- read.table("/Users/serachoi/Documents/EpiCompare_extra/picard/TIP/S_5_R1.target.markdup.MarkDuplicates.metrics.txt", header = TRUE, fill = TRUE)
tip_6_28_01_ac_picard <- read.table("/Users/serachoi/Documents/EpiCompare_extra/picard/TIP/S_6_R1.target.markdup.MarkDuplicates.metrics.txt", header = TRUE, fill = TRUE)

# batch 03/02/22
tip_4_03_02_ac <- ChIPseeker::readPeakFile("/Users/serachoi/Documents/EpiCompare_extra/peakfiles/TIP/S_4_03_02_22_R1.peaks.bed.stringent.bed", as="GRanges")
tip_4_03_02_ac_picard <- read.table("/Users/serachoi/Documents/EpiCompare_extra/picard/TIP/S_4_R1_03_03_2022.target.markdup.MarkDuplicates.metrics.txt", header = TRUE, fill = TRUE)

# batch 18/02/22
tip_8_18_02_me <- ChIPseeker::readPeakFile("/Users/serachoi/Documents/EpiCompare_extra/peakfiles/TIP/S_8_18_02_22_R1.peaks.bed.stringent.bed", as="GRanges")
tip_9_18_02_ac <- ChIPseeker::readPeakFile("/Users/serachoi/Documents/EpiCompare_extra/peakfiles/TIP/S_9_18_02_22_R1.peaks.bed.stringent.bed", as="GRanges")
tip_10_18_02_me <- ChIPseeker::readPeakFile("/Users/serachoi/Documents/EpiCompare_extra/peakfiles/TIP/S_10_18_02_22_R1.peaks.bed.stringent.bed", as="GRanges")
tip_11_18_02_ac <- ChIPseeker::readPeakFile("/Users/serachoi/Documents/EpiCompare_extra/peakfiles/TIP/S_11_18_02_22_R1.peaks.bed.stringent.bed", as="GRanges")

tip_8_18_02_me_picard <- read.table("/Users/serachoi/Documents/EpiCompare_extra/picard/TIP/S_8_R1.target.markdup.MarkDuplicates.metrics.txt", header = TRUE, fill = TRUE)
tip_9_18_02_ac_picard <- read.table("/Users/serachoi/Documents/EpiCompare_extra/picard/TIP/S_9_R1.target.markdup.MarkDuplicates.metrics.txt", header = TRUE, fill = TRUE)
tip_10_18_02_me_picard <- read.table("/Users/serachoi/Documents/EpiCompare_extra/picard/TIP/S_10_R1.target.markdup.MarkDuplicates.metrics.txt", header = TRUE, fill = TRUE)
tip_11_18_02_ac_picard <- read.table("/Users/serachoi/Documents/EpiCompare_extra/picard/TIP/S_11_R1.target.markdup.MarkDuplicates.metrics.txt", header = TRUE, fill = TRUE)

## create peaklist
peaklist <- list(tip_1_05_01_R1_ac, tip_1_05_01_R2_ac,
                 tip_2_28_01_ac, tip_4_28_01_ac, tip_6_28_01_ac,
                 tip_4_03_02_ac,
                 tip_9_18_02_ac, tip_11_18_02_ac,
                 tip_3_28_01_me, tip_5_28_01_me,
                 tip_8_18_02_me, tip_10_18_02_me)

## create picard list
picard_list <- list(tip_1_05_01_R1_ac_picard, tip_1_05_01_R2_ac_picard,
                    tip_2_28_01_ac_picard, tip_4_28_01_ac_picard, tip_6_28_01_ac_picard,
                    tip_4_03_02_ac_picard,
                    tip_9_18_02_ac_picard, tip_11_18_02_ac_picard,
                    tip_3_28_01_me_picard, tip_5_28_01_me_picard,
                    tip_8_18_02_me_picard, tip_10_18_02_me_picard)

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
names(picard_list) <- my_label

## reference files
encode_hg38 <- ChIPseeker::readPeakFile("/Users/serachoi/Documents/EpiCompare_extra/peakfiles/ENCODE/ENCODE_H3K27ac_hg38_ENCFF038DDS.bed", as="GRanges")
reference_H3K27ac <- list("H3K27ac_ENCODE_reference"=encode_hg38)
H3Kme3_ENCODE <- ChIPseeker::readPeakFile("/Users/serachoi/Documents/EpiCompare_extra/peakfiles/ENCODE/ENCODE_H3K27me3_ENCFF881ONN.bed", as="GRanges")
reference_H3K27me3 <- list("H3K27me3_ENCODE_reference"=H3Kme3_ENCODE)

## blacklist
data("hg38_blacklist")

## Run EpiCompare H3K27ac
EpiCompare(peakfiles = peaklist,
           genome_build = "hg38",
           blacklist = hg38_blacklist,
           reference = reference_H3K27ac,
           picard_files = picard_list,
           upset_plot = TRUE,
           stat_plot = TRUE,
           save_output = TRUE,
           chrmHMM_plot = TRUE,
           chrmHMM_annotation = "K562",
           chipseeker_plot = TRUE,
           enrichment_plot = TRUE,
           tss_plot = FALSE,
           interact = TRUE,
           output_filename = "TIPseq_H3K27ac_EpiCompare",
           output_dir = "/Users/serachoi/Documents/EpiCompare_extra/")

## Run EpiCompare H3K27me3
EpiCompare(peakfiles = peaklist,
           genome_build = "hg38",
           blacklist = hg38_blacklist,
           reference = reference_H3K27me3,
           picard_files = picard_list,
           upset_plot = TRUE,
           stat_plot = TRUE,
           save_output = TRUE,
           chrmHMM_plot = TRUE,
           chrmHMM_annotation = "K562",
           chipseeker_plot = TRUE,
           enrichment_plot = TRUE,
           tss_plot = FALSE,
           interact = TRUE,
           output_filename = "TIPseq_H3K27me3_EpiCompare",
           output_dir = "/Users/serachoi/Documents/EpiCompare_extra/")
