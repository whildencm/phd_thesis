library(riboWaltz)

# Set parameters
baseDir <- "/n/lab_storage/whipple_lab/share_root/Users/cwhilden/PROJECTS/manuscript_2026/riboseq"

gtf_path <- "/n/lab_storage/whipple_lab/share_root/Lab/Computing_Resources/Genome_mmu/gencode.vM30.annotation.gtf"
bam_dir <- file.path(baseDir,"results/alignments/sortedRpfTrs")
output_dir <- file.path(baseDir,"results/ribowaltz")

dir.create(output_dir, recursive = TRUE)

# Generate genome annotation and load bam files. Filter read lengths
annotation_dt <- create_annotation(gtfpath = gtf_path,dataSource = "gencode6", organism = "Mus musculus")
reads_list <- bamtolist(bamfolder = bam_dir, annotation = annotation_dt)
filtered_list <- length_filter(data = reads_list, length_filter_mode = "custom", length_range = 27:36)

# Define sample subsets
p0_WT <- c("EW6.RPF.Transcriptome.Sort",
"ER2.RPF.Transcriptome.Sort",
"ER7.RPF.Transcriptome.Sort",
"ER9.RPF.Transcriptome.Sort",
"EW1.RPF.Transcriptome.Sort",
"ER1.RPF.Transcriptome.Sort",
"ER8.RPF.Transcriptome.Sort",
"EW3.RPF.Transcriptome.Sort")

p0_MUT <- c("ER3.RPF.Transcriptome.Sort",
"ER4.RPF.Transcriptome.Sort",
"ER5.RPF.Transcriptome.Sort",
"EW5.RPF.Transcriptome.Sort",
"ER6.RPF.Transcriptome.Sort",
"EW2.RPF.Transcriptome.Sort",
"EW4.RPF.Transcriptome.Sort")

p7_WT <- c("EL3.RPF.Transcriptome.Sort",
"EP3.RPF.Transcriptome.Sort",
"ET2.RPF.Transcriptome.Sort",
"EL1.RPF.Transcriptome.Sort",
"EP4.RPF.Transcriptome.Sort",
"ET3.RPF.Transcriptome.Sort")

p7_MUT <- c("EL4.RPF.Transcriptome.Sort",
"EL5.RPF.Transcriptome.Sort",
"ET1.RPF.Transcriptome.Sort",
"EL2.RPF.Transcriptome.Sort",
"EP6.RPF.Transcriptome.Sort",
"ET6.RPF.Transcriptome.Sort")

# Generate heat maps of read ends near start/stop codons

pdf(file.path(output_dir,"P0_WT_PSiteOffset_Heatmap.pdf"),12,6)
rends_heat(filtered_list, annotation_dt, sample = p0_WT, multisamples="average", utr5l = 25, cdsl = 40, utr3l = 25, colour = "#333f50")
dev.off()

pdf(file.path(output_dir,"P0_MUT_PSiteOffset_Heatmap.pdf"),12,6)
rends_heat(filtered_list, annotation_dt, sample = p0_MUT, multisamples="average", utr5l = 25, cdsl = 40, utr3l = 25, colour = "#333f50")
dev.off()

pdf(file.path(output_dir,"P7_WT_PSiteOffset_Heatmap.pdf"),12,6)
rends_heat(filtered_list, annotation_dt, sample = p7_WT, multisamples="average", utr5l = 25, cdsl = 40, utr3l = 25, colour = "#333f50")
dev.off()

pdf(file.path(output_dir,"P7_MUT_PSiteOffset_Heatmap.pdf"),12,6)
rends_heat(filtered_list, annotation_dt, sample = p7_MUT, multisamples="average", utr5l = 25, cdsl = 40, utr3l = 25, colour = "#333f50")
dev.off()

# Calculate p-site offsets
psite_offset <- psite(reads_list, flanking = 6, extremity = "auto")
reads_psite_list <- psite_info(reads_list, psite_offset)

input_samples <- list("P0_WT" = p0_WT, "P0_MUT" = p0_MUT, "P7_WT" = p7_WT, "P7_MUT" = p7_MUT)

psites_per_region<- region_psite(reads_psite_list, annotation_dt,
					 sample = input_samples,
					 multisamples = "average",
					 plot_style = "stack",
					 cl = 85,
					 colour = c("#333f50", "gray70", "#39827c"))

pdf(file.path(output_dir,"PSitesPerRegion_Stacked_Barplot.pdf"),6,6)
psites_per_region[["plot"]]
dev.off()

write.csv(psites_per_region[["count_dt"]], file = file.path(output_dir,"PSitesPerRegion.csv"), row.names = FALSE)

# Periodicity per region heat maps
example_frames_stratified <- frame_psite_length(reads_psite_list, annotation_dt,
                                                 sample = input_samples,
                                                 multisamples = "average",
                                                 plot_style = "facet",
                                                 region = "all",
                                                 cl = 85, colour = "#333f50")

write.csv(example_frames_stratified[["count_dt"]], file = file.path(output_dir,"StratifiedFrames.csv"), row.names = FALSE)

pdf(file.path(output_dir,"StratifiedFrames.pdf"))
example_frames_stratified[["plot"]]
dev.off()

# P sites per frame in each region, all read lengths together

example_frames <- frame_psite(reads_psite_list, annotation_dt,
                               sample = input_samples,
                               multisamples = "average",
                               plot_style = "facet",
                               region = "all",
                               colour = c("#333f50", "#39827c"))

write.csv(example_frames[["count_dt"]], file = file.path(output_dir,"AllReadLengths_Frames.csv"), row.names = FALSE)

pdf(file.path(output_dir,"AllReadLengths_Frames.pdf"))
example_frames[["plot"]]
dev.off()

# Periodicity metaprofiles, facet
metaprofile_facet <- metaprofile_psite(reads_psite_list, annotation_dt,
                                          sample = input_samples,
                                          multisamples = "average",
                                          plot_style = "facet",
                                          utr5l = 20, cdsl = 40, utr3l = 20,
                                          colour = c("#333f50", "#39827c"))

pdf(file.path(output_dir,"Periodicity_Metaprofile_Facet.pdf"))
metaprofile_facet[["plot"]]
dev.off()

metaprofile_overlap <- metaprofile_psite(reads_psite_list, annotation_dt,
                                          sample = input_samples,
                                          multisamples = "average",
                                          plot_style = "overlap",
                                          utr5l = 20, cdsl = 40, utr3l = 20,
                                          colour = c("#333f50", "#39827c"))

pdf(file.path(output_dir,"Periodicity_Metaprofile_Overlap.pdf"))
metaprofile_overlap[["plot"]]
dev.off()

write.csv(metaprofile_facet[["count_dt"]], file = file.path(output_dir,"Metaprofile_Periodicity.csv"), row.names = FALSE)

# Periodicity heatmaps with offsites corrected
metaheatmap <- metaheatmap_psite(reads_psite_list, annotation_dt,
                                         sample = input_samples,
                                         multisamples = "average",
                                         utr5l = 20, cdsl = 40, utr3l = 20,
                                         colour = "#333f50")

pdf(file.path(output_dir,"HeatMaps_OffSetCorrected.pdf"))
metaheatmap[["plot"]]
dev.off()

write.csv(metaheatmap[["count_dt"]], file = file.path(output_dir,"HeatMaps_OffSetCorrected.csv"), row.names = FALSE)
