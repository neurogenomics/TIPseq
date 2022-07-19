make_fragments_file <- function(dir,
                                peak_files,
                                keep.chr=NULL,
                                style="UCSC",
                                save_path=tempfile(fileext = "fragments.tsv.gz")){
    # dir <- "/home/bms20/RDS/project/neurogenomics-lab/live/Data/tip_seq/processed_data/scTS_3_30_jun_2022/03_peak_calling/06_fragments/"
    # save_path <- "/home/bms20/RDS/project/neurogenomics-lab/live/Data/tip_seq/processed_data/scTS_3_30_jun_2022/fragments_merged.tsv"
    #### Find bin files ####
    ## Bin size set here in nf-core/cutandrun pipeline:
    ## https://github.com/nf-core/cutandrun/blob/971984a48ad4dc5b39fc20b56c5729f4ca20379a/conf/modules.config#L666
    bin_files <- list.files(dir, 
                            pattern = "\\.frags\\.bin[0-9]+\\.awk\\.bed", 
                            full.names = TRUE) 
    #### Find fragment files ####
    fragment_files <- list.files(dir, 
                                 pattern = "\\.frags\\.cut\\.bed", 
                                 full.names = TRUE) 
    #### Merge all fragment files into one ####
    gr.frag <- lapply(bin_files, function(x){ 
        message("Processing: ",basename(x))
        dt <- data.table::fread(x, col.names = c("seqnames","start",
                                                 "overlap","name")) 
        dt[,name:=(gsub("\\.frags\\.cut\\.bed","",name))]
        binsize <- as.integer(
            gsub(
                "bin","",
                grep("^bin",strsplit(basename(x),"\\.")[[1]],value = TRUE)
            )
        )
        dt[,end:=(start+binsize)]
        gr <- GenomicRanges::makeGRangesFromDataFrame(dt, 
                                                      keep.extra.columns = TRUE) 
        return(gr)
    }) |> 
        GenomicRanges::GRangesList() |> 
        unlist() |>
        GenomicRanges::sort() 
    GenomeInfoDb::seqlevelsStyle(gr.frag) <- style
    #### Remove non-standard chromosomes #### 
    if(!is.null(keep.chr)){
        gr.frag <- gr.frag[GenomicRanges::seqnames(gr.frag) %in% keep.chr,]
    }
    #### Save merged fragments file ####
    if(!is.null(save_path)){
        message("Saving merged fragments file ==> ",save_path)
        save_path_unzipped <- gsub("\\.gz","",save_path)
        df <- data.frame(gr.frag)[,c("seqnames","start","end","name","overlap")]
        data.table::fwrite(df, save_path_unzipped, 
                           col.names = FALSE,
                           row.names = FALSE, 
                           sep = "\t")
        message("Tabix-indexing merged fragments file.")
        save_path <- Rsamtools::bgzip(save_path_unzipped, 
                                      dest = save_path, 
                                      overwrite = TRUE)
        tbi <- Rsamtools::indexTabix(file = save_path, 
                                     format = "bed")
    }
    return(list(gr=gr.frag, 
                file=save_path,
                tbi=tbi))
    # getFragmentsOverlapsIDs = function(file_path){
    #     cell = ChIPseeker::readPeakFile(file_path)
    #     cell_name = basename(file_path)
    #     overlap = GenomicRanges:::countOverlaps_GenomicRanges(query = fragments, subject = cell)
    #     fragments$name = cell_name
    #     fragments$score = overlap
    #     return(fragments)
    # } 
     
}