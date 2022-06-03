prepare_barcodes <- function(metadata_excel,
                             barcode_col="barcode",
                             samplename_col="barcode.ID",
                             adapter_col="adapter.sequence",
                             i=NULL,
                             save_path = file.path(dirname(metadata_excel),
                                                   "barcodes.csv")){
    # metadata_excel = file.path(root,"raw_data/scTIP-seq/scTIP-seq K562.xlsx")
    # save_path <- file.path("/Volumes/bms20/projects/neurogenomics-lab/live/Data/tip_seq/raw_data/scTIP-seq/phase_1_06_apr_2022/X204SC21103786-Z01-F003/raw_data/scTS_1/barcodes.csv")
    # metadata_excel = file.path(root,"raw_data/scTIP-seq/scTIP-seq PFC_male.xlsx")
    Tn5 <- NULL;
    
    meta <- xlsx::read.xlsx(file = metadata_excel, 
                            startRow = 6,
                            sheetIndex = 1) |>
        data.table::data.table()
    #### Split our custom "adapter" sequence into its three components: ####
    ## Tn5 sequence + barcode + adapter sequence (in the traditional sense)
    split_seq <- stringr::str_split(meta$adapter.sequence, 
                                    pattern = meta$barcode, 
                                    simplify = TRUE)  
    ## get the sequence that is most likely to be the Tn5
    meta$Tn5 <- names(utils::tail(sort(table(split_seq[,1])), 1)) 
    meta <- dplyr::mutate(meta, 
                          adapter=stringr::str_replace(
                              string = get(adapter_col), 
                              pattern= paste0(Tn5,get(barcode_col)),""))
    
    #### Report forward adapter sequence ####
    forward_adapter <- unique(meta$adapter)
    message("Unique 'adapter' sequence(s): ",
            paste("\n -",forward_adapter, collapse = "")
    )
    #### Report reverse complement of adapter sequence #####
    reverse_adapter <- as.character(
        Biostrings::reverseComplement(
            Biostrings::DNAStringSet(unique(meta$adapter))
        )
    )
    message("Unique reverse complement 'adapter' sequence(s): ",
            paste("\n -",reverse_adapter, collapse = "")
    ) 
    #### Create (named) barcode sequences #####
    ## Include the Tn5 sequence as part of the "barcode" 
    if(is.null(samplename_col)){
        data <- data.table::copy(meta)[,id:=paste(paste0(Tn5,get(barcode_col)),
                                                  sep=":")][,list(id)] 
    } else {
        data <- data.table::copy(meta)[,id:=paste(get(barcode_col), #paste0(Tn5,get(barcode_col)),
                                                  get(samplename_col),
                                                  sep=":")][,list(id)] 
    } 
    #### Get specific rows ####
    if(!is.null(i)){
        data <- data[i,] 
    }
    #### Save ####
    if(!is.null(save_path)){
        message("Saving barcodes -->",save_path)
        data.table::fwrite(x = data, 
                           file = save_path,
                           sep=",",
                           col.names = FALSE, row.names = FALSE)
    } 
    #### Return ####
    return(list(data=data, 
                save_path=save_path))
}
