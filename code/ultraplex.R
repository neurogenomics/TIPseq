ultraplex <- function(ultraplex_path="ultraplex",
                      inputfastq,
                      inputfastq2=NULL, 
                      metadata_excel,
                      barcode_col="barcode",
                      samplename_col="barcode.ID",
                      adapter_col="adapter.sequence",
                      directory=file.path(dirname(inputfastq),
                                          "demultiplexed2"),
                      threads=4,
                      conda_env="ultra",
                      verbose=TRUE){  
    # inputfastq <- file.path(
    #     root,
    #     "raw_data/scTIP-seq/phase_1_06_apr_2022/X204SC21103786-Z01-F003",
    #     "raw_data/scTS_1/scTS_1_EKDL220003595-1a_HHN53DSX3_L1_1.fq.gz")
    # inputfastq2 <- file.path(
    #     root,
    #     "raw_data/scTIP-seq/phase_1_06_apr_2022/X204SC21103786-Z01-F003",
    #     "raw_data/scTS_1/scTS_1_EKDL220003595-1a_HHN53DSX3_L1_2.fq.gz")
    #### Get ultraplex executable ####
    env_url <- "https://github.com/ulelab/ultraplex/raw/master/environment.yml"
    if(!is.null(conda_env)){
        if(conda_env=="ultra"){
            conda_env <- echoconda::yaml_to_env(env_url)
        }
    } 
    ultraplex_path <- echoconda::find_packages(packages = "ultraplex",
                                               conda_env = "ultra",
                                               return_path = TRUE)   
    #### Import metadata ####
    meta <- xlsx::read.xlsx(file = metadata_excel, 
                            startRow = 6,
                            sheetIndex = 1)
    if(adapter_col %in% colnames(meta)){
        adapters <- unique(meta[,adapter_col])
    } 
    dir.create(directory, showWarnings = FALSE, recursive = TRUE)
    outputs <- lapply(seq_len(length(adapters)), function(i){
        adapter <- adapters[i]
        #### Prepare barcodes ####
        barcodes_path <- prepare_barcodes(metadata_excel = metadata_excel, 
                                          barcode_col = barcode_col, 
                                          samplename_col = samplename_col,
                                          i = i, 
                                          save_path = tempfile(
                                              fileext = paste0(
                                                  "_barcode_",i,".csv")
                                              )
                                          )
        cmd <- paste(ultraplex_path,
                     "-i",inputfastq,
                     if(is.null(inputfastq2)) {
                         NULL
                     } else {
                         paste("-i2",inputfastq2)
                     },
                     "-a",adapter,
                     "-b",barcodes_path$save_path,
                     "-d",directory,
                     "-t",threads
        ) 
        echoconda::cmd_print(cmd = cmd)
        system(cmd)
    })  
    out_files <- list.files(path = directory,
                            pattern = "*.fastq.gz$",
                            full.names = TRUE)
    return(out_files)
}
