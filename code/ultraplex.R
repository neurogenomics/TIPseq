
ultraplex <- function(ultraplex_path="ultraplex",
                      inputfastq,
                      inputfastq2=NULL,
                      barcodes,
                      directory=file.path(dirname(inputfastq),"demultiplexed"),
                      threads=4,
                      verbose=TRUE){
    # echoconda::import_cli()
    cmd <- paste(ultraplex_path,
                 "-i",inputfastq,
                 if(is.null(inputfastq2)) NULL else paste("-i2",inputfastq2),
                 "-b",barcodes,
                 "-d",directory,
                 "-t",threads
    ) 
    echoconda::cmd_print(cmd = cmd, verbose = verbose)
    echoconda::find_conda()
    out_files <- list.files(path = directory,
                            pattern = "*.fastq.gz$",
                            full.names = TRUE)
    return(out_files)
}
