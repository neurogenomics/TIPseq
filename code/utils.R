
gather_files <- function(path,
                         type = "peaks"){
    type <- tolower(type[1])
    pattern <- if(type=="peaks"){
        "*.peaks.bed.stringent.bed$"
    } else if (type=="picard") {
        "*.target.markdup.MarkDuplicates.metrics.txt$"
    } else {
        stop("Type must be 'peaks' or 'picard'.")
    }
    paths <- list.files(path = path,
                        pattern = pattern, 
                        recursive = TRUE, 
                        full.names = TRUE)
    if(type=="peaks"){
        paths <- grep("03_peak_calling", paths, value = TRUE)
        names <- paste(basename(dirname(dirname(dirname(paths)))),
                       stringr::str_split(basename(paths),"[.]",
                                          simplify = TRUE)[,1], 
                       sep=".")
    } else if(type=="picard"){
        names <- paste(basename(dirname(dirname(dirname(dirname(dirname(dirname(paths))))))),
                       stringr::str_split(basename(paths),"[.]",
                                          simplify = TRUE)[,1], 
                       sep=".")
    } 
    files <- lapply(paths, function(x){
        message(x)
        if(type=="peaks"){
            dat <- ChIPseeker::readPeakFile(x, as = "GRanges")
        } else if(type=="picard"){ 
            dat <- data.table::fread(x, skip = "LIBRARY",
                                     fill = TRUE,
                                     nrows = 1)
        }
        return(dat)
    }) %>% `names<-`(names)
    return(files)
}


#' Import narrow peaks
#' 
#'   @source \href{https://support.bioconductor.org/p/80535/#80537}{
#' Solution from Bioconductor forum}
import.narrowPeak <- function(con){
    extraCols <- c(signalValue = "numeric", 
                   pValue = "numeric",
                   qValue = "numeric",
                   peak = "integer")
    V2 <- NULL;
    #### Parse column metadat ####
    as_file <- file.path(
        "https://www.encodeproject.org/documents",
        "948203bb-8eb2-42a2-8b12-1c10f356c998/@@download/attachment",
        "narrowPeak.as"
    )
    f <- readLines(as_file)
    f <- f[!f %in% c("(",")")][-seq_len(2)]
    type_dict <- c("string"="character", 
                   "char[1]"="character",
                   "uint"="integer",
                   "int"="integer",
                   "float"="numeric")
    meta <- data.table::fread(text = f, 
                              header = FALSE) %>%
        tidyr::separate(col = "V1", 
                        into = c("type","name"), 
                        sep = " +") %>%
        dplyr::rename(info=V2) %>% 
        dplyr::mutate(typer=type_dict[type])
    extraCols <- stats::setNames(meta$typer, meta$name)
    reserved <- c("chrom","name","score","strand","chromStart","chromEnd")
    extraCols <- extraCols[!names(extraCols) %in% reserved]
    #### Import peak data ####
    peaks <- rtracklayer::import(con = con, 
                                 extraCols = extraCols)
    message(formatC(length(peaks), big.mark = ",")," peaks imported.")
    return(peaks)
}



rename_files <- function(files,
                         sample_dict,
                         max_char=NULL){
    split_names <- stringr::str_split(names(files),"[.]", simplify = TRUE)
    names(files) <- paste(
        sample_dict[split_names[,1]],
        split_names[,1],
        split_names[,2],
        sep = "."
    )
    if(!is.null(max_char)){
        names(files) <- stringr::str_trunc(names(files), width = max_char) 
    }
    return(files)
}
