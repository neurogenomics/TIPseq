
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
