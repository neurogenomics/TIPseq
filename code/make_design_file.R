make_design_file <- function(
        root="/Volumes/bms20/projects/neurogenomics-lab/live/Data/tip_seq",
        fastq_dir="raw_data/scTIP-seq/scTS_3_30_jun_2022/X204SC21103786-Z01-F007/",
        input_type=c("nfcore_cutandrun","ultraplex"),
        replace_dir=list("/rds/general/project/neurogenomics-lab/live/"=
                        c("/Volumes/bms20/projects/neurogenomics-lab/live/",
                          "/rds/general/user/bms20/home/neurogenomics/")),
        replicates=1,
        force_new=TRUE,
        exclude_processed=FALSE,
        save_dir=NULL,
        return_paths=FALSE){
    
    fq_num <- fastq_1 <- fastq_2 <- control <- group <- NULL;
    
    input_type <- tolower(input_type)[1]
    fq <- list.files(path = file.path(root, fastq_dir), 
                     pattern = "*.fq.gz|*.fastq.gz", 
                     recursive = TRUE,
                     full.names = TRUE)
    #### Replace dir root ####
    gsub_dir <- function(x,
                         replace_dir,
                         invert=FALSE){
        if(isTRUE(invert)){
            gsub(names(replace_dir),
                 replace_dir[[1]][1], x, fixed = FALSE) 
        } else {
            gsub(paste(replace_dir[[1]], collapse = "|"),
                 names(replace_dir), x, fixed = FALSE) 
        } 
    }
    if(!is.null(replace_dir)){
        fq <- gsub_dir(x=fq, replace_dir=replace_dir)
    }
    #### Create design matrix ####
    if(input_type=="nfcore_cutandrun"){
        d <- data.table::data.table(fq=fq) |>
            dplyr::mutate(batch=basename(dirname(dirname(dirname(dirname(fq))))),
                          batch_id=basename(dirname(dirname(dirname(fq)))),
                          group=basename(dirname(fq)),
                          id=gsub(".fq.gz$","",basename(fq)), 
                          control="") 
    } else if(input_type=="ultraplex"){
        d <- data.table::data.table(fq=fq) |>
            dplyr::mutate(batch=basename(dirname(dirname(dirname(dirname(fq))))), 
                          group=basename(dirname(dirname(fq))),
                          id=gsub(".fq.gz$|.fastq.gz$","",basename(fq)), 
                          control="") 
        d$batch_id <- sapply(stringr::str_split(d$id,"_"),
                             function(x)paste(x[-c(1,2,length(x))],collapse = "_"))
    }
    
    d$lane <- sapply(stringr::str_split(d$id,"_"), 
                     function(x){
                         ln <- toupper(rev(x)[2])
                         if(startsWith(ln,"L")) ln else NA
                     })
    d$fq_num <- sapply(stringr::str_split(d$id,"_"), 
                       function(x)gsub(".fastq.gz$|.fq.gz$","",rev(x)[1],
                                       ignore.case = TRUE))
    ### Ultraplex outputs forward/reverse reads as Fwd and Rev (instead of 1/2)
    d <- dplyr::mutate(d, 
                       fq_num=as.numeric(
                           ifelse(tolower(fq_num)=="fwd",1,
                                  ifelse(tolower(fq_num)=="rev",2,fq_num)))
    ) 
    #### Remove files with reads that couldn't be matched to any barcode ####
    d <- d[!grepl("ultraplex_demux_5bc_no_match", d$fq),]
    #### Omit any data that's already been processed ####
    if(isTRUE(exclude_processed)){
        processed <- basename(list.dirs(file.path(root,"processed_data"),
                                        recursive = FALSE))
        d <- subset(d, !batch %in% processed)
    } 
    #### Create design files ####
    batches <- unique(d$batch)
    design_files <- lapply(batches, function(b){
        message("Preparing: ",b)
        #### Spread data into wide format #### 
        ## With separate columns for fastq_1 (Fwd) and fastq_2 (Rev) 
        dat <- tidyr::pivot_wider(data = subset(d, batch==b), 
                                  id_cols = c("batch","batch_id",
                                              "control",
                                              "group","lane"), 
                                  names_from = "fq_num",
                                  names_glue="fastq_{fq_num}",
                                  values_from = "fq") |>
            ## Group changes depending on analysis goals 
            ## Required cols: "group","replicate","fastq_1","fastq_2","control"
            dplyr::select(group,fastq_1,fastq_2,control)
        #### Add replicates info ####
        if(length(replicates)>1 && length(replicates)!=nrow(dat)){
            stop("replicates must either be of length 1, ",
                 "or equal to the number of samples (",nrow(dat),").")
        }
        dat <- dplyr::mutate(dat,replicate=replicates, .after=1) 
        #### Construct save_dir ####
        if(is.null(save_dir)){
            if(input_type=="nfcore_cutandrun"){
                save_dir <- dirname(dirname(dirname(dirname(dat$fastq_1[1]))))
            } else if(input_type=="ultraplex"){
                save_dir <- dirname(dirname(dirname(dirname(dirname(dat$fastq_1[1])))))
            } 
        } 
        if(!is.null(replace_dir)){
            save_dir <- gsub_dir(x=save_dir, 
                                 replace_dir=replace_dir, 
                                 invert = TRUE)
        } 
        fpath <- file.path(save_dir,
                           "design.csv") 
        if(!file.exists(fpath) || force_new){
            message("Writing ==> ",fpath) 
            write.csv(x = dat, file = fpath, row.names = FALSE)
        }
        if(return_paths){
            return(fpath)
        } else {
            return(dat)
        } 
    }) |> `names<-`(batches)
    return(design_files)
}
