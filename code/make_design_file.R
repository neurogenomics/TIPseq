library(dplyr)
library(tidyr)

root <- "/Volumes/bms20/projects/neurogenomics-lab/live/Data/tip_seq" 
fastq_dir <- "raw_data/scTIP-seq/phase_1_06_apr_2022/X204SC21103786-Z01-F003/raw_data/scTS_1/demultiplexed"
fq <- list.files(path = file.path(root, fastq_dir), 
                 pattern = "*.fq.gz|*.fastq.gz", 
                 recursive = TRUE, full.names = TRUE)
fq <- gsub(paste(
    "/rds/general/user/bms20/home/neurogenomics/",
    "/Volumes/bms20/projects/neurogenomics-lab/live/",
    sep = "|"
),
     "/rds/general/project/neurogenomics-lab/live/", 
     fq)

input_type <- "ultraplex" # "nfcore_cutandrun"
if(input_type=="nfcore_cutandrun"){
    d <- data.table::data.table(fq=fq) %>% 
        dplyr::mutate(batch=basename(dirname(dirname(dirname(dirname(fq))))),
                      batch_id=basename(dirname(dirname(dirname(fq)))),
                      group=basename(dirname(fq)),
                      id=gsub(".fq.gz$","",basename(fq)), 
                      control_group=1) 
} else if(input_type=="ultraplex"){
    d <- data.table::data.table(fq=fq) %>% 
        dplyr::mutate(batch=basename(dirname(dirname(dirname(dirname(fq))))), 
                      group=basename(dirname(dirname(fq))),
                      id=gsub(".fq.gz$|.fastq.gz$","",basename(fq)), 
                      control_group=1) 
    d$batch_id <- sapply(stringr::str_split(d$id,"_"),
                         function(x)paste(x[-c(1,2,length(x))],collapse = "_"))
}

d$lane <- sapply(stringr::str_split(d$id,"_"), 
                 function(x){
                     ln <- toupper(rev(x)[2])
                     if(startsWith(ln,"L")) ln else NA
                 })
d$fq_num <- sapply(stringr::str_split(d$id,"_"), 
                   function(x)gsub(".fastq.gz$|.fq.gz$","",rev(x)[1]))
### Ultraplex outputs forward/reverse reads as Fwd and Rev (instead of 1/2)
d <- dplyr::mutate(d, 
    fq_num=as.numeric(ifelse(tolower(fq_num)=="fwd",1,
                      ifelse(tolower(fq_num)=="rev",2,fq_num)))
    )
d$replicate <- 1

#### Remove files with reads that couldn't be matched to any barcode ####
d <- d[!grepl("ultraplex_demux_5bc_no_match", d$fq),]

#### Omit any data that's already been processed ####
processed <- basename(list.dirs(file.path(root,"processed_data"),
                                recursive = FALSE))
# d <- subset(d, !batch %in% processed)


force_new <- TRUE
for(b in unique(d$batch)){
    message(b)
    #### Spread data into wide format #### 
    ## With separate columns for fastq_1 (Fwd) and fastq_2 (Rev) 
    dat <- tidyr::pivot_wider(data = subset(d, batch==b), 
                              id_cols = c("batch","batch_id",
                                          "replicate","control_group",
                                          "group","lane"), 
                              names_from = "fq_num",
                              names_glue="fastq_{fq_num}",
                              values_from = "fq") %>%
        ## Group change depending on analysis goals 
        ## Required cols: "group","replicate","control_group","fastq_1","fastq_2","sample"
        dplyr::select(group,
                      replicate,fastq_1,fastq_2)
    if(input_type=="nfcore_cutandrun"){
        save_dir <- dirname(dirname(dirname(dirname(dat$fastq_1[1]))))
    } else if(input_type=="ultraplex"){
        save_dir <- dirname(dirname(dirname(dirname(dirname(dat$fastq_1[1])))))
    }
   
    save_dir <-  gsub("/rds/general/project/neurogenomics-lab/live/", 
                      "/Volumes/bms20/projects/neurogenomics-lab/live/",
                      save_dir)
    fpath <- file.path(save_dir,
                       "design.csv") 
    if(!file.exists(fpath) || force_new){
        message("Writing ==> ",fpath)
        write.csv(x = dat, file = fpath, row.names = FALSE)
    }
}
