library(dplyr)
library(tidyr)

root <- "/Volumes/bms20/projects/neurogenomics-lab/live/Data/tip_seq" 

fq <- list.files(path = file.path(root, "raw_data"), 
                 pattern = "*.fq.gz", 
                 recursive = TRUE, full.names = TRUE)
fq <- gsub(paste(
    "/rds/general/user/bms20/home/neurogenomics/",
    "/Volumes/bms20/projects/neurogenomics-lab/live/",
    sep = "|"
),
     "/rds/general/project/neurogenomics-lab/live/", 
     fq)


d <- data.table::data.table(fq=fq) %>% 
    dplyr::mutate(batch=basename(dirname(dirname(dirname(dirname(fq))))),
                  batch_id=basename(dirname(dirname(dirname(fq)))),
                  sample=basename(dirname(fq)),
                  id=gsub(".fq.gz$","",basename(fq)), 
                  control_group=1) 
d$lane <- sapply(stringr::str_split(d$id,"_"), function(x)rev(x)[2])
d$fq_num <- sapply(stringr::str_split(d$id,"_"), function(x)rev(x)[1])
d$replicate <- 1

#### Omit any data that's already been processed ####
processed <- basename(list.dirs(file.path(root,"processed_data"),
                                recursive = FALSE))
# d <- subset(d, !batch %in% processed)


force_new <- TRUE
for(b in unique(d$batch)){
    message(b)
    #### Spread data into wide format #### 
    ## With separate columns for fastq_1 and fastq_2
    dat <- tidyr::pivot_wider(data = subset(d, batch==b), 
                              id_cols = c("batch","batch_id",
                                          "replicate","control_group",
                                          "sample","lane"), 
                              names_from = "fq_num",
                              names_glue="fastq_{fq_num}",
                              values_from = "fq") %>%
        ## Group change depending on analysis goals 
        ## Required cols: "group","replicate","control_group","fastq_1","fastq_2","sample"
        dplyr::select(group=sample,replicate,control_group,fastq_1,fastq_2,
                      sample,lane,batch)
    save_dir <- dirname(dirname(dirname(dirname(dat$fastq_1[1]))))
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
