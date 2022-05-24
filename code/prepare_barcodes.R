perpare_barcodes <- function(metadata_excel,
                             save_path = file.path(dirname(metadata_excel),
                                                   "barcodes.csv")){
    # metadata_excel = file.path(root,"raw_data/scTIP-seq/scTIP-seq K562.xlsx")
    # metadata_excel = file.path(root,"raw_data/scTIP-seq/scTIP-seq PFC_male.xlsx")
    barcode <- barcode.ID <- NULL;
    
    meta <- xlsx::read.xlsx(file = metadata_excel, 
                            startRow = 6,
                            sheetIndex = 1) |>
        data.table::data.table()
    data <- data.table::copy(meta)[,id:=paste(barcode,
                                              barcode.ID,
                                              sep=":")][,list(id)] 
    data.table::fwrite(x = data, 
                       file = save_path,
                       sep=",",
                       col.names = FALSE, row.names = FALSE)
    return(list(data=data, 
                save_path=save_path))
}



