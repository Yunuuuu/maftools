#' Merge multiple mafs into single MAF
#' @description Merges multiple maf files/objects/data.frames into a single MAF.
#' @param mafs a list of \code{\link{MAF}} objects or data.frames or paths to MAF files.
#' @param verbose Default TRUE
#' @param ... additional arguments passed \code{\link{read.maf}}
#' @param idcol Creates a column in the \code{clinicalData} showing which list
#' item those rows came from. \code{TRUE} names this column \code{".id"}.
#' \code{idcol="file"} names this column \code{"file"}. If the input list has
#' names, those names are the values placed in this id column, otherwise the
#' values are an integer vector \code{1:length(mafs)}. Only used when
#' \code{mafs} is a list of \code{\link{MAF}} objects.
#' @return \code{\link{MAF}} object
#' @export
#'
merge_mafs = function(mafs, verbose = TRUE, ..., idcol = TRUE){

  if(all(unlist(lapply(mafs, is, "MAF")))){

    if(verbose){
      cat(paste0("Merging ", length(mafs) ," MAF objects\n"))
    }
    mafs_dat = lapply(mafs, function(m) {
      data.table::rbindlist(l = list(m@data, m@maf.silent), use.names = TRUE, fill = TRUE)
    })
    mafs_dat = data.table::rbindlist(l = mafs_dat, fill = TRUE, use.names = TRUE)

    mafs_clin = lapply(mafs, function(m) {
      m@clinical.data
    })
    mafs_clin = data.table::rbindlist(l = mafs_clin, fill = TRUE, use.names = TRUE, idcol = idcol)
    maf = read.maf(maf = mafs_dat, clinicalData = mafs_clin, verbose = verbose, ...)
  }else if(all(unlist(lapply(mafs, is, "data.frame")))){
    if(verbose){
      cat(paste0("Merging ", length(mafs) ," MAF data.frames\n"))
    }
    mafs_dat = data.table::rbindlist(l = mafs, fill = TRUE, use.names = TRUE)
    maf = read.maf(maf = mafs_dat, verbose = verbose, ...)
  }else{
    if(verbose){
      cat(paste0("Merging ", length(mafs) ," MAF files\n"))
    }
    maf = lapply(mafs, function(x) {
      x = data.table::fread(x, stringsAsFactors = FALSE, fill = TRUE, showProgress = TRUE, header = TRUE, skip = "Hugo_Symbol")

      required.fields = c(
        'Hugo_Symbol',
        'Chromosome',
        'Start_Position',
        'End_Position',
        'Reference_Allele',
        'Tumor_Seq_Allele2',
        'Variant_Classification',
        'Variant_Type',
        'Tumor_Sample_Barcode'
      )

      #Change column names to standard names; i.e, camel case
      for (i in 1:length(required.fields)) {
        colId = suppressWarnings(grep(
          pattern = paste("^", required.fields[i], "$", sep = ""),
          x = colnames(x),
          ignore.case = TRUE
        ))
        if (length(colId) > 0) {
          colnames(x)[colId] = required.fields[i]
        }
      }
      x
    })
    #names(maf) = gsub(pattern = "\\.maf$", replacement = "", x = basename(path = unlist(mafs)), ignore.case = TRUE)
    names(maf) = basename(mafs)
    maf = data.table::rbindlist(l = maf, fill = TRUE, idcol = "Source_MAF", use.names = TRUE)

    maf = read.maf(maf = maf, verbose = verbose, ...)
  }

  maf
}
