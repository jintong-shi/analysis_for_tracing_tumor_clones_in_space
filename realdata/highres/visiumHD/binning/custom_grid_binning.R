#' @include load_data.R
#' @include get_hd_barcodes.R
#'
NULL

#' Perform custom grid binning for VisiumHD data
#' 
#' This function groups 2um squares into lower resolution squares. If one wants to
#' have count data corresponding to 4um squares instead of 2um squares, the function
#' identifies the \code{2*2 = 4} 2um squares inside each 4um square, adds up the 
#' count data and maps it to the histology image.
#' 
#' @details (
#' \code{CustomBinning} takes a \code{data.frame} or \code{tbl} as input that should
#' hold certain spaceranger output file paths for the 2um bin size. The table should
#' consist of four columns: "samples", "imgs", "spotfiles" and "json".
#'  \itemize{
#'    \item{"samples" : file paths to expression matrices, e.g. 
#'      \code{binned_outputs/square_002um/raw_bc_matrix.h5}}
#'    \item{"imgs" : file paths to images, e.g. \code{tissue_hires_image.jpg}}
#'    \item{"spotfiles" : file paths to spot coordinate \code{.parquet} files 
#'      \code{binned_outputs/square_002um/spatial/tissue_positions.parquet}}
#'    \item{"json" : file paths to scale factor JSON files, e.g. \code{scalefactors_json.json}}
#' }
#'
#' @param infoTable A \code{data.frame} or \code{tbl} with paths to spaceranger output files for the 2um bin
#' @param bin Integer specifying the desired bin size. Defaults to \code{48}. The bin size specified here must be divisible by 2 (i.e. \code{bin \link[base]{\%\%} 2 == 0})
#' @param remove_spots_outside_tissue Should spots outside the tissue be removed?
#' @param verbose Print messages 
#'
#' 
#' @import cli
#' @import rlang
#' @import dplyr
#' @importFrom glue glue
#' @importFrom Matrix sparseMatrix
#' @importFrom Seurat CreateSeuratObject
#' @importFrom glue glue
#' @importFrom dplyr case_when mutate
#' @importFrom tibble tibble
#' 
#' @return Returns a \code{Seurat} object compatible with \code{semla} and all its functions
#' 
#' @author Javier Escudero Morlanes
#' 
#' @family segment-visHD
#'
#' @export
#'
CustomBinning <- function(
    infoTable,
    bin = 48,
    verbose = TRUE
) {
  # Set global variables to NULL
  # bin <- bin_scale <- NULL
  library(Matrix)
  if (!requireNamespace("rhdf5", quietly = TRUE)) {
    abort(glue("Package {cli::col_br_magenta('rhdf5')} is required. Please install it with: \n",
               "BiocManager::install('rhdf5')"))
  }
  
  if (verbose) cli_h2("Reading and binning 10x VisiumHD data")
  
  # Check infoTable
  if (!all(c("samples", "imgs", "spotfiles") %in% colnames(infoTable)))
    abort("One or several of 'samples', 'imgs' and 'spotfiles' are missing from infoTable.")
  if (!any(c("json", "scalefactor") %in% colnames(infoTable)))
    abort("One of 'json' or 'scalefactor' columns needs to be provided")
  if (!any(class(infoTable) %in% c("tbl", "data.frame"))) abort(glue("Invalid class '{class(infoTable)}' of 'infoTable'."))
  if (nrow(infoTable) == 0) abort(glue("'infoTable' is empty."))
  if (!all(sapply(infoTable[, c("samples", "imgs", "spotfiles")], class) %in% "character")) abort(glue("Invalid column classes in 'infoTable'. Expected 'character' vectors."))
  if ("scalefactor" %in% colnames(infoTable)) {
    if (!inherits(infoTable[, "scalefactor", drop = TRUE], what = c("numeric", "integer")))
      abort(glue("Invalid column class for 'scalefactor'. Expected a 'numeric'."))
    if (!all(between(x = infoTable[, "scalefactor", drop = TRUE], left = 0.001, right = 1)))
      abort(glue("Scalefactors need to be in the range 0.001-1."))
  }
  missing_files <- infoTable |>
    select(all_of(c("samples", "spotfiles")), contains("json")) |>
    mutate(across(everything(), ~ file.exists(.x))) |>
    summarize(across(everything(), ~ any(!.x))) |>
    unlist()
  if (any(missing_files)) {
    checks <- missing_files[missing_files]
    sapply(names(checks), function(n) {
      abort(glue("Missing file(s) in the '{n}' column."))
    })
  }
  if ("json" %in% colnames(infoTable)) {
    if (!all(file.exists(infoTable |> pull(all_of("json"))))) {
      abort(glue("Missing json file(s)."))
    }
  }
  
  # Get bin_scale
  if(bin %% 2 != 0){
    abort(glue("The desired bin should be divisible by 2. Provided bin size is not: '{bin} %% 2 == {bin %% 2}'"))
  } 
  bin_scale <- bin / 2
  
  # Keep additional infoTable columns
  if (ncol(infoTable) > 4) {
    additionalMetaData <- infoTable |>
      select(-all_of(c("samples", "imgs", "spotfiles", "json")))
  } else {
    additionalMetaData <- NULL
  }
  
  # Pull information from infoTable
  h5_path <- infoTable$samples
  coord_path <- infoTable$spotfiles
  scalef_path <- infoTable$json
  img_path <- infoTable$imgs
  
  # Retrieve count data
  if (verbose) cli_alert_info("Loading matrix")
  ## get information from h5
  h5_list <- .get_from_h5(h5_file = h5_path)
  h5closeAll()
  ## assemble COO matrix from CSR 002um data
  cnt <- Matrix::sparseMatrix(i = h5_list$row_indices,
                              p = h5_list$col_ptrs,
                              x = h5_list$values,
                              dims = h5_list$shape,
                              repr = "T")
  
  # Spots to bin
  coord_df <- LoadSpatialCoordinates(coordinatefiles = coord_path,
                                     remove_spots_outside_tissue = FALSE,
                                     verbose = verbose) |> 
    filter(barcode %in% h5_list$barcodes)
  
  if (verbose) cli_alert_info("Keeping coordinates for {cli::col_br_magenta(nrow(coord_df))} spots with count data")
  
  coord_df <- coord_df[match(h5_list$barcodes, coord_df$barcode), ]
  coord_bin <- .get_bin_indices(bin_scale = bin_scale,
                                coord_df = coord_df)
  if (verbose) cli_alert_info("Identified indices of spots to bin together")
  ## assing new synthetic barcodes to the new spots
  array_bin <- paste0("s_", sprintf("%03d", bin), "um_", sprintf("%05d", coord_bin$y_bin), "_", sprintf("%05d", coord_bin$x_bin))
  
  if (verbose) cli_alert_info("Created synthetic barcode IDs.")
  if(max(table(array_bin)) > bin_scale^2){abort(glue("Expected a max of '{bin_scale^2}' hires spot per lowres spot. Got '{max(table(array_bin))}'"))}
  
  # Create the binned count data
  cnt_bin <- .bin_counts(new_barcodes = array_bin,
                         i = cnt@i,
                         j = cnt@j,
                         x = cnt@x,
                         shape = c(nrow(cnt), length(unique(array_bin))))
  rownames(cnt_bin) <- h5_list$nm # the genes should be in the same order as in the h5
  colnames(cnt_bin) <- unique(array_bin) # the new barcodes should be in the same order binned as the originals were on the h5
  if (verbose) cli_alert_info(glue("Assembled binned count matrix for {cli::col_br_magenta(ncol(cnt_bin))} spots"))
  
  # Spatial information for binned spots
  ## retrieve indices of unique lowres occurrences
  keep <- !duplicated(array_bin)
  ## get pixel coordinates of lowres spots
  pxl_coords <- .get_bin_pxl(array_bin = array_bin,
                             coord_df = coord_df,
                             keep = keep)
  
  # Assemble Seurat object with the count data
  # remove duplicateds added by jintong
  #duplicated_genenames <- rownames(cnt_bin)[duplicated(rownames(cnt_bin))]
  #cnt_bin_unique <- cnt_bin[setdiff(rownames(cnt_bin),unique_genenames),]

  if (verbose) cli_h3(text = "Creating `Seurat` object")
  object <- CreateSeuratObject(counts = cnt_bin,
                               assay = "Spatial",
                               meta.data = tibble(barcode = unique(array_bin),
                                                  resolution = paste0("0", bin)) |> 
                                 column_to_rownames(var = "barcode"))
  if (verbose) cli_alert_info("Created `Seurat` object")
  
  # Assemble Staffli
  sample_ID <- as.integer(1)
  ## coordiantes dataframe
  coordinates_bin <- tibble(barcode = array_bin[keep],
                            x = as.integer(coord_bin$x_bin[keep]),
                            y = as.integer(coord_bin$y_bin[keep]),
                            pxl_col_in_fullres = as.integer(pxl_coords$pxl_col_bin),
                            pxl_row_in_fullres = as.integer(pxl_coords$pxl_row_bin),
                            sampleID = sample_ID
  )
  ## read scalefactors
  scalef_info <- LoadScaleFactors(scalef_path) |>
    mutate(spot_diameter_fullres = spot_diameter_fullres * bin_scale,
           bin_size_um = bin,
           sampleID = sample_ID)
  ## read image and append scalefactors
  img_info <- LoadImageInfo(img_path) |> 
    mutate(sampleID = sample_ID)
  img_info <- UpdateImageInfo(image_info = img_info, 
                              scalefactors = scalef_info) |> 
    mutate(sampleID = sampleID)
  ## create staffli
  staffli_object <- CreateStaffliObject(imgs = img_path,
                                        meta_data = coordinates_bin,
                                        image_info = img_info,
                                        scalefactors = scalef_info)
  if (verbose) cli_alert_info("Created `Staffli` object")
  
  # Put together
  object@tools$Staffli <- staffli_object
  
  if (verbose) cli_alert_success(glue("Returning a `Seurat` object with {cli::col_br_blue(nrow(object))}",
                                      " features and {cli::col_br_magenta(ncol(object))} spots"))
  return(object)
}




#' Retrieve information from specified h5 file
#'
#' @param h5_file Full path to the 2um h5 file
#' @return A list containing row indices, column pointers, count values, barcode and gene names
#'
#' @noRd
.get_from_h5 <- function(h5_file) {
  # read h5
  data <- rhdf5::h5read(h5_file, "matrix")
  # # genes targetted
  # target_genes <- data$features$target_sets[[1]] + 1 # R indexing
  # gene indices, keep 
  row_indices <- data$indices + 1
  # spot indices
  col_ptrs <- data$indptr # no need to add +1 to row pointers, as they are supposed to start at 0
  # count data
  values <- data$data
  # barcode names
  barcodes <- data$barcodes
  # gene names remove gene names of genes that are not captured
  nm <- data$features$name
  # assemble in list
  out <- list("row_indices" = row_indices,
              "col_ptrs" = col_ptrs,
              "values" = values,
              "barcodes" = barcodes,
              "nm" = nm,
              "shape" = data$shape)
  
  return(out)
}

#' #' Read hires coordinates and filter spots
#' #' 
#' #' @param coord_path Full path to the hires coordinates file
#' #' @param barcodes Vector of strings indicating the barcodes to keep. Usually those in the count matrix
#' #' 
#' #' @return A tibble with the spatial coordinates of the spots of interest
#' #' 
#' #' @noRd
#' .read_coords <- function(
    #'     coord_path,
#'     barcodes
#' ){
#'   arrow::read_parquet(coord_path) |> 
#'     filter(barcode %in% barcodes)
#' }

#' Given the desired bin size, identify which spots should be binned
#' 
#' @param bin_scale Integer specifying the desired bin scale
#' @param coord_df Dataframe containing the spatial information of the highres
#' @param barcode Vector of barcode IDs, used to order the spatial dataframe
#' 
#' @return A list containing the \code{x} and \code{y} coordinates of the spots
#' that should be binned
#' 
#' @noRd
.get_bin_indices <- function(
    bin_scale,
    coord_df,
    barcodes){
  # order the spatial df to match the counts
  coord <- coord_df
  # coord <- coord_df[match(barcodes, coord_df$barcode), ]
  # modify indices
  ## retrieve OG coordinates
  x <- coord$x
  y <- coord$y
  ## which spots should be grouped
  x_bin <- floor(x / bin_scale)
  y_bin <- floor(y / bin_scale)
  # put together output
  out <- list("x_bin" = x_bin,
              "y_bin" = y_bin)
  
  return(out)
}

#' Given the desired bin size, identify which spots should be binned
#' 
#' @param new_barcodes Vector of strings with the new barcode IDs
#' @param i Vector of integers storing the gene indices in sparse format
#' @param j Vector of integers storing the spot indices in sparse format
#' @param x Vector of integers storing the count data in sparse format
#' @param shape Vector indicating nrows and ncols that the binned matrix should have
#' 
#' @return A sparse matrix with the binned count data
#' 
#' @noRd
.bin_counts <- function(
    new_barcodes,
    i,
    j,
    x,
    shape
) {
  # create map of unique array coordinates = numerical value. these integers should represent the unique values of array_bin
  unique_re <- unique(new_barcodes)
  mapping <- setNames(seq_along(unique_re), unique_re)  # map strings to numbers
  # generate numerical vector of highres spots to be binned
  j_re <- unname(mapping[new_barcodes[j + 1]])
  # create the new count matrix
  if(!identical(length(i), length(x)) && identical(length(j_re), length(i))){
    abort("Something weird happeniiiiiing")
  }
  cnt_bin <- Matrix::sparseMatrix(i = i + 1,
                                  j = j_re, # we already addressed R's indexing above
                                  x = x,
                                  dims = shape,
                                  repr = "T")
  
  return(cnt_bin)
}

#' Identify pixel coordinates of lowres spots by computing the centroid of all highres spots inside of it.
#' The new pixel coordinates should be in the same scale as the hires image, so don't need to divide the coordinates by the bin size if i want
#' the lowres coordinates to be mapped to the image.
#' 
#' @param array_bin Vector storing the names of the lowres spots. Should be the same lenght as length of highres spots
#' @param coord_df Dataframe containing the spatial information of the highres spots
#' @param keep Vector with the indices of unique lowres spots
#' 
#' @return A list storing the col and row highres coordinates of the new lowres spots
#' 
#' @noRd
.get_bin_pxl <- function(
    array_bin,
    coord_df,
    keep
){
  # compute centroid pixel coordinates for lowres spots
  ## row
  pxl_row_bin <- tibble(grp = array_bin,
                        row = pull(coord_df, pxl_row_in_fullres)) |> 
    group_by(grp) |> 
    summarise(mn = mean(row)
    ) 
  pxl_row_bin <- pxl_row_bin[match(array_bin[keep], pxl_row_bin$grp), ] |>
    pull(mn)
  ## col
  pxl_col_bin <- tibble(grp = array_bin,
                        col = pull(coord_df, pxl_col_in_fullres)) |> 
    group_by(grp) |> 
    summarise(mn = mean(col)
    )
  pxl_col_bin <- pxl_col_bin[match(array_bin[keep], pxl_col_bin$grp), ] |>
    pull(mn)
  # assemble in list
  out <- list("pxl_row_bin" = pxl_row_bin,
              "pxl_col_bin" = pxl_col_bin)
  
  return(out)
}