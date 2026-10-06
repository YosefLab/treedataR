#' Read a TreeData file as a TreeSummarizedExperiment
#'
#' Reads a `.h5td` or `.zarr` file written by the Python `treedata` package.
#' The AnnData part is read with anndataR (`X` becomes the `"X"` assay, `obs`
#' becomes `colData`, `obsm` becomes `reducedDims`, ...).
#'
#' Trees in `obst` become the `colTree`s and trees in `vart` the `rowTree`s.
#' Each cell is linked to the tree named in its `obs[[label]]` column, or else
#' to the first tree that contains it. Cells that are not in any tree are kept,
#' with `NA` links. The full trees, including all node and edge annotations, are
#' kept as [tidytree::treedata-class] objects and are available through [obst()]
#' and [vart()].
#'
#' @inheritParams read_trees
#' @param characters_key Name of the `obsm` element holding the character
#'   matrix. It is stored as `reducedDim(tse, characters_key)`. Use `NULL` to
#'   skip it.
#'
#' @return A [TreeSummarizedExperiment::TreeSummarizedExperiment-class] object.
#' @export
#' @examples
#' path <- system.file("extdata", "forest.h5td", package = "treedataR")
#' tse <- read_tse(path)
#' obst(tse)
read_tse <- function(path, format = c("auto", "h5td", "zarr"), characters_key = "characters",
                     length_key = "length", time_key = "time") {
    format <- .detect_format(path, match.arg(format))
    if (!file.exists(path)) stop("File '", path, "' does not exist.", call. = FALSE)
    sce <- .quietly(if (format == "zarr") {
        anndataR::read_zarr(path, as = "SingleCellExperiment")
    } else {
        anndataR::read_h5ad(path, as = "SingleCellExperiment")
    })
    trees <- read_trees(path, format, length_key = length_key, time_key = time_key)

    if (!is.null(characters_key)) {
        store <- .open_store(path, format)
        key <- paste0("obsm/", characters_key)
        if (.has_elem(store, format, key)) {
            chars <- as.matrix(.read_elem(store, format, key))
            rownames(chars) <- colnames(sce)
            SingleCellExperiment::reducedDim(sce, characters_key) <- chars
        }
    }

    tse <- methods::new("TreeSummarizedExperiment", sce)
    S4Vectors::metadata(tse)$treedata <- trees[c("label", "alignment", "allow_overlap")]
    obst(tse) <- trees$obst
    vart(tse) <- trees$vart
    tse
}
