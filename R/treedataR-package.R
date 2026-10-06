#' treedataR: read and write TreeData lineage tracing files
#'
#' Reads `.h5td` and `.zarr` files written by the Python `treedata` package
#' into [TreeSummarizedExperiment::TreeSummarizedExperiment-class] objects, and
#' writes them back.
#'
#' * [read_tse()] / [write_tse()]: whole files.
#' * [read_trees()]: only the trees.
#' * [obst()] / [vart()]: annotated trees as [tidytree::treedata-class] objects.
#' * [character_matrix()]: node x character matrix, including internal nodes.
#'
#' @keywords internal
#' @import methods
#' @importClassesFrom tidytree treedata
#' @importFrom rhdf5 H5Fopen
#' @importClassesFrom TreeSummarizedExperiment TreeSummarizedExperiment
"_PACKAGE"
