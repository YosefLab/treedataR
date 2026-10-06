#' Write a TreeSummarizedExperiment to a TreeData file
#'
#' Writes a `.h5td` or `.zarr` file that can be read with the Python
#' `treedata` package (`treedata.read_h5td()` / `treedata.read_zarr()`). The
#' AnnData part is written with anndataR. The `"X"` assay becomes `X` and the
#' other assays become `layers`.
#'
#' Trees come from [obst()] and [vart()], so node and edge annotations are
#' preserved. For a TSE not created by [read_tse()], the plain `colTree()` and
#' `rowTree()` are written instead. Branch lengths that were derived from node
#' times when reading are not written back as `length` edge attributes.
#'
#' If `colTree()`/`rowTree()` were changed after reading (e.g. with
#' `TreeSummarizedExperiment::changeTree()`) so that they no longer match
#' [obst()]/[vart()], a warning is given and the `colTree()`/`rowTree()`
#' topology and branch lengths are written. Annotations are kept for nodes
#' whose labels still exist. Trees that TreeSummarizedExperiment dropped from
#' `colTree()`/`rowTree()` (e.g. because subsetting removed all their cells)
#' are not written.
#'
#' If the trees contain observations that are no longer in the object (e.g.
#' after subsetting), `alignment` is relaxed to `"subset"`. If an observation
#' appears in several trees, `allow_overlap` is set to `TRUE`.
#'
#' @param x A [TreeSummarizedExperiment::TreeSummarizedExperiment-class].
#' @param path Output path. A path ending in `.zarr` is written as zarr unless
#'   `format` says otherwise.
#' @param format One of `"auto"`, `"h5td"` or `"zarr"`.
#' @param characters_key Name of the `reducedDim` holding the character
#'   matrix. It is written to `obsm` as a dataframe so it can hold strings.
#' @param overwrite Whether to replace an existing file.
#'
#' @return `path`, invisibly.
#' @export
#' @examples
#' path <- system.file("extdata", "forest.h5td", package = "treedataR")
#' tse <- read_tse(path)
#' out <- tempfile(fileext = ".h5td")
#' write_tse(tse, out)
write_tse <- function(x, path, format = c("auto", "h5td", "zarr"), characters_key = "characters",
                      overwrite = FALSE) {
    if (!methods::is(x, "TreeSummarizedExperiment")) {
        stop("`x` must be a TreeSummarizedExperiment.", call. = FALSE)
    }
    format <- match.arg(format)
    if (format == "auto") format <- if (grepl("\\.zarr/?$", path)) "zarr" else "h5td"
    if (file.exists(path)) {
        if (!overwrite) stop("'", path, "' already exists. Use `overwrite = TRUE` to replace it.", call. = FALSE)
        unlink(path, recursive = TRUE)
    }

    info <- S4Vectors::metadata(x)$treedata %||% list()
    obst <- .reconcile_trees(obst(x), x@colTree, "obst", "colTree")
    vart <- .reconcile_trees(vart(x), x@rowTree, "vart", "rowTree")
    alignment <- info$alignment %||% "leaves"
    allow_overlap <- isTRUE(info$allow_overlap)
    label <- info$label %||% "tree"

    sce <- methods::as(x, "SingleCellExperiment")
    S4Vectors::metadata(sce)$treedata <- NULL
    chars <- NULL
    if (!is.null(characters_key) && characters_key %in% SingleCellExperiment::reducedDimNames(sce)) {
        chars <- SingleCellExperiment::reducedDim(sce, characters_key)
        SingleCellExperiment::reducedDim(sce, characters_key) <- NULL
    }
    assays <- SummarizedExperiment::assayNames(sce)
    layers <- stats::setNames(setdiff(assays, "X"), setdiff(assays, "X"))
    adata <- anndataR::as_AnnData(
        sce,
        x_mapping = if ("X" %in% assays) "X" else NULL,
        layers_mapping = if (length(layers)) layers else FALSE
    )
    .quietly(if (format == "zarr") anndataR::write_zarr(adata, path) else anndataR::write_h5ad(adata, path))

    alignment <- .check_alignment(c(obst, vart), c(colnames(x), rownames(x)), alignment)
    if (!allow_overlap && (.has_overlap(obst) || .has_overlap(vart))) {
        message("Trees share nodes; setting `allow_overlap = TRUE`.")
        allow_overlap <- TRUE
    }

    store <- .open_store(path, format)
    if (!is.null(chars)) {
        chars <- as.data.frame(chars, stringsAsFactors = FALSE)
        rownames(chars) <- colnames(x)
        .write_elem(chars, store, format, paste0("obsm/", characters_key))
    }
    .write_scalar(label, store, format, "label")
    .write_scalar(allow_overlap, store, format, "allow_overlap")
    .write_scalar(alignment, store, format, "alignment")
    .write_trees(store, format, "obst", obst)
    .write_trees(store, format, "vart", vart)
    .write_root_encoding(store, format)
    invisible(path)
}

# The annotated trees in metadata and the phylo trees in colTree()/rowTree()
# can drift apart if the latter are edited with TreeSummarizedExperiment
# functions. colTree()/rowTree() win: annotations are carried over to nodes
# whose labels still exist. Raw (degenerate) trees only live in metadata and are
# kept as is.
.reconcile_trees <- function(annotated, phylos, key, slot) {
    raw <- Filter(function(t) inherits(t, "treedataR_raw_tree"), annotated)
    annotated <- Filter(function(t) !inherits(t, "treedataR_raw_tree"), annotated)
    if (is.null(phylos)) phylos <- .named_list()
    # TreeSummarizedExperiment drops trees that no longer have linked columns
    # or rows (e.g. after subsetting), so a missing tree was removed on purpose.
    dropped <- setdiff(names(annotated), names(phylos))
    if (length(dropped)) {
        message(sprintf("Not writing %s tree(s) %s, which are no longer in %s().",
            key, paste0("'", dropped, "'", collapse = ", "), slot))
        annotated <- annotated[setdiff(names(annotated), dropped)]
    }
    differ <- Filter(function(n) !.same_topology(annotated[[n]], phylos[[n]]), names(phylos))
    if (!length(differ)) {
        return(c(annotated, raw))
    }
    warning(sprintf(
        "%s() and %s() disagree for tree(s) %s; writing the %s() topology and keeping annotations of matching nodes.",
        key, slot, paste0("'", differ, "'", collapse = ", "), slot
    ), call. = FALSE)
    trees <- lapply(stats::setNames(names(phylos), names(phylos)), function(n) {
        if (n %in% differ) .transfer_annotations(annotated[[n]], phylos[[n]]) else annotated[[n]]
    })
    c(trees, raw)
}

.node_labels <- function(phy) c(phy$tip.label, phy$node.label %||% rep(NA_character_, phy$Nnode))

# Edges keyed by "parent->child" label, with their lengths.
.edge_table <- function(phy) {
    labels <- .node_labels(phy)
    lens <- phy$edge.length %||% rep(NA_real_, nrow(phy$edge))
    stats::setNames(lens, paste(labels[phy$edge[, 1]], labels[phy$edge[, 2]], sep = "->"))
}

.same_topology <- function(td, phy) {
    if (is.null(td) || is.null(phy)) {
        return(FALSE)
    }
    a <- .edge_table(td@phylo)
    b <- .edge_table(phy)
    setequal(names(a), names(b)) && isTRUE(all.equal(a[names(b)], b))
}

.transfer_annotations <- function(td, phy) {
    new <- methods::new("treedata", phylo = phy)
    if (is.null(td)) {
        return(new)
    }
    old_labels <- .node_labels(td@phylo)
    new_labels <- .node_labels(phy)
    data <- td@data
    data$node <- match(old_labels[data$node], new_labels)
    new@data <- data[!is.na(data$node), ]
    # Only lengths that are still the time-derived values stay unwritten.
    old_lens <- .edge_table(td@phylo)
    new_lens <- .edge_table(phy)
    old_child <- sub(".*->", "", names(old_lens))
    new_child <- sub(".*->", "", names(new_lens))
    from_time <- td@info$length_from_time
    from_time <- from_time[vapply(from_time, function(l) {
        l %in% new_child && isTRUE(all.equal(old_lens[[match(l, old_child)]], new_lens[[match(l, new_child)]]))
    }, logical(1))]
    new@info <- list(edge_attrs = td@info$edge_attrs, length_from_time = from_time)
    new
}

.tree_nodes <- function(tree, leaves_only = FALSE) {
    ct <- .treedata_to_columnar(tree)
    if (leaves_only) setdiff(ct$nodes, ct$edges[, 1]) else ct$nodes
}

.check_alignment <- function(trees, names, alignment) {
    if (alignment == "subset" || !length(trees)) {
        return(alignment)
    }
    nodes <- unlist(lapply(trees, .tree_nodes, leaves_only = alignment == "leaves"), use.names = FALSE)
    if (!all(nodes %in% names)) {
        message("Some tree ", if (alignment == "leaves") "leaves" else "nodes",
            " are not in the object; setting `alignment = \"subset\"`.")
        return("subset")
    }
    alignment
}

.has_overlap <- function(trees) {
    if (length(trees) < 2) {
        return(FALSE)
    }
    anyDuplicated(unlist(lapply(trees, .tree_nodes), use.names = FALSE)) > 0
}
