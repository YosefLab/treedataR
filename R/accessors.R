#' Annotated trees of a TreeSummarizedExperiment
#'
#' `obst()` and `vart()` return the observation (column) and variable (row)
#' trees as a named list of [tidytree::treedata-class] objects. Unlike `colTree()`
#' and `rowTree()`, these keep every node and edge annotation, for example
#' ancestral `characters` and `time` on internal nodes.
#'
#' Assigning trees also updates `colTree()`/`colLinks()` (or `rowTree()`/
#' `rowLinks()`). Each column or row is linked to the tree named in its
#' `label` column if there is one (see [read_tse()]), otherwise to the first
#' tree that contains it.
#'
#' @param x A [TreeSummarizedExperiment::TreeSummarizedExperiment-class].
#' @param value A named list of [tidytree::treedata-class] or [ape::phylo]
#'   objects.
#'
#' @return A named list of trees.
#' @export
obst <- function(x) .get_trees(x, "obst")

#' @rdname obst
#' @export
`obst<-` <- function(x, value) .set_trees(x, "obst", value)

#' @rdname obst
#' @export
vart <- function(x) .get_trees(x, "vart")

#' @rdname obst
#' @export
`vart<-` <- function(x, value) .set_trees(x, "vart", value)

.get_trees <- function(x, key) {
    trees <- S4Vectors::metadata(x)$treedata[[key]]
    if (!is.null(trees)) {
        return(trees)
    }
    # Fall back to the plain phylo trees of a TSE not created by read_tse().
    phylos <- if (key == "obst") x@colTree else x@rowTree
    if (is.null(phylos)) .named_list() else lapply(phylos, function(p) methods::new("treedata", phylo = p))
}

.set_trees <- function(x, key, value) {
    if (is.null(value)) value <- .named_list()
    if (length(value) && (is.null(names(value)) || any(names(value) == ""))) {
        stop("Trees must be a named list.", call. = FALSE)
    }
    value <- lapply(value, function(t) if (inherits(t, "phylo")) methods::new("treedata", phylo = t) else t)
    md <- S4Vectors::metadata(x)
    if (is.null(md$treedata)) md$treedata <- list()
    md$treedata[[key]] <- value
    S4Vectors::metadata(x) <- md

    phylos <- Filter(Negate(is.null), lapply(value, function(t) if (methods::is(t, "treedata")) t@phylo))
    on_col <- key == "obst"
    ids <- if (on_col) colnames(x) else rownames(x)
    label <- md$treedata$label
    meta <- if (on_col) SummarizedExperiment::colData(x) else SummarizedExperiment::rowData(x)
    assigned <- if (!is.null(label) && label %in% names(meta)) as.character(meta[[label]]) else NULL
    links <- if (length(phylos)) .build_links(ids, phylos, assigned) else NULL
    if (on_col) {
        x@colTree <- if (length(phylos)) phylos else NULL
        x@colLinks <- links
    } else {
        x@rowTree <- if (length(phylos)) phylos else NULL
        x@rowLinks <- links
    }
    methods::validObject(x)
    x
}

# Link each id to a node in one of `phylos`, preferring the tree named in
# `assigned`. Ids found in no tree get NA links.
.build_links <- function(ids, phylos, assigned = NULL) {
    n <- length(ids)
    node_lab <- rep(NA_character_, n)
    node_num <- rep(NA_integer_, n)
    is_leaf <- rep(NA, n)
    which_tree <- rep(NA_character_, n)
    for (tree in names(phylos)) {
        phy <- phylos[[tree]]
        hit <- match(ids, c(phy$tip.label, phy$node.label))
        preferred <- if (is.null(assigned)) FALSE else !is.na(assigned) & assigned == tree
        take <- !is.na(hit) & (is.na(which_tree) | preferred)
        node_lab[take] <- ids[take]
        node_num[take] <- hit[take]
        is_leaf[take] <- hit[take] <= length(phy$tip.label)
        which_tree[take] <- tree
    }
    TreeSummarizedExperiment::LinkDataFrame(
        nodeLab = node_lab,
        nodeLab_alias = ifelse(is.na(node_num), NA_character_, paste0("alias_", node_num)),
        nodeNum = node_num,
        isLeaf = is_leaf,
        whichTree = which_tree
    )
}

#' Character matrix of a tree
#'
#' Builds a node x character matrix from a list-valued node attribute (by
#' default `characters`). Internal nodes are included, which gives access to
#' ancestral states reconstructed in Python.
#'
#' @param tree A [tidytree::treedata-class] object, e.g. an element of [obst()].
#' @param nodes Which nodes to return: `"all"`, `"leaves"` or `"internal"`.
#' @param key Name of the node attribute holding the characters.
#'
#' @return A matrix with one row per node (named by node label) and one column
#'   per character. Nodes without the attribute have `NA` rows.
#' @export
#' @examples
#' path <- system.file("extdata", "forest.h5td", package = "treedataR")
#' tse <- read_tse(path)
#' character_matrix(obst(tse)$t1, nodes = "internal")
character_matrix <- function(tree, nodes = c("all", "leaves", "internal"), key = "characters") {
    nodes <- match.arg(nodes)
    phy <- tree@phylo
    data <- tree@data
    if (!key %in% names(data)) stop("Tree has no '", key, "' node attribute.", call. = FALSE)
    ntip <- length(phy$tip.label)
    keep <- switch(nodes,
        all = seq_len(ntip + phy$Nnode),
        leaves = seq_len(ntip),
        internal = ntip + seq_len(phy$Nnode)
    )
    labels <- c(phy$tip.label, phy$node.label %||% rep(NA_character_, phy$Nnode))
    col <- data[[key]][match(keep, data$node)]
    if (!is.list(col)) col <- as.list(col)
    width <- max(0L, lengths(col))
    rows <- lapply(col, function(v) {
        v <- unlist(v, use.names = FALSE)
        if (is.null(v)) v <- NA
        length(v) <- width
        v
    })
    mat <- do.call(rbind, rows)
    if (is.null(mat)) mat <- matrix(NA, 0, width)
    rownames(mat) <- labels[keep]
    mat
}
