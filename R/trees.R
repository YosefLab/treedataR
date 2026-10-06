# Conversion between the on-disk tree representation and tidytree::treedata.
#
# Both on-disk formats are first normalised to a "columnar" tree, mirroring the
# zarr layout:
#   list(nodes      = <character node ids>,
#        edges      = <character matrix, n_edges x 2 (source, target)>,
#        node_attrs = <named list of columns; each a list, one entry per node>,
#        edge_attrs = <named list of columns; each a list, one entry per edge>)
# Within attribute columns `NULL` marks an absent attribute. Values that must
# stay JSON arrays (list-column entries) are wrapped in `I()`.

# ---------------------------------------------------------------------------
# networkx dict (HDF5 JSON) <-> columnar
# ---------------------------------------------------------------------------

.dict_to_columnar <- function(d) {
    node_dict <- d$nodes
    edge_dict <- d$edges
    if (is.null(node_dict)) node_dict <- list()
    if (is.null(edge_dict)) edge_dict <- list()
    src <- rep(names(edge_dict), lengths(edge_dict))
    dst <- unlist(lapply(edge_dict, names), use.names = FALSE)
    if (is.null(dst)) dst <- character(0)
    edges <- matrix(c(src, dst), ncol = 2)
    ids <- unique(c(names(node_dict), src, dst))
    node_dict <- node_dict[ids]
    edge_vals <- do.call(c, unname(lapply(edge_dict, unname)))
    list(
        nodes = ids,
        edges = edges,
        node_attrs = .dicts_to_columns(node_dict),
        edge_attrs = .dicts_to_columns(edge_vals)
    )
}

.dicts_to_columns <- function(dicts) {
    keys <- unique(unlist(lapply(dicts, names), use.names = FALSE))
    stats::setNames(lapply(keys, function(k) lapply(dicts, function(a) a[[k]])), keys)
}

.columnar_to_dict <- function(ct) {
    nodes <- lapply(seq_along(ct$nodes), function(i) .row_dict(ct$node_attrs, i))
    names(nodes) <- ct$nodes
    edges <- lapply(ct$nodes, function(x) .empty_dict())
    names(edges) <- ct$nodes
    for (j in seq_len(nrow(ct$edges))) {
        edges[[ct$edges[j, 1]]][[ct$edges[j, 2]]] <- .row_dict(ct$edge_attrs, j)
    }
    list(edges = edges, nodes = nodes)
}

.row_dict <- function(cols, i) {
    vals <- lapply(cols, function(col) col[[i]])
    vals <- vals[!vapply(vals, is.null, logical(1))]
    if (length(vals) == 0) .empty_dict() else vals
}

# A named empty list serialises to `{}` (an unnamed one would become `[]`).
.empty_dict <- function() stats::setNames(list(), character(0))

# ---------------------------------------------------------------------------
# columnar -> tidytree::treedata
# ---------------------------------------------------------------------------

.columnar_to_treedata <- function(ct, name = "", length_key = "length", time_key = "time") {
    ids <- ct$nodes
    E <- ct$edges
    n <- length(ids)
    if (n < 2) {
        warning("Tree '", name, "' has fewer than two nodes and cannot be represented as a phylo; ",
            "it is kept as a raw tree.",
            call. = FALSE
        )
        return(structure(ct, class = "treedataR_raw_tree"))
    }
    roots <- setdiff(ids, E[, 2])
    if (length(roots) != 1 || nrow(E) != n - 1 || anyDuplicated(E[, 2])) {
        stop("Tree '", name, "' is not a rooted tree.", call. = FALSE)
    }
    children <- split(E[, 2], factor(E[, 1], levels = ids))

    # Iterative preorder traversal (children in stored order).
    pre <- character(n)
    stack <- character(n)
    stack[1] <- roots
    top <- 1L
    k <- 0L
    while (top > 0L) {
        v <- stack[top]
        top <- top - 1L
        k <- k + 1L
        pre[k] <- v
        ch <- children[[v]]
        if (length(ch)) {
            stack[top + seq_along(ch)] <- rev(ch)
            top <- top + length(ch)
        }
    }
    if (k != n) stop("Tree '", name, "' is not connected.", call. = FALSE)

    # phylo numbering: tips 1..ntip, then internal nodes with the root first.
    is_tip <- lengths(children[pre]) == 0
    ntip <- sum(is_tip)
    num <- integer(n)
    names(num) <- pre
    num[is_tip] <- seq_len(ntip)
    num[!is_tip] <- ntip + seq_len(n - ntip)

    # Edges in cladewise (preorder) order, indexed by child.
    child <- pre[-1]
    eidx <- match(child, E[, 2])
    edge <- cbind(num[E[eidx, 1]], num[child])
    dimnames(edge) <- NULL

    phy <- list(
        edge = edge,
        Nnode = n - ntip,
        tip.label = pre[is_tip],
        node.label = pre[!is_tip]
    )
    class(phy) <- "phylo"

    # Edge lengths: `length` edge attribute, falling back to node time differences.
    lens <- rep(NA_real_, length(child))
    if (!is.null(length_key) && length_key %in% names(ct$edge_attrs)) {
        stored <- ct$edge_attrs[[length_key]][eidx]
        has <- vapply(stored, function(v) is.numeric(v) && length(v) == 1, logical(1))
        lens[has] <- as.numeric(unlist(stored[has]))
    }
    from_time <- character(0)
    if (!is.null(time_key) && time_key %in% names(ct$node_attrs) && anyNA(lens)) {
        times <- ct$node_attrs[[time_key]]
        names(times) <- ids
        tnum <- vapply(times, function(v) if (is.numeric(v) && length(v) == 1) as.numeric(v) else NA_real_, 0)
        derived <- tnum[child] - tnum[E[eidx, 1]]
        fill <- is.na(lens) & !is.na(derived)
        lens[fill] <- derived[fill]
        from_time <- child[fill]
    }
    if (!all(is.na(lens))) phy$edge.length <- lens

    # Node data: node attributes plus non-length edge attributes on the child node.
    node_order <- order(num[ids])
    data <- list(node = num[ids][node_order])
    for (k in names(ct$node_attrs)) {
        data[[k]] <- .simplify_column(ct$node_attrs[[k]][node_order])
    }
    edge_cols <- character(0)
    in_edge <- match(ids, E[, 2])[node_order]
    for (k in setdiff(names(ct$edge_attrs), length_key)) {
        col_name <- if (k %in% names(data)) paste0("edge_", k) else k
        vals <- lapply(in_edge, function(e) if (is.na(e)) NULL else ct$edge_attrs[[k]][[e]])
        data[[col_name]] <- .simplify_column(vals)
        edge_cols[col_name] <- k
    }
    data$node <- unname(data$node)
    data <- tibble::as_tibble(data)

    td <- methods::new("treedata", phylo = phy, data = data)
    td@info <- list(edge_attrs = edge_cols, length_from_time = from_time)
    td
}

# Turn a list column (NULL = absent) into an atomic vector with NAs when every
# present value is a scalar of a compatible type; otherwise keep a list column.
.simplify_column <- function(col) {
    present <- !vapply(col, is.null, logical(1))
    vals <- col[present]
    is_scalar <- vapply(vals, function(v) is.atomic(v) && length(v) == 1 && !inherits(v, "AsIs"), logical(1))
    if (length(vals) && all(is_scalar)) {
        types <- unique(vapply(vals, typeof, character(1)))
        if (all(types %in% c("integer", "double"))) types <- if ("double" %in% types) "double" else "integer"
        if (length(types) == 1) {
            out <- vector(types, length(col))
            out[] <- NA
            out[present] <- unlist(vals, use.names = FALSE)
            return(out)
        }
    }
    lapply(unname(col), .simplify_value)
}

.simplify_value <- function(v) {
    if (is.list(v) && length(v) && is.null(names(v))) {
        scalar <- vapply(v, function(x) is.atomic(x) && length(x) == 1, logical(1))
        if (all(scalar) && length(unique(vapply(v, typeof, character(1)))) == 1) {
            return(unlist(v, use.names = FALSE))
        }
    }
    if (inherits(v, "AsIs")) class(v) <- setdiff(class(v), "AsIs")
    v
}

# ---------------------------------------------------------------------------
# tidytree::treedata / phylo -> columnar
# ---------------------------------------------------------------------------

.treedata_to_columnar <- function(x) {
    if (inherits(x, "treedataR_raw_tree")) {
        return(unclass(x))
    }
    if (inherits(x, "phylo")) {
        phy <- x
        data <- NULL
        info <- list()
    } else if (methods::is(x, "treedata")) {
        phy <- x@phylo
        data <- x@data
        info <- x@info
    } else {
        stop("Trees must be tidytree::treedata or ape::phylo objects.", call. = FALSE)
    }
    ntip <- length(phy$tip.label)
    nn <- ntip + phy$Nnode
    node_label <- phy$node.label
    if (is.null(node_label)) node_label <- rep(NA_character_, phy$Nnode)
    labels <- c(phy$tip.label, node_label)
    missing <- is.na(labels) | labels == ""
    labels[missing] <- paste0("node", seq_len(nn)[missing])
    if (anyDuplicated(labels)) {
        stop("Node labels must be unique; duplicated: ",
            paste(unique(labels[duplicated(labels)]), collapse = ", "),
            call. = FALSE
        )
    }

    parent <- phy$edge[, 1]
    child <- phy$edge[, 2]
    edge_cols <- info$edge_attrs
    if (is.null(edge_cols)) edge_cols <- character(0)

    node_attrs <- list()
    edge_attrs <- list()
    if (!is.null(data) && nrow(data)) {
        rows <- match(seq_len(nn), data$node)
        for (k in setdiff(names(data), "node")) {
            col <- .column_to_values(data[[k]])
            vals <- lapply(rows, function(r) if (is.na(r)) NULL else col[[r]])
            if (k %in% names(edge_cols)) {
                edge_attrs[[edge_cols[[k]]]] <- vals[child]
            } else {
                node_attrs[[k]] <- vals
            }
        }
    }
    if (!is.null(phy$edge.length)) {
        from_time <- info$length_from_time
        lens <- lapply(seq_along(child), function(e) {
            l <- phy$edge.length[e]
            if (is.na(l) || labels[child[e]] %in% from_time) NULL else l
        })
        if (any(!vapply(lens, is.null, logical(1)))) edge_attrs[["length"]] <- lens
    }
    edge_attrs <- edge_attrs[vapply(edge_attrs, function(c) any(!vapply(c, is.null, logical(1))), logical(1))]
    node_attrs <- node_attrs[vapply(node_attrs, function(c) any(!vapply(c, is.null, logical(1))), logical(1))]

    list(
        nodes = labels,
        edges = matrix(c(labels[parent], labels[child]), ncol = 2),
        node_attrs = node_attrs,
        edge_attrs = edge_attrs
    )
}

# Convert a data column into per-row values (NULL = absent). List-column
# entries are wrapped in I() so they are written as JSON arrays.
.column_to_values <- function(col) {
    if (is.factor(col)) col <- as.character(col)
    if (is.list(col)) {
        return(lapply(col, function(v) {
            if (is.null(v) || (is.atomic(v) && length(v) == 1 && is.na(v))) NULL else I(v)
        }))
    }
    lapply(col, function(v) if (is.na(v)) NULL else v)
}
