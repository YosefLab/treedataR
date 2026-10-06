#' Read the trees from a TreeData file
#'
#' Reads only the tree-specific elements of a `.h5td` or `.zarr` file written
#' by the Python `treedata` package, without loading the expression data.
#'
#' Branch lengths are taken from the `length_key` edge attribute. Edges without
#' it fall back to the difference in the `time_key` node attribute between the
#' child and parent (set `time_key = NULL` to disable this).
#'
#' @param path Path to a `.h5td` file or `.zarr` store.
#' @param format One of `"auto"`, `"h5td"` or `"zarr"`. `"auto"` treats
#'   directories and paths ending in `.zarr` as zarr and everything else as HDF5.
#' @param length_key Edge attribute holding branch lengths.
#' @param time_key Node attribute used to derive branch lengths when
#'   `length_key` is absent.
#'
#' @return A list with elements `obst` and `vart` (named lists of
#'   [tidytree::treedata-class] objects), `label`, `alignment` and
#'   `allow_overlap`.
#' @export
read_trees <- function(path, format = c("auto", "h5td", "zarr"), length_key = "length", time_key = "time") {
    format <- .detect_format(path, match.arg(format))
    store <- .open_store(path, format)
    out <- list(obst = .named_list(), vart = .named_list(), label = NULL, alignment = "leaves", allow_overlap = FALSE)
    for (key in c("obst", "vart")) {
        if (.has_elem(store, format, key)) {
            cts <- if (format == "zarr" && .zarr_is_group(store, key)) {
                .read_columnar_zarr(store, key)
            } else {
                .parse_tree_json(.read_elem(store, format, key))
            }
            out[[key]] <- .to_treedata_list(cts, length_key, time_key)
        }
    }
    for (key in c("label", "alignment", "allow_overlap")) {
        if (.has_elem(store, format, key)) {
            value <- .read_elem(store, format, key)
            if (!is.null(value)) out[[key]] <- value
        }
    }
    out$allow_overlap <- as.logical(out$allow_overlap)
    # Files from treedata < 0.1 kept the trees as JSON in `raw.treedata`.
    if (.has_elem(store, format, "raw.treedata")) {
        legacy <- jsonlite::fromJSON(.read_elem(store, format, "raw.treedata"), simplifyVector = FALSE)
        for (key in c("obst", "vart")) {
            if (!is.null(legacy[[key]])) {
                out[[key]] <- .to_treedata_list(lapply(legacy[[key]], .dict_to_columnar), length_key, time_key)
            }
        }
        out$allow_overlap <- isTRUE(as.logical(legacy$allow_overlap))
        out$label <- legacy$label
    }
    out
}

.named_list <- function() stats::setNames(list(), character(0))

.to_treedata_list <- function(cts, length_key, time_key) {
    trees <- .named_list()
    for (name in names(cts)) {
        trees[[name]] <- .columnar_to_treedata(cts[[name]], name, length_key, time_key)
    }
    trees
}

.parse_tree_json <- function(json) {
    if (length(json) == 0 || !nzchar(json)) {
        return(.named_list())
    }
    lapply(jsonlite::fromJSON(json, simplifyVector = FALSE), .dict_to_columnar)
}

# ---------------------------------------------------------------------------
# zarr columnar layout
# ---------------------------------------------------------------------------

.read_columnar_zarr <- function(store, key) {
    cts <- .named_list()
    for (name in .zarr_group_keys(store, key)) {
        base <- paste(key, name, sep = "/")
        nodes <- as.character(.read_elem(store, "zarr", paste0(base, "/nodes")))
        edges <- .read_elem(store, "zarr", paste0(base, "/edges"))
        edges <- if (length(edges)) matrix(as.character(edges), ncol = 2) else matrix(character(0), ncol = 2)
        attrs <- lapply(c(node = "node_attrs", edge = "edge_attrs"), function(k) {
            path <- paste(base, k, sep = "/")
            if (!.has_elem(store, "zarr", path)) {
                return(list())
            }
            lapply(.read_elem(store, "zarr", path), .decode_zarr_column)
        })
        cts[[name]] <- list(nodes = nodes, edges = edges, node_attrs = attrs$node, edge_attrs = attrs$edge)
    }
    cts
}

# Native numeric/bool columns are dense; string columns hold per-element JSON
# with "" marking an absent attribute.
.decode_zarr_column <- function(values) {
    values <- as.vector(values)
    if (!is.character(values)) {
        return(as.list(values))
    }
    absent <- is.na(values) | values == ""
    parsed <- jsonlite::fromJSON(
        paste0("[", paste(ifelse(absent, "null", values), collapse = ","), "]"),
        simplifyVector = FALSE
    )
    parsed[absent] <- list(NULL)
    parsed
}

.encode_zarr_column <- function(col) {
    present <- !vapply(col, is.null, logical(1))
    native <- all(present) && length(col) > 0 && all(vapply(col, function(v) {
        (is.numeric(v) || is.logical(v)) && length(v) == 1 && !inherits(v, "AsIs") && !is.na(v)
    }, logical(1)))
    if (native) {
        types <- unique(vapply(col, typeof, character(1)))
        if (length(types) == 1 || all(types %in% c("integer", "double"))) {
            return(unlist(col, use.names = FALSE))
        }
    }
    vapply(col, function(v) if (is.null(v)) "" else .to_json(v), character(1))
}

# 17 significant digits so that doubles survive the round trip exactly.
.to_json <- function(x) {
    as.character(jsonlite::toJSON(x, auto_unbox = TRUE, digits = I(17), null = "null", always_decimal = TRUE))
}

# ---------------------------------------------------------------------------
# Writing
# ---------------------------------------------------------------------------

.write_trees <- function(store, format, key, trees) {
    cts <- lapply(trees, .treedata_to_columnar)
    if (format == "h5td") {
        json <- if (length(cts)) .to_json(lapply(cts, .columnar_to_dict)) else "{}"
        .write_elem(json, store, format, key)
        return(invisible())
    }
    groups <- lapply(cts, function(ct) {
        g <- list(nodes = ct$nodes, edges = ct$edges)
        if (length(ct$node_attrs)) g$node_attrs <- lapply(ct$node_attrs, .encode_zarr_column)
        if (length(ct$edge_attrs)) g$edge_attrs <- lapply(ct$edge_attrs, .encode_zarr_column)
        g
    })
    .write_elem(if (length(groups)) groups else .named_list(), store, format, key)
    invisible()
}
