# Thin wrappers around anndataR's element-level readers and writers.
#
# anndataR does not export these, but they implement the AnnData on-disk
# encodings for both HDF5 and zarr. Every use goes through this file so there is
# a single place to update if the internals change.

.anndataR <- function(name) get(name, envir = asNamespace("anndataR"), inherits = FALSE)

`%||%` <- function(x, y) if (is.null(x)) y else x

.detect_format <- function(path, format = c("auto", "h5td", "zarr")) {
    format <- match.arg(format)
    if (format != "auto") {
        return(format)
    }
    if (dir.exists(path) || grepl("\\.zarr/?$", path)) "zarr" else "h5td"
}

# A handle for element-level I/O: an anndataR HDF5File for HDF5, the store
# path for zarr.
.open_store <- function(path, format) {
    if (format == "zarr") normalizePath(path, mustWork = TRUE) else .anndataR("HDF5File")$new(path)
}

.has_elem <- function(store, format, name) {
    if (format == "zarr") .anndataR("zarr_path_exists")(store, name) else .anndataR("hdf5_path_exists")(store, name)
}

.read_elem <- function(store, format, name) {
    reader <- if (format == "zarr") "read_zarr_element" else "read_h5ad_element"
    .quietly(.anndataR(reader)(store, name, stop_on_error = TRUE))
}

.write_elem <- function(value, store, format, name) {
    if (format == "zarr") {
        .quietly(.anndataR("write_zarr_element")(value, store, name, zarr_format = 3, stop_on_error = TRUE))
    } else {
        .anndataR("write_h5ad_element")(value, store, name, stop_on_error = TRUE)
    }
}

.write_scalar <- function(value, store, format, name) {
    .write_elem(value, store, format, name)
    if (format == "zarr") .zarr_make_0d(store, name)
}

# anndataR writes zarr scalars as shape-[1] arrays, but the AnnData spec (and
# Python anndata) expect 0-d arrays; otherwise Python reads e.g. "['leaves']".
# The single chunk holds the same bytes either way, so only the metadata and
# the chunk key change.
.zarr_make_0d <- function(store, name) {
    meta_path <- file.path(store, name, "zarr.json")
    meta <- jsonlite::read_json(meta_path)
    if (length(meta$shape) != 1) {
        return(invisible())
    }
    meta$shape <- list()
    meta$chunk_grid$configuration$chunk_shape <- list()
    jsonlite::write_json(meta, meta_path, auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA)
    chunk <- file.path(store, name, "c", "0")
    if (file.exists(chunk)) {
        tmp <- tempfile(tmpdir = file.path(store, name))
        file.rename(chunk, tmp)
        unlink(file.path(store, name, "c"), recursive = TRUE)
        file.rename(tmp, file.path(store, name, "c"))
    }
    invisible()
}

.write_root_encoding <- function(store, format) {
    if (format == "zarr") {
        .anndataR("write_zarr_encoding")(store, "", "treedata", "0.1.0", zarr_format = 3)
    } else {
        .anndataR("write_h5ad_encoding")(store, "/", "treedata", "0.1.0")
    }
}

.zarr_group_keys <- function(store, name) .anndataR("read_zarr_mapping_keys")(store, name)

.zarr_is_group <- function(store, name) {
    meta <- file.path(store, name, "zarr.json")
    if (file.exists(meta)) {
        return(identical(jsonlite::read_json(meta)$node_type, "group"))
    }
    file.exists(file.path(store, name, ".zgroup"))
}

# Rarr warns about zstd support whenever extSoftVersion() does not list zstd,
# even when memDecompress() handles it fine. Drop only that warning.
.quietly <- function(expr) {
    withCallingHandlers(expr, warning = function(w) {
        if (grepl("memCompress() and memDecompress()", conditionMessage(w), fixed = TRUE)) {
            invokeRestart("muffleWarning")
        }
    })
}
