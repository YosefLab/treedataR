fixture <- function(name) {
    path <- system.file("extdata", name, package = "treedataR")
    if (!nzchar(path)) path <- test_path("..", "..", "inst", "extdata", name)
    path
}

formats <- c("h5td", "zarr")

# Trees reduced to comparable plain objects (column order is not meaningful).
tree_summary <- function(trees) {
    lapply(trees, function(t) {
        d <- as.data.frame(t@data)
        list(phylo = t@phylo, data = d[order(names(d))])
    })
}

read_trees_json <- function(path) {
    rhdf5::h5read(path, "obst")
}
