dict <- function(json) .dict_to_columnar(jsonlite::fromJSON(json, simplifyVector = FALSE))

test_that("networkx dicts become phylo trees with node data", {
    ct <- dict('{"nodes": {"r": {"time": 0}, "a": {"time": 1.5, "characters": ["1", "0"]}, "b": {"time": 2}},
                 "edges": {"r": {"a": {"length": 1.5, "n_edits": 2}, "b": {}}, "a": {}, "b": {}}}')
    td <- .columnar_to_treedata(ct, "t")
    phy <- td@phylo
    expect_s3_class(phy, "phylo")
    expect_equal(phy$tip.label, c("a", "b"))
    expect_equal(phy$node.label, "r")
    expect_equal(phy$edge, rbind(c(3L, 1L), c(3L, 2L)))
    # r->a has a stored length; r->b falls back to time[b] - time[r].
    expect_equal(phy$edge.length, c(1.5, 2))
    expect_equal(td@info$length_from_time, "b")
    expect_equal(td@data$time, c(1.5, 2, 0))
    expect_equal(td@data$characters, list(c("1", "0"), NULL, NULL))
    expect_equal(td@data$n_edits, c(2L, NA, NA))
    expect_equal(td@info$edge_attrs, c(n_edits = "n_edits"))
})

test_that("time fallback can be disabled", {
    ct <- dict('{"nodes": {"r": {"time": 0}, "a": {"time": 1}, "b": {"time": 2}},
                 "edges": {"r": {"a": {}, "b": {}}}}')
    expect_null(.columnar_to_treedata(ct, "t", time_key = NULL)@phylo$edge.length)
    expect_equal(.columnar_to_treedata(ct, "t")@phylo$edge.length, c(1, 2))
})

test_that("edge attributes clashing with node attributes are prefixed", {
    ct <- dict('{"nodes": {"r": {"w": 1}, "a": {"w": 2}, "b": {}},
                 "edges": {"r": {"a": {"w": 9}, "b": {}}}}')
    td <- .columnar_to_treedata(ct, "t")
    expect_equal(td@info$edge_attrs, c(edge_w = "w"))
    expect_equal(td@data$edge_w, c(9L, NA, NA))
})

test_that("unary nodes and deep trees are supported", {
    n <- 2000
    ids <- paste0("n", seq_len(n))
    ct <- list(nodes = ids, edges = cbind(ids[-n], ids[-1]), node_attrs = list(), edge_attrs = list())
    td <- .columnar_to_treedata(ct, "chain")
    expect_equal(td@phylo$Nnode, n - 1)
    expect_equal(td@phylo$tip.label, ids[n])
})

test_that("degenerate and invalid trees are handled", {
    single <- list(nodes = "x", edges = matrix(character(0), ncol = 2), node_attrs = list(), edge_attrs = list())
    expect_warning(raw <- .columnar_to_treedata(single, "s"), "fewer than two nodes")
    expect_s3_class(raw, "treedataR_raw_tree")
    expect_equal(.treedata_to_columnar(raw)$nodes, "x")

    two_roots <- dict('{"nodes": {"a": {}, "b": {}, "c": {}}, "edges": {"a": {"c": {}}}}')
    expect_error(.columnar_to_treedata(two_roots, "bad"), "not a rooted tree")
})

test_that("treedata converts back to the same columnar tree", {
    ct <- dict('{"nodes": {"r": {"time": 0.0, "characters": ["0"]}, "a": {"time": 1.0, "flag": true}, "b": {"time": 3.0}},
                 "edges": {"r": {"a": {"length": 1.0, "n": 1}, "b": {}}}}')
    back <- .treedata_to_columnar(.columnar_to_treedata(ct, "t"))
    expect_setequal(back$nodes, ct$nodes)
    d1 <- .columnar_to_dict(ct)
    d2 <- .columnar_to_dict(back)
    expect_equal(.to_json(d2$nodes[names(d1$nodes)]), .to_json(d1$nodes))
    # The time-derived length of r->b is not written back.
    expect_equal(d2$edges$r$b, .empty_dict())
    expect_equal(d2$edges$r$a, list(n = 1L, length = 1))
})

test_that("plain phylo objects are converted with generated labels", {
    phy <- ape::read.tree(text = "((a:1,b:2):1,c:3);")
    ct <- .treedata_to_columnar(phy)
    expect_equal(ct$nodes, c("a", "b", "c", "node4", "node5"))
    expect_equal(unlist(ct$edge_attrs$length), phy$edge.length)
})

test_that("JSON keeps list values as arrays and doubles as floats", {
    ct <- list(
        nodes = c("r", "a"), edges = matrix(c("r", "a"), ncol = 2),
        node_attrs = list(characters = list(I("1"), NULL), time = list(0, 1)), edge_attrs = list()
    )
    json <- .to_json(.columnar_to_dict(ct))
    expect_match(json, '"characters":["1"]', fixed = TRUE)
    expect_match(json, '"time":0.0', fixed = TRUE)
    x <- 0.6043894636567326
    expect_identical(jsonlite::fromJSON(.to_json(x)), x)
    expect_match(json, '"a":{}', fixed = TRUE)
})

test_that("zarr attribute columns are encoded like Python treedata", {
    expect_equal(.encode_zarr_column(list(1, 2.5)), c(1, 2.5))
    expect_equal(.encode_zarr_column(list(TRUE, FALSE)), c(TRUE, FALSE))
    expect_equal(.encode_zarr_column(list(1, NULL)), c("1.0", ""))
    expect_equal(.encode_zarr_column(list("x", I(c("1", "2")))), c('"x"', '["1","2"]'))
    expect_equal(.decode_zarr_column(c('"x"', "", '["1","2"]')), list("x", NULL, list("1", "2")))
})
