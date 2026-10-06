library(TreeSummarizedExperiment)

for (fmt in formats) {
    for (name in c("forest", "nodes")) {
        test_that(paste("round trip through R:", name, fmt), {
            tse <- read_tse(fixture(paste0(name, ".", fmt)))
            out <- file.path(withr::local_tempdir(), paste0("out.", fmt))
            write_tse(tse, out)
            back <- read_tse(out)
            expect_equal(tree_summary(obst(back)), tree_summary(obst(tse)))
            expect_equal(tree_summary(vart(back)), tree_summary(vart(tse)))
            expect_equal(colLinks(back), colLinks(tse))
            expect_equal(as.matrix(assay(back, "X")), as.matrix(assay(tse, "X")))
            expect_equal(as.data.frame(colData(back)), as.data.frame(colData(tse)))
            expect_equal(reducedDimNames(back), reducedDimNames(tse))
            if ("characters" %in% reducedDimNames(tse)) {
                expect_equal(reducedDim(back, "characters"), reducedDim(tse, "characters"))
            }
            expect_equal(S4Vectors::metadata(back)$treedata[c("label", "alignment", "allow_overlap")],
                S4Vectors::metadata(tse)$treedata[c("label", "alignment", "allow_overlap")])
        })
    }
}

test_that("write_tse refuses to overwrite unless asked", {
    tse <- read_tse(fixture("forest.h5td"))
    out <- file.path(withr::local_tempdir(), "out.h5td")
    write_tse(tse, out)
    expect_error(write_tse(tse, out), "already exists")
    expect_silent(write_tse(tse, out, overwrite = TRUE))
})

test_that("subsetting relaxes the alignment", {
    tse <- read_tse(fixture("forest.h5td"))[, 1:4]
    out <- file.path(withr::local_tempdir(), "out.h5td")
    expect_message(write_tse(tse, out), "subset")
    expect_equal(read_trees(out)$alignment, "subset")
})

test_that("trees dropped from colTree by subsetting are not written", {
    tse <- read_tse(fixture("forest.h5td"))[, 1:5]
    out <- file.path(withr::local_tempdir(), "out.h5td")
    expect_no_warning(expect_message(write_tse(tse, out), "'t3'"))
    trees <- read_trees(out)
    expect_named(trees$obst, c("t1", "t2"))
    expect_equal(trees$alignment, "leaves")
})

test_that("a TSE built from scratch can be written", {
    phy <- ape::read.tree(text = "((c1:1,c2:2)a:1,c3:1)r;")
    m <- matrix(1:6, 2, 3, dimnames = list(c("g1", "g2"), c("c1", "c2", "c3")))
    tse <- TreeSummarizedExperiment(assays = list(X = m), colTree = phy)
    out <- file.path(withr::local_tempdir(), "out.h5td")
    write_tse(tse, out)
    back <- read_tse(out)
    expect_equal(back@colTree[[1]]$edge.length, phy$edge.length)
    expect_equal(obst(back)[[1]]@phylo$tip.label, phy$tip.label)
})

test_that("unchanged trees write without a consistency warning", {
    tse <- read_tse(fixture("forest.h5td"))
    out <- file.path(withr::local_tempdir(), "out.h5td")
    expect_no_warning(write_tse(tse, out))
})

test_that("edited colTree lengths win over obst() and keep annotations", {
    tse <- read_tse(fixture("forest.h5td"))
    phy <- tse@colTree$t1
    phy$edge.length <- phy$edge.length * 10
    tse@colTree$t1 <- phy
    out <- file.path(withr::local_tempdir(), "out.h5td")
    expect_warning(write_tse(tse, out), "'t1'")
    t1 <- read_trees(out)$obst$t1
    expect_equal(t1@phylo$edge.length, phy$edge.length)
    a1 <- t1@data[t1@data$node == match("a1", c(t1@phylo$tip.label, t1@phylo$node.label)), ]
    expect_equal(a1$time, 1)
    expect_equal(a1$characters[[1]], c("1", "0", "0"))
})

test_that("time-derived lengths stay unwritten only where unchanged", {
    tse <- read_tse(fixture("forest.h5td"))
    phy <- tse@colTree$t3
    child <- c(phy$tip.label, phy$node.label)[phy$edge[, 2]]
    phy$edge.length[child == "c6"] <- 5
    tse@colTree$t3 <- phy
    out <- file.path(withr::local_tempdir(), "out.h5td")
    expect_warning(write_tse(tse, out), "'t3'")
    lens <- jsonlite::fromJSON(read_trees_json(out))$t3$edges$a3
    expect_equal(lens$c6$length, 5)
    expect_null(lens$c7$length)
})

test_that("dropped tips in colTree are dropped from the written tree", {
    tse <- read_tse(fixture("forest.h5td"))
    tse@colTree$t1 <- ape::drop.tip(tse@colTree$t1, "c3", collapse.singles = FALSE)
    out <- file.path(withr::local_tempdir(), "out.h5td")
    expect_warning(write_tse(tse, out), "'t1'")
    t1 <- read_trees(out)$obst$t1
    expect_false("c3" %in% t1@phylo$tip.label)
    labels <- c(t1@phylo$tip.label, t1@phylo$node.label)
    expect_equal(t1@data$depth[t1@data$node == match("a1", labels)], 1L)
})
