library(TreeSummarizedExperiment)

for (fmt in formats) {
    test_that(paste("forest fixture reads as a TSE:", fmt), {
        tse <- read_tse(fixture(paste0("forest.", fmt)))
        expect_s4_class(tse, "TreeSummarizedExperiment")
        expect_equal(dim(tse), c(4L, 7L))
        expect_equal(colTreeNames(tse), c("t1", "t2", "t3"))
        expect_equal(colLinks(tse)$whichTree, c("t1", "t1", "t1", "t2", "t2", "t3", "t3"))
        expect_true(all(colLinks(tse)$isLeaf))
        expect_null(rowTree(tse))

        trees <- obst(tse)
        expect_named(trees, c("t1", "t2", "t3"))
        t1 <- trees$t1
        expect_equal(t1@phylo$edge.length, c(1, 2, 2, 3))
        expect_equal(trees$t2@phylo$edge.length, c(1.5, 2.5))
        expect_equal(trees$t3@phylo$edge.length, c(1.25, 1, 2.5))
        expect_equal(trees$t3@info$length_from_time, c("c6", "c7"))

        # Internal-node annotations.
        internal <- character_matrix(t1, "internal")
        expect_equal(rownames(internal), c("r1", "a1"))
        expect_equal(unname(internal["a1", ]), c("1", "0", "0"))
        data <- as.data.frame(t1@data)
        a1 <- data[data$node == match("a1", c(t1@phylo$tip.label, t1@phylo$node.label)), ]
        expect_equal(a1$time, 1)
        expect_equal(a1$depth, 1L)
        expect_true(a1$stump)
        expect_equal(a1$n_edits, 1L)

        chars <- reducedDim(tse, "characters")
        expect_equal(dim(chars), c(7L, 3L))
        expect_equal(colnames(chars), c("site1", "site2", "site3"))
        expect_equal(unname(chars["c2", ]), c("1", "3", "-"))
        expect_equal(chars[rownames(chars) %in% t1@phylo$tip.label, ],
            character_matrix(t1, "leaves")[c("c1", "c2", "c3"), ],
            ignore_attr = TRUE
        )

        info <- S4Vectors::metadata(tse)$treedata
        expect_equal(info$alignment, "leaves")
        expect_equal(info$label, "tree")
        expect_true(is.logical(info$allow_overlap))
    })

    test_that(paste("nodes-aligned fixture links internal nodes and vart:", fmt), {
        tse <- read_tse(fixture(paste0("nodes.", fmt)))
        links <- colLinks(tse)
        expect_equal(links$nodeLab, colnames(tse))
        expect_equal(links$isLeaf[match(c("root", "m1"), links$nodeLab)], c(FALSE, FALSE))
        expect_equal(rowTreeNames(tse), "genes")
        expect_false(rowLinks(tse)$isLeaf[rownames(tse) == "gr"])
        expect_equal(S4Vectors::metadata(tse)$treedata$alignment, "nodes")
    })
}

test_that("read_trees returns only the trees", {
    trees <- read_trees(fixture("forest.h5td"))
    expect_named(trees, c("obst", "vart", "label", "alignment", "allow_overlap"))
    expect_length(trees$vart, 0)
    expect_s4_class(trees$obst$t1, "treedata")
})

test_that("cells outside every tree are kept with NA links", {
    tse <- read_tse(fixture("forest.h5td"))
    obst(tse) <- obst(tse)[c("t1", "t3")]
    expect_equal(ncol(tse), 7L)
    expect_equal(colLinks(tse)$whichTree, c("t1", "t1", "t1", NA, NA, "t3", "t3"))
})

test_that("plain phylo trees can be assigned", {
    tse <- read_tse(fixture("forest.h5td"))
    obst(tse) <- list(t = ape::read.tree(text = "((c1,c2),c3);"))
    expect_equal(colTreeNames(tse), "t")
    expect_s4_class(obst(tse)$t, "treedata")
})
