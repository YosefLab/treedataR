# treedataR

Read and write [treedata](https://github.com/YosefLab/treedata) lineage tracing
files (`.h5td` and `.zarr`) in R as
[TreeSummarizedExperiment](https://bioconductor.org/packages/TreeSummarizedExperiment)
objects.

```r
library(treedataR)

tse <- read_tse("data.h5td")          # or "data.zarr"

obst(tse)                              # annotated trees (tidytree::treedata)
TreeSummarizedExperiment::colTree(tse) # plain phylo trees
reducedDim(tse, "characters")          # cell x site character matrix

# ancestral states on internal nodes
character_matrix(obst(tse)[[1]], nodes = "internal")

write_tse(tse, "out.h5td")            # readable by treedata.read_h5td()
```

## How the data are mapped

| treedata (Python)                 | TreeSummarizedExperiment                                  |
| --------------------------------- | --------------------------------------------------------- |
| `X`, `layers`                     | assays (`"X"` + one per layer)                            |
| `obs`, `var`                      | `colData`, `rowData`                                      |
| `obsm`, `varm`, `obsp`, `varp`, `uns` | `reducedDims`, ..., `metadata` (via anndataR)          |
| `obsm["characters"]` (dataframe)  | `reducedDim(tse, "characters")` (character matrix)        |
| `obst`, `vart`                    | `colTree`/`colLinks`, `rowTree`/`rowLinks` (plain `phylo`) |
| tree node / edge attributes       | `obst(tse)`, `vart(tse)`: `tidytree::treedata` node data  |
| `label`, `alignment`, `allow_overlap` | `metadata(tse)$treedata`                              |

- **Forests.** Each cell is linked to the tree named in its `obs[[label]]`
  column (default `"tree"`). Cells not in any tree are kept, with `NA` links.
- **Internal nodes.** Node attributes such as `time`, `depth` or ancestral
  `characters` are columns of the `treedata` tibble. List-valued attributes
  become list-columns. With `alignment = "nodes"`, cells mapped to internal
  nodes are linked to those nodes.
- **Edge attributes.** Each edge attribute is stored on the edge's child node.
  If its name clashes with a node attribute, it gets an `edge_` prefix.
- **Branch lengths** come from the `length` edge attribute. Where it is absent,
  they are `time[child] - time[parent]`. Lengths derived this way are not
  written back as `length` attributes.
- **Limitation.** An attribute that is explicitly `None` in Python and one that
  is absent both read as `NA`.

## Installation

```r
BiocManager::install(c("TreeSummarizedExperiment", "tidytree", "rhdf5"))
remotes::install_github("scverse/anndataR")   # >= 1.3.2
remotes::install_github("Huber-group-EMBL/Rarr") # only needed for .zarr
remotes::install_local("path/to/treedataR")
```

The development version of anndataR (>= 1.3.2) is required. It is the first
version that reads the `nullable-string-array` encoding written by
anndata >= 0.13.

## Development

Tests use small fixtures in `inst/extdata`, regenerated with
`python inst/scripts/make_fixtures.py` (needs Python `treedata` >= 0.3).

To also check that Python reads the files written by R, set `TREEDATAR_PYTHON`
to a Python interpreter that has `treedata` installed:

```sh
TREEDATAR_PYTHON=/path/to/python Rscript -e 'devtools::test()'
```
