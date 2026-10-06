# Checks that files written by treedataR load in Python treedata with the same
# trees and data as the original fixtures. Set TREEDATAR_PYTHON to a Python
# with `treedata` installed to enable.

python <- Sys.getenv("TREEDATAR_PYTHON")

check_py <- '
import sys, pandas as pd, treedata as td
orig, new, fmt = sys.argv[1], sys.argv[2], sys.argv[3]
a = td.read_h5td(orig)
b = td.read_h5td(new) if fmt == "h5td" else td.read_zarr(new)
assert a.alignment == b.alignment, (a.alignment, b.alignment)
assert a.allow_overlap == b.allow_overlap
assert a.label == b.label
for k in ["obst", "vart"]:
    A, B = getattr(a, k), getattr(b, k)
    assert set(A.keys()) == set(B.keys()), (k, list(A.keys()), list(B.keys()))
    for n in A.keys():
        ga, gb = A[n], B[n]
        assert set(ga.nodes) == set(gb.nodes) and set(ga.edges) == set(gb.edges), n
        for v in ga.nodes:
            assert ga.nodes[v] == gb.nodes[v], (n, v, ga.nodes[v], gb.nodes[v])
        for e in ga.edges:
            assert ga.edges[e] == gb.edges[e], (n, e, ga.edges[e], gb.edges[e])
if "characters" in a.obsm:
    # zarr string columns written by anndataR come back as object dtype
    pd.testing.assert_frame_equal(a.obsm["characters"], b.obsm["characters"], check_dtype=False)
assert list(a.obs_names) == list(b.obs_names)
assert abs(a.X - b.X).sum() == 0
print("OK")
'

for (fmt in formats) {
    for (name in c("forest", "nodes")) {
        test_that(paste("Python treedata reads R-written files:", name, fmt), {
            skip_if(python == "", "TREEDATAR_PYTHON not set")
            tse <- read_tse(fixture(paste0(name, ".", fmt)))
            out <- file.path(withr::local_tempdir(), paste0("out.", fmt))
            write_tse(tse, out)
            script <- withr::local_tempfile(fileext = ".py")
            writeLines(check_py, script)
            res <- suppressWarnings(system2(python, c(script, fixture(paste0(name, ".h5td")), out, fmt),
                stdout = TRUE, stderr = TRUE
            ))
            expect_equal(tail(res, 1), "OK", info = paste(res, collapse = "\n"))
        })
    }
}
