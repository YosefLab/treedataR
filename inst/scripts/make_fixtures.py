"""Build small TreeData test fixtures for treedataR.

Usage: python make_fixtures.py [outdir]   (default: inst/extdata)

Requires the Python `treedata` package (>= 0.3, columnar zarr trees).
"""

import sys
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
import scipy.sparse as sp
import treedata as td


def _tree(edges, times, characters, lengths=None, n_edits=None):
    t = nx.DiGraph()
    t.add_edges_from(edges)
    for n in t.nodes:
        t.nodes[n]["time"] = float(times[n])
        t.nodes[n]["characters"] = characters[n]
    if lengths is not None:
        for (u, v), length in lengths.items():
            t[u][v]["length"] = length
    if n_edits is not None:
        for (u, v), k in n_edits.items():
            t[u][v]["n_edits"] = k
    return t


def forest() -> td.TreeData:
    """Three-tree forest aligned to leaves.

    t1: every edge has `length` (+ int `n_edits`), internal node `depth`, sparse bool `stump`
    t2: no edge lengths; lengths must come from node `time`
    t3: mixed; only some edges have `length`
    """
    # t1: r1 -> (a1 -> (c1, c2), c3)
    t1 = _tree(
        [("r1", "a1"), ("a1", "c1"), ("a1", "c2"), ("r1", "c3")],
        times={"r1": 0, "a1": 1, "c1": 3, "c2": 3, "c3": 3},
        characters={
            "r1": ["0", "0", "0"],
            "a1": ["1", "0", "0"],
            "c1": ["1", "2", "0"],
            "c2": ["1", "3", "-"],
            "c3": ["0", "0", "4"],
        },
        lengths={("r1", "a1"): 1.0, ("a1", "c1"): 2.0, ("a1", "c2"): 2.0, ("r1", "c3"): 3.0},
        n_edits={("r1", "a1"): 1, ("a1", "c1"): 1, ("a1", "c2"): 2, ("r1", "c3"): 1},
    )
    t1.nodes["r1"]["depth"] = 0
    t1.nodes["a1"]["depth"] = 1
    t1.nodes["a1"]["stump"] = True
    # t2: r2 -> (c4, c5)
    t2 = _tree(
        [("r2", "c4"), ("r2", "c5")],
        times={"r2": 0.5, "c4": 2.0, "c5": 3.0},
        characters={"r2": ["0", "0", "0"], "c4": ["5", "0", "0"], "c5": ["0", "6", "*"]},
    )
    # t3: r3 -> (a3 -> (c6, c7)); only r3->a3 has a length
    t3 = _tree(
        [("r3", "a3"), ("a3", "c6"), ("a3", "c7")],
        times={"r3": 0, "a3": 1.5, "c6": 2.5, "c7": 4.0},
        characters={
            "r3": ["0", "0", "0"],
            "a3": ["0", "0", "7"],
            "c6": ["0", "1", "7"],
            "c7": ["2", "0", "7"],
        },
        lengths={("r3", "a3"): 1.25},
    )
    cells = ["c1", "c2", "c3", "c4", "c5", "c6", "c7"]
    rng = np.random.default_rng(0)
    X = sp.csr_matrix(rng.poisson(1.0, size=(len(cells), 4)).astype(np.float32))
    obs = pd.DataFrame(
        {"cluster": pd.Categorical(["A", "A", "B", "B", "C", "C", "A"]), "score": np.arange(7, dtype=float)},
        index=cells,
    )
    var = pd.DataFrame(index=[f"g{i}" for i in range(4)])
    chars = pd.DataFrame(
        [t.nodes[c]["characters"] for t in (t1, t2, t3) for c in t.nodes if c in cells],
        index=[c for t in (t1, t2, t3) for c in t.nodes if c in cells],
        columns=["site1", "site2", "site3"],
    ).loc[cells]
    tdata = td.TreeData(X=X, obs=obs, var=var, obst={"t1": t1, "t2": t2, "t3": t3})
    tdata.obsm["characters"] = chars
    tdata.obsm["X_umap"] = rng.normal(size=(len(cells), 2))
    return tdata


def nodes_aligned() -> td.TreeData:
    """Single tree where internal nodes are also observations (alignment='nodes'), plus a vart tree."""
    t = nx.DiGraph([("root", "m1"), ("m1", "l1"), ("m1", "l2"), ("root", "l3")])
    for i, n in enumerate(t.nodes):
        t.nodes[n]["time"] = float(i)
    for u, v in t.edges:
        t[u][v]["length"] = 1.0
    cells = list(t.nodes)
    vt = nx.DiGraph([("gr", "g0"), ("gr", "g1"), ("gr", "g2")])
    var = pd.DataFrame(index=["g0", "g1", "g2", "gr"])
    X = np.arange(len(cells) * len(var), dtype=float).reshape(len(cells), len(var))
    return td.TreeData(
        X=X, obs=pd.DataFrame(index=cells), var=var, obst={"lineage": t}, vart={"genes": vt}, alignment="nodes"
    )


def main() -> None:
    outdir = Path(sys.argv[1] if len(sys.argv) > 1 else "inst/extdata")
    outdir.mkdir(parents=True, exist_ok=True)
    for name, make in [("forest", forest), ("nodes", nodes_aligned)]:
        make().write_h5td(outdir / f"{name}.h5td")
        make().write_zarr(outdir / f"{name}.zarr")


if __name__ == "__main__":
    main()
