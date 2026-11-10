"""The Kang 2018 interferon-beta PBMC dataset.

Peripheral blood mononuclear cells from 8 lupus patients, each split into a
control aliquot and an aliquot stimulated with interferon-beta for 6 hours
(Kang et al. 2018, Nat Biotechnol 36:89-94, GEO GSE96583).

This dataset is the right one for the question this project asks because it
has genuine biological replication: 8 donors, each contributing to both
conditions. That paired design is what makes it possible to measure how
badly a per-cell test behaves, because donors can be shuffled to build a
null where the true answer is known to be "no differences".

The interferon response is also one of the strongest and best characterised
transcriptional programmes in immunology, so the genes that should come out
(ISG15, IFI6, IFIT1/3, MX1, OAS1) are known in advance. A method that misses
those is broken; a method that returns them plus three thousand others is
not obviously working either.
"""

from pathlib import Path

import requests

DATA_DIR = Path(__file__).resolve().parents[2] / "data"

# The scverse example-data host, which is what pertpy.data.kang_2018()
# fetches. 38.4 MB, md5 adb2246232e8493031c576982c0c02a3.
#
# Not figshare.com/ndownloader/files/34464122, which looks like a download
# link and is not: it is an async job-dispatch endpoint that answers 202
# with Content-Length 0, so a naive download writes an empty file and the
# failure only shows up later as "file signature not found" from h5py. The
# working figshare mirror is ndownloader.figshare.com/files/34464122, which
# serves a byte-identical copy but signs its redirect for GET only and
# expires it after ten seconds.
KANG_URL = "https://exampledata.scverse.org/pertpy/kang_2018.h5ad"
KANG_FILE = "kang_2018.h5ad"
KANG_MD5 = "adb2246232e8493031c576982c0c02a3"

HDF5_MAGIC = b"\x89HDF\r\n\x1a\n"

# canonical interferon-stimulated genes, used as a positive control
ISG_MARKERS = [
    "ISG15", "IFI6", "IFIT1", "IFIT3", "MX1", "OAS1", "IFI44L",
    "RSAD2", "IFITM3", "STAT1", "IRF7", "CXCL10",
]


def download_kang(force=False):
    """Download and verify, rather than download and hope.

    The check is here because the first URL tried for this file silently
    produced a zero-byte download; validating the HDF5 magic bytes turns
    that into an immediate, legible error instead of a confusing h5py
    failure several steps later.
    """
    import hashlib

    dest = DATA_DIR / KANG_FILE
    if dest.exists() and not force and _looks_like_hdf5(dest):
        return dest

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    print(f"downloading {KANG_FILE} ...", flush=True)
    r = requests.get(KANG_URL, timeout=1800, stream=True, allow_redirects=True)
    r.raise_for_status()

    digest = hashlib.md5()
    with open(dest, "wb") as fh:
        for chunk in r.iter_content(chunk_size=1 << 20):
            fh.write(chunk)
            digest.update(chunk)

    size = dest.stat().st_size
    if size == 0:
        raise RuntimeError(f"{KANG_URL} returned an empty file "
                           f"(HTTP {r.status_code})")
    if not _looks_like_hdf5(dest):
        raise RuntimeError(f"{dest} is not an HDF5 file; the download "
                           "probably returned an error page")
    got = digest.hexdigest()
    if got != KANG_MD5:
        print(f"  warning: md5 {got} != expected {KANG_MD5}")

    print(f"  -> {dest} ({size / 1e6:.1f} MB, md5 {got})", flush=True)
    return dest


def _looks_like_hdf5(path):
    try:
        with open(path, "rb") as fh:
            return fh.read(8) == HDF5_MAGIC
    except OSError:
        return False


def load_kang():
    """Returns an AnnData with raw counts in .X.

    In this file the donor column is called 'replicate' and the condition
    column is called 'label', neither of which is guessable. They are
    renamed to 'donor' and 'condition' so the rest of the package does not
    have to know, and the candidate lists cover the other namings this
    dataset circulates under.
    """
    import anndata as ad

    adata = ad.read_h5ad(download_kang())

    rename = {}
    for cand in ("label", "stim", "condition"):
        if cand in adata.obs.columns:
            rename[cand] = "condition"
            break
    for cand in ("replicate", "ind", "patient", "donor", "individual"):
        if cand in adata.obs.columns:
            rename[cand] = "donor"
            break
    for cand in ("cell_type", "cell", "celltype", "cell_abbr"):
        if cand in adata.obs.columns:
            rename[cand] = "cell_type"
            break
    adata.obs = adata.obs.rename(columns=rename)

    for required in ("condition", "donor"):
        if required not in adata.obs.columns:
            raise KeyError(
                f"expected a {required!r} column; found {list(adata.obs.columns)}")

    adata.obs["donor"] = adata.obs["donor"].astype(str)
    adata.obs["condition"] = adata.obs["condition"].astype(str)
    return adata


def describe(adata):
    """What is actually in the object, printed rather than assumed."""
    print(f"{adata.n_obs:,} cells x {adata.n_vars:,} genes")
    for col in ("condition", "donor", "cell_type"):
        if col in adata.obs.columns:
            counts = adata.obs[col].value_counts()
            print(f"  {col}: {len(counts)} levels -> "
                  f"{dict(list(counts.items())[:6])}"
                  f"{' ...' if len(counts) > 6 else ''}")
    import numpy as np
    X = adata.X
    sample = X[:100].toarray() if hasattr(X, "toarray") else X[:100]
    is_int = np.allclose(sample, np.round(sample))
    print(f"  .X looks like {'raw counts' if is_int else 'transformed values'}"
          f" (max {sample.max():.2f})")
