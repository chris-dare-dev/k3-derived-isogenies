"""Helpers shared by the notebooks: locate the repository and make sure the stored results are present.

The stored results are not in git. They ship as one archive (see DATA.md). ensure_data() looks for them
in experiments/; if they are missing it unpacks a local copy of the archive, or downloads it from
DATA_URL, checks its SHA-256 and unpacks it.
"""
import hashlib
import os
import pathlib
import subprocess
import sys
import tarfile
import urllib.request

DATA_URL = "https://github.com/chris-dare-dev/k3-derived-isogenies/releases/download/v1/k3-derived-isogenies-data-v1.tar.xz"
DATA_SHA256 = "7a16e7d759533558269147fef4b931c3cdaae66eae78fd41463cf2c9977b3714"
ARCHIVE = "k3-derived-isogenies-data-v1.tar.xz"


def repo_root():
    here = pathlib.Path.cwd().resolve()
    for p in (here, *here.parents):
        if (p / "experiments" / "sixlines_ledger.py").exists():
            return p
    raise FileNotFoundError("run the notebooks from inside the repository")


ROOT = repo_root()
EXP = ROOT / "experiments"


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def _unpack(path):
    if _sha256(path) != DATA_SHA256:
        raise ValueError(f"{path}: SHA-256 does not match DATA.md")
    with tarfile.open(path, "r:xz") as tar:
        names = tar.getnames()
        if not all(n == "experiments" or n.startswith("experiments/") for n in names):
            raise ValueError("unexpected paths in the archive")
        if sys.version_info >= (3, 12):
            tar.extractall(ROOT, filter="data")
        else:
            tar.extractall(ROOT)


def ensure_data():
    """Return True when the stored results are in experiments/, fetching them if possible."""
    if (EXP / "sixlines_groupsize.jsonl").exists():
        return True
    for cand in (os.environ.get("K3DI_DATA_ARCHIVE"), ROOT / "dist" / ARCHIVE, ROOT / ARCHIVE):
        if cand and pathlib.Path(cand).exists():
            _unpack(cand)
            return True
    if DATA_URL:
        target = ROOT / ARCHIVE
        urllib.request.urlretrieve(DATA_URL, target)
        _unpack(target)
        return True
    print("The stored results are not present and no download URL is set; see DATA.md.")
    return False


def run(*args, timeout=900):
    """Run one of the scripts in experiments/ with this kernel's Python, as the notes do (python -I)."""
    out = subprocess.run([sys.executable, "-I", *map(str, args)], cwd=EXP, capture_output=True,
                         text=True, timeout=timeout)
    if out.returncode != 0:
        raise RuntimeError(out.stderr[-2000:])
    return out.stdout
