"""List and fetch only the tables we need from the Jackson et al. 2020 Zenodo record.

The record (DOI 10.5281/zenodo.4607374) contains a ~37 GB archive
(SingleCell_and_Metadata.zip) that holds the clinical metadata. We never download it:
instead an HTTP-range-backed file object lets ``zipfile`` read the central directory
and extract individual members.

Usage:
    python scripts/download_data.py --list               # list record files + archive members
    python scripts/download_data.py                      # fetch needed tables into data/raw
    python scripts/download_data.py --metabric           # also METABRIC IMC + clinical (part II)
    python scripts/download_data.py --genomics           # also METABRIC mutations + gene panel (part IV)
    python scripts/download_data.py --basel-expression   # also Basel per-cell marker intensities (part IV)
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import shutil
import sys
import time
import urllib.request
import zipfile
from pathlib import Path

RECORD_ID = "4607374"
API = f"https://zenodo.org/api/records/{RECORD_ID}"
BIG_ARCHIVE = "SingleCell_and_Metadata.zip"
SMALL_FILES = ["singlecell_cluster_labels.zip", "singlecell_locations.zip"]
MAX_MEMBER_BYTES = 2_000_000_000  # refuse to pull any single member over ~2 GB

# Members of the big archive we need. Only Basel has a standalone patient metadata table;
# Zurich metadata exists only inside a 3.7 GB histoCAT .mat session, so we do not use it.
NEEDED_MEMBERS = [
    "Data_publication/BaselTMA/Basel_PatientMetadata.csv",
]

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"


class HTTPRangeFile(io.RawIOBase):
    """Read-only seekable file over HTTP using Range requests."""

    def __init__(self, url: str, size: int):
        self.url, self.size, self.pos = url, size, 0
        self.bytes_fetched = 0

    def seekable(self):
        return True

    def readable(self):
        return True

    def tell(self):
        return self.pos

    def seek(self, offset, whence=0):
        if whence == 0:
            self.pos = offset
        elif whence == 1:
            self.pos += offset
        else:
            self.pos = self.size + offset
        return self.pos

    def read(self, n=-1):
        if n is None or n < 0:
            n = self.size - self.pos
        n = min(n, self.size - self.pos)
        if n <= 0:
            return b""
        req = urllib.request.Request(
            self.url, headers={"Range": f"bytes={self.pos}-{self.pos + n - 1}"}
        )
        with urllib.request.urlopen(req, timeout=120) as r:
            if r.status != 206:
                raise RuntimeError(f"server ignored Range request (status {r.status})")
            data = r.read()
        self.pos += len(data)
        self.bytes_fetched += len(data)
        return data

    def readinto(self, b):
        data = self.read(len(b))
        b[: len(data)] = data
        return len(data)


def record_files() -> dict[str, dict]:
    with urllib.request.urlopen(API, timeout=60) as r:
        meta = json.load(r)
    return {
        f["key"]: {"size": f["size"], "checksum": f["checksum"], "url": f["links"]["self"]}
        for f in meta["files"]
    }


def open_remote_zip(info: dict, buffer_size: int = 1 << 20) -> tuple[zipfile.ZipFile, HTTPRangeFile]:
    raw = HTTPRangeFile(info["url"], info["size"])
    buf = io.BufferedReader(raw, buffer_size=buffer_size)
    return zipfile.ZipFile(buf), raw


def md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url: str, dest: Path) -> None:
    with urllib.request.urlopen(url, timeout=300) as r, open(dest, "wb") as f:
        shutil.copyfileobj(r, f, 1 << 20)


def list_all(files: dict[str, dict]) -> None:
    print(f"Zenodo record {RECORD_ID}")
    for k, v in files.items():
        print(f"  {k:40s} {v['size'] / 1e6:12.1f} MB  {v['checksum']}")
    zf, raw = open_remote_zip(files[BIG_ARCHIVE])
    print(f"\nMembers of {BIG_ARCHIVE} (central directory read with {raw.bytes_fetched / 1e6:.1f} MB of range requests):")
    for zi in zf.infolist():
        print(f"  {zi.file_size / 1e6:12.2f} MB  {zi.filename}")


def fetch(files: dict[str, dict]) -> None:
    RAW.mkdir(parents=True, exist_ok=True)
    manifest = {}
    for name in SMALL_FILES:
        dest = RAW / name
        expected = files[name]["checksum"].removeprefix("md5:")
        if not dest.exists() or md5(dest) != expected:
            print(f"downloading {name} ({files[name]['size'] / 1e6:.1f} MB)")
            download(files[name]["url"], dest)
        got = md5(dest)
        if got != expected:
            sys.exit(f"checksum mismatch for {name}: {got} != {expected}")
        manifest[name] = {"md5": got, "source": files[name]["url"]}
        out = RAW / name.removesuffix(".zip")
        with zipfile.ZipFile(dest) as z:
            z.extractall(out)

    zf, raw = open_remote_zip(files[BIG_ARCHIVE])
    by_name = {zi.filename: zi for zi in zf.infolist()}
    for member in NEEDED_MEMBERS:
        zi = by_name.get(member)
        if zi is None:
            sys.exit(f"member not found in archive: {member}")
        if zi.file_size > MAX_MEMBER_BYTES:
            sys.exit(f"{member} is {zi.file_size / 1e9:.1f} GB; refusing (ask first)")
        dest = RAW / Path(member).name
        if not dest.exists():
            print(f"range-extracting {member} ({zi.file_size / 1e6:.2f} MB)")
            with zf.open(zi) as src, open(dest, "wb") as f:  # zipfile verifies CRC32
                shutil.copyfileobj(src, f)
        manifest[member] = {"md5": md5(dest), "crc32": f"{zi.CRC:08x}", "source": f"{BIG_ARCHIVE}::{member}"}
    print(f"range requests fetched {raw.bytes_fetched / 1e6:.1f} MB from {BIG_ARCHIVE}")
    (RAW / "MANIFEST.json").write_text(json.dumps(manifest, indent=2))
    print(f"wrote {RAW / 'MANIFEST.json'}")


# --------------------------------------------------------------------------- METABRIC (part II)

METABRIC_RECORD = "5850952"  # Danenberg et al. 2022, CC-BY-4.0
METABRIC_MEMBER = "toPublicRepository/SingleCells.csv"  # ~383 MB compressed, ~850 MB extracted
RAW_MB = ROOT / "data" / "raw_metabric"
CBIO = "https://www.cbioportal.org/api/studies/brca_metabric/clinical-data"


def fetch_metabric() -> None:
    RAW_MB.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(f"https://zenodo.org/api/records/{METABRIC_RECORD}", timeout=60) as r:
        f = json.load(r)["files"][0]
    zf, raw = open_remote_zip({"url": f["links"]["self"], "size": f["size"]}, buffer_size=32 << 20)
    zi = zf.getinfo(METABRIC_MEMBER)
    dest = RAW_MB / "SingleCells.csv"
    if not dest.exists() or dest.stat().st_size != zi.file_size:
        print(f"range-extracting {METABRIC_MEMBER} ({zi.compress_size / 1e6:.0f} MB compressed)")
        with zf.open(zi) as src, open(dest, "wb") as out:  # zipfile verifies CRC32 at EOF
            shutil.copyfileobj(src, out, 1 << 22)
    # Public cBioPortal clinical data: patient- and sample-level attributes, long format.
    rows = []
    for kind in ["PATIENT", "SAMPLE"]:
        req = urllib.request.Request(  # cBioPortal rejects urllib's default User-Agent (403)
            f"{CBIO}?clinicalDataType={kind}&projection=SUMMARY",
            headers={"User-Agent": "spatialsurv/0.1 (research script)", "Accept": "application/json"},
        )
        with urllib.request.urlopen(req, timeout=120) as r:
            rows += json.load(r)
    clin = {}
    for d in rows:
        clin.setdefault(d["patientId"], {})[d["clinicalAttributeId"]] = d["value"]
    (RAW_MB / "metabric_clinical.json").write_text(json.dumps(clin, indent=1))
    manifest = {
        "SingleCells.csv": {"md5": md5(dest), "crc32": f"{zi.CRC:08x}", "source": f"zenodo {METABRIC_RECORD}::{METABRIC_MEMBER}"},
        "metabric_clinical.json": {"md5": md5(RAW_MB / "metabric_clinical.json"), "source": CBIO,
                                   "n_patients": len(clin)},
    }
    (RAW_MB / "MANIFEST.json").write_text(json.dumps(manifest, indent=2))
    print(f"METABRIC: {len(clin)} clinical patients; wrote {RAW_MB / 'MANIFEST.json'}")


# --------------------------------------------------------------------------- METABRIC genomics (part IV)

CBIO_API = "https://www.cbioportal.org/api"
MUT_PROFILE = "brca_metabric_mutations"
SEQUENCED = "brca_metabric_sequenced"


def cbio(path: str, body: dict | None = None):
    req = urllib.request.Request(  # cBioPortal rejects urllib's default User-Agent (403)
        f"{CBIO_API}/{path}", data=None if body is None else json.dumps(body).encode(),
        headers={"User-Agent": "spatialsurv/0.1 (research script)", "Accept": "application/json",
                 "Content-Type": "application/json"},
        method="GET" if body is None else "POST",
    )
    with urllib.request.urlopen(req, timeout=300) as r:
        return json.load(r)


def fetch_genomics() -> None:
    """Somatic mutations (targeted panel), the sequenced-sample list and the gene panel(s).

    Subtype (CLAUDIN_SUBTYPE), IntClust and ER/PR/HER2 are clinical attributes, already in
    metabric_clinical.json (--metabric). The API is live, so the snapshot date is recorded.
    """
    RAW_MB.mkdir(parents=True, exist_ok=True)
    muts = cbio(f"molecular-profiles/{MUT_PROFILE}/mutations/fetch?projection=DETAILED", {"sampleListId": SEQUENCED})
    keep = ["sampleId", "patientId", "entrezGeneId", "mutationType", "variantType", "proteinChange", "mutationStatus"]
    rows = [{k: m.get(k) for k in keep} | {"hugoGeneSymbol": m["gene"]["hugoGeneSymbol"]} for m in muts]
    sequenced = cbio(f"sample-lists/{SEQUENCED}")["sampleIds"]
    panel_data = cbio(f"molecular-profiles/{MUT_PROFILE}/gene-panel-data/fetch", {"sampleListId": SEQUENCED})
    panel_of = {d["sampleId"]: d.get("genePanelId") for d in panel_data}
    panels = {}
    for pid in sorted({p for p in panel_of.values() if p}):
        panels[pid] = sorted(g["hugoGeneSymbol"] for g in cbio(f"gene-panels/{pid}")["genes"])
    out = {
        "mutations": RAW_MB / "metabric_mutations.json",
        "sequenced": RAW_MB / "metabric_sequenced_samples.json",
        "panels": RAW_MB / "metabric_gene_panels.json",
    }
    out["mutations"].write_text(json.dumps(rows, indent=0))
    out["sequenced"].write_text(json.dumps({"sampleIds": sorted(sequenced), "genePanelId": panel_of}, indent=0))
    out["panels"].write_text(json.dumps(panels, indent=1))
    snapshot = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    manifest_path = ROOT / "results" / "genotype" / "data_manifest.json"
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    manifest["cbioportal_genomics"] = {
        "snapshot": snapshot, "api": CBIO_API, "study": "brca_metabric", "mutation_profile": MUT_PROFILE,
        "sample_list": SEQUENCED, "n_mutations": len(rows), "n_sequenced_samples": len(sequenced),
        "gene_panels": {k: len(v) for k, v in panels.items()},
        "files": {p.name: {"md5": md5(p), "bytes": p.stat().st_size} for p in out.values()},
    }
    clin = RAW_MB / "metabric_clinical.json"
    if clin.exists():  # subtype / IntClust / receptor status come from the --metabric clinical snapshot
        manifest["cbioportal_clinical"] = {
            "file": clin.name, "md5": md5(clin), "source": CBIO,
            "snapshot_file_mtime": time.strftime("%Y-%m-%dT%H:%M:%S%z", time.localtime(clin.stat().st_mtime)),
            "note": "snapshot taken by --metabric (part II); not re-fetched, so part II numbers are unaffected"}
    manifest_path.write_text(json.dumps(manifest, indent=2))
    print(f"genomics: {len(rows)} mutations, {len(sequenced)} sequenced samples, panels {list(panels)}; "
          f"wrote {manifest_path}")


# --------------------------------------------------------------------------- Basel expression (part IV)

SC_DAT = "Data_publication/BaselTMA/SC_dat.csv"  # long format: core, CellId, id, channel, mc_counts (~4.9 GB)
# Metal tags of the channels kept (Basel_Zuri_StainingPanel.csv). Lineage markers shared with the
# METABRIC panel, plus a few Basel-only ones used only for descriptive heatmaps.
BASEL_METALS = {
    "Lu175": "panCK", "Nd144": "CK8-18", "Pr141": "CK5", "Sm152": "CD3", "Dy164": "CD20", "Nd146": "CD68",
    "Yb172": "CD31-vWF", "Nd148": "SMA", "Dy162": "CD45", "Sm149": "Vimentin", "Nd142": "Fibronectin",
    "Eu151": "HER2", "Er168": "Ki-67", "In113": "Histone H3",
}


def fetch_basel_expression(files: dict[str, dict]) -> None:
    """Stream SC_dat.csv out of the 37 GB archive by HTTP range requests (~0.6 GB compressed),
    keep only BASEL_METALS channels and pivot to one row per cell. The full member never touches
    disk; zipfile verifies its CRC32 once the stream is read to the end."""
    import pandas as pd

    dest = RAW / "basel_expression.csv.gz"
    zf, raw = open_remote_zip(files[BIG_ARCHIVE], buffer_size=32 << 20)
    zi = zf.getinfo(SC_DAT)
    print(f"streaming {SC_DAT} ({zi.compress_size / 1e6:.0f} MB compressed, {zi.file_size / 1e9:.1f} GB)")
    pattern = "|".join(f"{m}Di" for m in BASEL_METALS)
    parts, seen = [], set()
    with zf.open(zi) as src:
        for chunk in pd.read_csv(src, chunksize=5_000_000, usecols=["core", "CellId", "channel", "mc_counts"],
                                 dtype={"core": "category", "CellId": "int32", "channel": "category",
                                        "mc_counts": "float32"}):
            seen |= set(chunk["channel"].cat.categories)
            chunk = chunk[chunk["channel"].astype(str).str.contains(pattern)]
            parts.append(chunk.astype({"core": str, "channel": str}))
    long = pd.concat(parts, ignore_index=True)
    metal = long["channel"].str.extract(r"([A-Z][a-z]?\d{2,3})Di", expand=False)
    long["marker"] = metal.map(BASEL_METALS)
    wide = long.pivot_table(index=["core", "CellId"], columns="marker", values="mc_counts", aggfunc="first")
    wide = wide.reset_index()
    wide.columns.name = None
    wide.to_csv(dest, index=False, compression="gzip", float_format="%.6g")
    kept = sorted({c for c in seen if any(f"{m}Di" in c for m in BASEL_METALS)})
    manifest = json.loads((RAW / "MANIFEST.json").read_text()) if (RAW / "MANIFEST.json").exists() else {}
    manifest[SC_DAT] = {"crc32": f"{zi.CRC:08x}", "source": f"{BIG_ARCHIVE}::{SC_DAT}", "channels_kept": kept,
                        "derived_file": dest.name, "derived_md5": md5(dest), "n_cells": int(len(wide)),
                        "range_bytes_fetched": raw.bytes_fetched}
    (RAW / "MANIFEST.json").write_text(json.dumps(manifest, indent=2))
    print(f"Basel expression: {len(wide)} cells x {len(kept)} channels; {raw.bytes_fetched / 1e6:.0f} MB fetched")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", action="store_true", help="only list files and archive members")
    ap.add_argument("--metabric", action="store_true", help="also fetch the METABRIC IMC cells + clinical (part II)")
    ap.add_argument("--genomics", action="store_true", help="also fetch METABRIC mutations + gene panel (part IV)")
    ap.add_argument("--basel-expression", action="store_true",
                    help="also stream Basel per-cell marker intensities out of the big archive (part IV)")
    args = ap.parse_args()
    files = record_files()
    if args.list:
        list_all(files)
        return
    fetch(files)
    if args.metabric:
        fetch_metabric()
    if args.genomics:
        fetch_genomics()
    if args.basel_expression:
        fetch_basel_expression(files)


if __name__ == "__main__":
    main()
