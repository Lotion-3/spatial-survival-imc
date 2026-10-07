"""List and fetch only the tables we need from the Jackson et al. 2020 Zenodo record.

The record (DOI 10.5281/zenodo.4607374) contains a ~37 GB archive
(SingleCell_and_Metadata.zip) that holds the clinical metadata. We never download it:
instead an HTTP-range-backed file object lets ``zipfile`` read the central directory
and extract individual members.

Usage:
    python scripts/download_data.py --list               # list record files + archive members
    python scripts/download_data.py                      # fetch needed tables into data/raw
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import shutil
import sys
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


def open_remote_zip(info: dict) -> tuple[zipfile.ZipFile, HTTPRangeFile]:
    raw = HTTPRangeFile(info["url"], info["size"])
    buf = io.BufferedReader(raw, buffer_size=1 << 20)
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


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", action="store_true", help="only list files and archive members")
    args = ap.parse_args()
    files = record_files()
    list_all(files) if args.list else fetch(files)


if __name__ == "__main__":
    main()
