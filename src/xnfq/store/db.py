from __future__ import annotations

import os
import sqlite3
import sys
from pathlib import Path
from typing import Iterable
import urllib.request


def _packaged_path() -> Path:
    """Return the packaged path `.../xnfq/store/classnb.db` inside the project/package."""
    return Path(__file__).resolve().parent / "classnb.db"


def find_db() -> Path | None:
    """Return the packaged `classnb.db` path if it exists, otherwise None.

    This package is configured to prefer the DB file located exactly in the
    package's `xnfq/store` directory (e.g. `src/xnfq/store/classnb.db`).
    """
    p = _packaged_path()
    try:
        return p if p.exists() else None
    except Exception:
        return None


def query_classnb(D: int) -> int | None:
    """Query the precomputed DB for class number of discriminant `D`.

    Returns the class number as `int` if present, otherwise `None`.
    """
    db = find_db()
    if db is None:
        return None
    try:
        conn = sqlite3.connect(str(db))
        cur = conn.execute("SELECT h FROM classnb WHERE D = ?", (D,))
        row = cur.fetchone()
        conn.close()
        if row is None:
            return None
        return int(row[0])
    except Exception:
        return None


def _ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def download_db(url: str, dest: Path | None = None, *, chunk_size: int = 8192) -> Path:
    """Download a DB from `url` and save it to `dest`.

    By default (when `dest` is None) the file is saved to the packaged store
    path: `.../xnfq/store/classnb.db` so the repository's `xnfq/store` location
    is used as the canonical place for the DB.

    Returns the final `Path` to the downloaded file.
    """
    if dest is None:
        dest = _packaged_path()
    dest = Path(dest)
    _ensure_parent(dest)

    tmp = dest.with_suffix(dest.suffix + ".tmp")
    with urllib.request.urlopen(url) as resp, open(tmp, "wb") as out:
        total = resp.getheader("Content-Length")
        if total is not None:
            total = int(total)
        downloaded = 0
        while True:
            chunk = resp.read(chunk_size)
            if not chunk:
                break
            out.write(chunk)
            downloaded += len(chunk)
            if total:
                pct = downloaded * 100 // total
                print(f"\rDownloading: {pct}% ({downloaded}/{total} bytes)", end="", flush=True)
    tmp.replace(dest)
    print(f"\nSaved DB to: {dest}")
    return dest
