#!/usr/bin/env python3
"""Download a precomputed `classnb.db` into the user data dir or provided path.

This is intentionally separate from installation so callers can opt-in to
downloading a large binary file from a remote host.

Usage:
    python scripts/get_classnb_db.py --url https://example.com/classnb.db

If `--dest` is omitted the file is saved to `$XDG_DATA_HOME/xnfq/classnb.db`
or `~/.local/share/xnfq/classnb.db`.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import os
import sys


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--url", help="URL to the precomputed classnb.db."
                   " If omitted, reads from XNFQ_CLASSNB_DB_URL env var.")
    p.add_argument("--dest", type=Path, help="Destination path for the DB file")
    p.add_argument("--force", action="store_true", help="Overwrite destination if exists")
    args = p.parse_args()

    url = args.url or os.getenv("XNFQ_CLASSNB_DB_URL")
    if not url:
        print("Error: no URL provided. Pass --url or set XNFQ_CLASSNB_DB_URL.")
        raise SystemExit(2)

    # Delegate to package CLI (ensures consistent behavior when installed).
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
    from xnfq.cli import get_classnb_db_main

    raise SystemExit(get_classnb_db_main([arg for arg in sys.argv[1:]]))


if __name__ == "__main__":
    main()
