from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path
from typing import Sequence


def get_classnb_db_main(argv: Sequence[str] | None = None) -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--url", help="URL to the precomputed classnb.db."
                   " If omitted, reads from XNFQ_CLASSNB_DB_URL env var.")
    p.add_argument("--dest", type=Path, help="Destination path for the DB file")
    p.add_argument("--force", action="store_true", help="Overwrite destination if exists")
    args = p.parse_args(argv)

    url = args.url or os.getenv("XNFQ_CLASSNB_DB_URL")
    if not url:
        print("Error: no URL provided. Pass --url or set XNFQ_CLASSNB_DB_URL.")
        return 2

    from xnfq.store.db import download_db, find_db

    dest = args.dest if args.dest else None

    existing = find_db()
    if existing and not args.force:
        print(f"Existing DB found at: {existing}. Use --force to overwrite or --dest to pick another location.")
        return 0

    download_db(url, dest)
    return 0


def main() -> None:
    raise SystemExit(get_classnb_db_main(sys.argv[1:]))


if __name__ == "__main__":
    main()
