"""Populate store/classnb.db with class numbers for fundamental discriminants
-1 > D >= -max-abs-d, D % 4 in (0, 1), using the pure-Python classnumber().

Usage:
    python scripts/build_db.py --max 1000000
"""

from __future__ import annotations
import argparse
import sqlite3
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from xnfq.arithmetic.quadratic import cln
from xnfq.config import add_config_args, apply_config_from_args

DB_PATH = (
    Path(__file__).resolve().parents[1] / "src" / "xnfq" / "store" / "classnb.db"
)

BATCH_SIZE = 100_000

def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--start-d", type=int, default=1)
    p.add_argument("--max-d", type=int, default=100_000)
    p.add_argument("--batch-size", type=int, default=BATCH_SIZE)
    p.add_argument("--db", type=Path, default=DB_PATH)
    # add session configuration flags (e.g. --use-pari)
    add_config_args(p)
    args = p.parse_args()
    # apply configuration (will enable PARI if --use-pari was passed)
    apply_config_from_args(args)
    # Ensure the parent directory for the DB exists so sqlite can create the file
    try:
        args.db.parent.mkdir(parents=True, exist_ok=True)
    except Exception:
        # If args.db has no parent or cannot be created, let sqlite raise a clearer error
        pass

    conn = sqlite3.connect(args.db)
    conn.execute(
        "CREATE TABLE IF NOT EXISTS classnb (D INTEGER PRIMARY KEY, h INTEGER)"
    )

    start = time.time()
    batch: list[tuple[int, int]] = []
    cur = conn.cursor()
    added = 0
    skipped = 0
    for D in range(-args.start_d, -args.max_d - 1, -1):
        if D % 4 not in (0, 1):
            continue
        # skip if already present in DB
        cur.execute("SELECT 1 FROM classnb WHERE D = ?", (D,))
        if cur.fetchone():
            skipped += 1
            if skipped % args.batch_size == 0:
                print(f"\rSkipped {skipped} existing entries...", end="", flush=True)
            continue
        h = cln(D)
        batch.append((D, h))
        if len(batch) >= args.batch_size:
            conn.executemany("INSERT OR REPLACE INTO classnb VALUES (?,?)", batch)
            conn.commit()
            added += len(batch)
            print(f"\rD={D}  ({time.time() - start:.1f}s)", end="", flush=True)
            batch.clear()

    if batch:
        conn.executemany("INSERT OR REPLACE INTO classnb VALUES (?,?)", batch)
        conn.commit()
        added += len(batch)

    conn.execute("PRAGMA wal_checkpoint(FULL)")
    conn.close()
    print(f"\nDone in {time.time() - start:.1f}s -> {args.db}")
    print(f"Added {added} new entries, skipped {skipped} existing entries")


if __name__ == "__main__":
    main()
