from __future__ import annotations

from typing import Callable

from .data import ResultData

_tqdm = None

class _DummyTqdm:
    def __init__(self, iterable=None, total=None, desc=None, **kwargs):
        self._iterable = iterable
        self.total = total
        self.desc = desc or ""
        self.count = 0
        self._last_len = 0
        self._closed = False

    def __iter__(self):
        if self._iterable is None:
            return iter(())
        for item in self._iterable:
            # update per yielded item
            self.update(1)
            yield item
        # ensure we print final line after iteration
        self.close()

    def update(self, n=1):
        self.count += n
        if self.total:
            pct = 100.0 * self.count / float(self.total)
            msg = f"{self.desc} {pct:5.1f}% ({self.count}/{self.total})"
        else:
            msg = f"{self.desc} {self.count} items"
        pad = max(0, self._last_len - len(msg))
        print(msg + (" " * pad), end="\r", flush=True)
        self._last_len = len(msg)

    def close(self):
        if not self._closed:
            if self.total:
                pct = 100.0 * self.count / float(self.total)
                final = f"{self.desc} {pct:5.1f}% ({self.count}/{self.total})"
            else:
                final = f"{self.desc} {self.count} items"
            print(final + (" " * max(0, self._last_len - len(final))))
            self._closed = True

    def __enter__(self):
        if self.total:
            msg = f"{self.desc} 0.0% (0/{self.total})"
        else:
            msg = f"{self.desc} 0 (0/?)"
        print(msg, end="\r", flush=True)
        self._last_len = len(msg)
        return self

    def __exit__(self, exc_type, exc, tb):
        self.close()
        return False

    @staticmethod
    def write(msg):
        print(msg)


class Colors:
    HEADER = "\033[95m"
    BRIGHT_MAGENTA = "\033[95m"
    BRIGHT_CYAN = "\033[96m"
    BRIGHT_BLUE = "\033[94m"
    BRIGHT_GREEN = "\033[92m"
    BRIGHT_YELLOW = "\033[93m"
    BRIGHT_RED = "\033[91m"
    BRIGHT_WHITE = "\033[97m"

    # Neon / lime green (truecolor). Hex #39FF14 -> RGB (57,255,20)
    NEON_LIME = "\033[38;2;57;255;20m"
    # Neon pink and purple (truecolor)
    # Neon pink: Hex #FF6EC7 -> RGB (255,110,199)
    NEON_PINK = "\033[38;2;255;110;199m"
    # Neon purple: Hex #9B30FF -> RGB (155,48,255)
    NEON_PURPLE = "\033[38;2;155;48;255m"

    RED = "\033[31m"
    GREEN = "\033[32m"
    YELLOW = "\033[33m"
    BLUE = "\033[34m"
    MAGENTA = "\033[35m"
    CYAN = "\033[36m"
    TURQUOISE = "\033[36m"
    WHITE = "\033[37m"

    WARNING = "\033[93m"
    FAIL = "\033[91m"
    SUCCESS = "\033[92m"
    INFO = "\033[94m"

    BOLD = "\033[1m"
    DIM = "\033[2m"
    ITALIC = "\033[3m"
    UNDERLINE = "\033[4m"

    ENDC = "\033[0m"


class Logger:
    HEADLINE = Colors.BRIGHT_CYAN + Colors.BOLD
    INFO = Colors.BRIGHT_BLUE
    SUCCESS = Colors.BRIGHT_GREEN
    WARNING = Colors.BRIGHT_YELLOW
    ERROR = Colors.BRIGHT_RED
    NOTICE = Colors.BRIGHT_MAGENTA
    NEON = Colors.NEON_LIME
    NEON_PINK = Colors.NEON_PINK
    NEON_PURPLE = Colors.NEON_PURPLE
    NORMAL = Colors.ENDC

    @staticmethod
    def cprint(message: str, color: str = "", bold: bool = False) -> None:
        style = color
        if bold:
            style += Colors.BOLD
        print(f"{style}{message}{Colors.ENDC}")

    @staticmethod
    def progress(iterable=None, **kwargs):
        """Return a tqdm progress iterator when available, otherwise a
        lightweight fallback that yields the iterable and implements the
        minimal tqdm API (`update`, `close`, context manager, and `write`).

        Usage:
            for item in Logger.progress(items, desc="..."):
                ...

            with Logger.progress(total=100) as p:
                p.update(10)
        """
        # Always use the built-in minimal progress implementation.
        return _DummyTqdm(iterable=iterable, **kwargs)

    @staticmethod
    def progress_write(msg: str) -> None:
        """Write a message that won't break a running progress bar."""
        print(msg)

    @staticmethod
    def header(message: str, color: str = "") -> None:
        Logger.cprint("-" * 50, color)
        Logger.cprint(message, color, True)
        Logger.cprint("-" * 50, color)

    @staticmethod
    def print_results(
        rows: list[list[ResultData]],
        color_fn: Callable[[list[ResultData]], str] | None = None,
    ) -> None:
        if not rows:
            return

        visible_indices = [i for i, rd in enumerate(rows[0]) if rd.show]
        if not visible_indices:
            return

        headers = ["|  " + rows[0][i].label + " " for i in visible_indices]
        widths = [len(h) for h in headers]

        all_formatted: list[list[str]] = []
        for row in rows:
            fmted = []
            for col, i in enumerate(visible_indices):
                rd = row[i]
                s = rd.formatted()
                fmted.append(s)
                display_len = rd.width if rd.width > 0 else len(s)
                widths[col] = max(widths[col], display_len)
            all_formatted.append(fmted)

        header_str = " ".join(h.rjust(w) for h, w in zip(headers, widths))
        print(header_str)
        print("-" * len(header_str))

        for row, fmted in zip(rows, all_formatted):
            # If any cell provides an ANSI color in `fmt`, use that as the
            # row color (first match). Otherwise fall back to flag/color_fn/default.
            row_color = None
            for rd in row:
                if rd.fmt and "\033" in str(rd.fmt):
                    row_color = rd.fmt
                    break
            clr = (
                row_color
                if row_color is not None
                else (
                    Colors.FAIL
                    if any(rd.is_flagged() for rd in row)
                    else (color_fn(row) if color_fn else Colors.BOLD)
                )
            )
            parts = []
            for col, (i, s) in enumerate(zip(visible_indices, fmted)):
                rd = row[i]
                if rd.width > 0:
                    parts.append(" " * max(0, widths[col] - rd.width) + s)
                else:
                    parts.append(s.rjust(widths[col]))
            print(f"{clr}{' '.join(parts)}{Colors.ENDC}")
