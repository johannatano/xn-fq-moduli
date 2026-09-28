from dataclasses import dataclass
from typing import Any


@dataclass
class Data:
    label: str
    value: Any


@dataclass
class ResultData(Data):
    show: bool = True
    fmt: str = ""
    width: int = 0  # 0 = auto; set explicitly when value contains ANSI codes
    flag: str = ""  # "non-zero" | "zero" | "negative" — triggers error color on the row

    def formatted(self) -> str:
        # Support special formatter names
        if self.fmt == "factors":
            from utils.fmt import fmt_factored
            return fmt_factored(int(self.value))

        # If fmt looks like an ANSI color escape, treat it as a color wrapper
        if self.fmt and "\033" in self.fmt:
            try:
                s = str(self.value)
            except Exception:
                s = repr(self.value)
            return f"{self.fmt}{s}\033[0m"

        # Otherwise treat fmt as a normal format specifier
        if self.fmt:
            try:
                return format(self.value, self.fmt)
            except Exception:
                return str(self.value)
        return str(self.value)

    def is_flagged(self) -> bool:
        if self.flag == "non-zero":
            return self.value != 0
        if self.flag == "zero":
            return self.value == 0
        if self.flag == "negative":
            return self.value < 0
        return False

    def __str__(self) -> str:
        return f"{self.label}: {self.formatted()}"

    def __repr__(self) -> str:
        return str(self)
