from __future__ import annotations

from dataclasses import dataclass


@dataclass
class SessionConfig:
    use_pari: bool = False
    fast_trace: bool = True


_CONFIG = SessionConfig()
_PARI = None


def get_config() -> SessionConfig:
    return _CONFIG


def get_pari():
    global _PARI

    if not _CONFIG.use_pari:
        return None
    if _PARI is None:
        from cypari2 import Pari

        _PARI = Pari()
    return _PARI


def configure_session(
    *,
    use_pari: bool | None = None,
    fast_trace: bool | None = None,
) -> SessionConfig:
    if use_pari is not None:
        _CONFIG.use_pari = use_pari
    if fast_trace is not None:
        _CONFIG.fast_trace = fast_trace
    return _CONFIG


def add_config_args(parser) -> None:
    group = parser.add_argument_group("session configuration")
    group.add_argument(
        "--use-pari",
        action="store_true",
        help="Use cypari2 for fundamental discriminant computation.",
    )
    group.add_argument(
        "--fast-trace",
        action="store_true",
        default=True,
        help="Iterate only the trace residue classes compatible with N-torsion.",
    )


def apply_config_from_args(args) -> SessionConfig:
    return configure_session(
        use_pari=getattr(args, "use_pari", None),
        fast_trace=getattr(args, "fast_trace", None),
    )
