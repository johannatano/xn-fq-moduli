from __future__ import annotations

import argparse

from xnfq.config import add_config_args


def _weight_type(v: str) -> int:
    try:
        k = int(v)
    except Exception:
        raise argparse.ArgumentTypeError("weight k must be an integer >= 2")
    if k < 2:
        raise argparse.ArgumentTypeError("weight k must be >= 2")
    return k


def parse_example_args(
    description: str,
    *,
    include_prime_range: bool = False,
    include_weight: bool = False,
    include_sage: bool = False,
    overrides: dict | None = None,
):
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-p", "--p", type=int, default=5, help="Field characteristic")
    parser.add_argument("-n", "--n", type=int, default=1, help="Extension degree")
    parser.add_argument("-N", "--N", type=int, default=11, help="Level N")
    parser.add_argument("-type", "--type", type=int, default=1, help="Level type")

    if include_prime_range:
        parser.add_argument(
            "--range",
            type=int,
            default=None,
            help="Interval width: include primes in [p, p+range) starting at baseline p",
        )

    if include_weight:
        parser.add_argument("-k", type=_weight_type, default=2, help="Weight k (integer >= 2)")

    if include_sage:
        parser.add_argument(
            "--sage", action="store_true", help="Compare against Sage's Hecke operator"
        )

    add_config_args(parser)
    args = parser.parse_args()

    # If overrides supplied, apply them only when the corresponding
    # option was not provided on the command line.
    if overrides:
        import sys

        provided = set(sys.argv[1:])
        # mapping from attribute name (dest) to option tokens to check
        opt_map = {
            "p": ("-p", "--p"),
            "n": ("-n", "--n"),
            "N": ("-N", "--N"),
            "type": ("-type", "--type"),
            "range": ("--range",),
            "k": ("-k",),
            "sage": ("--sage",),
        }

        for key, val in overrides.items():
            toks = opt_map.get(key, (f"--{key}",))
            if not any(tok in provided for tok in toks):
                setattr(args, key, val)
    
        # Determine which option flags were present on the command line so
        # callers can distinguish "not provided" from the parser's defaults.
        import sys
    
        provided_tokens = set(sys.argv[1:])
        # mapping from attribute name (dest) to option tokens
        opt_map = {
            "p": ("-p", "--p"),
            "n": ("-n", "--n"),
            "N": ("-N", "--N"),
            "type": ("-type", "--type"),
            "range": ("--range",),
            "k": ("-k",),
            "sage": ("--sage",),
        }
        provided_dest = {key for key, toks in opt_map.items() if any(tok in provided_tokens for tok in toks)}
    
        # attach the set of provided dest names so callers may check presence
        setattr(args, "_provided", provided_dest)
    
        # If overrides supplied, apply them only when the corresponding
        # option was not provided on the command line.
        if overrides:
            for key, val in overrides.items():
                if key not in provided_dest:
                    setattr(args, key, val)

    if include_prime_range:
        if args.range is not None:
            from utils.primes import next_primes
            args.primes = next_primes(int(args.p), max(1, int(args.range)))
        else:
            args.primes = [int(args.p)]

    return args
