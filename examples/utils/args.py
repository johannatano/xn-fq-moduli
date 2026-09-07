from __future__ import annotations

import argparse

from xnfq.config import add_config_args


def parse_example_args(
    description: str,
    *,
    include_prime_range: bool = False,
    include_sym_power: bool = False,
    include_sage: bool = False,
):
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-p", "--p", type=int, default=5, help="Field characteristic")
    parser.add_argument("-n", "--n", type=int, default=1, help="Extension degree")
    parser.add_argument("-N", "--N", type=int, default=11, help="Level N")
    parser.add_argument("-type", "--type", type=int, default=1, help="Level type")
    
    if include_prime_range:
        parser.add_argument("--pmin", type=int, default=2, help="Lower prime bound")
        parser.add_argument("--pmax", type=int, default=100, help="Upper prime bound")
        parser.add_argument("--plist", type=int, nargs="*", help="Explicit list of primes")
        parser.add_argument(
            "--random", action="store_true", help="Pick one random prime from prange"
        )

    if include_sym_power:
        parser.add_argument("-k", type=int, default=2, help="Symmetric power k")

    if include_sage:
        parser.add_argument(
            "--sage", action="store_true", help="Compare against Sage's Hecke operator"
        )

    add_config_args(parser)
    return parser.parse_args()
