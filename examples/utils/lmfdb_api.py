from __future__ import annotations

"""Simple LMFDB trace fetcher for the full Gamma_1(N) trace sequence."""

import os
import pickle
import time

import requests


_API_URL = "https://www.lmfdb.org/api/mf_gamma1/"
_HEADERS = {"User-Agent": "x1fq-moduli/0.1", "Accept": "application/json"}
_CACHE_DIR = os.path.join(os.path.dirname(__file__), ".lmfdb_cache")
_LAST_REQUEST_TIME = 0.0
_MIN_REQUEST_DELAY = 1.0

os.makedirs(_CACHE_DIR, exist_ok=True)

def fetch_traces(
    level: int,
    weight: int,
    *,
    q: int | None = None,
    p: int | None = None,
    pmin: int | None = None,
    pmax: int | None = None,
    verbose: bool = False,
    use_cache: bool = True,
) -> int | dict[int, int]:
    """Fetch full-space Hecke traces for one index or for primes in a range.

    Exactly one selection mode must be used:
    - `q=<int>` for one exact trace index
    - `p=<prime>` for one prime trace
    - `pmin=<int>, pmax=<int>` for all primes in that interval
    """

    indices, return_scalar = _requested_indices(q=q, p=p, pmin=pmin, pmax=pmax)
    trace_sums = _fetch_trace_sums(
        level,
        weight,
        max_n=max(indices),
        verbose=verbose,
        use_cache=use_cache,
    )

    if return_scalar:
        return trace_sums.get(indices[0], 0)
    return {index: trace_sums.get(index, 0) for index in indices}


def _requested_indices(
    *,
    q: int | None,
    p: int | None,
    pmin: int | None,
    pmax: int | None,
) -> tuple[list[int], bool]:
    modes = [q is not None, p is not None, pmin is not None or pmax is not None]
    if sum(modes) != 1:
        raise ValueError("choose exactly one of q, p, or pmin/pmax")

    if q is not None:
        if q < 1:
            raise ValueError("q must be positive")
        return [q], True

    if p is not None:
        if not _is_prime(p):
            raise ValueError("p must be prime")
        return [p], True

    if pmin is None or pmax is None:
        raise ValueError("both pmin and pmax are required for a prime range")
    if pmin > pmax:
        raise ValueError("pmin must be <= pmax")

    primes = [candidate for candidate in range(max(2, pmin), pmax + 1) if _is_prime(candidate)]
    return primes, False


def _fetch_trace_sums(
    level: int,
    weight: int,
    *,
    max_n: int,
    verbose: bool,
    use_cache: bool,
) -> dict[int, int]:
    cache_file = os.path.join(_CACHE_DIR, f"traces_l{level}_w{weight}_n{max_n}.pkl")
    if use_cache and os.path.exists(cache_file):
        if verbose:
            print(f"Loading cached traces from {cache_file}")
        with open(cache_file, "rb") as handle:
            return pickle.load(handle)

    _rate_limit()

    if verbose:
        print(f"Fetching LMFDB Gamma_1 traces for level={level}, weight={weight}")

    response = requests.get(
        _API_URL,
        params={
            "level": level,
            "weight": weight,
            "_fields": "traces,label",
            "_format": "json",
        },
        headers=_HEADERS,
        timeout=30,
    )
    response.raise_for_status()

    content_type = response.headers.get("content-type", "")
    if "json" not in content_type:
        raise RuntimeError(
            "LMFDB returned a non-JSON response, likely a reCAPTCHA or rate-limit page"
        )

    data = response.json()
    spaces = data.get("data", [])
    if verbose:
        print(f"Fetched {len(spaces)} Gamma_1 space row(s)")

    trace_sums = {n: 0 for n in range(1, max_n + 1)}
    if not spaces:
        return trace_sums

    traces = spaces[0].get("traces", [])
    upper = min(len(traces), max_n)
    for n in range(upper):
        trace_sums[n + 1] = traces[n]

    if use_cache:
        with open(cache_file, "wb") as handle:
            pickle.dump(trace_sums, handle)

    return trace_sums


def _rate_limit() -> None:
    global _LAST_REQUEST_TIME

    elapsed = time.time() - _LAST_REQUEST_TIME
    if elapsed < _MIN_REQUEST_DELAY:
        time.sleep(_MIN_REQUEST_DELAY - elapsed)
    _LAST_REQUEST_TIME = time.time()


def _is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n == 2:
        return True
    if n % 2 == 0:
        return False

    factor = 3
    while factor * factor <= n:
        if n % factor == 0:
            return False
        factor += 2
    return True


def main() -> None:
    print(fetch_traces(level=11, weight=2, p=5, verbose=True))
    print(fetch_traces(level=11, weight=2, pmin=2, pmax=11, verbose=True))


if __name__ == "__main__":
    main()


