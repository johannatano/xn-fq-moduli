from __future__ import annotations

from typing import List


def next_primes(start: int, count: int) -> List[int]:
    """Return `count` primes >= start.

    Prefer SymPy's `nextprime` when available for speed and correctness.
    Falls back to a simple primality test when SymPy is not installed.
    """
    start = max(2, int(start))
    cnt = max(0, int(count))
    if cnt == 0:
        return []

    try:
        from sympy import nextprime

        res: List[int] = []
        cur = nextprime(start - 1)
        for _ in range(cnt):
            res.append(int(cur))
            cur = nextprime(cur)
        return res
    except Exception:
        # simple primality check fallback
        def is_prime(x: int) -> bool:
            if x < 2:
                return False
            if x % 2 == 0:
                return x == 2
            r = int(x ** 0.5)
            i = 3
            while i <= r:
                if x % i == 0:
                    return False
                i += 2
            return True

        res = []
        cand = start
        while len(res) < cnt:
            if is_prime(cand):
                res.append(cand)
            cand += 1
        return res
