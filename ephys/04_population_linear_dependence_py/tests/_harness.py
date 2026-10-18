"""Minimal test runner (no pytest in the neuro_mi env): runs every test_* function of
the calling module, prints PASS/FAIL with the evidence each test prints, and exits
non-zero if anything failed. Usage at the bottom of a test file:
    if __name__ == "__main__": run(globals())
"""
import sys, time, traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def run(ns):
    tests = [(k, v) for k, v in ns.items() if k.startswith("test_") and callable(v)]
    n_fail = 0
    for name, fn in tests:
        t0 = time.time()
        try:
            fn()
            print(f"PASS  {name}  ({time.time()-t0:.1f}s)")
        except Exception:
            n_fail += 1
            print(f"FAIL  {name}")
            traceback.print_exc()
    print(f"\n{len(tests) - n_fail}/{len(tests)} passed")
    sys.exit(1 if n_fail else 0)
