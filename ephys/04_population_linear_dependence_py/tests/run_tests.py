"""Run every tests/test_*.py in its own process; non-zero exit if any fails.
    python tests/run_tests.py"""
import subprocess, sys
from pathlib import Path

here = Path(__file__).resolve().parent
failed = []
for f in sorted(here.glob("test_*.py")):
    print(f"\n==== {f.name} " + "=" * (60 - len(f.name)), flush=True)
    r = subprocess.run([sys.executable, str(f)], cwd=str(here))
    if r.returncode:
        failed.append(f.name)
print("\n" + ("ALL PASSED" if not failed else f"FAILED: {failed}"))
sys.exit(1 if failed else 0)
