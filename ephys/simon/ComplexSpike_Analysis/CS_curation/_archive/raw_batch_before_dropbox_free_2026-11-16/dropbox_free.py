# -*- coding: utf-8 -*-
"""
dropbox_free.py — give back the D: space of a raw .ap.bin that Dropbox downloaded, after we copied/used it.

Simon's explicit permission (2026-11-16), for the raw .ap.bin files of the step0 v7 sessions ONLY:
  1. Dropbox's own right-click action "Make online-only" (shell verb, matched by its exact name; nothing else on the
     menu is ever invoked) → the file goes back to a 0-byte placeholder; content in the cloud is unchanged.
  2. Dropbox then moves the local bytes to D:\\Learning Lab Dropbox\\.dropbox.cache\\old_files\\<hash> (it would clear them
     after ~3 days). We delete exactly that cache copy: same byte size AND same modified time as the raw file, exactly
     one match, inside .dropbox.cache\\old_files only. Nothing else under Dropbox is changed or deleted.
Usage: python dropbox_free.py <raw .ap.bin under D:\\Learning Lab Dropbox\\> [...]
"""
from pathlib import Path
import os
import shutil
import subprocess
import sys
import time

DROPBOX = "D:\\Learning Lab Dropbox\\"
CACHE_OLD = Path(DROPBOX) / ".dropbox.cache" / "old_files"
OFFLINE, REPARSE = 0x1000, 0x400
VERB = "Make online-only"
# PowerShell 7 (MTA): the verb runs; Windows PowerShell 5.1 (STA) queues it on a message loop that never runs
PWSH = shutil.which("pwsh") or r"C:\Program Files\PowerShell\pwsh.exe"

PS = r"""
$d = $env:DBX_DIR; $n = $env:DBX_NAME
$sh = New-Object -ComObject Shell.Application
$item = $sh.Namespace($d).ParseName($n)
if ($null -eq $item) { Write-Output "NOITEM"; exit 2 }
$v = @($item.Verbs() | Where-Object { $_.Name -ceq "Make online-only" })
if ($v.Count -ne 1) { Write-Output "VERB_COUNT $($v.Count)"; exit 3 }
$v[0].DoIt()
# the verb runs asynchronously inside this process: stay alive until Dropbox has made the file a placeholder
$f = Join-Path $d $n
for ($i = 0; $i -lt 60; $i++) { Start-Sleep 2; if ([int64](Get-Item -LiteralPath $f -Force).Attributes.value__ -band 0x1000) { break } }
Write-Output "OK"
"""


def attrs(p):
    return os.stat(p).st_file_attributes


def online_only(p):
    return bool(attrs(p) & OFFLINE)


def free(src, log=print, wait_s=300):
    """Make `src` online-only and delete its cache copy. Returns GB given back on D:."""
    src = Path(src)
    assert str(src).startswith(DROPBOX) and src.name.endswith(".ap.bin") and ".dropbox.cache" not in str(src), src
    st = src.stat()
    size, mtime = st.st_size, st.st_mtime
    free0 = shutil.disk_usage("D:/").free
    if not online_only(src):
        r = subprocess.run([PWSH, "-NoProfile", "-NonInteractive", "-Command", PS], capture_output=True, text=True,
                           timeout=120, env=dict(os.environ, DBX_DIR=str(src.parent), DBX_NAME=src.name))
        if "OK" not in r.stdout:
            raise RuntimeError(f"'{VERB}' failed for {src.name}: {r.stdout.strip()} {r.stderr.strip()[:300]}")
        t0 = time.time()
        while not online_only(src):
            if time.time() - t0 > wait_s:
                raise RuntimeError(f"{src.name} not online-only {wait_s} s after '{VERB}'")
            time.sleep(3)
        log(f"  made online-only (Dropbox '{VERB}'): {src.name}")
    # its bytes now sit in .dropbox.cache/old_files: find exactly that copy
    t0 = time.time()
    while True:
        hits = [f for f in CACHE_OLD.iterdir() if f.is_file() and f.stat().st_size == size and abs(f.stat().st_mtime - mtime) < 2]
        if len(hits) == 1 or time.time() - t0 > 120:
            break
        time.sleep(5)
    if len(hits) != 1:
        log(f"  cache copy of {src.name}: {len(hits)} matches (size {size}, mtime) — nothing deleted")
        return (shutil.disk_usage("D:/").free - free0) / 1e9
    h = hits[0]
    assert str(h).startswith(str(CACHE_OLD) + os.sep) and h.parent == CACHE_OLD
    h.unlink()
    time.sleep(2)
    gb = (shutil.disk_usage("D:/").free - free0) / 1e9
    log(f"  deleted its Dropbox cache copy {h.name} ({size / 1e9:.1f} GB); D: +{gb:.1f} GB, free {shutil.disk_usage('D:/').free / 1e9:.0f} GB")
    return gb


if __name__ == "__main__":
    for a in sys.argv[1:]:
        free(a)
