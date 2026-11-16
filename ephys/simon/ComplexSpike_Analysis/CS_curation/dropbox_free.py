# -*- coding: utf-8 -*-
r"""
dropbox_free.py — set a raw .ap.bin that Dropbox downloaded back to online-only, after we copied/used it.

Uses Dropbox's own right-click action "Make online-only" (shell verb, matched by its exact name; nothing else on the
menu is ever invoked). The file becomes a 0-byte placeholder again; its content in the cloud is unchanged.
Dropbox moves the local bytes to .dropbox.cache/old_files and clears that itself (~3 days, Dropbox Help Center).
NOTHING under the Dropbox folder is ever deleted by this script (Simon, 2026-11-16: "Never ever delete anything from
dropbox directly"). Before the action we wait until Dropbox is not downloading anything (no file > 100 MB in
.dropbox.cache/new_files): while a download ran, the action was silently ignored (2026-11-16, FerreroRocher R4).
Usage: python dropbox_free.py <raw .ap.bin under D:\Learning Lab Dropbox\> [...]
"""
from pathlib import Path
import os
import shutil
import subprocess
import sys
import time

DROPBOX = "D:\\Learning Lab Dropbox\\"
CACHE_NEW = Path(DROPBOX) / ".dropbox.cache" / "new_files"
IDLE_MIN_BYTES = 100e6
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


def downloading():
    """Large files Dropbox is downloading right now (it stages them in .dropbox.cache/new_files)."""
    try:
        return [f for f in CACHE_NEW.iterdir() if f.is_file() and f.stat().st_size > IDLE_MIN_BYTES]
    except OSError:
        return []


def wait_idle(log=print, poll_s=60, max_s=24 * 3600):
    t0, said = time.time(), False
    while downloading():
        if not said:
            log("  Dropbox is downloading; waiting until it is done before 'Make online-only'"); said = True
        if time.time() - t0 > max_s:
            raise RuntimeError("Dropbox still downloading after 24 h")
        time.sleep(poll_s)
    return True


def free(src, log=print, wait_s=300):
    """Make `src` online-only (no deletion of anything). Returns the change of D: free space in GB (usually ~0: the
    bytes stay in Dropbox's cache until Dropbox clears it)."""
    src = Path(src)
    assert str(src).startswith(DROPBOX) and src.name.endswith(".ap.bin") and ".dropbox.cache" not in str(src), src
    st = src.stat()
    size, mtime = st.st_size, st.st_mtime
    free0 = shutil.disk_usage("D:/").free
    if not online_only(src):
        wait_idle(log)
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
    return (shutil.disk_usage("D:/").free - free0) / 1e9


if __name__ == "__main__":
    for a in sys.argv[1:]:
        free(a)
