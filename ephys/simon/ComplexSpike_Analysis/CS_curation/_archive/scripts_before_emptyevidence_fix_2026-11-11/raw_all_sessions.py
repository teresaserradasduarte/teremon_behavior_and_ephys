# -*- coding: utf-8 -*-
"""
raw_all_sessions.py — build the raw-data layer (Raw data tab, single duplicate events, SS same-cell check) for every
step0 v7 session, one session at a time: hydrate the raw .ap.bin in Dropbox → copy it to C: → process → delete the
C: copy → set the Dropbox file back to online-only → next session.

Dropbox: read-only EXCEPT the sync state of the session's raw .ap.bin (Simon's explicit, narrow exception of
2026-11-09): `attrib +P` (keep on this device → Dropbox downloads it) before the copy and `attrib -P +U`
(online-only) right after it. Content is never changed, moved, renamed or deleted; no other Dropbox file is touched.
Files that were already available offline before the batch (e.g. CoteDor R5) are left exactly as they are.

Per session, in order (env RAW_SESSION / RAW_BIN / RAW_META, logs in CS_curation/_raw_batch_logs/<session>/):
  raw_waveforms_session.py → raw_waveforms_figures.py → raw_dup_events.py → ss_same_cell_check.py
Resumable: a session with <session>/raw_waveforms/_raw_batch_done.json is skipped. On a failing script the batch
stops and the local copy is KEPT (so the fix can be rerun without copying again).
Usage: python raw_all_sessions.py [session ...]    (default: all v7 sessions except 1_CoteDor_R5)
"""
from pathlib import Path
import json
import os
import re
import shutil
import stat
import subprocess
import sys
import time

CODE = Path(__file__).resolve().parent
LOCAL_ROOT = Path(r"C:\Users\simon\Desktop\Data\ephys_and_behavior")
CUR = LOCAL_ROOT / "CS_curation"
RAW_COPIES = LOCAL_ROOT / "raw_copies"
LOG_ROOT = CUR / "_raw_batch_logs"
STEPS = ["raw_waveforms_session.py", "raw_waveforms_figures.py", "raw_dup_events.py", "ss_same_cell_check.py"]
OUTPUTS = ["raw_waveforms/raw_data.js", "raw_waveforms/raw_dup_events.js", "ss_same_cell/ss_same_cell.js"]
SKIP_DEFAULT = {"1_CoteDor_R5"}             # done interactively; its local copy is kept (Simon, 2026-11-09)
DROPBOX = "D:\\Learning Lab Dropbox\\"
OFFLINE, PINNED, RECALL = 0x1000, 0x80000, 0x400000
HYDRATE_TIMEOUT = 6 * 3600
MIN_FREE_C = 20e9                            # keep this much free on C: after the copy


def log(msg, f=None):
    line = f"{time.strftime('%Y-%m-%d %H:%M:%S')}  {msg}"
    print(line, flush=True)
    with open(LOG_ROOT / "batch.log", "a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def attrs(p):
    return os.stat(p).st_file_attributes


def meta_of(b):
    return dict(l.strip().split("=", 1) for l in b.with_suffix(".meta").read_text().splitlines() if "=" in l)


def resolve_bin(sess):
    """The raw .ap.bin the sorter ran on (dat_path name; size must match its .meta; spikes must fit)."""
    import numpy as np
    p0 = json.loads((CUR / sess / "step0_explore_bhvspan_v7batch" / "params_used.json").read_text())
    sd = Path(p0["sorter_dir"])
    txt = (sd / "params.py").read_text()
    dat = re.search(r"dat_path\s*=\s*r?['\"](.+?)['\"]", txt)
    name = Path(dat.group(1).replace("\\\\", "\\")).name if dat else None
    smax = int(np.load(sd / "spike_times.npy", mmap_mode="r").max())
    cands = [b for d in (sd.parent, sd) for b in sorted(d.glob("*.ap.bin")) if b.with_suffix(".meta").exists()]
    ok = []
    for b in cands:
        m = meta_of(b)
        ns = int(m["fileSizeBytes"]) // (2 * int(m["nSavedChans"]))
        if int(m["fileSizeBytes"]) == b.stat().st_size and smax < ns:
            ok.append(b)
    ok.sort(key=lambda b: (b.name != name, b.parent != sd.parent))       # dat_path name first, then the imec folder
    assert ok, f"{sess}: no raw .ap.bin matching the sorter in {sd.parent}"
    return ok[0]


def is_local(p):
    a = attrs(p)
    return not (a & OFFLINE) and not (a & RECALL)


TRIGGER = """import sys
try:
    open(sys.argv[1], 'rb').read(1)
except OSError:
    pass
"""
MAX_TRIGGERS = 5                             # Simon 2026-11-09: re-trigger a dropped download up to 5x
STALL_S = 15 * 60                            # no D: space used for this long while not on device → re-trigger


def trigger(src):
    """Dropbox ignores the pin flag alone; one read of the file starts the download (the read itself times out)."""
    return subprocess.Popen([sys.executable, "-c", TRIGGER, str(src)], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def hydrate(src):
    """Pin (keep on device), trigger the download, wait until the file is local. Returns True if WE pinned it."""
    assert str(src).startswith(DROPBOX) and src.name.endswith(".ap.bin")
    if is_local(src):
        log(f"  already on this device (left as it is): {src.name}")
        return False
    subprocess.run(["attrib", "+P", str(src)], check=True)
    log(f"  pinned (attrib +P) {src.name}, triggering the Dropbox download of {src.stat().st_size / 1e9:.1f} GB")
    t0, last, n_trig = time.time(), 0, 1
    trigger(src)
    free_min, t_prog = shutil.disk_usage("D:/").free, time.time()
    while not is_local(src):
        if time.time() - t0 > HYDRATE_TIMEOUT:
            raise RuntimeError(f"hydration timeout for {src}")
        f = shutil.disk_usage("D:/").free
        if f < free_min - 1e9:
            free_min, t_prog = f, time.time()
        if time.time() - t_prog > STALL_S:
            try:                                     # fallback: downloaded but the attribute did not change?
                subprocess.run([sys.executable, "-c", TRIGGER.replace("pass", "sys.exit(1)"), str(src)], timeout=60, check=True)
                log("    file readable although still flagged online-only: treating it as on device"); break
            except (subprocess.TimeoutExpired, subprocess.CalledProcessError):
                pass
            if n_trig >= MAX_TRIGGERS:
                raise RuntimeError(f"download stalled {MAX_TRIGGERS} times: {src}")
            n_trig += 1; trigger(src); free_min, t_prog = f, time.time()
            log(f"    download stalled, re-triggered ({n_trig}/{MAX_TRIGGERS})")
        if time.time() - last > 600:
            log(f"    … downloading ({(time.time() - t0) / 60:.0f} min, D: free {f / 1e9:.0f} GB)")
            last = time.time()
        time.sleep(30)
    log(f"  on device after {(time.time() - t0) / 60:.0f} min")
    return True


def release(src):
    """Back to online-only (frees D:). Only for files this batch pinned."""
    assert str(src).startswith(DROPBOX) and src.name.endswith(".ap.bin")
    subprocess.run(["attrib", "-P", "+U", str(src)], check=True)
    log(f"  set back to online-only (attrib -P +U) {src.name}")


def free_c_space(sess, keep=True):
    """After a failure: keep the session's local copy only if C: still has room for the largest next copy."""
    d = RAW_COPIES / sess
    if not d.exists():
        return
    if keep and shutil.disk_usage("C:/").free > 100e9 + MIN_FREE_C:
        log(f"  local copy of {sess} kept for the retry"); return
    assert str(d).startswith(str(RAW_COPIES)) and not str(d).startswith(DROPBOX)
    shutil.rmtree(d); log(f"  local copy of {sess} deleted (space)")


def run_session(sess):
    ld = LOG_ROOT / sess; ld.mkdir(parents=True, exist_ok=True)
    src = resolve_bin(sess)
    dst_dir = RAW_COPIES / sess
    dst = dst_dir / src.name
    log(f"== {sess}: {src}")
    if not (dst.exists() and dst.stat().st_size == src.stat().st_size):
        need = src.stat().st_size - (dst_dir / (src.name + ".part")).stat().st_size if (dst_dir / (src.name + ".part")).exists() else src.stat().st_size
        assert shutil.disk_usage("C:/").free - need > MIN_FREE_C, "not enough space on C:"
        pinned = hydrate(src)
        t0 = time.time()
        r = subprocess.run([sys.executable, str(CODE / "copy_raw_local.py"), str(src), str(dst_dir)],
                           stdout=open(ld / "copy.log", "w", encoding="utf-8"), stderr=subprocess.STDOUT)
        if pinned:
            release(src)          # also after a failed copy
        if r.returncode:
            raise RuntimeError(f"{sess}: copy failed, see {ld / 'copy.log'}")
        log(f"  copied to C: in {(time.time() - t0) / 60:.0f} min")
    else:
        log("  local copy already complete")
    env = dict(os.environ, RAW_SESSION=sess, RAW_BIN=str(dst), RAW_META=str(dst.with_suffix(".meta")),
               MPLBACKEND="Agg", PYTHONIOENCODING="utf-8")
    for s in STEPS:
        t0 = time.time()
        r = subprocess.run([sys.executable, str(CODE / s)], env=env, cwd=str(CODE),
                           stdout=open(ld / f"{s[:-3]}.log", "w", encoding="utf-8"), stderr=subprocess.STDOUT)
        if r.returncode:
            raise RuntimeError(f"{sess}: {s} failed (exit {r.returncode}), see {ld / (s[:-3] + '.log')}; local copy kept")
        log(f"  {s} ok ({(time.time() - t0) / 60:.1f} min)")
    missing = [o for o in OUTPUTS if not (CUR / sess / o).exists()]
    if missing:
        raise RuntimeError(f"{sess}: outputs missing {missing}; local copy kept")
    # delete the LOCAL copy only (C:), never anything on Dropbox
    assert str(dst).startswith(str(RAW_COPIES)) and not str(dst).startswith(DROPBOX)
    dst.unlink(); dst.with_suffix(".meta").unlink(missing_ok=True)
    try:
        dst_dir.rmdir()
    except OSError:
        pass
    (CUR / sess / "raw_waveforms" / "_raw_batch_done.json").write_text(json.dumps(
        dict(session=sess, raw_bin=str(src), finished=time.strftime("%Y-%m-%d %H:%M"), steps=STEPS), indent=1))
    log(f"  done, local copy deleted (C: free {shutil.disk_usage('C:/').free / 1e9:.0f} GB)")


if __name__ == "__main__":
    LOG_ROOT.mkdir(parents=True, exist_ok=True)
    import ctypes                                    # keep the PC from idle-sleeping while this process runs
    ctypes.windll.kernel32.SetThreadExecutionState(0x80000000 | 0x00000001)   # ES_CONTINUOUS | ES_SYSTEM_REQUIRED
    sessions = sys.argv[1:] or [d.parent.name for d in sorted(CUR.glob("*/step0_explore_bhvspan_v7batch"))
                                if d.parent.name not in SKIP_DEFAULT]
    todo = [s for s in sessions if not (CUR / s / "raw_waveforms" / "_raw_batch_done.json").exists()]
    log(f"batch: {len(todo)} sessions to do ({len(sessions) - len(todo)} already done): {', '.join(todo)}")
    failed = []
    for s in todo:                                   # Simon 2026-11-09: a failing session is retried at the end
        try:
            run_session(s)
        except Exception as e:
            log(f"  FAILED {s}: {e}"); failed.append(s); free_c_space(s)
    if failed:
        log(f"retrying at the end: {', '.join(failed)}")
        still = []
        for s in failed:
            try:
                run_session(s)
            except Exception as e:
                log(f"  FAILED again {s}: {e}"); still.append(s); free_c_space(s, keep=False)
        failed = still
    log(f"batch finished; failed: {', '.join(failed) if failed else 'none'}")
