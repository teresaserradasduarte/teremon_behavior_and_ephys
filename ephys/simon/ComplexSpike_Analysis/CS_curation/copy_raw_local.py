# -*- coding: utf-8 -*-
"""
copy_raw_local.py — copy one raw .ap.bin (+ .ap.meta) from Dropbox to a LOCAL folder, read-only on Dropbox.

Plain buffered sequential read (Dropbox hydrates the online-only file as it is read) → write to a .part file on C:,
renamed when the byte count matches the meta's fileSizeBytes. Re-running resumes from the .part file's size.
Usage: python copy_raw_local.py <src .ap.bin> <dst folder>
"""
import sys
import shutil
import time
from pathlib import Path

src, dst_dir = Path(sys.argv[1]), Path(sys.argv[2])
dst_dir.mkdir(parents=True, exist_ok=True)
meta_src = src.with_suffix(".meta")
shutil.copyfile(meta_src, dst_dir / meta_src.name)                       # small; read from Dropbox, written locally
size = int(dict(l.strip().split("=", 1) for l in meta_src.read_text().splitlines() if "=" in l)["fileSizeBytes"])
dst, part = dst_dir / src.name, dst_dir / (src.name + ".part")
if dst.exists() and dst.stat().st_size == size:
    print("already complete:", dst); sys.exit(0)
done = part.stat().st_size if part.exists() else 0
CH = 64 * 1024 * 1024
t0, last = time.time(), 0.0
with open(src, "rb") as fi, open(part, "ab") as fo:
    fi.seek(done)
    while True:
        for attempt in range(20):                                         # hydration hiccups: wait and retry
            try:
                b = fi.read(CH); break
            except OSError as e:
                print(f"read error at {done / 1e9:.2f} GB ({e}); retry {attempt + 1} in 30 s", flush=True)
                time.sleep(30); fi.seek(done)
        else:
            raise SystemExit("giving up")
        if not b:
            break
        fo.write(b); done += len(b)
        if time.time() - last > 30:
            el = time.time() - t0; last = time.time()
            print(f"{done / 1e9:6.1f} / {size / 1e9:.1f} GB  ({done / max(el, 1e-9) / 1e6:.0f} MB/s avg)", flush=True)
assert done == size, f"size mismatch {done} vs {size}"
part.rename(dst)
print(f"done: {dst} ({size / 1e9:.1f} GB, {(time.time() - t0) / 60:.1f} min)")
