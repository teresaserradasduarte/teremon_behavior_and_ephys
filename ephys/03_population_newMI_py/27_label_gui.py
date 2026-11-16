"""27_label_gui.py - label the 2D MI maps by hand, in the browser.

Run from a terminal (not the VS Code interactive window), in the neuro_mi env:

    python 27_label_gui.py              first pass (200 maps)
    python 27_label_gui.py --repeat     repeat pass (50 maps, new order)
    python 27_label_gui.py --relabel "Ferrero R1 224" "Lindt R6 996"
                                        correct the label of chosen maps (any map of 12,
                                        also ones outside the label set: those are ADDED
                                        as new labelled maps). Or give a .txt file with
                                        one map per line, e.g. relabel_<tag>.txt from 28.
                                        In this mode the map's name is shown.

This starts a small local server (only on this computer, 127.0.0.1) and opens the
labelling page in your browser. Stop it with the "Done" button or Ctrl+C.

The maps come from 26_label_set.py (label_set_<tag>.pkl under the data root). The
labels go to labels_<tag>.csv next to it, written after EVERY label, so you can stop
at any time; reopening resumes at the first unlabelled map.

Nothing that identifies a map is shown (no uid, region, mouse, or 23/25 class).

Keys
    1-7        label (column, diagonal, tiling, 2 tiles, broad, nothing, mixed/unsure)
    Left/Right previous / next map        (or the slider / the strip under it)
    v          cycle the view: excess MI / z / smoothed z
    Enter      (in the note box) back to the map
"""
import os, re, sys, csv, time, json, argparse, pickle, threading, webbrowser
from http.server import ThreadingHTTPServer, BaseHTTPRequestHandler
from urllib.parse import urlparse, parse_qs
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path: sys.path.insert(0, _HERE)
os.chdir(_HERE)
import config

DEFAULT_TAG = "x-y-z-speed-dist_water_nw70_cat1"
FIELDS = ("pass", "pos", "key", "uid", "kind", "stratum", "label", "view", "note", "time", "prev")
# prev: the label before the last change (relabelling), empty otherwise
# set here, not taken from the label set, so a set made before "2 tiles" existed still works
LABELS = ("column", "diagonal", "tiling", "2 tiles", "broad", "nothing", "mixed/unsure")
LABEL_CLR = ("#3b82f6", "#f59e0b", "#10b981", "#ec4899", "#8b5cf6", "#9ca3af", "#a16207")


# ── data and labels ─────────────────────────────────────────────────────────
def load_set(path):
    with open(path, "rb") as fh:
        return pickle.load(fh)


def labels_path(set_path):
    p = str(set_path)
    return os.path.join(os.path.dirname(p),
                        os.path.basename(p).replace("label_set_", "labels_").replace(".pkl", ".csv"))


def read_labels(path):
    """-> {(pass, key): row dict}"""
    out = {}
    if os.path.isfile(path):
        with open(path, newline="", encoding="utf-8") as fh:
            for row in csv.DictReader(fh):
                out[(row["pass"], row["key"])] = row
    return out


def write_labels(path, rows, tries=8):
    """whole file to a temporary file, then swapped in (Dropbox may lock it: retry)"""
    tmp = path + ".tmp"
    with open(tmp, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, FIELDS)
        w.writeheader()
        for r in sorted(rows.values(), key=lambda r: (r["pass"], int(r["pos"]))):
            w.writerow(r)
    for t in range(tries):
        try:
            os.replace(tmp, path); return True
        except OSError:
            time.sleep(0.3 * (t + 1))
    return False


def _viridis_lut(n=256):
    from matplotlib import colormaps
    return (colormaps["viridis"](np.linspace(0, 1, n))[:, :3] * 255).round().astype(int).tolist()


def _map_json(M):
    """one view -> {data (rows, bottom row first, NaN -> null), vmin, vmax} (1st-99th pct)"""
    ok = np.isfinite(M)
    lo, hi = (float(np.nanpercentile(M, 1)), float(np.nanpercentile(M, 99))) if ok.any() else (0., 1.)
    if hi <= lo: hi = lo + 1e-12
    data = [[float(v) if f else None for v, f in zip(r, fr)] for r, fr in zip(M, ok)]
    return {"data": data, "vmin": lo, "vmax": hi}


# ── which maps: a pass of the label set, or a relabel list ──────────────────
def pass_entries(S, repeat):
    """the first pass (or the repeat pass): entries in the order of 26"""
    seq = S["repeat"] if repeat else S["order"]
    return [dict(item=S["items"][int(j)], V=S["views"][int(j)], pos=p) for p, j in enumerate(seq)]


def parse_spec(line):
    """'CB 5_FerreroRocher R4 278', 'Ferrero R1 224 (shuffle)', 'real:CB|...|R4|278' or a uid
    -> (region or None, mouse part, session, unit, kind); text after '#' is ignored"""
    line = line.split("#")[0].strip()
    if not line: return None
    kind = "shuffle" if "shuffle" in line.lower() else "real"
    for pre in ("real:", "shuffle:"):
        if line.startswith(pre): line = line[len(pre):]
    tok = [t for t in re.split(r"[|\s,]+", line.replace("(shuffle)", "")) if t and t.lower() != "shuffle"]
    if len(tok) == 4: return tok[0], tok[1], tok[2], tok[3], kind
    if len(tok) == 3: return None, tok[0], tok[1], tok[2], kind
    raise ValueError(f"cannot read the map {line!r}: give 'region mouse session unit' or 'mouse session unit'")


def _match(uid, spec):
    reg, mo, se, un, _ = spec
    r, m, s, u = uid.split("|")
    return (reg is None or reg == r) and mo.lower() in m.lower() and s == se and u == un


def relabel_entries(S, specs):
    """maps to (re)label. Maps of the label set keep their views; any other map of 12 is
    built here (same three views as 26) and saved as a new first-pass label ('added')."""
    from scipy.ndimage import gaussian_filter
    out, C12, C23 = [], None, None
    pos_of = {int(j): p for p, j in enumerate(S["order"])}
    for spec, group in specs:
        hit = [j for j, it in enumerate(S["items"]) if it["kind"] == spec[4] and _match(it["uid"], spec)]
        if hit:
            out += [dict(item=S["items"][j], V=S["views"][j], pos=pos_of[j], group=group) for j in hit]
            continue
        if C12 is None:
            print("loading 12's maps for maps outside the label set ...")
            with open(config.var_path(S["params"]["SRC_12"]), "rb") as fh: C12 = pickle.load(fh)
        rows = [i for i, u in enumerate(C12["uids"]) if _match(u, spec)]
        if not rows: print(f"  NOT FOUND: {' '.join(x for x in spec[:4] if x)} -- skipped"); continue
        if len(rows) > 1: print(f"  {spec}: {len(rows)} maps match, all added")
        for i in rows:
            uid = C12["uids"][i]
            if spec[4] == "shuffle":
                if C23 is None:
                    with open(config.var_path(S["params"]["SRC_23"]), "rb") as fh: C23 = pickle.load(fh)
                M = C23["sh"][list(C23["key"]["uids"]).index(uid)]
            else:
                M = C12["mi"][i]
            with np.errstate(all="ignore"):
                exc = M - C12["null"][i]
                z = exc / np.where(C12["null_std"][i] > 0, C12["null_std"][i], np.nan)
                ok = np.isfinite(z)
                num = gaussian_filter(np.where(ok, z, 0.0), S["params"]["SMOOTH_SIG"])
                den = gaussian_filter(ok.astype(float), S["params"]["SMOOTH_SIG"])
                zs = np.where(ok, num / den, np.nan)
            item = dict(uid=uid, region=uid.split("|")[0], kind=spec[4], stratum="added",
                        key=f"{spec[4]}:{uid}")
            out.append(dict(item=item, V=np.stack([exc, z, zs]).astype(np.float32), pos=10000 + i,
                            group=group))
    seen, uniq = set(), []
    for e in out:
        if e["item"]["key"] not in seen: seen.add(e["item"]["key"]); uniq.append(e)
    return uniq


def read_specs(args):
    """--relabel arguments: map names, or files with one map per line. In a file, a line
    '## title' starts a group (shown as a section of the gallery) -> [(spec, group)]"""
    specs = []
    for a in args:
        if not os.path.isfile(a):
            specs.append((parse_spec(a), "")); continue
        group = os.path.basename(a)
        for line in open(a, encoding="utf-8").read().splitlines():
            if line.startswith("##"): group = line.lstrip("#").strip(); continue
            s = parse_spec(line)
            if s: specs.append((s, group))
    return [(s, g) for s, g in specs if s]


# ── the labelling session (no web code) ─────────────────────────────────────
class Session:
    def __init__(self, S, entries, lpath, PASS="1", mode="first pass"):
        self.S, self.E, self.lpath, self.lock = S, entries, lpath, threading.Lock()
        self.PASS, self.mode = PASS, mode
        self.views = list(S["view_names"])
        self.rows = read_labels(lpath)

    def _row(self, pos):
        return self.rows.get((self.PASS, self.E[pos]["item"]["key"]))

    def state(self):
        rs = [self._row(p) for p in range(len(self.E))]
        return {"labels": [r["label"] if r else "" for r in rs],
                "notes": [r["note"] if r else "" for r in rs]}

    def meta(self):
        st = self.state()
        start = 0 if self.mode == "relabel" else next((p for p, l in enumerate(st["labels"]) if not l), 0)
        names = ([e["item"]["uid"].replace("|", " ") + (" (shuffle)" if e["item"]["kind"] == "shuffle" else "")
                  for e in self.E] if self.mode == "relabel" else None)   # blind unless relabelling
        groups = [e.get("group", "") for e in self.E] if self.mode == "relabel" else None
        return dict(st, n=len(self.E), start=start, mode=self.mode, names=names, groups=groups,
                    label_names=list(LABELS), label_clr=list(LABEL_CLR), views=self.views,
                    lut=_viridis_lut(), wc=[float(v) for v in self.S["win_centers"]],
                    kc=[float(v) for v in self.S["kin_centers"]], csv=os.path.basename(self.lpath))

    def map(self, pos):
        V = self.E[pos]["V"]
        return {"pos": pos, "views": [_map_json(V[v]) for v in range(len(self.views))]}

    def label(self, pos, label, view, note):
        if label not in LABELS: raise ValueError(f"unknown label {label!r}")
        e = self.E[pos]; it = e["item"]
        with self.lock:
            old = self.rows.get((self.PASS, it["key"]))
            prev = old.get("prev", "") if old else ""
            if old and old["label"] != label: prev = old["label"]
            self.rows[(self.PASS, it["key"])] = {
                "pass": self.PASS, "pos": str(e["pos"]), "key": it["key"], "uid": it["uid"],
                "kind": it["kind"], "stratum": it["stratum"], "label": label, "view": view,
                "note": (note or "").strip(), "time": time.strftime("%Y-%m-%d %H:%M:%S"),
                "prev": prev}
            ok = write_labels(self.lpath, self.rows)
        return dict(self.state(), ok=ok)


# ── the local web server ────────────────────────────────────────────────────
def make_server(sess, port=0):
    class H(BaseHTTPRequestHandler):
        def log_message(self, *a): pass

        def _send(self, body, ctype="application/json", code=200):
            b = body.encode("utf-8") if isinstance(body, str) else body
            self.send_response(code)
            self.send_header("Content-Type", ctype + "; charset=utf-8")
            self.send_header("Content-Length", str(len(b)))
            self.send_header("Cache-Control", "no-store")
            self.end_headers(); self.wfile.write(b)

        def _json(self, obj, code=200):
            self._send(json.dumps(obj), code=code)

        def do_GET(self):
            u = urlparse(self.path)
            try:
                if u.path == "/": self._send(PAGE, "text/html")
                elif u.path == "/api/meta": self._json(sess.meta())
                elif u.path == "/api/map": self._json(sess.map(int(parse_qs(u.query)["pos"][0])))
                else: self._json({"error": "not found"}, 404)
            except Exception as e:
                self._json({"error": repr(e)}, 500)

        def do_POST(self):
            u = urlparse(self.path)
            try:
                body = json.loads(self.rfile.read(int(self.headers.get("Content-Length", 0))) or b"{}")
                if u.path == "/api/label":
                    self._json(sess.label(int(body["pos"]), body["label"], body.get("view", ""),
                                          body.get("note", "")))
                elif u.path == "/api/quit":
                    self._json({"ok": True})
                    threading.Thread(target=self.server.shutdown, daemon=True).start()
                else: self._json({"error": "not found"}, 404)
            except Exception as e:
                self._json({"error": repr(e)}, 500)

    return ThreadingHTTPServer(("127.0.0.1", port), H)


def run_gui(set_path, repeat, port=0, open_browser=True, lpath=None, relabel=None):
    S = load_set(set_path)
    if relabel:
        E = relabel_entries(S, read_specs(relabel))
        if not E: sys.exit("none of the maps to relabel were found")
        sess = Session(S, E, lpath or labels_path(set_path), "1", "relabel")
    else:
        sess = Session(S, pass_entries(S, repeat), lpath or labels_path(set_path),
                       "repeat" if repeat else "1", "repeat pass" if repeat else "first pass")
    srv = make_server(sess, port)
    url = f"http://127.0.0.1:{srv.server_address[1]}/"
    print(f"{sess.mode}: {len(sess.E)} maps, "
          f"{sum(bool(l) for l in sess.state()['labels'])} already labelled")
    if relabel:
        for e in sess.E:
            r = sess.rows.get(("1", e["item"]["key"]))
            print(f"  {e['item']['key']:<45} {'now: ' + r['label'] if r else 'new (not labelled yet)'}")
    print(f"labels -> {sess.lpath}")
    print(f"open {url}   (stop: 'Done' button or Ctrl+C)")
    if open_browser: threading.Timer(0.5, lambda: webbrowser.open(url)).start()
    try: srv.serve_forever()
    except KeyboardInterrupt: pass
    srv.server_close(); print("stopped")


# ── self-test: no browser ───────────────────────────────────────────────────
def selftest(set_path, out_dir):
    from urllib.request import urlopen, Request
    os.makedirs(out_dir, exist_ok=True)
    lp = os.path.join(out_dir, "labels_selftest.csv")
    if os.path.isfile(lp): os.replace(lp, lp + ".old")
    S = load_set(set_path)
    sess = Session(S, pass_entries(S, False), lp)
    srv = make_server(sess); base = f"http://127.0.0.1:{srv.server_address[1]}"
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    get = lambda p: json.loads(urlopen(base + p).read())
    post = lambda p, o: json.loads(urlopen(Request(base + p, json.dumps(o).encode(),
                                                   {"Content-Type": "application/json"})).read())
    assert b"<canvas" in urlopen(base + "/").read()
    m = get("/api/meta"); print(f"meta: n={m['n']}, start={m['start']}, labels={m['label_names']}")
    mp = get("/api/map?pos=0")
    d = np.array(mp["views"][1]["data"], dtype=float)
    assert d.shape == S["views"].shape[2:], d.shape
    assert np.allclose(d, S["views"][S["order"][0], 1], equal_nan=True, atol=1e-6)
    for pos, lb in enumerate(["2 tiles", "column", "mixed/unsure"]):
        r = post("/api/label", {"pos": pos, "label": lb, "view": "z", "note": "test, comma"})
        assert r["ok"] and r["labels"][pos] == lb
    r = post("/api/label", {"pos": 1, "label": "diagonal", "view": "smoothed z", "note": ""})
    back = read_labels(lp)
    rk = lambda p: back[("1", S["items"][S["order"][p]]["key"])]
    assert [rk(p)["label"] for p in range(3)] == ["2 tiles", "diagonal", "mixed/unsure"]
    assert rk(1)["prev"] == "column" and rk(0)["prev"] == ""
    assert Session(S, pass_entries(S, False), lp).meta()["start"] == 3
    post("/api/quit", {})
    # relabel: one map of the set (by name) + one map outside it (built from 12)
    it0 = S["items"][S["order"][0]]
    reg, mo, se, un = it0["uid"].split("|")
    inset = {it["uid"] for it in S["items"]}
    with open(config.var_path(S["params"]["SRC_12"]), "rb") as fh: _u = pickle.load(fh)["uids"]
    out_uid = next(u for u in _u if u not in inset)
    E = relabel_entries(S, read_specs([f"{mo.split('_')[-1]} {se} {un}  # a comment",
                                       out_uid.replace("|", " ")]))
    assert [e["item"]["uid"] for e in E] == [it0["uid"], out_uid], [e["item"]["uid"] for e in E]
    rs = Session(S, E, lp, "1", "relabel")
    assert rs.meta()["labels"][0] == "2 tiles" and rs.meta()["names"][1] == out_uid.replace("|", " ")
    assert E[1]["V"].shape == S["views"].shape[1:]
    rs.label(1, "diagonal", "z", ""); rs.label(0, "broad", "z", "")
    back = read_labels(lp)
    assert back[("1", "real:" + out_uid)]["stratum"] == "added"
    assert back[("1", it0["key"])]["label"] == "broad" and back[("1", it0["key"])]["prev"] == "2 tiles"
    print(f"relabel OK: {it0['uid']} 2 tiles -> broad (prev kept); added {out_uid} from 12")
    print(f"CSV round trip OK ({len(back)} rows, relabel overwrote); would resume at map 4")
    print("SELFTEST OK")


# ── the page ────────────────────────────────────────────────────────────────
PAGE = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Map labelling</title>
<style>
:root{--bg:#f4f5f7;--card:#fff;--ink:#1d2330;--mute:#6b7280;--line:#e4e7ec;--soft:#eef0f4;--acc:#2563eb}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);font:14px/1.45 system-ui,-apple-system,"Segoe UI",Roboto,sans-serif}
header{display:flex;align-items:center;gap:14px;padding:14px 24px 10px}
h1{font-size:16px;font-weight:650;margin:0}
.pill{font-size:12px;padding:2px 10px;border-radius:999px;background:var(--soft);color:var(--mute)}
.spacer{flex:1}
.progress{width:180px;height:6px;border-radius:3px;background:var(--line);overflow:hidden}
.progress div{height:100%;background:var(--acc);transition:width .2s}
main{display:grid;grid-template-columns:minmax(0,1fr) 270px;gap:18px;padding:0 24px 24px;max-width:1280px}
.card{background:var(--card);border:1px solid var(--line);border-radius:14px;padding:16px;box-shadow:0 1px 2px rgba(16,24,40,.04)}
.plotwrap{position:relative}
canvas#plot{width:100%;height:min(66vh,620px);display:block}
.maphead{display:flex;align-items:baseline;gap:10px;margin-bottom:6px}
.maphead b{font-size:15px}
.nav{display:flex;align-items:center;gap:10px;margin-top:10px}
.nav input[type=range]{flex:1;accent-color:var(--acc)}
.btn{border:1px solid var(--line);background:#fff;border-radius:8px;padding:6px 12px;cursor:pointer;font:inherit;color:var(--ink)}
.btn:hover{background:var(--soft)}
.btn.primary{background:var(--acc);border-color:var(--acc);color:#fff}
#strip{display:flex;gap:1px;height:12px;margin-top:8px;cursor:pointer}
#strip div{flex:1;background:var(--line);border-radius:2px;min-width:1px}
#strip div.cur{outline:2px solid var(--ink);outline-offset:1px}
.sec{font-size:11px;letter-spacing:.06em;text-transform:uppercase;color:var(--mute);margin:14px 0 6px}
.sec:first-child{margin-top:0}
.seg{display:flex;background:var(--soft);border-radius:9px;padding:3px}
.seg button{flex:1;border:0;background:transparent;padding:6px 4px;border-radius:7px;cursor:pointer;color:var(--mute);font:inherit;font-size:13px}
.seg button.on{background:#fff;color:var(--ink);box-shadow:0 1px 3px rgba(16,24,40,.12)}
.lab{display:flex;align-items:center;gap:10px;width:100%;padding:8px 10px;margin:5px 0;border:1px solid var(--line);
     border-radius:9px;background:#fff;cursor:pointer;font:inherit;text-align:left;color:var(--ink);transition:background .12s}
.lab:hover{background:var(--soft)}
.lab .dot{width:10px;height:10px;border-radius:50%;background:var(--c);flex:none}
.lab kbd{margin-left:auto;font:12px ui-monospace,Consolas,monospace;color:var(--mute);border:1px solid var(--line);border-radius:5px;padding:0 6px}
.lab.on{background:var(--c);border-color:var(--c);color:#fff;font-weight:600;box-shadow:inset 0 2px 4px rgba(0,0,0,.18)}
.lab.on .dot{background:#fff}
.lab.on kbd{color:#fff;border-color:rgba(255,255,255,.6)}
input#note{width:100%;padding:7px 10px;border:1px solid var(--line);border-radius:8px;font:inherit}
input#note:focus{outline:2px solid var(--acc);border-color:transparent}
label.chk{display:flex;align-items:center;gap:8px;margin-top:12px;color:var(--mute);cursor:pointer}
#status{margin-top:12px;font-size:12px;color:var(--mute);min-height:2.6em}
#status.err{color:#b42318;font-weight:600}
.keys{font-size:12px;color:var(--mute);margin-top:10px}
.keys kbd{font:11px ui-monospace,Consolas,monospace;border:1px solid var(--line);border-radius:4px;padding:0 4px;background:#fff}
#closed{display:none;padding:60px;text-align:center;color:var(--mute);font-size:16px}
#gallery{display:none;margin:0 24px 24px;max-width:1280px}
#gallery h3{font-size:13px;font-weight:600;margin:18px 0 8px;color:var(--ink)}
#gallery h3:first-child{margin-top:0}
#gallery h3 span{color:var(--mute);font-weight:400}
.grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(150px,1fr));gap:10px}
.thumb{border:1px solid var(--line);border-radius:10px;padding:6px;cursor:pointer;background:#fff;transition:border-color .12s,box-shadow .12s}
.thumb:hover{border-color:var(--acc);box-shadow:0 2px 8px rgba(37,99,235,.15)}
.thumb canvas{width:100%;height:105px;display:block;border-radius:6px;background:var(--soft)}
.thumb .nm{font-size:11px;color:var(--mute);margin-top:5px;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}
.chip{display:inline-block;font-size:11px;padding:1px 8px;border-radius:999px;margin-top:4px;color:#fff;background:var(--c)}
.chip.none{background:var(--soft);color:var(--mute)}
.btn.on{background:var(--ink);border-color:var(--ink);color:#fff}
@media (max-width:860px){main{grid-template-columns:1fr}}
</style></head>
<body>
<header>
  <h1>Map labelling</h1><span class="pill" id="pass"></span>
  <span class="spacer"></span>
  <span id="count" style="color:var(--mute)"></span>
  <div class="progress"><div id="bar" style="width:0"></div></div>
  <button class="btn" id="galbtn" style="display:none" title="all maps as thumbnails (g)">Gallery</button>
  <button class="btn" id="done" title="stop the server">Done</button>
</header>
<main id="app">
  <section class="card">
    <div class="maphead"><b id="title">map</b><span id="cur" class="pill"></span></div>
    <div class="plotwrap"><canvas id="plot"></canvas></div>
    <div class="nav">
      <button class="btn" id="prev">&larr;</button>
      <input type="range" id="slider" min="1" value="1">
      <button class="btn" id="next">&rarr;</button>
      <button class="btn" id="first" title="first unlabelled map">first unlabelled</button>
    </div>
    <div id="strip" title="each cell is one map, coloured by its label; click to jump"></div>
  </section>
  <aside class="card">
    <div class="sec">View <span style="text-transform:none;letter-spacing:0">(v)</span></div>
    <div class="seg" id="views"></div>
    <div class="sec">Label</div>
    <div id="labels"></div>
    <div class="sec">Note (optional)</div>
    <input id="note" placeholder="e.g. stripe + blob" autocomplete="off">
    <label class="chk"><input type="checkbox" id="auto"> go to the next map after labelling</label>
    <div id="status"></div>
    <div class="keys"><kbd>1</kbd>&ndash;<kbd>7</kbd> label &middot; <kbd>&larr;</kbd><kbd>&rarr;</kbd> move &middot; <kbd>v</kbd> view</div>
  </aside>
</main>
<section id="gallery" class="card"></section>
<div id="closed">Server stopped &mdash; you can close this tab.</div>
<script>
const $ = id => document.getElementById(id);
let M = null, pos = 0, view = 1, labs = [], notes = [], tok = 0;
const cache = new Map();

async function api(path, body) {
  const r = await fetch(path, body ? {method: "POST", headers: {"Content-Type": "application/json"},
                                      body: JSON.stringify(body)} : {});
  const j = await r.json(); if (j.error) throw new Error(j.error); return j;
}
function getMap(p) {
  if (!cache.has(p)) cache.set(p, api("/api/map?pos=" + p).catch(e => { cache.delete(p); throw e; }));
  return cache.get(p);
}
function status(t, err) { const s = $("status"); s.textContent = t; s.className = err ? "err" : ""; }

async function init() {
  M = await api("/api/meta");
  labs = M.labels; notes = M.notes;
  $("pass").textContent = M.mode;
  $("slider").max = M.n;
  try { $("auto").checked = localStorage.getItem("mapAuto") === "1"; } catch (e) {}
  $("auto").onchange = () => { try { localStorage.setItem("mapAuto", $("auto").checked ? "1" : "0"); } catch (e) {} };
  M.views.forEach((v, k) => { const b = document.createElement("button"); b.textContent = v;
    b.onclick = () => setView(k); $("views").appendChild(b); });
  M.label_names.forEach((l, k) => { const b = document.createElement("button"); b.className = "lab";
    b.style.setProperty("--c", M.label_clr[k]);
    b.innerHTML = `<span class="dot"></span><span>${l}</span><kbd>${k + 1}</kbd>`;
    b.onclick = () => setLabel(k); $("labels").appendChild(b); });
  for (let p = 0; p < M.n; p++) { const d = document.createElement("div"); d.onclick = () => go(p);
    $("strip").appendChild(d); }
  $("prev").onclick = () => go(pos - 1); $("next").onclick = () => go(pos + 1);
  $("first").onclick = () => { const p = labs.findIndex(l => !l);
    if (p < 0) status("all maps labelled"); else go(p); };
  $("slider").oninput = e => go(+e.target.value - 1);
  $("note").addEventListener("keydown", e => { if (e.key === "Enter") e.target.blur(); });
  $("note").addEventListener("change", () => { if (labs[pos]) setLabel(M.label_names.indexOf(labs[pos]), true); });
  $("done").onclick = async () => { try { await api("/api/quit", {}); } catch (e) {}
    $("app").style.display = "none"; document.querySelector("header").style.display = "none";
    $("closed").style.display = "block"; };
  window.addEventListener("resize", draw);
  document.addEventListener("keydown", onKey);
  status("labels save to " + M.csv);
  if (M.mode === "relabel") {
    $("galbtn").style.display = ""; $("galbtn").onclick = () => showGallery(!galOn);
    $("first").style.display = "none";
    go(M.start); showGallery(true);
  } else go(M.start);
}

// ── gallery (relabel mode): every map as a thumbnail; click one to (re)label it ──
let galOn = false;
function showGallery(on) {
  galOn = on;
  $("gallery").style.display = on ? "block" : "none";
  $("app").style.display = on ? "none" : "";
  $("galbtn").classList.toggle("on", on);
  $("galbtn").textContent = on ? "Back to map" : "Gallery";
  if (on) buildGallery(); else draw();
}
function buildGallery() {
  const G = $("gallery"); G.innerHTML = "";
  const counts = {}; M.groups.forEach(g => counts[g] = (counts[g] || 0) + 1);
  let grp = null, grid = null;
  M.groups.forEach((g, p) => {
    if (g !== grp) {
      grp = g; const h = document.createElement("h3");
      h.innerHTML = `${g || "maps"} <span>(${counts[g]})</span>`; G.appendChild(h);
      grid = document.createElement("div"); grid.className = "grid"; G.appendChild(grid);
    }
    const t = document.createElement("div"); t.className = "thumb";
    const cv = document.createElement("canvas"); t.appendChild(cv);
    const nm = document.createElement("div"); nm.className = "nm"; nm.textContent = M.names[p];
    nm.title = M.names[p]; t.appendChild(nm);
    const i = M.label_names.indexOf(labs[p]), ch = document.createElement("span");
    ch.className = "chip" + (i < 0 ? " none" : ""); ch.textContent = i < 0 ? "not labelled" : labs[p];
    if (i >= 0) ch.style.setProperty("--c", M.label_clr[i]);
    t.appendChild(ch);
    t.onclick = () => { showGallery(false); go(p); };
    grid.appendChild(t);
    getMap(p).then(m => drawThumb(cv, m.views[view])).catch(() => {});
  });
  window.scrollTo(0, 0);
}
function heat(V) {
  // the map as a small image (row 0 of the data is the bottom row)
  const D = V.data, nr = D.length, nc = D[0].length, lut = M.lut;
  const off = document.createElement("canvas"); off.width = nc; off.height = nr;
  const og = off.getContext("2d"), im = og.createImageData(nc, nr);
  for (let r = 0; r < nr; r++) for (let c = 0; c < nc; c++) {
    const v = D[r][c], i = 4 * ((nr - 1 - r) * nc + c);
    let rgb = [215, 218, 224];
    if (v !== null) { let u = (v - V.vmin) / (V.vmax - V.vmin); u = Math.max(0, Math.min(1, u));
      rgb = lut[Math.round(u * (lut.length - 1))]; }
    im.data[i] = rgb[0]; im.data[i + 1] = rgb[1]; im.data[i + 2] = rgb[2]; im.data[i + 3] = 255;
  }
  og.putImageData(im, 0, 0);
  return off;
}
function drawThumb(cv, V) {
  const dpr = window.devicePixelRatio || 1, W = cv.clientWidth, H = cv.clientHeight;
  cv.width = W * dpr; cv.height = H * dpr;
  const g = cv.getContext("2d"); g.setTransform(dpr, 0, 0, dpr, 0, 0);
  g.imageSmoothingEnabled = false; g.drawImage(heat(V), 0, 0, W, H);
  const x0 = M.wc[0], x1 = M.wc[M.wc.length - 1], y0 = M.kc[0], y1 = M.kc[M.kc.length - 1];
  const X = x => (x - x0) / (x1 - x0) * W, Y = y => H - (y - y0) / (y1 - y0) * H;
  const lo = Math.max(x0, y0), hi = Math.min(x1, y1);
  g.strokeStyle = "rgba(255,255,255,.8)"; g.setLineDash([4, 3]); g.lineWidth = 1;
  g.beginPath(); g.moveTo(X(lo), Y(lo)); g.lineTo(X(hi), Y(hi)); g.stroke();
}

async function go(p) {
  p = Math.max(0, Math.min(M.n - 1, p)); pos = p; const t = ++tok;
  sidebar();
  try { const m = await getMap(p); if (t !== tok) return; cur = m; draw(); }
  catch (e) { status("could not load map: " + e.message, true); }
  getMap(Math.min(M.n - 1, p + 1)); getMap(Math.max(0, p - 1));
}
let cur = null;
function setView(k) { view = (k + M.views.length) % M.views.length; sidebar(); draw(); }

async function setLabel(k, quiet) {
  const p = pos, lb = M.label_names[k];
  labs[p] = lb; notes[p] = $("note").value; sidebar();
  try {
    const r = await api("/api/label", {pos: p, label: lb, view: M.views[view], note: $("note").value});
    labs = r.labels; notes = r.notes;
    if (r.ok) status((quiet ? "note saved" : `map ${p + 1}: ${lb}`) + " — saved");
    else status("COULD NOT SAVE (file locked?). The label is kept and the file is rewritten with the next one.", true);
  } catch (e) { status("save failed: " + e.message, true); }
  if (!quiet && $("auto").checked && p === pos && p < M.n - 1) go(p + 1); else sidebar();
}

function onKey(e) {
  if (e.target.tagName === "INPUT" && e.target.type !== "range" && e.target.type !== "checkbox") return;
  const n = M.label_names.length;
  if (M.mode === "relabel" && (e.key === "g" || e.key === "G" || (e.key === "Escape" && !galOn))) {
    showGallery(!galOn); return; }
  if (galOn) { if (e.key === "Escape") showGallery(false); return; }
  if (e.key >= "1" && e.key <= String(n)) { setLabel(+e.key - 1); e.preventDefault(); }
  else if (e.key === "ArrowLeft") { go(pos - 1); e.preventDefault(); }
  else if (e.key === "ArrowRight") { go(pos + 1); e.preventDefault(); }
  else if (e.key === "v" || e.key === "V") setView(view + 1);
}

function sidebar() {
  const done = labs.filter(Boolean).length;
  $("count").textContent = `${done} / ${M.n} labelled`;
  $("bar").style.width = (100 * done / M.n) + "%";
  $("title").textContent = `map ${pos + 1} / ${M.n}` + (M.names ? "  —  " + M.names[pos] : "");
  $("cur").textContent = labs[pos] || "unlabelled";
  $("slider").value = pos + 1;
  [...$("views").children].forEach((b, k) => b.classList.toggle("on", k === view));
  [...$("labels").children].forEach((b, k) => b.classList.toggle("on", labs[pos] === M.label_names[k]));
  if (document.activeElement !== $("note")) $("note").value = notes[pos] || "";
  [...$("strip").children].forEach((d, k) => {
    const i = M.label_names.indexOf(labs[k]);
    d.style.background = i >= 0 ? M.label_clr[i] : "";
    d.classList.toggle("cur", k === pos);
  });
}

function ticks(a, b, n) {
  const s = (b - a) / n, p = Math.pow(10, Math.floor(Math.log10(s))), m = s / p;
  const st = (m < 1.5 ? 1 : m < 3 ? 2 : m < 7 ? 5 : 10) * p, out = [];
  for (let t = Math.ceil(a / st - 1e-9) * st; t <= b + 1e-9 * st; t += st) out.push(+t.toFixed(12));
  return out;
}
const fmt = t => Math.abs(t) < 1e-12 ? "0" : String(+t.toPrecision(3));

function draw() {
  if (!cur) return;
  const cv = $("plot"), dpr = window.devicePixelRatio || 1, W = cv.clientWidth, H = cv.clientHeight;
  cv.width = W * dpr; cv.height = H * dpr;
  const g = cv.getContext("2d"); g.setTransform(dpr, 0, 0, dpr, 0, 0); g.clearRect(0, 0, W, H);
  const V = cur.views[view], lut = M.lut;
  const L = 74, R = 96, T = 10, B = 50, pw = W - L - R, ph = H - T - B;
  g.imageSmoothingEnabled = false; g.drawImage(heat(V), L, T, pw, ph);
  const x0 = M.wc[0], x1 = M.wc[M.wc.length - 1], y0 = M.kc[0], y1 = M.kc[M.kc.length - 1];
  const X = x => L + (x - x0) / (x1 - x0) * pw, Y = y => T + ph - (y - y0) / (y1 - y0) * ph;
  // zero-lag diagonal and zero lines
  g.save(); g.beginPath(); g.rect(L, T, pw, ph); g.clip();
  g.strokeStyle = "rgba(255,255,255,.8)"; g.lineWidth = 1.2; g.setLineDash([7, 5]);
  const lo = Math.max(x0, y0), hi = Math.min(x1, y1);
  g.beginPath(); g.moveTo(X(lo), Y(lo)); g.lineTo(X(hi), Y(hi)); g.stroke();
  g.setLineDash([2, 4]); g.lineWidth = 1;
  g.beginPath(); g.moveTo(X(0), T); g.lineTo(X(0), T + ph); g.moveTo(L, Y(0)); g.lineTo(L + pw, Y(0)); g.stroke();
  g.restore();
  // axes
  g.strokeStyle = "#c4c9d2"; g.setLineDash([]); g.lineWidth = 1; g.strokeRect(L + .5, T + .5, pw - 1, ph - 1);
  g.fillStyle = "#4b5563"; g.font = "12px system-ui, Segoe UI, sans-serif";
  g.textAlign = "center"; g.textBaseline = "top";
  ticks(x0, x1, 6).forEach(t => { g.fillRect(X(t), T + ph, 1, 4); g.fillText(fmt(t), X(t), T + ph + 7); });
  g.textAlign = "right"; g.textBaseline = "middle";
  ticks(y0, y1, 6).forEach(t => { g.fillRect(L - 4, Y(t), 4, 1); g.fillText(fmt(t), L - 8, Y(t)); });
  g.fillStyle = "#1d2330"; g.font = "13px system-ui, Segoe UI, sans-serif";
  g.textAlign = "center"; g.textBaseline = "bottom";
  g.fillText("neural window centre re reach endpoint (s)", L + pw / 2, H - 4);
  g.save(); g.translate(16, T + ph / 2); g.rotate(-Math.PI / 2); g.textBaseline = "top";
  g.fillText("kinematic window centre (s)", 0, 0); g.restore();
  // colour bar
  const cx = L + pw + 18, cw = 14;
  for (let k = 0; k < ph; k++) { const rgb = lut[Math.round((1 - k / ph) * (lut.length - 1))];
    g.fillStyle = `rgb(${rgb})`; g.fillRect(cx, T + k, cw, 1.5); }
  g.strokeStyle = "#c4c9d2"; g.strokeRect(cx + .5, T + .5, cw - 1, ph - 1);
  g.fillStyle = "#4b5563"; g.font = "11px system-ui, Segoe UI, sans-serif"; g.textAlign = "left"; g.textBaseline = "middle";
  ticks(V.vmin, V.vmax, 5).forEach(t => { const yy = T + ph - (t - V.vmin) / (V.vmax - V.vmin) * ph;
    g.fillRect(cx + cw, yy, 3, 1); g.fillText(fmt(t), cx + cw + 6, yy); });
  g.save(); g.translate(W - 6, T + ph / 2); g.rotate(-Math.PI / 2); g.textAlign = "center"; g.textBaseline = "bottom";
  g.fillStyle = "#1d2330"; g.font = "12px system-ui, Segoe UI, sans-serif";
  g.fillText(M.views[view] + (view === 0 ? " (bits)" : ""), 0, 0); g.restore();
}
init().catch(e => status("could not start: " + e.message, true));
</script>
</body></html>
"""


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--set", default=None, help="label set .pkl (default: under the data root)")
    ap.add_argument("--repeat", action="store_true", help="the repeat pass")
    ap.add_argument("--port", type=int, default=0, help="server port (default: any free one)")
    ap.add_argument("--no-browser", action="store_true", help="do not open the browser")
    ap.add_argument("--relabel", nargs="+", default=None, metavar="MAP",
                    help="maps to (re)label: 'Ferrero R1 224', 'CB 5_FerreroRocher R1 224 (shuffle)', "
                         "or a .txt file with one map per line (28 writes relabel_<tag>.txt)")
    ap.add_argument("--selftest", default=None, metavar="DIR", help="no browser: test into DIR")
    a = ap.parse_args()
    sp = a.set or str(config.var_path(f"label_set_{DEFAULT_TAG}.pkl"))
    if not os.path.isfile(sp):
        sys.exit(f"label set not found: {sp}\nrun 26_label_set.py first")
    if a.selftest: selftest(sp, a.selftest)
    else: run_gui(sp, a.repeat, a.port, not a.no_browser, relabel=a.relabel)
