# STEP 07 -- interactive review of the cell types, in the browser.
#
#   conda activate celltype
#   python 07_review_gui.py            (opens http://127.0.0.1:<port>/ ; stop with Ctrl+C)
#
# LEFT   every neuron of the chosen region as a dot. Choose what the axes and the
#        colour show (any feature, the embedding, C4's confidence ...), and filter
#        by type, mouse, session, flag, confidence. Click a dot to open the neuron.
# RIGHT  the neuron: waveform (with the mean of its type in grey), the waveform
#        on the neighbouring channels, the autocorrelogram, C4's class
#        probabilities and CS-triggered histograms (cerebellum), the criteria
#        values (striatum), and every label it has: s04, C4, automatic (+ the
#        rule and flags), manual.
# RELABEL  click a type button (or keys 1-6), optionally with a note. "clear"
#        removes your manual label. <- / -> step through the filtered neurons.
#
# What it writes: ONLY overrides.csv (one line appended per click; the last line
# per neuron counts, an empty label = cleared). It never changes the automatic
# tables or any .mat. Run 06_build_cell_type_table.py afterwards to put the manual
# labels into cell_types.csv.
import json, math, os, sys, threading, time, webbrowser
from http.server import ThreadingHTTPServer, BaseHTTPRequestHandler
from urllib.parse import urlparse, parse_qs
import numpy as np
import ct_config as C, ct_io as io

PORT = 0                # 0 = any free port
OPEN_BROWSER = True
SCATTER_FEATURES = ["emb_x", "emb_y", "pc1", "pc2", "fr", "wf_duration_us", "wf_halfwidth_us",
                    "wf_pt_ratio", "wf_repol_slope", "cv2", "isi_median_ms", "prop_long_isi",
                    "post_spike_supp_ms", "refractory_ms", "burst_index", "c4_ratio", "c4_votes",
                    "c4_prob", "clf_prob", "spread", "pre_peak", "peak_ch", "n_spikes"]
LOG_FEATURES = ["fr", "isi_median_ms", "post_spike_supp_ms", "burst_index", "c4_ratio", "n_spikes"]
OVR_COLS = ["uid", "label", "note", "time", "prev"]


def _clean(o):
    """NaN / numpy -> plain JSON."""
    if isinstance(o, dict): return {k: _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)): return [_clean(v) for v in o]
    if isinstance(o, np.ndarray): return _clean(o.tolist())
    if isinstance(o, (np.floating, float)):
        o = float(o); return o if math.isfinite(o) else None
    if isinstance(o, np.integer): return int(o)
    if isinstance(o, np.bool_): return bool(o)
    return o


class State:
    """All the data of the GUI. Takes the tables as arguments so it can be tested
    without touching the disk; State.from_disk() loads the real ones."""

    def __init__(self, feat, c4_rows, cb_auto, cp_auto, pairs, overrides_path):
        self.F, self.ovr_path, self.lock = feat, str(overrides_path), threading.Lock()
        self.idx = {r["uid"]: i for i, r in enumerate(feat["rows"])}
        self.c4 = {r["uid"]: r for r in c4_rows}
        self.auto = {r["uid"]: r for r in list(cb_auto) + list(cp_auto)}
        self.pairs = pairs or dict(pairs=[], lags=None)
        self.by_unit = {}
        for p in self.pairs["pairs"]:
            self.by_unit.setdefault(p["cs"], []).append(p); self.by_unit.setdefault(p["tgt"], []).append(p)
        self.manual = {}
        for r in io.read_csv(self.ovr_path): self.manual[r["uid"]] = r
        self.pc = self._pcs()

    @classmethod
    def from_disk(cls):
        return cls(io.load_features(), io.read_csv(C.CB_C4_CSV), io.read_csv(C.CB_AUTO_CSV),
                   io.read_csv(C.CP_AUTO_CSV),
                   io.load_pkl(C.CB_PAIRS_PKL) if C.CB_PAIRS_PKL.is_file() else None, C.OVERRIDES_CSV)

    def _pcs(self):
        """A quick 2-D view per region (PCA of waveform + log ACG) for the scatter."""
        from sklearn.decomposition import PCA
        out = {}
        A = self.F["arrays"]
        for reg in ("CB", "CP"):
            ii = [i for i, r in enumerate(self.F["rows"]) if r["region"] == reg]
            if len(ii) < 5: continue
            W = np.nan_to_num(A["wf_norm"][ii].astype(float))
            G = A["acg"][ii, :200].astype(float)
            G = np.log1p(np.nan_to_num(G / np.maximum(np.nanmean(A["acg"][ii, 100:400], axis=1, keepdims=True), 1e-9)))
            X = np.c_[W / max(W.std(), 1e-9), G / max(G.std(), 1e-9)]
            P = PCA(2, random_state=0).fit_transform(X)
            for i, p in zip(ii, P): out[self.F["rows"][i]["uid"]] = p
        return out

    def _labels(self, u):
        a, c, m = self.auto.get(u, {}), self.c4.get(u, {}), self.manual.get(u, {})
        auto = a.get("auto_cell_type", "") or C.UNLAB
        man = m.get("label", "")
        return dict(old=a.get("old_cell_type", "") or c.get("old_cell_type", ""),
                    c4=c.get("c4_pred", ""), c4_status=c.get("c4_status", ""), auto=auto,
                    source=a.get("auto_source", ""), conf=a.get("auto_conf", ""),
                    flags=a.get("flags", ""), partner=a.get("partner", ""), manual=man,
                    note=m.get("note", "") if man else "", final=man or auto,
                    crit=a.get("crit_cell_type", ""), cluster=a.get("cluster", ""),
                    cluster_label=a.get("cluster_label", ""), clf=a.get("clf_pred", ""),
                    clf_prob=io.fnum(a.get("clf_prob")), suggested=a.get("suggested_cell_type", ""))

    def meta(self):
        neurons = []
        for r in self.F["rows"]:
            u = r["uid"]; a, c = self.auto.get(u, {}), self.c4.get(u, {})
            f = {k: r.get(k) for k in SCATTER_FEATURES if k in r}
            for k in ("c4_ratio", "c4_votes", "c4_prob"): f[k] = io.fnum(c.get(k))
            for k in ("emb_x", "emb_y", "clf_prob", "spread", "pre_peak"): f[k] = io.fnum(a.get(k))
            p = self.pc.get(u, (np.nan, np.nan)); f["pc1"], f["pc2"] = float(p[0]), float(p[1])
            neurons.append(dict(uid=u, region=r["region"], mouse=r["mouse"], sess=r["sess"],
                                phyID=r["phyID"], f=f, **self._labels(u)))
        # mean waveform / ACG per region and FINAL type, drawn behind each neuron
        A = self.F["arrays"]; means = {}
        for reg in ("CB", "CP"):
            for t in C.TYPES[reg]:
                ii = [self.idx[n["uid"]] for n in neurons if n["region"] == reg and n["final"] == t]
                if len(ii) < 3: continue
                g = A["acg"][ii, :100].astype(float)
                g = g / np.maximum(np.nanmean(A["acg"][ii, 100:400], axis=1, keepdims=True), 1e-9)
                means[f"{reg}|{t}"] = dict(wf=np.nanmean(A["wf_norm"][ii], axis=0), acg=np.nanmean(g, axis=0))
        return _clean(dict(neurons=neurons, types=C.TYPES, colors=C.COLORS, features=SCATTER_FEATURES,
                           log=LOG_FEATURES, means=means, wf_t=self.F["meta"]["wf_norm_t_ms"],
                           unlab=C.UNLAB, csv=os.path.basename(self.ovr_path)))

    def neuron(self, u):
        i = self.idx[u]; r = self.F["rows"][i]; A = self.F["arrays"]
        acg = A["acg"][i].astype(float)
        base = float(np.nanmean(acg[100:400])) if np.isfinite(acg).any() else np.nan
        c = self.c4.get(u, {})
        probs = {k[2:]: io.fnum(v) for k, v in c.items() if k.startswith("p_")}
        prs = []
        for p in sorted(self.by_unit.get(u, []), key=lambda p: (not p["passed"], p["ratio"]))[:4]:
            other = p["tgt"] if p["cs"] == u else p["cs"]
            prs.append(dict(cs=p["cs"], tgt=p["tgt"], role="CS" if p["cs"] == u else "paused unit",
                            other=other, other_label=self._labels(other)["final"], rate=p["rate"],
                            **{k: p.get(k) for k in ("ratio", "asym", "refr", "share", "z", "base_hz",
                                                     "n_cs", "ch_dist", "passed")},
                            used=bool(p.get("used", False))))
        return _clean(dict(uid=u, feat={k: v for k, v in r.items() if not isinstance(v, str) or k == "kind"},
                           wf=A["wf_norm"][i], wf_multi=A["wf_multi"][i], acg=acg, acg_base=base,
                           probs=probs, pairs=prs, lags=self.pairs.get("lags"), **self._labels(u)))

    def label(self, u, label, note):
        if u not in self.idx: raise ValueError(f"unknown neuron {u!r}")
        reg = self.F["rows"][self.idx[u]]["region"]
        if label and label not in C.TYPES[reg]: raise ValueError(f"{label!r} is not a {reg} type")
        with self.lock:
            prev = self.manual.get(u, {}).get("label", "")
            row = dict(uid=u, label=label, note=(note or "").strip(),
                       time=time.strftime("%Y-%m-%d %H:%M:%S"), prev=prev)
            import csv
            new = not os.path.isfile(self.ovr_path)
            C.out_path(self.ovr_path)
            with open(self.ovr_path, "a", newline="", encoding="utf-8") as fh:
                w = csv.DictWriter(fh, OVR_COLS)
                if new: w.writeheader()
                w.writerow(row)
            self.manual[u] = row
        return dict(uid=u, n_manual=sum(bool(r.get("label")) for r in self.manual.values()),
                    **self._labels(u))


def _plotly_js():
    try:
        from plotly.offline import get_plotlyjs
        return get_plotlyjs()
    except ImportError:
        return None


def make_server(S, port=0):
    js = _plotly_js()
    page = PAGE.replace("/*PLOTLY*/", '<script src="/plotly.js"></script>' if js else
                        '<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>')

    class H(BaseHTTPRequestHandler):
        def log_message(self, *a): pass

        def _send(self, body, ctype="application/json", code=200):
            b = body.encode("utf-8") if isinstance(body, str) else body
            self.send_response(code); self.send_header("Content-Type", ctype + "; charset=utf-8")
            self.send_header("Content-Length", str(len(b))); self.send_header("Cache-Control", "no-store")
            self.end_headers(); self.wfile.write(b)

        def do_GET(self):
            u = urlparse(self.path)
            try:
                if u.path == "/": self._send(page, "text/html")
                elif u.path == "/plotly.js" and js: self._send(js, "application/javascript")
                elif u.path == "/api/meta": self._send(json.dumps(S.meta()))
                elif u.path == "/api/neuron": self._send(json.dumps(S.neuron(parse_qs(u.query)["uid"][0])))
                else: self._send(json.dumps({"error": "not found"}), code=404)
            except Exception as e:
                self._send(json.dumps({"error": repr(e)}), code=500)

        def do_POST(self):
            try:
                body = json.loads(self.rfile.read(int(self.headers.get("Content-Length", 0))) or b"{}")
                if urlparse(self.path).path == "/api/label":
                    self._send(json.dumps(_clean(S.label(body["uid"], body.get("label", ""), body.get("note", "")))))
                else: self._send(json.dumps({"error": "not found"}), code=404)
            except Exception as e:
                self._send(json.dumps({"error": repr(e)}), code=500)

    return ThreadingHTTPServer(("127.0.0.1", port), H)


def run():
    S = State.from_disk()
    srv = make_server(S, PORT)
    url = f"http://127.0.0.1:{srv.server_address[1]}/"
    print(f"{len(S.idx)} neurons, {sum(bool(r.get('label')) for r in S.manual.values())} manual labels")
    print(f"manual labels -> {S.ovr_path}\nopen {url}   (stop: Ctrl+C)")
    if OPEN_BROWSER: threading.Timer(0.5, lambda: webbrowser.open(url)).start()
    try: srv.serve_forever()
    except KeyboardInterrupt: pass
    srv.server_close(); print("stopped -- run 06_build_cell_type_table.py to update cell_types.csv")


PAGE = r"""<!doctype html><html><head><meta charset="utf-8"><title>Cell types</title>
/*PLOTLY*/
<style>
body{margin:0;font:13px/1.35 system-ui,Segoe UI,sans-serif;color:#222;background:#fafafa}
#bar{display:flex;flex-wrap:wrap;gap:8px 14px;align-items:center;padding:7px 10px;background:#fff;border-bottom:1px solid #ddd}
#bar label{display:flex;gap:4px;align-items:center;white-space:nowrap}
select,input[type=text]{font:inherit;padding:2px 4px}
#main{display:grid;grid-template-columns:minmax(380px,44%) 1fr;height:calc(100vh - 78px)}
#left{border-right:1px solid #ddd;display:flex;flex-direction:column;min-width:0}
#scatter{flex:1;min-height:0}
#right{overflow:auto;padding:8px 12px}
#chips{display:flex;flex-wrap:wrap;gap:4px;padding:5px 10px;background:#fff;border-bottom:1px solid #eee}
.chip{padding:1px 8px;border-radius:10px;border:1.5px solid;cursor:pointer;user-select:none;font-size:12px}
.chip.off{opacity:.3}
.lab{display:inline-block;padding:1px 8px;border-radius:4px;color:#fff;font-weight:600;margin-right:4px}
.row{display:grid;grid-template-columns:repeat(3,1fr);gap:6px}
.plot{height:215px;background:#fff;border:1px solid #e5e5e5;border-radius:4px}
#btns button{font:inherit;font-weight:600;padding:5px 11px;margin:2px;border-radius:5px;border:2px solid;background:#fff;cursor:pointer}
#btns button.on{color:#fff}
table{border-collapse:collapse;font-size:12px}td{padding:1px 8px 1px 0;white-space:nowrap}
.dim{color:#777}.flag{color:#b00020;font-weight:600}h3{margin:6px 0 4px;font-size:15px}
#status{margin-left:auto;color:#555}
</style></head><body>
<div id="bar">
 <label>region <select id="region"><option>CB</option><option>CP</option></select></label>
 <label>x <select id="fx"></select></label><label>y <select id="fy"></select></label>
 <label>colour <select id="col"><option value="final">final label</option><option value="auto">automatic</option>
  <option value="old">s04</option><option value="c4">C4 best guess</option><option value="crit">criteria (CP)</option>
  <option value="cluster_label">cluster (CP)</option><option value="clf">classifier (CP)</option>
  <option value="suggested">suggested relabel (CP)</option><option value="conf">confidence</option>
  <option value="source">rule</option><option value="mouse">mouse</option></select></label>
 <label>mouse <select id="mouse"></select></label><label>session <select id="sess"></select></label>
 <label><input type="checkbox" id="fflag">flagged</label>
 <label><input type="checkbox" id="flow">low conf.</label>
 <label><input type="checkbox" id="fsug">suggested</label>
 <label><input type="checkbox" id="fchg">changed vs s04</label>
 <label><input type="checkbox" id="fman">manual</label>
 <label>find <input type="text" id="find" size="9" placeholder="phyID"></label>
 <span id="status"></span>
</div>
<div id="chips"></div>
<div id="main"><div id="left"><div id="scatter"></div></div>
<div id="right"><div id="info" class="dim">Click a neuron. Keys: &larr; &rarr; previous / next in the filtered set, 1-6 set label.</div>
 <div id="btns"></div>
 <div style="margin:4px 0"><input type="text" id="note" size="50" placeholder="note (saved with the next label click)"></div>
 <div class="row"><div id="p_wf" class="plot"></div><div id="p_multi" class="plot"></div><div id="p_acg" class="plot"></div></div>
 <div class="row" style="margin-top:6px"><div id="p_acgl" class="plot"></div><div id="p_a" class="plot"></div><div id="p_b" class="plot"></div></div>
 <div class="row" style="margin-top:6px"><div id="p_c" class="plot"></div><div id="p_d" class="plot"></div><div id="feat"></div></div>
</div></div>
<script>
let M=null, cur=null, shown=[], off=new Set();
const $=id=>document.getElementById(id), CFG={displayModeBar:false,responsive:true};
const L=(t,x,y,extra)=>Object.assign({margin:{l:42,r:8,t:26,b:32},title:{text:t,font:{size:12}},xaxis:{title:{text:x,font:{size:11}}},yaxis:{title:{text:y,font:{size:11}}},showlegend:false,font:{size:10}},extra||{});
// Colours: cell types keep their fixed colours (ct_config.COLORS). Anything else
// (mice, rules, confidence) gets a colour from OTHER_PAL in sorted order, so two
// values on the screen never share a colour, and none reuses a cell-type colour.
// CS-SS pair histograms are dark blue (SS blue + CS black).
const PAIR_COL='#0b2a6f';
const OTHER_PAL=['#8c564b','#17becf','#bcbd22','#ff7f0e','#7b4173','#8c6d31','#5254a3','#637939','#ad494a','#3182bd'];
let CM={};
function setOther(cats){CM={};let i=0;for(const c of cats){if(c==='-'||M.colors[c])continue;CM[c]=OTHER_PAL[i%OTHER_PAL.length];i++}}
const colOf=v=>M.colors[v]||CM[v]||'#bbb';
async function api(p,body){const r=await fetch(p,body?{method:'POST',body:JSON.stringify(body)}:{});const j=await r.json();if(j.error)throw new Error(j.error);return j}
function opt(sel,vals,keep){const v=keep?sel.value:null;sel.innerHTML=vals.map(x=>`<option>${x}</option>`).join('');if(v&&vals.includes(v))sel.value=v}
function region(){return $('region').value}
function setup(){
  const reg=region(), N=M.neurons.filter(n=>n.region===reg);
  const feats=M.features.filter(k=>N.some(n=>n.f[k]!=null));
  opt($('fx'),feats,true);opt($('fy'),feats,true);
  if(!$('fx').dataset.set||$('fx').dataset.set!==reg){$('fx').value=feats.includes('emb_x')?'emb_x':'pc1';$('fy').value=feats.includes('emb_y')?'emb_y':'pc2';$('fx').dataset.set=reg}
  opt($('mouse'),['all',...new Set(N.map(n=>n.mouse))]);opt($('sess'),['all',...[...new Set(N.map(n=>n.sess))].sort()]);
  off.clear();draw()}
function filt(){
  const reg=region(),mo=$('mouse').value,se=$('sess').value,q=$('find').value.trim(),k=$('col').value;
  return M.neurons.filter(n=>n.region===reg&&(mo==='all'||n.mouse===mo)&&(se==='all'||n.sess===se)
    &&(!$('fflag').checked||n.flags)&&(!$('flow').checked||n.conf==='low')&&(!$('fsug').checked||n.suggested)&&(!$('fman').checked||n.manual)
    &&(!$('fchg').checked||(n.old&&n.old!==n.final))&&(!q||String(n.phyID)===q)&&!off.has(n[k]||'-'))}
function draw(){
  const k=$('col').value,fx=$('fx').value,fy=$('fy').value,reg=region();
  const all=M.neurons.filter(n=>n.region===reg), cats=[...new Set(all.map(n=>n[k]||'-'))].sort();
  setOther(cats);
  $('chips').innerHTML=cats.map(c=>`<span class="chip ${off.has(c)?'off':''}" data-c="${c}" style="border-color:${colOf(c)};color:${colOf(c)}">${c} ${all.filter(n=>(n[k]||'-')===c).length}</span>`).join('');
  document.querySelectorAll('.chip').forEach(e=>e.onclick=()=>{const c=e.dataset.c;off.has(c)?off.delete(c):off.add(c);draw()});
  shown=filt();
  const tr=cats.filter(c=>!off.has(c)).map(c=>{const s=shown.filter(n=>(n[k]||'-')===c);return{type:'scattergl',mode:'markers',name:c,
    x:s.map(n=>n.f[fx]),y:s.map(n=>n.f[fy]),customdata:s.map(n=>n.uid),text:s.map(n=>`${n.mouse} ${n.sess} #${n.phyID}<br>final ${n.final} | auto ${n.auto} (${n.source||'-'}) | s04 ${n.old||'-'}${n.flags?'<br>'+n.flags:''}`),
    hovertemplate:'%{text}<extra></extra>',marker:{size:s.map(n=>n.manual?9:6),color:colOf(c),opacity:.75,symbol:s.map(n=>n.manual?'diamond':'circle'),line:{width:s.map(n=>n.flags?1:0),color:'#000'}}}});
  if(cur){const n=M.neurons.find(n=>n.uid===cur);if(n&&n.region===reg)tr.push({type:'scattergl',mode:'markers',x:[n.f[fx]],y:[n.f[fy]],hoverinfo:'skip',marker:{size:17,color:'rgba(0,0,0,0)',line:{width:2.5,color:'#000'}}})}
  const ax=f=>({title:{text:f},type:M.log.includes(f)?'log':'linear'});
  Plotly.react('scatter',tr,{margin:{l:55,r:10,t:10,b:45},xaxis:ax(fx),yaxis:ax(fy),showlegend:false,hovermode:'closest',uirevision:reg+fx+fy},{responsive:true,scrollZoom:true});
  $('scatter').removeAllListeners&&$('scatter').removeAllListeners('plotly_click');
  $('scatter').on('plotly_click',e=>{const u=e.points[0].customdata;if(u)show(u)});
  $('status').textContent=`${shown.length} / ${all.length} shown · ${M.neurons.filter(n=>n.manual).length} manual · diamonds = manual, black edge = flagged · ${M.csv}`}
const fmt=v=>v==null?'–':(typeof v==='number'?(Math.abs(v)>=100?v.toFixed(0):v.toPrecision(3)):v);
async function show(u){
  cur=u;const d=await api('/api/neuron?uid='+encodeURIComponent(u));const n=M.neurons.find(n=>n.uid===u),reg=n.region;
  const tag=(t,v)=>v?`<span class="dim">${t}</span> <span class="lab" style="background:${colOf(v)}">${v}</span>`:'';
  $('info').className='';$('info').innerHTML=`<h3>${u.split('|').join(' · ')}</h3>${tag('final',d.final)}${tag('manual',d.manual)}${tag('automatic',d.auto)}
   <span class="dim">${d.source?('rule: '+d.source+(d.conf?', '+d.conf+' confidence':'')):''}</span> ${tag('s04',d.old)}${tag('C4',d.c4)}<span class="dim">${d.c4_status||''}</span>
   ${tag('criteria',d.crit)}${d.cluster_label?`<span class="dim">cluster ${d.cluster}: ${d.cluster_label}</span>`:''}
   ${d.clf?`${tag('classifier',d.clf)}<span class="dim">p ${fmt(d.clf_prob)}</span>`:''}
   ${d.suggested?`<div class="flag">classifier suggests <b>${d.suggested}</b> (p ${fmt(d.clf_prob)}): press its button to accept, or leave it</div>`:''}
   ${d.flags?`<div class="flag">flags: ${d.flags.split(';').join(', ')}</div>`:''}${d.partner?`<div>partner: <a href="#" id="plink">${d.partner}</a></div>`:''}${d.note?`<div class="dim">note: ${d.note}</div>`:''}`;
  if(d.partner)$('plink').onclick=e=>{e.preventDefault();show(d.partner)};
  $('btns').innerHTML=M.types[reg].map((t,i)=>`<button data-t="${t}" class="${d.manual===t?'on':''}" style="border-color:${colOf(t)};${d.manual===t?'background:'+colOf(t):'color:'+colOf(t)}">${i+1} ${t}</button>`).join('')+`<button data-t="" style="border-color:#999;color:#555">clear manual</button>`;
  document.querySelectorAll('#btns button').forEach(b=>b.onclick=()=>setLabel(b.dataset.t));
  const c=colOf(d.final), mean=M.means[reg+'|'+d.final];
  Plotly.react('p_wf',[...(mean?[{x:M.wf_t,y:mean.wf,line:{color:'#bbb',width:4},name:'type mean'}]:[]),{x:M.wf_t,y:d.wf,line:{color:c,width:2}}],L('waveform (grey: mean of '+d.final+')','ms','norm.'),CFG);
  const nch=d.wf_multi.length,amp=Math.max(...d.wf_multi.flat().filter(v=>v!=null).map(Math.abs))||1;
  Plotly.react('p_multi',d.wf_multi.map((w,i)=>({y:w.map(v=>v==null?null:v/amp+ (i-(nch-1)/2)*0.6),line:{color:i===(nch-1)/2?c:'#888',width:i===(nch-1)/2?2:1}})),L('peak channel ± '+((nch-1)/2)+' channels','sample','',{yaxis:{showticklabels:false}}),CFG);
  const lag=d.acg.map((_,i)=>i+0.5);
  Plotly.react('p_acg',[{type:'bar',x:lag.slice(0,100),y:d.acg.slice(0,100),marker:{color:c}},...(mean&&d.acg_base?[{x:lag.slice(0,100),y:mean.acg.map(v=>v*d.acg_base),line:{color:'#999',width:1.5}}]:[])],L('ACG 0–100 ms (grey: type mean, scaled)','lag (ms)','sp/s',{bargap:0}),CFG);
  const pss=d.feat.post_spike_supp_ms;
  Plotly.react('p_acgl',[{type:'bar',x:lag,y:d.acg,marker:{color:c}}],L('ACG 0–500 ms','lag (ms)','sp/s',{bargap:0,shapes:[...(d.acg_base?[{type:'line',x0:0,x1:500,y0:d.acg_base,y1:d.acg_base,line:{dash:'dot',width:1}}]:[]),...(pss!=null?[{type:'line',x0:pss,x1:pss,y0:0,y1:1,yref:'paper',line:{color:'red',width:1}}]:[])]}),CFG);
  const slots=['p_a','p_b','p_c','p_d'];let s=0;
  if(Object.keys(d.probs).length&&Object.values(d.probs).some(v=>v!=null)){const ks=Object.keys(d.probs);
    Plotly.react(slots[s++],[{type:'bar',x:ks,y:ks.map(k=>d.probs[k]),marker:{color:ks.map(colOf)}}],L(`C4 probabilities (ratio ${fmt(n.f.c4_ratio)}, votes ${fmt(n.f.c4_votes)})`,'','p',{yaxis:{range:[0,1]}}),CFG)}
  if(reg==='CP'){const f=d.feat;
    Plotly.react(slots[s++],[{type:'bar',orientation:'h',y:['duration / 400 us','post-spike supp. / 40 ms','long-ISI prop. / 0.1'],x:[f.wf_duration_us/400,f.post_spike_supp_ms/40,f.prop_long_isi/0.1],marker:{color:c}}],L('criteria values ÷ threshold (line = threshold)','','',{margin:{l:150,r:8,t:26,b:30},xaxis:{type:'log'},shapes:[{type:'line',x0:1,x1:1,y0:0,y1:1,yref:'paper',line:{dash:'dash'}}]}),CFG)}
  for(const p of d.pairs){if(s>=slots.length)break;
    Plotly.react(slots[s++],[{type:'bar',x:d.lags.map(v=>v*1000),y:p.rate,marker:{color:p.passed?PAIR_COL:'#888'}}],
      L(`${p.role==='CS'?'this CS → #'+p.other.split('|')[3]:'CS #'+p.other.split('|')[3]+' → this unit'} (${p.other_label}) r=${fmt(p.ratio)} asym=${fmt(p.asym)} refr=${fmt(p.refr)} ${p.passed?(p.used?'PASS, used':'pass'):'fail'}`,'ms from CS','sp/s',{bargap:0,shapes:[{type:'rect',x0:3,x1:8,y0:0,y1:1,yref:'paper',fillcolor:'rgba(0,0,0,.08)',line:{width:0}}]}),CFG)}
  for(;s<slots.length;s++)Plotly.purge(slots[s]);
  const keys=['fr','n_spikes','cv2','isi_median_ms','refractory_ms','burst_index','post_spike_supp_ms','prop_long_isi','wf_duration_us','wf_halfwidth_us','wf_pt_ratio','wf_polarity','peak_ch','depth'];
  $('feat').innerHTML='<table>'+keys.filter(k=>k in d.feat).map(k=>`<tr><td class="dim">${k}</td><td>${fmt(d.feat[k])}</td></tr>`).join('')+'</table>';
  draw()}
async function setLabel(t){
  if(!cur)return;const r=await api('/api/label',{uid:cur,label:t,note:$('note').value});
  Object.assign(M.neurons.find(n=>n.uid===cur),r);$('note').value='';show(cur)}
function step(k){if(!shown.length)return;let i=shown.findIndex(n=>n.uid===cur);i=(i+k+shown.length)%shown.length;show(shown[i<0?0:i].uid)}
document.addEventListener('keydown',e=>{if(e.target.tagName==='INPUT')return;
  if(e.key==='ArrowRight')step(1);else if(e.key==='ArrowLeft')step(-1);
  else if(cur&&/^[1-9]$/.test(e.key)){const t=M.types[M.neurons.find(n=>n.uid===cur).region][+e.key-1];if(t)setLabel(t)}});
['fx','fy','col','mouse','sess','fflag','flow','fsug','fchg','fman'].forEach(id=>$(id).onchange=()=>{if(id==='col')off.clear();draw()});
$('find').oninput=draw;$('region').onchange=()=>{cur=null;setup()};
api('/api/meta').then(m=>{M=m;setup()}).catch(e=>{$('info').textContent='could not load: '+e.message});
</script></body></html>"""


if __name__ == "__main__":
    run()
