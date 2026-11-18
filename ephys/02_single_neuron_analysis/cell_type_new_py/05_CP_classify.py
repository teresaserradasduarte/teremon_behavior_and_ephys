# %% ── 0. parameters ─────────────────────────────────────────────────────────
# STEP 05 -- striatum (CP) cell types, in two independent ways that are then
# compared.
#
# A. CRITERIA (Peters et al. 2021; bombcell +bc/+clsfy/classifyCells.m, values
#    from +bc/+ep/ephysPropValues.m):
#      waveform duration (trough -> following peak) > 400 us   "wide"
#          post-spike suppression <  40 ms  -> MSN
#          post-spike suppression >= 40 ms  -> TAN
#      waveform duration <= 400 us                              "narrow"
#          proportion of time in ISIs > 2 s <= 0.1 -> FSI
#          proportion of time in ISIs > 2 s >  0.1 -> UIN
#    post-spike suppression = first lag at which the ACG is back at its
#    100-400 ms level. The features come from step 01.
#    STRICTER FSI (not in bombcell; its criteria gave far too many FSI here):
#    a narrow unit with few long ISIs is FSI only if it ALSO has
#      duration <= FSI_DUR_US        clearly narrow, not just under the 400 us split
#      firing rate >= FSI_FR_MIN     FSIs fire fast and tonically
#      a negative main deflection    (FSI_NEG_ONLY) positive spikes are mostly axonal
#    A narrow unit that fails any of these becomes UIN, with a flag saying which
#    one (fsi_not_narrow, fsi_fr_low, fsi_positive). Set a parameter to None /
#    False to switch that condition off; all off = bombcell's rule.
# B. UNSUPERVISED: every neuron is described by its normalised waveform (PCA),
#    its ACG shape (PCA of the log ACG) and a few scalars; a Gaussian mixture
#    finds clusters WITHOUT using the criteria, and a 2-D embedding (UMAP if
#    installed, else t-SNE) shows them. Each cluster is then named after the
#    criteria label most of its neurons carry.
# C. CLASSIFIER CHECK (to find mislabelled neurons): a random forest learns the
#    criteria labels from the features of B, but WITHOUT the three numbers the
#    criteria are made of (duration, post-spike suppression, prop. long ISI) --
#    otherwise it would only re-learn the thresholds. It still sees the whole
#    waveform and ACG shape, so it judges "does this neuron look like the others
#    that carry its label". Every neuron is predicted by a model that never saw
#    it (N_FOLDS-fold cross-validation); unlabelled neurons by a model trained
#    on all labelled ones. When it confidently (p >= RELABEL_P) says another
#    type, the neuron gets the flag clf_suggests_<type> and suggested_cell_type.
#    It NEVER changes the label: accept or reject the suggestion in the GUI.
# The automatic label is the CRITERIA label.
#   auto_conf  high = the classifier agrees and no threshold is close
#                     (+ the cluster agrees, if CONF_USES_CLUSTER)
#              low  = otherwise
#   flags      clf_suggests_<type>, cluster_disagrees, near_threshold, few_spikes,
#              msn_fr_high, fsi_not_narrow, fsi_fr_low, fsi_positive
import collections, warnings
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
warnings.filterwarnings("ignore", message="KMeans is known to have a memory leak")   # sklearn/MKL noise
from sklearn.mixture import GaussianMixture
import ct_config as C, ct_io as io

DUR_US, PSS_MS, PROP_ISI = 400.0, 40.0, 0.1      # bombcell's thresholds
FSI_DUR_US   = 350.0    # stricter FSI: max duration (us);       None = off
FSI_FR_MIN   = 5.0      # stricter FSI: min firing rate (Hz);    None = off
FSI_NEG_ONLY = True     # stricter FSI: negative-going spikes only
MIN_SPIKES = 300        # fewer -> unlabelled (ACG too noisy for the criteria)
NEAR_FRAC  = 0.10       # "near a threshold" = within 10% of it
MSN_FR_MAX = 10.0       # Hz; an MSN above this is flagged
N_WF_PC, N_ACG_PC = 5, 3
N_CLUST    = 6          # mixture components; see the BIC panel before changing
K_RANGE    = range(2, 11)
EMBED      = "auto"     # "auto" (UMAP if installed, else t-SNE) | "umap" | "tsne"
SPREAD_FRAC = 0.3       # spatial spread = channels (of peak +-4) above this fraction of the peak amplitude
RELABEL_P  = 0.8        # classifier probability needed to suggest another type
N_FOLDS    = 5          # cross-validation folds
CONF_USES_CLUSTER = False   # True: high confidence also needs the cluster to agree
SEED       = 0

SHOW_FIGS = True
SAVE_FIGS = True
SAVE_VARS = True        # writes cp_auto.csv and cp_embedding.pkl
TYPES = ["MSN", "FSI", "TAN", "UIN"]


# %% ── 1. A: the criteria ────────────────────────────────────────────────────
FEAT = io.load_features()
I = [i for i, r in enumerate(FEAT["rows"]) if r["region"] == "CP"]
R = [FEAT["rows"][i] for i in I]
f = lambda k: np.array([r[k] for r in R], float)
dur, pss, pli, fr, nsp = f("wf_duration_us"), f("post_spike_supp_ms"), f("prop_long_isi"), f("fr"), f("n_spikes")
pol = f("wf_polarity")

crit = np.full(len(R), C.UNLAB, object)
okk = (nsp >= MIN_SPIKES) & np.isfinite(dur)
wide = dur > DUR_US
crit[okk & wide & (pss < PSS_MS)] = "MSN"
crit[okk & wide & (pss >= PSS_MS)] = "TAN"
crit[okk & ~wide & (pli <= PROP_ISI)] = "FSI"
crit[okk & ~wide & (pli > PROP_ISI)] = "UIN"
near = (np.abs(dur - DUR_US) <= NEAR_FRAC * DUR_US) \
       | (wide & (np.abs(pss - PSS_MS) <= NEAR_FRAC * PSS_MS)) \
       | (~wide & (np.abs(pli - PROP_ISI) <= NEAR_FRAC * PROP_ISI))
# stricter FSI: each failed condition is kept, so the flag can say which one
n_bomb = int(np.sum(crit == "FSI")); isf = crit == "FSI"
FSI_FAIL = {
    "fsi_not_narrow": isf & (dur > FSI_DUR_US) if FSI_DUR_US is not None else isf & False,
    "fsi_fr_low":     isf & ~(fr >= FSI_FR_MIN) if FSI_FR_MIN is not None else isf & False,
    "fsi_positive":   isf & (pol > 0) if FSI_NEG_ONLY else isf & False}
for m in FSI_FAIL.values(): crit[m] = "UIN"
if FSI_DUR_US is not None: near |= ~wide & (np.abs(dur - FSI_DUR_US) <= NEAR_FRAC * FSI_DUR_US)
if FSI_FR_MIN is not None: near |= ~wide & (np.abs(fr - FSI_FR_MIN) <= NEAR_FRAC * FSI_FR_MIN)
print(f"FSI: {n_bomb} by bombcell's rule -> {int(np.sum(crit == 'FSI'))} after the stricter conditions;"
      f" moved to UIN:", {k: int(m.sum()) for k, m in FSI_FAIL.items()}, "(a unit can fail several)")
print(f"{len(R)} CP neurons; criteria:", dict(collections.Counter(crit)))


# %% ── 2. B: features, clusters, embedding ───────────────────────────────────
def _z(x):
    x = np.asarray(x, float)
    m = np.nanmedian(x, axis=0)
    x = np.where(np.isfinite(x), x, m)                     # missing -> median
    s = x.std(axis=0); s[s == 0] = 1
    return (x - x.mean(axis=0)) / s

A = FEAT["arrays"]
W = A["wf_norm"][I].astype(float); W = np.where(np.isfinite(W), W, 0.0)
G = A["acg"][I, :200].astype(float)                        # 0-200 ms
G = G / np.maximum(np.nanmean(A["acg"][I, 100:400], axis=1, keepdims=True), 1e-9)
G = np.log1p(np.where(np.isfinite(G), G, 0.0))
G = G.reshape(len(I), 40, 5).mean(axis=2)                  # 5 ms bins
wf_pc = PCA(N_WF_PC, random_state=SEED).fit_transform(W)
acg_pc = PCA(N_ACG_PC, random_state=SEED).fit_transform(G)
# two shape features that separate cell bodies from axons:
#   spread    = on how many of the 9 channels (peak +-4) the spike is >= SPREAD_FRAC of its peak
#   pre-peak  = the positive bump BEFORE the trough, relative to the trough (axons: large)
amp9 = np.nanmax(A["wf_multi"][I], axis=2) - np.nanmin(A["wf_multi"][I], axis=2)
spread = np.sum(amp9 >= SPREAD_FRAC * np.nanmax(amp9, axis=1, keepdims=True), axis=1).astype(float)
t_wf = np.asarray(FEAT["meta"]["wf_norm_t_ms"])
prepk = np.nanmax(W[:, t_wf < -0.05], axis=1)               # trough = -1, so this is the ratio
scal = np.c_[dur, f("wf_halfwidth_us"), f("wf_pt_ratio"), f("wf_repol_slope"),
             np.log10(np.clip(fr, 1e-3, None)), np.log10(np.clip(pss, 0.5, None)), pli, f("cv2"),
             spread, prepk]
SCAL_NAMES = ["duration", "half-width", "peak/trough", "repol. slope", "log FR",
              "log post-spike supp.", "prop. long ISI", "CV2", "spatial spread", "pre-peak"]
CRIT_FEATS = ["duration", "log post-spike supp.", "prop. long ISI"]   # left out of the classifier
# three blocks, each scaled to the same total weight
X = np.c_[_z(wf_pc) / np.sqrt(N_WF_PC), _z(acg_pc) / np.sqrt(N_ACG_PC), _z(scal) / np.sqrt(scal.shape[1])]

BIC = [GaussianMixture(k, covariance_type="full", n_init=3, random_state=SEED).fit(X).bic(X)
       for k in K_RANGE]
gm = GaussianMixture(N_CLUST, covariance_type="full", n_init=5, random_state=SEED).fit(X)
clu = gm.predict(X)
clu_name, clu_purity = {}, {}
for k in range(N_CLUST):
    c = collections.Counter(t for t in crit[clu == k] if t != C.UNLAB)
    clu_name[k] = c.most_common(1)[0][0] if c else C.UNLAB
    clu_purity[k] = c.most_common(1)[0][1] / max(sum(c.values()), 1) if c else 0.0
    print(f"  cluster {k}: n={np.sum(clu == k):4d} -> {clu_name[k]:5s} "
          f"(purity {clu_purity[k]:.2f})  {dict(c)}")

method = EMBED
if EMBED in ("auto", "umap"):
    try:
        import umap
        XY = umap.UMAP(n_neighbors=30, min_dist=0.2, random_state=SEED).fit_transform(X); method = "umap"
    except ImportError:
        if EMBED == "umap": raise
        method = "tsne"
if method == "tsne":
    from sklearn.manifold import TSNE
    XY = TSNE(2, perplexity=40, init="pca", random_state=SEED).fit_transform(X)
print(f"embedding: {method}")


# %% ── 3. C: classifier check ────────────────────────────────────────────────
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_predict
keep = [j for j, n in enumerate(SCAL_NAMES) if n not in CRIT_FEATS]
Xc = np.c_[_z(wf_pc), _z(acg_pc), _z(scal[:, keep])]
XC_NAMES = [f"wf PC{i + 1}" for i in range(N_WF_PC)] + [f"ACG PC{i + 1}" for i in range(N_ACG_PC)] \
           + [SCAL_NAMES[j] for j in keep]
cnt = collections.Counter(t for t in crit if t != C.UNLAB)
CLF_TYPES = sorted(t for t, n in cnt.items() if n >= 2)          # a class needs >= 2 neurons to be learned
lab_m = np.isin(crit, CLF_TYPES)
rf = lambda: RandomForestClassifier(500, min_samples_leaf=2, class_weight="balanced", n_jobs=-1,
                                    random_state=SEED)
P = np.full((len(R), len(CLF_TYPES)), np.nan)
if len(CLF_TYPES) >= 2:
    k = int(min(N_FOLDS, min(cnt[t] for t in CLF_TYPES)))
    P[lab_m] = cross_val_predict(rf(), Xc[lab_m], crit[lab_m].astype(str), method="predict_proba",
                                 cv=StratifiedKFold(k, shuffle=True, random_state=SEED))
    model = rf().fit(Xc[lab_m], crit[lab_m].astype(str))         # classes_ are sorted = CLF_TYPES
    if (~lab_m).any(): P[~lab_m] = model.predict_proba(Xc[~lab_m])
    IMPORT = model.feature_importances_
    print(f"classifier: {k}-fold, trained on", dict(cnt))
else:
    IMPORT = np.zeros(Xc.shape[1]); print("!! classifier skipped: fewer than two labelled types")
ok_p = np.isfinite(P).all(axis=1)
clf_pred = np.array([CLF_TYPES[i] if o else "" for i, o in zip(np.nan_to_num(P).argmax(1), ok_p)], object)
clf_prob = np.where(ok_p, np.nan_to_num(P).max(1), np.nan)
suggest = (clf_pred != "") & (clf_pred != crit) & (clf_prob >= RELABEL_P)
print(f"classifier agrees with the criteria for {np.mean(clf_pred[lab_m] == crit[lab_m]):.0%} of labelled neurons;"
      f" suggests another type for {suggest.sum()}:",
      dict(collections.Counter(f"{a}->{b}" for a, b in zip(crit[suggest], clf_pred[suggest]))))


# %% ── 4. the automatic label ────────────────────────────────────────────────
rows = []
for j, r in enumerate(R):
    lab = str(crit[j]); flags = []
    if nsp[j] < MIN_SPIKES: flags.append("few_spikes")
    if suggest[j]: flags.append(f"clf_suggests_{clf_pred[j]}")
    clu_ok = lab != C.UNLAB and clu_name[clu[j]] == lab
    if lab != C.UNLAB and not clu_ok: flags.append("cluster_disagrees")
    if lab != C.UNLAB and near[j]: flags.append("near_threshold")
    if lab == "MSN" and fr[j] > MSN_FR_MAX: flags.append("msn_fr_high")
    flags += [k for k, m in FSI_FAIL.items() if m[j]]
    high = clf_pred[j] == lab and not near[j] and (clu_ok or not CONF_USES_CLUSTER)
    rows.append(dict(uid=r["uid"], mouse=r["mouse"], sess=r["sess"], phyID=r["phyID"],
                     crit_cell_type=lab, cluster=int(clu[j]), cluster_label=clu_name[clu[j]],
                     clf_pred=clf_pred[j], clf_prob=float(clf_prob[j]),
                     suggested_cell_type=clf_pred[j] if suggest[j] else "",
                     **{f"clf_p_{t}": float(P[j, i]) for i, t in enumerate(CLF_TYPES)},
                     auto_cell_type=lab, auto_source="criteria" if lab != C.UNLAB else "",
                     auto_conf=("" if lab == C.UNLAB else "high" if high else "low"),
                     flags=";".join(flags), emb_x=float(XY[j, 0]), emb_y=float(XY[j, 1]),
                     spread=float(spread[j]), pre_peak=float(prepk[j])))
print("confidence:", dict(collections.Counter((r["auto_cell_type"], r["auto_conf"]) for r in rows)))
print("flags     :", dict(collections.Counter(x for r in rows for x in r["flags"].split(";") if x)))
if SAVE_VARS:
    io.write_csv(C.CP_AUTO_CSV, rows)
    io.save_pkl(C.CP_EMBED_PKL, dict(uids=[r["uid"] for r in R], xy=XY.astype(np.float32),
                                     cluster=clu, method=method, cluster_label=clu_name,
                                     X=X.astype(np.float32)))
    print(f"saved -> {C.CP_AUTO_CSV}\n         {C.CP_EMBED_PKL}")


# %% ── 5. figure: the criteria planes ────────────────────────────────────────
fdir = C.fig_dir("05_CP_classify") if SAVE_FIGS else None
col = np.array([C.COLORS[t] for t in crit])
fig, axes = plt.subplots(1, 3, figsize=(18, 5))
ax = axes[0]
ax.hist(dur[np.isfinite(dur)], np.linspace(0, 1000, 81), color="0.5")
ax.axvline(DUR_US, color="k", ls="--"); ax.set(xlabel="waveform duration (us)", ylabel="neurons",
                                               title="narrow | wide  (should be bimodal)")
ax = axes[1]
ax.scatter(dur, np.clip(pss, 0.5, 500), s=7, c=col, alpha=0.6, lw=0)
ax.axvline(DUR_US, color="k", ls="--", lw=0.8); ax.plot([DUR_US, 1000], [PSS_MS] * 2, "k--", lw=0.8)
ax.set_yscale("log"); ax.set(xlabel="waveform duration (us)", ylabel="post-spike suppression (ms)",
                             xlim=(0, 1000), title="wide side: MSN (below) vs TAN (above)")
ax = axes[2]
ax.scatter(dur, pli, s=7, c=col, alpha=0.6, lw=0)
ax.axvline(DUR_US, color="k", ls="--", lw=0.8); ax.plot([0, DUR_US], [PROP_ISI] * 2, "k--", lw=0.8)
if FSI_DUR_US is not None: ax.plot([FSI_DUR_US] * 2, [0, PROP_ISI], "k:", lw=1)
fail = np.any(list(FSI_FAIL.values()), axis=0)
ax.scatter(dur[fail], pli[fail], s=22, facecolors="none", edgecolors="k", lw=0.6)
ax.set(xlabel="waveform duration (us)", ylabel="proportion of time in ISIs > 2 s", xlim=(0, 1000),
       title=f"narrow side: FSI (below) vs UIN (above)\ncircled: failed stricter FSI (n={fail.sum()}); dotted: FSI max duration")
for t in TYPES + [C.UNLAB]:
    ax.scatter([], [], color=C.COLORS[t], label=f"{t} ({np.sum(crit == t)})")
ax.legend(frameon=False, fontsize=8)
for a in axes: a.spines[["top", "right"]].set_visible(False)
fig.suptitle("05 — A. criteria (Peters 2021 / bombcell)"); fig.tight_layout()
if SAVE_FIGS: fig.savefig(fdir / "criteria_planes.png", dpi=150)
if SHOW_FIGS: plt.show()
plt.close(fig)


# %% ── 6. figure: embedding, clusters, and how they agree with the criteria ──
fig, axes = plt.subplots(2, 3, figsize=(18, 11))
ax = axes[0, 0]; ax.scatter(XY[:, 0], XY[:, 1], s=7, c=col, alpha=0.7, lw=0)
ax.set_title(f"{method}: criteria label")
ax = axes[0, 1]; cm = plt.cm.tab10(np.arange(N_CLUST) % 10)
ax.scatter(XY[:, 0], XY[:, 1], s=7, c=cm[clu], alpha=0.7, lw=0)
for k in range(N_CLUST):
    if np.any(clu == k):
        ax.text(*np.median(XY[clu == k], axis=0), f"{k}:{clu_name[k]}", fontsize=10, fontweight="bold",
                ha="center", bbox=dict(fc="w", ec=cm[k], alpha=0.8, pad=1.5))
ax.set_title("mixture clusters (named after their majority criteria label)")
ax = axes[0, 2]
dis = np.array(["cluster_disagrees" in r["flags"] for r in rows])
ax.scatter(XY[~dis, 0], XY[~dis, 1], s=6, color="0.8", lw=0)
ax.scatter(XY[dis, 0], XY[dis, 1], s=10, c=col[dis], lw=0)
ax.set_title(f"criteria and cluster disagree (n={dis.sum()}), coloured by criteria label")
ax = axes[1, 0]; s = ax.scatter(XY[:, 0], XY[:, 1], s=7, c=np.log10(np.clip(fr, 0.05, None)), cmap="viridis", lw=0)
plt.colorbar(s, ax=ax, label="log10 firing rate (Hz)"); ax.set_title("firing rate")
ax = axes[1, 1]; ax.plot(list(K_RANGE), BIC, "ko-"); ax.axvline(N_CLUST, color="r", ls="--")
ax.set(xlabel="number of mixture components", ylabel="BIC (lower = better)",
       title="how many clusters do the data support?")
ax = axes[1, 2]
M = np.array([[np.sum((clu == k) & (crit == t)) for t in TYPES + [C.UNLAB]] for k in range(N_CLUST)])
ax.imshow(M / np.maximum(M.sum(1, keepdims=True), 1), cmap="Blues", vmin=0, vmax=1, aspect="auto")
for k in range(N_CLUST):
    for j in range(M.shape[1]):
        if M[k, j]: ax.text(j, k, M[k, j], ha="center", va="center", fontsize=9)
ax.set_xticks(range(M.shape[1]), TYPES + [C.UNLAB]); ax.set_yticks(range(N_CLUST), [f"{k}:{clu_name[k]}" for k in range(N_CLUST)])
ax.set(xlabel="criteria label", ylabel="cluster", title="neurons per cluster and criteria label")
for a in axes.ravel()[[0, 1, 2, 3]]: a.set_xticks([]); a.set_yticks([])
fig.suptitle("05 — B. unsupervised view vs criteria"); fig.tight_layout()
if SAVE_FIGS: fig.savefig(fdir / "embedding_clusters.png", dpi=150)
if SHOW_FIGS: plt.show()
plt.close(fig)


# %% ── 7. figure: the classifier check ───────────────────────────────────────
# a) criteria label (rows) vs the cross-validated prediction (columns): the
#    off-diagonal cells are the candidates for relabelling
# b) how sure the classifier is, when it agrees and when it does not
# c) the embedding with the suggestions: filled = criteria colour, edge = suggested type
# d) which features it relied on (the criteria features are not among them)
fig, axes = plt.subplots(1, 4, figsize=(22, 5))
ax = axes[0]
Mc = np.array([[np.sum(lab_m & (crit == a) & (clf_pred == b)) for b in CLF_TYPES] for a in CLF_TYPES])
ax.imshow(Mc / np.maximum(Mc.sum(1, keepdims=True), 1), cmap="Blues", vmin=0, vmax=1)
for a in range(len(CLF_TYPES)):
    for b in range(len(CLF_TYPES)):
        ax.text(b, a, Mc[a, b], ha="center", va="center", fontsize=10)
ax.set_xticks(range(len(CLF_TYPES)), CLF_TYPES); ax.set_yticks(range(len(CLF_TYPES)), CLF_TYPES)
ax.set(xlabel="classifier (cross-validated)", ylabel="criteria label",
       title=f"agreement {np.mean(clf_pred[lab_m] == crit[lab_m]):.0%}")
ax = axes[1]; ag = lab_m & (clf_pred == crit); dg = lab_m & (clf_pred != crit)
ax.hist([clf_prob[ag], clf_prob[dg]], np.linspace(0, 1, 21), color=["0.6", "#e76f51"],
        label=[f"agrees ({ag.sum()})", f"disagrees ({dg.sum()})"], stacked=True)
ax.axvline(RELABEL_P, color="k", ls="--"); ax.legend(frameon=False, fontsize=8)
ax.set(xlabel="classifier probability of its prediction", ylabel="neurons",
       title=f"right of the line = suggested ({suggest.sum()})")
ax = axes[2]; ax.scatter(XY[:, 0], XY[:, 1], s=6, color="0.85", lw=0)
for t in CLF_TYPES:
    s = suggest & (clf_pred == t)
    ax.scatter(XY[s, 0], XY[s, 1], s=30, c=col[s], edgecolors=C.COLORS[t], lw=1.6)
    ax.scatter([], [], s=30, color="w", edgecolors=C.COLORS[t], lw=1.6, label=f"-> {t} ({s.sum()})")
ax.legend(frameon=False, fontsize=8); ax.set_xticks([]); ax.set_yticks([])
ax.set_title(f"{method}: suggestions (fill = criteria, edge = suggested)")
ax = axes[3]; o = np.argsort(IMPORT)
ax.barh(np.array(XC_NAMES)[o], IMPORT[o], color="0.5"); ax.set(xlabel="importance", title="what the classifier used")
for a in axes[[1, 3]]: a.spines[["top", "right"]].set_visible(False)
fig.suptitle("05 — C. classifier check (never changes a label; suggestions are reviewed in the GUI)")
fig.tight_layout()
if SAVE_FIGS: fig.savefig(fdir / "classifier_check.png", dpi=150)
if SHOW_FIGS: plt.show()
plt.close(fig)


# %% ── 8. figure: does each type look like itself? ───────────────────────────
# Per type: waveforms, ACGs (scaled to the 100-400 ms level) and firing rates.
# Bold = mean of the high-confidence neurons, thin grey = low-confidence ones.
t_ms = FEAT["meta"]["wf_norm_t_ms"]; lag = (np.arange(200) + 0.5)
Gr = A["acg"][I, :200].astype(float) / np.maximum(np.nanmean(A["acg"][I, 100:400], axis=1, keepdims=True), 1e-9)
conf = np.array([r["auto_conf"] for r in rows])
fig, axes = plt.subplots(3, len(TYPES), figsize=(4 * len(TYPES), 9))
for j, t in enumerate(TYPES):
    hi, lo = (crit == t) & (conf == "high"), (crit == t) & (conf == "low")
    for k in np.flatnonzero(lo)[:60]:
        axes[0, j].plot(t_ms, W[k], color="0.6", lw=0.4); axes[1, j].plot(lag, Gr[k], color="0.6", lw=0.4)
    if hi.any():
        axes[0, j].plot(t_ms, W[hi].mean(0), color=C.COLORS[t], lw=2.5)
        axes[1, j].plot(lag, np.nanmean(Gr[hi], 0), color=C.COLORS[t], lw=2.5)
    v = fr[(crit == t) & np.isfinite(fr) & (fr > 0)]
    axes[2, j].hist(v, np.geomspace(0.01, 100, 40), color=C.COLORS[t]); axes[2, j].set_xscale("log")
    axes[0, j].set_title(f"{t}: {hi.sum()} high, {lo.sum()} low confidence", color=C.COLORS[t], fontsize=10)
    axes[1, j].set_ylim(0, 4)
    axes[0, j].set_xlabel("ms"); axes[1, j].set_xlabel("ACG lag (ms)"); axes[2, j].set_xlabel("firing rate (Hz)")
    for a in axes[:, j]: a.spines[["top", "right"]].set_visible(False)
axes[0, 0].set_ylabel("normalised waveform"); axes[1, 0].set_ylabel("ACG / baseline"); axes[2, 0].set_ylabel("neurons")
fig.suptitle("05 — each striatal type: waveform, ACG, firing rate"); fig.tight_layout()
if SAVE_FIGS: fig.savefig(fdir / "types_overview.png", dpi=150)
if SHOW_FIGS: plt.show()
plt.close(fig)


# %% ── 9. figure: proportions per mouse and session ──────────────────────────
# Expected in dorsal striatum: mostly MSN, a few % FSI, ~1-2% TAN. A session far
# from the others points to a problem in that recording (or its region borders).
keys = sorted({(r["mouse"], r["sess"]) for r in rows})
fig, ax = plt.subplots(figsize=(15, 4.5)); bottom = np.zeros(len(keys))
for t in TYPES + [C.UNLAB]:
    n = np.array([sum(r["auto_cell_type"] == t and (r["mouse"], r["sess"]) == k for r in rows) for k in keys], float)
    tot = np.array([sum((r["mouse"], r["sess"]) == k for r in rows) for k in keys], float)
    ax.bar(range(len(keys)), n / tot, bottom=bottom, color=C.COLORS[t], label=t, edgecolor="w"); bottom += n / tot
ax.set_xticks(range(len(keys)), [f"{m[2:6]} {s}" for m, s in keys], rotation=60, fontsize=8)
ax.set_ylabel("fraction of CP neurons"); ax.legend(frameon=False, fontsize=8, ncol=5, loc="lower right")
ax.spines[["top", "right"]].set_visible(False); fig.tight_layout()
if SAVE_FIGS: fig.savefig(fdir / "proportions_per_session.png", dpi=150)
if SHOW_FIGS: plt.show()
plt.close(fig)
# %%
