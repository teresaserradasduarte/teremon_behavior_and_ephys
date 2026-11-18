# %% ── 0. parameters ─────────────────────────────────────────────────────────
# STEP 01 -- features of every cerebellar-cortex (CB) and striatum (CP) neuron.
#
# Reads, per session, eg_neurons_CB.mat / eg_neurons_BG.mat (read only): the
# whole-session spike times `st` and the template `templateWaveforms`.
# Per neuron it computes
#   waveform   peak channel, peak-to-trough duration, half-width, peak/trough
#              ratio, polarity, repolarisation slope       (ct_features.waveform_features)
#   train      firing rate, CV2, median ISI, bombcell's propLongISI and
#              postSpikeSuppression, refractory period, burst index
#   arrays     normalised peak-channel waveform, the waveform on the peak channel
#              +-N_NEIGH channels, the autocorrelogram (1 ms bins, 0-500 ms)
# and writes ONE file, features.pkl. Every later step starts from it.
# Figures: the distribution of each feature per region, and example neurons.
import time
import numpy as np
import matplotlib.pyplot as plt
import ct_config as C, ct_io as io, ct_features as F

MICE     = None        # None = all; or e.g. ["1_CoteDor"]
SESSIONS = None        # None = all; or e.g. ["R5"]
N_NEIGH  = 4           # channels kept on each side of the peak channel (for the GUI)
ACG_BIN, ACG_WIN = 0.001, 0.5      # s (bombcell: 1 ms bins, +-500 ms)
MIN_SPIKES = 100       # fewer spikes -> features left NaN

SHOW_FIGS = True
SAVE_FIGS = True
SAVE_VARS = True       # writes features.pkl

SCALARS = ["fr", "cv2", "isi_median_ms", "prop_long_isi", "post_spike_supp_ms",
           "refractory_ms", "burst_index", "acg_peak_ms", "wf_duration_us",
           "wf_halfwidth_us", "wf_pt_ratio", "wf_repol_slope", "wf_amp"]


# %% ── 1. extract ────────────────────────────────────────────────────────────
rows, WFN, WFM, ACG = [], [], [], []
t00 = time.time()
for mouse, sess, kind, path in io.sessions(mice=MICE, sess=SESSIONS):
    region = C.KINDS[kind]
    t0 = time.time()
    neus = io.load_eg_neurons(path, regions=(region,))
    for d in neus:
        W, st = d["W"], d["st"]
        pk = F.peak_channel(W)
        w = W[pk]
        r = dict(uid=C.uid(region, mouse, sess, d["phyID"]), region=region, mouse=mouse,
                 sess=sess, phyID=d["phyID"], kind=kind, depth=d["depth"], peak_ch=pk,
                 n_wf_samples=W.shape[1], t_first=float(st[0]) if st.size else np.nan,
                 t_last=float(st[-1]) if st.size else np.nan)
        r.update(F.waveform_features(w, C.FS))
        if st.size >= MIN_SPIKES:
            a, lags = F.acg(st, ACG_BIN, ACG_WIN)
        else:
            a = np.full(int(round(ACG_WIN / ACG_BIN)), np.nan)
            lags = (np.arange(a.size) + 0.5) * ACG_BIN
        r.update(F.train_features(st if st.size >= MIN_SPIKES else st[:0], a, lags))
        r["n_spikes"] = int(st.size)
        rows.append(r)
        WFN.append(F.norm_waveform(w))
        ch = np.clip(np.arange(pk - N_NEIGH, pk + N_NEIGH + 1), 0, W.shape[0] - 1)
        m = np.full((2 * N_NEIGH + 1, 82), np.nan)
        m[:, :min(82, W.shape[1])] = W[ch][:, :82]
        WFM.append(m); ACG.append(a)
    print(f"{mouse} {sess} {kind}: {len(neus)} {region} neurons  ({time.time() - t0:.0f}s)")

FEAT = dict(rows=rows,
            arrays=dict(wf_norm=np.asarray(WFN, np.float32), wf_multi=np.asarray(WFM, np.float32),
                        acg=np.asarray(ACG, np.float32)),
            meta=dict(acg_bin=ACG_BIN, acg_win=ACG_WIN, n_neigh=N_NEIGH, fs=C.FS,
                      wf_norm_t_ms=np.linspace(-1.0, 1.7, 82), scalars=SCALARS,
                      made=time.strftime("%Y-%m-%d %H:%M")))
print(f"\n{len(rows)} neurons: " + ", ".join(
    f"{reg} {sum(r['region'] == reg for r in rows)}" for reg in ("CB", "CP"))
    + f"   ({time.time() - t00:.0f}s)")
if SAVE_VARS:
    io.save_pkl(C.FEATURES_PKL, FEAT); print(f"saved -> {C.FEATURES_PKL}")


# %% ── 2. figure: feature distributions per region ──────────────────────────
# One histogram per feature, CB and CP overlaid. Sanity check of the extraction:
# e.g. waveform duration should be bimodal in CP (narrow FSI vs wide MSN/TAN).
LOGX = {"fr", "isi_median_ms", "prop_long_isi", "post_spike_supp_ms", "burst_index", "wf_amp"}
fig, axes = plt.subplots(3, 5, figsize=(19, 10)); axes = axes.ravel()
for ax, k in zip(axes, SCALARS):
    for reg, c in (("CB", "#1c78b5"), ("CP", "#d1495b")):
        v = np.array([r[k] for r in rows if r["region"] == reg], float)
        v = v[np.isfinite(v)]
        if k in LOGX: v = v[v > 0]
        if v.size < 2: continue
        lo, hi = np.percentile(v, [0.5, 99.5])
        if hi <= lo: continue
        bins = np.geomspace(lo, hi, 50) if k in LOGX else np.linspace(lo, hi, 50)
        ax.hist(v, bins, histtype="stepfilled", alpha=0.45, color=c, label=f"{reg} (n={v.size})")
    if k in LOGX: ax.set_xscale("log")
    ax.set_title(k, fontsize=10); ax.spines[["top", "right"]].set_visible(False)
axes[0].legend(fontsize=8, frameon=False)
for ax in axes[len(SCALARS):]: ax.set_visible(False)
fig.suptitle("01 — feature distributions (all neurons)", fontsize=12); fig.tight_layout()
if SAVE_FIGS: fig.savefig(C.fig_dir("01_features") / "feature_distributions.png", dpi=150)
if SHOW_FIGS: plt.show()
plt.close(fig)


# %% ── 3. figure: example neurons (waveform + ACG) ───────────────────────────
# A random sample per region: what the stored waveform and ACG look like, with
# the measured duration / suppression marked, so the features can be checked by eye.
rs = np.random.default_rng(0)
t_ms = FEAT["meta"]["wf_norm_t_ms"]; lag_ms = (np.arange(ACG[0].size) + 0.5) * ACG_BIN * 1e3
for reg in ("CB", "CP"):
    idx = [i for i, r in enumerate(rows) if r["region"] == reg and r["n_spikes"] >= MIN_SPIKES]
    if not idx: continue
    pick = rs.choice(idx, min(12, len(idx)), replace=False)
    fig, axes = plt.subplots(4, 6, figsize=(18, 10))
    for j, i in enumerate(pick):
        r = rows[i]
        a1, a2 = axes[(j // 6) * 2, j % 6], axes[(j // 6) * 2 + 1, j % 6]
        a1.plot(t_ms, WFN[i], "k"); a1.axvline(0, color="r", lw=0.6)
        a1.axvline(r["wf_duration_us"] / 1e3, color="r", lw=0.6, ls="--")
        a1.set_title(f"{r['mouse'][2:]} {r['sess']} #{r['phyID']}\n"
                     f"dur {r['wf_duration_us']:.0f} us, {r['fr']:.1f} Hz", fontsize=8)
        a2.bar(lag_ms[:150], ACG[i][:150], width=1.0, color="0.4")
        if np.isfinite(r["post_spike_supp_ms"]) and r["post_spike_supp_ms"] <= 150:
            a2.axvline(r["post_spike_supp_ms"], color="r", lw=0.8)
        a2.set_xlabel("lag (ms)", fontsize=7)
        for a in (a1, a2): a.tick_params(labelsize=7); a.spines[["top", "right"]].set_visible(False)
    fig.suptitle(f"01 — {reg} examples: normalised waveform (red: trough, dashed: following "
                 f"peak) and ACG (red: post-spike suppression)", fontsize=11)
    fig.tight_layout()
    if SAVE_FIGS: fig.savefig(C.fig_dir("01_features") / f"examples_{reg}.png", dpi=150)
    if SHOW_FIGS: plt.show()
    plt.close(fig)
# %%
