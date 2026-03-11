#%%
import os
import subprocess
import gc
from time import time, sleep
import numpy as np
from collections import Counter
import json
import tempfile
import csv
import shutil


def ensure_file_available(path, timeout=60):
    """Ensure a Dropbox (or OneDrive) cloud-only file is downloaded and available.

    On Windows, uses `attrib` to pin the file (make available offline) and waits
    until it can be opened. Returns True if file is ready, False on timeout.
    """
    if not os.path.exists(path):
        return False  # file doesn't exist at all

    # Quick check: try opening the file directly
    try:
        with open(path, 'rb') as f:
            f.read(1)
        return True
    except OSError:
        pass  # likely cloud placeholder, proceed to hydrate

    if os.name != 'nt':
        return False  # only Windows cloud placeholders handled here

    # Pin the file to make it available offline (works for Dropbox & OneDrive)
    # attrib -U +P removes Unpinned attribute and sets Pinned
    try:
        subprocess.run(['attrib', '-U', '+P', path], check=False,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except Exception:
        pass

    # Wait until the file becomes readable
    start = time()
    while time() - start < timeout:
        try:
            with open(path, 'rb') as f:
                f.read(1)
            return True
        except OSError:
            sleep(60)
    return False


def load_npy(path):
    """Robust np.load wrapper.

    1. Ensures cloud-only files (Dropbox/OneDrive) are downloaded first.
    2. Tries a normal load, then Windows long-path prefix, then temp-copy fallback.
    Raises the original exception if all strategies fail.
    """
    # Make sure cloud file is hydrated/available
    ensure_file_available(path, timeout=3600)
    print(f"{path} File now available")
    try:
        return np.load(path, allow_pickle=True)
    except OSError as first_exc:
        # Try Windows long-path prefix if applicable
        if os.name == 'nt':
            try:
                p2 = "\\\\?\\" + os.path.abspath(path)
                return np.load(p2, allow_pickle=True)
            except Exception:
                pass
        # Try copying to a temp file (handles some path/provider edge cases)
        try:
            tmp_fd, tmp_path = tempfile.mkstemp(suffix=os.path.splitext(path)[1])
            os.close(tmp_fd)
            shutil.copy2(path, tmp_path)
            try:
                arr = np.load(tmp_path, allow_pickle=True)
            finally:
                try:
                    os.remove(tmp_path)
                except Exception:
                    pass
            return arr
        except Exception:
            pass
        # re-raise the original error for visibility
        raise first_exc


def compute_positions_from_pc(pc_feat_path, pc_feat_ind_path, spike_templates, chan_x, chan_y):
    """
    Compute spike x/y positions from pc_features and pc_feature_ind (matches MATLAB logic).
    pc_features: nSpikes x nFeatures x nLocalChannels
    pc_feature_ind: nTemplates x nLocalChannels (indices of channels for each local channel)
    spike_templates: per-spike template ids (0-based)
    chan_x, chan_y: arrays of channel coordinates

    Returns: spike_x, spike_y (floats, length nSpikes)
    """
    pc = load_npy(pc_feat_path)
    pc_ind = load_npy(pc_feat_ind_path)

    # Ensure pc has axes (nSpikes, nFeatures, nLocalChannels).
    # Heuristics: we know number of channels from chan_x and number of spikes from spike_templates.
    n_channels = chan_x.shape[0]
    spikes_len = int(spike_templates.shape[0])
    if pc.ndim != 3:
        raise ValueError(f"Unexpected shape for pc_features.npy: {pc.shape}. Expected (nSpikes, nFeatures, nLocalChannels)")
    # If axes are in a different order, try to detect and permute to (spikes, features, channels)
    if not (pc.shape[0] == spikes_len and pc.shape[2] == n_channels):
        # attempt to locate axes
        axis_spikes = None
        axis_channels = None
        for i, s in enumerate(pc.shape):
            if s == spikes_len and axis_spikes is None:
                axis_spikes = i
            elif s == n_channels and axis_channels is None:
                axis_channels = i
        if axis_spikes is not None and axis_channels is not None:
            axis_features = [0,1,2]
            axis_features.remove(axis_spikes)
            axis_features.remove(axis_channels)
            axis_features = axis_features[0]
            perm = (axis_spikes, axis_features, axis_channels)
            if perm != (0,1,2):
                pc = np.transpose(pc, perm)
                # now pc should be (nSpikes, nFeatures, nLocalChannels)
        # else: leave as-is and let the subsequent shape check raise if still wrong

    # take first PC (feature 0) like MATLAB's pcFeat(:,1,:)
    pc1 = pc[:, 0, :].astype(float)  # (nSpikes, nLocalChannels)
    pc1[pc1 < 0] = 0.0

    if pc_ind.ndim != 2:
        raise ValueError(f"Unexpected shape for pc_feature_ind.npy: {pc_ind.shape}. Expected (nTemplates, nLocalChannels)")

    spike_templates = np.asarray(spike_templates).astype(int)
    n_spikes = pc1.shape[0]
    if spike_templates.shape[0] != n_spikes:
        raise ValueError(f"Mismatch: pc_features has {n_spikes} spikes but spike_templates has {spike_templates.shape[0]} entries")

    # get local channel indices per spike from templates
    # assume pc_feature_ind is 0-based; if values exceed channel count, handle 1-based
    spike_feat_ind = pc_ind[spike_templates, :].astype(int)  # (nSpikes, nLocalChannels)

    n_channels = chan_x.shape[0]
    max_ind = spike_feat_ind.max()
    min_ind = spike_feat_ind.min()
    if min_ind < 0:
        raise ValueError("pc_feature_ind contains negative indices")
    if max_ind >= n_channels:
        # maybe 1-based indexing from MATLAB; subtract 1
        if max_ind == n_channels:
            spike_feat_ind = spike_feat_ind - 1
        else:
            raise ValueError(f"pc_feature_ind max index {max_ind} >= number of channels {n_channels}")

    spike_feat_xcoords = chan_x[spike_feat_ind]
    spike_feat_ycoords = chan_y[spike_feat_ind]

    w = pc1 ** 2
    denom = w.sum(axis=1)
    denom_safe = denom.copy()
    denom_safe[denom_safe == 0] = 1.0

    spike_x = (spike_feat_xcoords * w).sum(axis=1) / denom_safe
    spike_y = (spike_feat_ycoords * w).sum(axis=1) / denom_safe

    return spike_x, spike_y

def main(ks_dir, out_txt):
    start_time = time()
    # Load channel positions
    ch_pos_path = os.path.join(ks_dir, 'channel_positions.npy')
    if not os.path.exists(ch_pos_path):
        raise FileNotFoundError(f"Missing {ch_pos_path}")
    ch_pos = load_npy(ch_pos_path)  # expected shape (nChannels, 2) with x,y columns
    if ch_pos.ndim != 2 or ch_pos.shape[1] < 2:
        raise ValueError(f"Unexpected shape for channel_positions.npy: {ch_pos.shape}")
    chan_x = ch_pos[:,0]
    chan_y = ch_pos[:,1]

    # load spike_templates early (used for pc_features fallback)
    templates_path = os.path.join(ks_dir, 'spike_templates.npy')
    spike_templates = None
    if os.path.exists(templates_path):
        spike_templates = load_npy(templates_path).astype(int)

    # Load spike positions (if present)
    sp_pos_path = os.path.join(ks_dir, 'spike_positions.npy')
    assigned_chan = None
    if os.path.exists(sp_pos_path):
        sp_pos = load_npy(sp_pos_path)
        if sp_pos.ndim == 2 and sp_pos.shape[1] >= 2:
            spike_x = sp_pos[:,0].astype(float)
            spike_y = sp_pos[:,1].astype(float)
        elif sp_pos.ndim == 1:
            # single column: likely channel indices (0-based)
            spike_chan = sp_pos.astype(int)
            # validate indices
            if spike_chan.max() >= ch_pos.shape[0] or spike_chan.min() < 0:
                raise ValueError("spike_positions.npy appears to contain integer values outside channel index range.")
            assigned_chan = spike_chan  # directly use these channel indices (0-based)
            spike_x = chan_x[assigned_chan]
            spike_y = chan_y[assigned_chan]
        else:
            raise ValueError(f"Unexpected spike_positions.npy shape: {sp_pos.shape}")
    else:
        # Try fallback: compute positions from pc_features.npy + pc_feature_ind.npy and spike_templates.npy
        pc_feat_path = os.path.join(ks_dir, 'pc_features.npy')
        pc_feat_ind_path = os.path.join(ks_dir, 'pc_feature_ind.npy')


        if os.path.exists(pc_feat_path) and os.path.exists(pc_feat_ind_path) and spike_templates is not None:
            spike_x, spike_y = compute_positions_from_pc(pc_feat_path, pc_feat_ind_path, spike_templates, chan_x, chan_y)
        else:
            raise FileNotFoundError("No spike_positions.npy found and missing pc_features/pc_feature_ind/spike_templates fallback files. "
                                    "Provide `spike_positions.npy` or the pc feature files `pc_features.npy` + `pc_feature_ind.npy` plus `spike_templates.npy`.")
    print(f"Loaded and computed spike positions for {spike_x.shape[0]} spikes")
    
    stop_time = time()
    print(f"Elapsed time: {stop_time - start_time:.2f} seconds")
    
    n_spikes = spike_x.shape[0]

    # We DO NOT build a full (n_spikes x n_channels) distance matrix (would OOM).
    # Instead, for each cluster we compute the median spike depth and assign the cluster
    # to the channel whose y-coordinate is closest to that median (matches requested approach).

    # Load cluster labels (prefer spike_clusters.npy)
    clusters_path = os.path.join(ks_dir, 'spike_clusters.npy')
    if os.path.exists(clusters_path):
        spike_clusters = load_npy(clusters_path).astype(int)
    else:
        templates_path = os.path.join(ks_dir, 'spike_templates.npy')
        if os.path.exists(templates_path):
            spike_clusters = load_npy(templates_path).astype(int)
        else:
            raise FileNotFoundError("Neither spike_clusters.npy nor spike_templates.npy found in ksDir.")

    if spike_clusters.shape[0] != n_spikes:
        raise ValueError(f"Mismatch: {n_spikes} spikes (from positions) but {spike_clusters.shape[0]} labels (from clusters/templates).")

    # Aggregate per-cluster: compute median spike depth per cluster and assign nearest channel
    unique_clusters = np.unique(spike_clusters)
    mapping = []
    for c in np.sort(unique_clusters):
        idx = np.where(spike_clusters == c)[0]
        if idx.size == 0:
            continue
        # median x,y of spikes in this cluster (use nanmedian to be robust)
        median_x = float(np.nanmedian(spike_x[idx]))
        median_y = float(np.nanmedian(spike_y[idx]))
        # compute squared Euclidean distance from cluster median to each channel
        dx = chan_x - median_x
        dy = chan_y - median_y
        d2 = dx*dx + dy*dy
        ch_nearest = int(np.nanargmin(d2))
        mapping.append((int(c), ch_nearest))

    # Save mapping to text file: two columns cluster_id channel (tab-separated)
    with open(out_txt, 'w') as fh:
        for c, ch in mapping:
            fh.write(f"{c}\t{ch}\n")

    print(f"Wrote {len(mapping)} cluster->channel mappings to {out_txt}")
    print("Note: channel indices are 0-based (matching typical Kilosort outputs).")


def create_cluster_region_file(cluster_map_path, channel_locations_path, out_path):
    """Read a cluster->channel map (two-column txt) and a channel locations file
    (JSON/.npy) or simplified TSV and write a three-column file: cluster_id, channel, brain_region
    """
    # Load cluster->channel mapping
    mapping = []
    with open(cluster_map_path, 'r', encoding='utf8') as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                c = int(parts[0])
                ch = int(parts[1])
            except ValueError:
                # skip header or malformed
                continue
            mapping.append((c, ch))

    if len(mapping) == 0:
        raise ValueError(f"No valid cluster->channel entries found in {cluster_map_path}")

    # load channel->region map (supports json, npy, or simplified TSV)
    chan_map = load_channel_map(channel_locations_path)

    # write output
    with open(out_path, 'w', encoding='utf8') as outfh:
        outfh.write('cluster_id\tchannel\tbrain_region\n')
        for c, ch in mapping:
            region = chan_map.get(ch)
            if region is None:
                region = chan_map.get(ch + 1)
            if region is None:
                region = ''
            outfh.write(f"{c}\t{ch}\t{region}\n")

    print(f"Wrote cluster->channel->region mapping to {out_path}")


def load_channel_map(channel_locations_path):
    """Load channel->brain_region mapping from JSON or .npy file into dict {channel_idx: region}.
    Reuses logic from previous implementation but factored out for reuse.
    """
    chan_map = {}
    path = channel_locations_path
    lower = path.lower()
    # support simplified TSV/CSV files too
    if lower.endswith('.tsv') or lower.endswith('.txt') or lower.endswith('.csv'):
        # try TSV with header: channel, brain_region, brain_region_id, axial, lateral, x, y, z
        with open(path, 'r', encoding='utf8') as fh:
            hdr = fh.readline()
            if not hdr:
                return chan_map
            headers = [h.strip().lower() for h in hdr.strip().split('\t')]
            # find channel and brain_region columns
            try:
                ch_i = headers.index('channel')
            except ValueError:
                ch_i = None
            try:
                region_i = headers.index('brain_region')
            except ValueError:
                region_i = None
            # fallback: try comma-split if no tabs
            delim = '\t' if '\t' in hdr else ','
            fh.seek(0)
            for line in fh:
                parts = [p.strip() for p in line.strip().split(delim)]
                if len(parts) == 0:
                    continue
                if ch_i is None:
                    # try first column as channel index or 'channel_0' string
                    ch_val = parts[0]
                    if ch_val.startswith('channel_'):
                        try:
                            idx = int(ch_val.split('_',1)[1])
                        except Exception:
                            continue
                    else:
                        try:
                            idx = int(ch_val)
                        except Exception:
                            continue
                else:
                    ch_val = parts[ch_i]
                    if ch_val.startswith('channel_'):
                        try:
                            idx = int(ch_val.split('_',1)[1])
                        except Exception:
                            continue
                    else:
                        try:
                            idx = int(ch_val)
                        except Exception:
                            continue
                region = ''
                if region_i is not None and region_i < len(parts):
                    region = parts[region_i]
                else:
                    # try column 1
                    if len(parts) > 1:
                        region = parts[1]
                chan_map[idx] = region
        return chan_map

    if lower.endswith('.json'):
        with open(channel_locations_path, 'r', encoding='utf8') as jh:
            data = json.load(jh)
        for k, v in data.items():
            if k.startswith('channel_'):
                try:
                    idx = int(k.split('_', 1)[1])
                except Exception:
                    continue
                region = v.get('brain_region') if isinstance(v, dict) else None
                if region is None and isinstance(v, dict):
                    region = str(v.get('brain_region_id', ''))
                chan_map[idx] = region
        return chan_map

    # else try .npy
    arr = np.load(channel_locations_path, allow_pickle=True)
    if isinstance(arr, dict):
        for k, v in arr.items():
            if k.startswith('channel_'):
                try:
                    idx = int(k.split('_', 1)[1])
                except Exception:
                    continue
                region = v.get('brain_region') if isinstance(v, dict) else None
                if region is None and isinstance(v, dict):
                    region = str(v.get('brain_region_id', ''))
                chan_map[idx] = region
        return chan_map

    if isinstance(arr, np.ndarray) and arr.dtype == object:
        if arr.size == 1 and isinstance(arr.item(), dict):
            for k, v in arr.item().items():
                if k.startswith('channel_'):
                    try:
                        idx = int(k.split('_', 1)[1])
                    except Exception:
                        continue
                    region = v.get('brain_region') if isinstance(v, dict) else None
                    if region is None and isinstance(v, dict):
                        region = str(v.get('brain_region_id', ''))
                    chan_map[idx] = region
            return chan_map
        else:
            for elem in arr:
                if isinstance(elem, dict):
                    idx = elem.get('original_channel_idx')
                    if idx is None:
                        continue
                    region = elem.get('brain_region') or str(elem.get('brain_region_id', ''))
                    chan_map[int(idx)] = region
            if chan_map:
                return chan_map

    try:
        for elem in arr:
            idx = None
            if hasattr(arr.dtype, 'names') and 'original_channel_idx' in arr.dtype.names:
                idx = int(elem['original_channel_idx'])
            elif hasattr(arr.dtype, 'names') and 'channel' in arr.dtype.names:
                idx = int(elem['channel'])
            if idx is None:
                continue
            region = None
            if hasattr(arr.dtype, 'names') and 'brain_region' in arr.dtype.names:
                region = elem['brain_region'].decode() if isinstance(elem['brain_region'], bytes) else str(elem['brain_region'])
            elif hasattr(arr.dtype, 'names') and 'brain_region_id' in arr.dtype.names:
                region = str(elem['brain_region_id'])
            chan_map[idx] = region
    except Exception:
        pass

    return chan_map


def create_cluster_region_from_clusterinfo(cluster_info_path, channel_locations_path, out_path):
    """Read a TSV/CSV cluster info file with headers (must include 'cluster_id' and 'ch')
    and produce an output file with columns: cluster_id, channel, brain_region.
    """
    # Detect delimiter
    with open(cluster_info_path, 'r', encoding='utf8') as fh:
        sample = fh.read(4096)
        sniffer = csv.Sniffer()
        dialect = sniffer.sniff(sample) if sample else csv.excel_tab
        fh.seek(0)
        reader = csv.DictReader(fh, dialect=dialect)

        # Normalize header names
        headers = [h.strip() for h in reader.fieldnames]
        lower = [h.lower() for h in headers]
        if 'cluster_id' in lower:
            cid_key = headers[lower.index('cluster_id')]
        elif 'cluster' in lower:
            cid_key = headers[lower.index('cluster')]
        else:
            raise ValueError('cluster info file missing "cluster_id" or "cluster" column')

        if 'ch' in lower:
            # load channel map via helper to support JSON/.npy/TSV
            chan_map = load_channel_map(channel_locations_path)

            # write output
            with open(out_path, 'w', encoding='utf8') as outfh:
                outfh.write('cluster_id\tchannel\tbrain_region\n')
                for c, ch in mapping:
                    region = chan_map.get(ch)
                    if region is None:
                        region = chan_map.get(ch + 1)
                    if region is None:
                        region = ''
                    outfh.write(f"{c}\t{ch}\t{region}\n")

            print(f'Wrote cluster->channel->region mapping to {out_path}')

    # load mapping
    cluster_region_dict = load_cluster_region_dict(cluster_region_map_path)

    # read extragood file and collect cluster ids
    extragood_ids = []
    with open(extragood_path, 'r', encoding='utf8') as fh:
        sample = fh.read(4096)
        fh.seek(0)
        try:
            sniffer = csv.Sniffer()
            dialect = sniffer.sniff(sample) if sample else csv.excel_tab
        except Exception:
            dialect = csv.excel_tab
        fh.seek(0)
        reader = csv.DictReader(fh, dialect=dialect)
        if reader.fieldnames:
            headers = [h.strip() for h in reader.fieldnames]
            lower = [h.lower() for h in headers]
            if 'cluster_id' in lower:
                cid_key = headers[lower.index('cluster_id')]
            elif 'cluster' in lower:
                cid_key = headers[lower.index('cluster')]
            else:
                cid_key = None
            if cid_key is not None:
                for row in reader:
                    try:
                        cid = int(row[cid_key])
                    except Exception:
                        continue
                    extragood_ids.append(cid)
        else:
            # fallback: file maybe a single-column list
            fh.seek(0)
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                try:
                    extragood_ids.append(int(line.split()[0]))
                except Exception:
                    continue

    # write output file
    with open(out_path, 'w', encoding='utf8') as outfh:
        outfh.write('cluster_id\tchannel\tbrain_region\n')
        for cid in extragood_ids:
            ch_region = cluster_region_dict.get(cid)
            if ch_region is None:
                outfh.write(f"{cid}\t\t\n")
            else:
                ch, region = ch_region
                outfh.write(f"{cid}\t{ch}\t{region}\n")

    # cleanup tmp if created
    if tmp_created is not None:
        try:
            os.remove(tmp_created)
        except Exception:
            pass

    print(f'Wrote extragood cluster->region list to {out_path}')


def find_session_paths(session_name, ephys_root=None, hist_root=None, prefer_ibl_subdir=True):
    """
    Try to locate session files given a session name.
    Returns a dict with possible keys: ks_dir, cluster_info, channel_locations, cluster_extragood

    ephys_root and hist_root default to environment variables EPHYS_ROOT and HIST_ROOT or cwd.
    The function does a shallow recursive search - stop on first reasonable match.
    """
    if ephys_root is None:
        ephys_root = os.environ.get('EPHYS_ROOT', os.getcwd())
    if hist_root is None:
        hist_root = os.environ.get('HIST_ROOT', os.getcwd())

    result = {'ks_dir': None, 'cluster_info': None, 'channel_locations': None, 'cluster_extragood': None}

    # find ks_dir: look for directory containing session_name and 'ibl_sorter_results'
    for root, dirs, files in os.walk(ephys_root):
        if session_name in root:
            # prefer explicitly named ibl_sorter_results subdir
            if 'ibl_sorter_results' in dirs and prefer_ibl_subdir:
                candidate = os.path.join(root, 'ibl_sorter_results')
                result['ks_dir'] = candidate
                break
            # otherwise, if current dir contains expected files, accept
            if any(fn in files for fn in ('channel_positions.npy', 'spike_clusters.npy', 'cluster_info.tsv')):
                result['ks_dir'] = root
                break

    # if not found yet, try a looser search for a folder name match
    if result['ks_dir'] is None:
        for root, dirs, files in os.walk(ephys_root):
            if session_name in os.path.basename(root):
                if 'ibl_sorter_results' in dirs:
                    result['ks_dir'] = os.path.join(root, 'ibl_sorter_results')
                    break

    # cluster_info and cluster_extragood expected under ks_dir
    if result['ks_dir']:
        ci = os.path.join(result['ks_dir'], 'cluster_info.tsv')
        if os.path.exists(ci):
            result['cluster_info'] = ci
        ce = os.path.join(result['ks_dir'], 'cluster_extragood.tsv')
        if os.path.exists(ce):
            result['cluster_extragood'] = ce


    # find channel_locations under hist_root
    for root, dirs, files in os.walk(hist_root):
        if session_name in root and 'channel_locations.json' in files:
            result['channel_locations'] = os.path.join(root, 'channel_locations.json')
            break
    # fallback: also check hist_root for any json named channel_locations.json containing session_name in path
    if result['channel_locations'] is None:
        for root, dirs, files in os.walk(hist_root):
            if 'channel_locations.json' in files and session_name in root:
                result['channel_locations'] = os.path.join(root, 'channel_locations.json')
                break
    

    return result
#
# batch automation helper
def automate_all_sessions(ephys_root: str, hist_root: str) -> None:
    """Walk ephys_root per-animal and per-session, compute cluster->channel and cluster->region files.

    For each session folder under each animal, finds the first subdirectory that contains
    `neurons_session.mat` and uses that as `ks_dir`. Then finds the corresponding
    histology folder under `hist_root` with the same session name and uses
    `channel_locations_simplified.tsv` (preferred) or `channel_locations.json`.
    Writes `cluster_channel_map.txt` and `cluster_channel_region_map.txt` into the `ks_dir`.
    """
    if not os.path.isdir(ephys_root):
        raise ValueError(f'ephys_root not found: {ephys_root}')
    if not os.path.isdir(hist_root):
        raise ValueError(f'hist_root not found: {hist_root}')

    animals = [d for d in sorted(os.listdir(ephys_root)) if os.path.isdir(os.path.join(ephys_root, d))]
    print(f'Found {len(animals)} animals under {ephys_root}')
    for a in animals:
        a_path = os.path.join(ephys_root, a)
        sessions = [s for s in sorted(os.listdir(a_path)) if os.path.isdir(os.path.join(a_path, s))]
        print(f'Animal {a}: {len(sessions)} sessions')
        print(sessions)
        for s in sessions:
            session_path = os.path.join(a_path, s)
            for fo in sorted(os.listdir(session_path)):
                imec_session_path = os.path.join(session_path, fo)
                if not os.path.isdir(imec_session_path):
                    continue
                print(f' Processing session {s}...')
                # find ks_dir by locating neurons_session.mat under session_path
                ks_dir = None

                for root, dirs, files in os.walk(imec_session_path):
                    if 'neurons_session.mat' in files:
                        ks_dir = root
                        break
                if ks_dir is None:
                    print(f'  Skipped {s}: neurons_session.mat not found under {imec_session_path}')
                    continue

                # find histology folder with same session name under hist_root
                hist_session_dir = None
                for root, dirs, files in os.walk(hist_root):
                    if os.path.basename(root) == fo:
                        hist_session_dir = root
                        break
                if hist_session_dir is None:
                    # fallback: look for any path containing session name with channel_locations files
                    for root, dirs, files in os.walk(hist_root):
                        if fo in root and any(fn.startswith('channel_locations') for fn in files):
                            hist_session_dir = root
                            break
                if hist_session_dir is None:
                    print(f'  Skipped {fo}: histology folder for session not found under {hist_root}')
                    continue

                # choose channel locations file
                tsv_path = os.path.join(hist_session_dir, 'channel_locations_simplified.tsv')
                json_path = os.path.join(hist_session_dir, 'channel_locations.json')
                channel_locations_path = None
                if os.path.exists(tsv_path):
                    channel_locations_path = tsv_path
                elif os.path.exists(json_path):
                    channel_locations_path = json_path
                else:
                    # pick any channel_locations* file
                    for fn in os.listdir(hist_session_dir):
                        if fn.startswith('channel_locations'):
                            channel_locations_path = os.path.join(hist_session_dir, fn)
                            break
                if channel_locations_path is None:
                    print(f'  Skipped {s}: no channel_locations file in {hist_session_dir}')
                    continue
                print(ks_dir)
                # prepare output paths and skip if final regions file already exists
                out_map = os.path.join(ks_dir, 'cluster_channel_map.txt')
                out_regions = os.path.join(ks_dir, 'cluster_channel_region_map.txt')
                if os.path.exists(out_regions):
                    print(f'  Skipping {fo}: {out_regions} already exists')
                    continue
                try:
                    print(f'  Computing cluster->channel mapping in {ks_dir} -> {out_map}')
                    main(ks_dir, out_map)
                except Exception as exc:
                    print(f'  Error computing mapping for {fo}: {exc}')
                    continue

                try:
                    out_regions = os.path.join(ks_dir, 'cluster_channel_region_map.txt')
                    print(f'  Mapping clusters to regions using {channel_locations_path} -> {out_regions}')
                    create_cluster_region_file(out_map, channel_locations_path, out_regions)
                except Exception as exc:
                    print(f'  Error mapping clusters->regions for {fo}: {exc}')
                    continue

                print(f'  Completed session {fo}: outputs in {ks_dir}')

def summarize_cell_counts(ephys_root: str, unknown_label: str = 'Unknown', output_filename: str = 'cell_counts.tsv') -> None:
    """Create EXTRAGOOD cell count summaries by brain region at four levels:
    - per imec within each session (writes session_path/cell_count_{fo}.tsv)
    - per session (writes output into each session folder)
    - per animal (writes output into each animal folder)
    - overall across animals (writes output into ephys_root)

    Filters clusters using cluster_extragood.tsv (expects columns: cluster_id, extragood='y').
    Uses existing cluster->channel->region outputs: 'cluster_channel_region_map.txt'.
    """
    if not os.path.isdir(ephys_root):
        raise ValueError(f'ephys_root not found: {ephys_root}')

    all_counts = Counter()

    animals = [d for d in sorted(os.listdir(ephys_root)) if os.path.isdir(os.path.join(ephys_root, d))]
    print(f'Found {len(animals)} animals under {ephys_root}')
    for a in animals:
        a_path = os.path.join(ephys_root, a)
        animal_counts = Counter()

        sessions = [s for s in sorted(os.listdir(a_path)) if os.path.isdir(os.path.join(a_path, s))]
        print(f'Animal {a}: {len(sessions)} sessions')
        for s in sessions:
            session_path = os.path.join(a_path, s)
            session_counts = Counter()

            # iterate possible imec subfolders within the session folder
            for fo in sorted(os.listdir(session_path)):
                imec_session_path = os.path.join(session_path, fo)
                if not os.path.isdir(imec_session_path):
                    continue

                # find ks_dir by locating neurons_session.mat
                ks_dir = None
                for root, dirs, files in os.walk(imec_session_path):
                    if 'neurons_session.mat' in files:
                        ks_dir = root
                        break
                if ks_dir is None:
                    continue

                # expected output file from automate_all_sessions
                ccrm = os.path.join(ks_dir, 'cluster_channel_region_map.txt')
                if not os.path.exists(ccrm):
                    # skip if mapping not available for this imec
                    continue

                # extragood list
                extragood_path = os.path.join(ks_dir, 'cluster_extragood.tsv')
                extragood_ids = set()
                if os.path.exists(extragood_path):
                    try:
                        with open(extragood_path, 'r', encoding='utf8') as efh:
                            sample = efh.read(4096)
                            efh.seek(0)
                            try:
                                sniffer = csv.Sniffer()
                                dialect = sniffer.sniff(sample) if sample else csv.excel_tab
                            except Exception:
                                dialect = csv.excel_tab
                            reader = csv.DictReader(efh, dialect=dialect)
                            if reader.fieldnames:
                                headers = [h.strip() for h in reader.fieldnames]
                                lower = [h.lower() for h in headers]
                                cid_key = headers[lower.index('cluster_id')] if 'cluster_id' in lower else None
                                good_key = headers[lower.index('extragood')] if 'extragood' in lower else None
                                if cid_key and good_key:
                                    for row in reader:
                                        try:
                                            if str(row.get(good_key, '')).strip().lower() == 'y':
                                                extragood_ids.add(int(str(row.get(cid_key)).strip()))
                                        except Exception:
                                            continue
                    except Exception:
                        pass

                if not extragood_ids:
                    # No extragood list found; skip this imec from counts
                    print(f"  No extragood list in {ks_dir}; skipping imec for summaries")
                    continue

                # counts for this imec only
                imec_counts = Counter()

                # read counts for this ks_dir and accumulate into session
                with open(ccrm, 'r', encoding='utf8') as fh:
                    reader = csv.DictReader(fh, delimiter='\t')
                    # if header is missing or malformed, fallback to manual parse
                    if reader.fieldnames and all(h in [hh.strip().lower() for hh in reader.fieldnames] for h in ['cluster_id','brain_region']):
                        headers = [h.strip() for h in reader.fieldnames]
                        lower = [h.lower() for h in headers]
                        cid_key = headers[lower.index('cluster_id')]
                        region_key = headers[lower.index('brain_region')]
                        for row in reader:
                            try:
                                cid = int(str(row.get(cid_key)).strip())
                            except Exception:
                                continue
                            if cid not in extragood_ids:
                                continue
                            region = (row.get(region_key) or '').strip() or unknown_label
                            imec_counts[region] += 1
                    else:
                        fh.seek(0)
                        # skip header
                        header = fh.readline()
                        for line in fh:
                            parts = line.strip().split('\t')
                            if len(parts) < 3:
                                continue
                            try:
                                cid = int(parts[0].strip())
                            except Exception:
                                continue
                            if cid not in extragood_ids:
                                continue
                            region = parts[2].strip() or unknown_label
                            imec_counts[region] += 1

                # write per-imec summary if counted
                if imec_counts:
                    per_imec_out = os.path.join(session_path, f'cell_count_{fo}.tsv')
                    with open(per_imec_out, 'w', encoding='utf8', newline='') as outfh:
                        w = csv.writer(outfh, delimiter='\t')
                        w.writerow(['brain_region', 'count'])
                        for region, count in sorted(imec_counts.items()):
                            w.writerow([region, count])
                    print(f'Wrote imec summary: {per_imec_out}')

                # add to session totals
                session_counts.update(imec_counts)

            # write per-session summary (if any counts)
            if session_counts:
                out_path = os.path.join(session_path, output_filename)
                with open(out_path, 'w', encoding='utf8', newline='') as outfh:
                    w = csv.writer(outfh, delimiter='\t')
                    w.writerow(['brain_region', 'count'])
                    for region, count in sorted(session_counts.items()):
                        w.writerow([region, count])
                print(f'Wrote session summary: {out_path}')

            # accumulate to animal/global totals
            animal_counts.update(session_counts)
            all_counts.update(session_counts)

        # write per-animal summary (if any counts)
        if animal_counts:
            out_path = os.path.join(a_path, output_filename)
            with open(out_path, 'w', encoding='utf8', newline='') as outfh:
                w = csv.writer(outfh, delimiter='\t')
                w.writerow(['brain_region', 'count'])
                for region, count in sorted(animal_counts.items()):
                    w.writerow([region, count])
            print(f'Wrote animal summary: {out_path}')

    # write overall summary (if any counts)
    if all_counts:
        out_path = os.path.join(ephys_root, output_filename)
        with open(out_path, 'w', encoding='utf8', newline='') as outfh:
            w = csv.writer(outfh, delimiter='\t')
            w.writerow(['brain_region', 'count'])
            for region, count in sorted(all_counts.items()):
                w.writerow([region, count])
        print(f'Wrote overall summary: {out_path}')

#%%

automate_all_sessions(ephys_root=r'D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\Data\TD\ephys_data\20230801_ChocolateGroup',
                    hist_root=r'D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\Data\TD\histology\20230801_ChocolateGroup')



#%% Summarize cell counts per region at session/animal/all levels
# Reuse the same roots as above; edit as needed.
summarize_cell_counts(ephys_root=r'D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\Data\TD\ephys_data\20230801_ChocolateGroup')
                        

# %%
