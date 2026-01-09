#%%
"""Simplify region map from a `channel_locations.json`.

Usage examples:
  python simplify_region_map.py --input path/to/channel_locations.json
  python simplify_region_map.py --input ... --tsv output.tsv

This prints channel->region mappings (first N), the traversal order of regions,
and the first/last channel for each region. It can also write a TSV assigning
each channel its region and metadata.
"""
from __future__ import annotations
import argparse
import json
from collections import defaultdict
import os
from typing import Dict, List, Tuple, Optional


def load_channel_locations(path: str) -> Dict[str, dict]:
    with open(path, 'r', encoding='utf8') as fh:
        return json.load(fh)


def get_channels_sorted(data: Dict[str, dict]) -> List[Tuple[int, str, Optional[str]]]:
    items = []
    for key, value in data.items():
        if not key.startswith('channel_'):
            continue
        try:
            idx = int(key.split('_', 1)[1])
        except Exception:
            continue
        region = value.get('brain_region') if isinstance(value, dict) else None
        items.append((idx, key, region))
    items.sort(key=lambda x: x[0])
    return items


def summarize_regions(channels: List[Tuple[int, str, Optional[str]]]) -> Tuple[List[str], Dict[str, Tuple[int,int]], Dict[str,int]]:
    region_order: List[str] = []
    reg_indices: Dict[str, List[int]] = defaultdict(list)
    for idx, _, region in channels:
        region_name = region if region is not None else 'UNKNOWN'
        if not region_order or region_order[-1] != region_name:
            region_order.append(region_name)
        reg_indices[region_name].append(idx)
    summary = {region: (min(idxs), max(idxs)) for region, idxs in reg_indices.items()}
    counts = {region: len(idxs) for region, idxs in reg_indices.items()}
    return region_order, summary, counts


def apply_region_assignment(channels: List[List], data: Dict[str, dict], first: int, last: int, new_region: str, new_region_id: Optional[int] = None) -> int:
    """Assign `new_region` (and optional `new_region_id`) to channels in [first,last].

    Returns number of channels updated.
    """
    idx_to_pos = {entry[0]: i for i, entry in enumerate(channels)}
    changed = 0
    for idx in range(first, last + 1):
        pos = idx_to_pos.get(idx)
        if pos is None:
            continue
        channels[pos][2] = new_region
        key = channels[pos][1]
        if key in data and isinstance(data[key], dict):
            data[key]['brain_region'] = new_region
            if new_region_id is not None:
                data[key]['brain_region_id'] = new_region_id
        changed += 1
    return changed


def interactive_edit(channels: List[List], data: Dict[str, dict], commands: Optional[List[str]] = None, source_json_path: Optional[str] = None) -> None:
    """Interactive loop to reassign brain regions for channel index ranges.

    Commands:
      set FIRST LAST REGION [REGION_ID]   e.g. set 100 120 SSp-ll1 1030
      show                                reprint summary
      save [path]                         write TSV now (path optional)
      done/quit                           exit interactive mode
    """
    print('\nEntering interactive edit mode. Type "help" for commands.')
    idx_to_pos = {entry[0]: i for i, entry in enumerate(channels)}

    def refresh_index():
        nonlocal idx_to_pos
        idx_to_pos = {entry[0]: i for i, entry in enumerate(channels)}

    def print_summary():
        region_order, summary, counts = summarize_regions([(e[0], e[1], e[2]) for e in channels])
        print('\nRegion traversal order:')
        print(' -> '.join(region_order))
        print('\nFirst and last channel per region:')
        for region in region_order:
            first, last = summary[region]
            print(f'{region}: first=channel_{first}, last=channel_{last} ({counts[region]} channels)')

    print_summary()

    cmd_iter = iter(commands) if commands is not None else None
    while True:
        try:
            if cmd_iter is not None:
                try:
                    cmd = next(cmd_iter).strip()
                    print(f'\nedit> {cmd}')
                except StopIteration:
                    print('\nNo more commands; exiting interactive mode.')
                    break
            else:
                cmd = input('\nedit> ').strip()
        except (EOFError, KeyboardInterrupt):
            print('\nExiting interactive mode.')
            break
        if not cmd:
            continue
        parts = cmd.split()
        if parts[0].lower() in ('done', 'quit', 'exit'):
            break
        if parts[0].lower() == 'help':
            print('Commands:')
            print('  set FIRST LAST REGION [REGION_ID]  - assign region to channel range')
            print('  show                               - show current region summary')
            print('  save [path]                        - write TSV (path optional)')
            print('  done | quit                        - exit interactive mode')
            continue
        if parts[0].lower() == 'show':
            print_summary()
            continue
        if parts[0].lower() == 'set' and len(parts) >= 4:
            try:
                first = int(parts[1])
                last = int(parts[2])
            except ValueError:
                print('Invalid indices. Use integers for FIRST and LAST.')
                continue
            new_region = parts[3]
            new_region_id = None
            if len(parts) >= 5:
                try:
                    new_region_id = int(parts[4])
                except ValueError:
                    new_region_id = parts[4]
            # apply changes via helper
            changed = apply_region_assignment(channels, data, first, last, new_region, new_region_id)
            refresh_index()
            print(f'Assigned region "{new_region}" to channel_{first}..channel_{last} ({changed} channels updated)')
            print_summary()
            continue
        if parts[0].lower() == 'save':
            outpath = None
            if len(parts) >= 2:
                outpath = parts[1]
            if not outpath:
                # if source JSON path provided, default to same folder
                if source_json_path:
                    outdir = os.path.dirname(os.path.abspath(source_json_path))
                    outpath = os.path.join(outdir, 'channel_locations_simplified.tsv')
                    write_tsv(channels, data, outpath)
                    print(f'Wrote TSV to: {outpath}')
                    print('Auto-saving completed — exiting interactive mode.')
                    break
                # otherwise ask the user (fallback)
                outpath = input('TSV path to write (leave blank to cancel): ').strip()
                if not outpath:
                    print('Save cancelled.')
                    continue
            write_tsv(channels, data, outpath)
            print(f'Wrote TSV to: {outpath}')
            continue
        print('Unknown command. Type "help" for usage.')


def write_tsv(channels: List[Tuple[int, str, Optional[str]]], data: Dict[str, dict], outpath: str) -> None:
    with open(outpath, 'w', encoding='utf8') as outfh:
        header = ['channel', 'brain_region', 'brain_region_id', 'axial', 'lateral', 'x', 'y', 'z']
        outfh.write('\t'.join(header) + '\n')
        for entry in channels:
            # support either tuple (idx,key,region) or mutable list [idx,key,region]
            if isinstance(entry, (list, tuple)):
                idx, key, region = entry[0], entry[1], entry[2]
            else:
                continue
            row = data.get(key, {})
            out = [f'channel_{idx}', region if region is not None else 'UNKNOWN', str(row.get('brain_region_id', '')),
                   str(row.get('axial', '')), str(row.get('lateral', '')),
                   str(row.get('x', '')), str(row.get('y', '')), str(row.get('z', ''))]
            outfh.write('\t'.join(out) + '\n')


def process_channel_locations(input_path: str, tsv_path: Optional[str] = None, print_limit: int = 0, interactive: bool = False) -> int:
    data = load_channel_locations(input_path)
    channels = get_channels_sorted(data)
    # make channels mutable lists for editing: [idx, key, region]
    channels = [[i, k, r] for (i, k, r) in channels]
    region_order, summary, counts = summarize_regions(channels)

    print(f'Total channels: {len(channels)}')
    print('\nChannel -> Region mapping (first {0} shown):'.format(print_limit))
    for idx, key, region in channels[:print_limit]:
        print(f'channel_{idx}: {region}')
    if len(channels) > print_limit:
        print('...')

    print('\nRegion traversal order:')
    print(' -> '.join(region_order))

    print('\nFirst and last channel per region:')
    for region in region_order:
        first, last = summary[region]
        print(f'{region}: first=channel_{first}, last=channel_{last} ({counts[region]} channels)')

    if interactive:
        interactive_edit(channels, data, source_json_path=input_path)

    if tsv_path:
        write_tsv(channels, data, tsv_path)
        print(f'\nTSV written to: {tsv_path}')

    return 0



# %% Run the session one by one
path = r"D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\Data\TD\histology\20230801_ChocolateGroup\5_FerreroRocher\LineFit\CB\outputForIBL\20082023_Ferrero_StrCer_S6_g0_imec0"
path_channel_locations = path + r"\channel_locations.json"
data = load_channel_locations(path_channel_locations)
channels = get_channels_sorted(data)
process_channel_locations(path_channel_locations, interactive=True, print_limit=0)
