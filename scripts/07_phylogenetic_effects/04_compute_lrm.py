"""Matching-split distance for binary and unresolved unrooted consensus trees.

"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from cogent3 import load_tree
from scipy.optimize import linear_sum_assignment

DEFAULT_DIR = Path('Gene_Trees_50')
SUFFIX = '.filter2.fasta.contree'


def read_splits(path, expected_names=None):
    tree = load_tree(str(path)).unrooted()
    names = sorted(tree.get_tip_names())
    if len(names) != len(set(names)):
        raise ValueError(f'Duplicate tip names: {path}')
    if expected_names is not None and names != expected_names:
        raise ValueError(f'Taxa differ from original: {path}')
    if len(names) < 4:
        raise ValueError(f'Fewer than four taxa: {path}')
    indices = {name: i for i, name in enumerate(names)}
    full = (1 << len(names)) - 1
    splits = set()
    for subset in tree.subsets():
        mask = sum(1 << indices[name] for name in subset)
        if 1 < mask.bit_count() < len(names) - 1:
            splits.add(min(mask, full ^ mask))
    return names, sorted(splits)


def matching_split_distance(a, b, n_tips):
    size = max(len(a), len(b))
    if size == 0:
        return 0
    costs = np.zeros((size, size), dtype=np.int64)
    for i, split_a in enumerate(a):
        for j, split_b in enumerate(b):
            difference = (split_a ^ split_b).bit_count()
            costs[i, j] = min(difference, n_tips - difference)
    if len(a) < size:
        costs[len(a):, :len(b)] = [
            min(s.bit_count(), n_tips - s.bit_count()) for s in b
        ]
    if len(b) < size:
        costs[:len(a), len(b):] = np.array([
            min(s.bit_count(), n_tips - s.bit_count()) for s in a
        ])[:, None]
    rows, cols = linear_sum_assignment(costs)
    return int(costs[rows, cols].sum())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tree-dir', type=Path, default=DEFAULT_DIR)
    parser.add_argument('--out-dir', type=Path, default=Path('.'))
    args = parser.parse_args()
    originals = [p for p in sorted(args.tree_dir.glob('*' + SUFFIX))
                 if '.rand' not in p.name and '.nors' not in p.name]
    if not originals:
        parser.error(f'No original trees found in {args.tree_dir}')
    records, missing = [], []
    for original in originals:
        base = original.name.removesuffix(SUFFIX)
        og = base.split('.')[0]
        names, reference = read_splits(original)
        targets = [('nors', None)] + [('rand', i) for i in range(1, 101)]
        for kind, replicate in targets:
            tag = kind if replicate is None else f'rand{replicate}'
            path = args.tree_dir / f'{base}.filter2.{tag}.fasta.contree'
            if not path.exists():
                missing.append({'file': str(path), 'reason': 'missing'})
                continue
            _, comparison = read_splits(path, names)
            distance = matching_split_distance(reference, comparison, len(names))
            records.append([og, base, kind, replicate, distance, len(names),
                            len(reference), len(comparison)])
        print(f'Finished {base}', flush=True)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(records, columns=[
        'Orthogroup', 'Family', 'type', 'replicate', 'MS', 'n_tips',
        'original_splits', 'comparison_splits'])
    df['replicate'] = df['replicate'].astype('Int64')
    df.to_csv(args.out_dir / 'matching_split_results.tsv', sep='\t', index=False, na_rep='NA')
    pd.DataFrame(missing, columns=['file', 'reason']).to_csv(
        args.out_dir / 'missing_trees.tsv', sep='\t', index=False)
    nors = df.loc[df['type'] == 'nors', ['Family', 'MS']].rename(columns={'MS': 'MS_nors'})
    delta = df.loc[df['type'] == 'rand'].merge(nors, on='Family', validate='many_to_one')
    delta['delta_MS'] = delta['MS_nors'] - delta['MS']
    delta.to_csv(args.out_dir / 'delta_matching_split_results.tsv', sep='\t', index=False)
    print(f'Wrote {len(df)} distances and {len(delta)} contrasts; {len(missing)} missing files.')


if __name__ == '__main__':
    main()
