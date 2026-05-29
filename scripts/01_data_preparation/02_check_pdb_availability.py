#!/usr/bin/env python3
from __future__ import annotations

import argparse
import time
from pathlib import Path
from typing import Dict, List, Set, Optional

import requests


"""
Identify orthogroups containing proteins with experimentally
resolved structures in the Protein Data Bank (PDB).

The script queries the PDBe API using UniProt identifiers
mapped from OMA orthogroups and records whether at least one
experimental structure exists for each orthogroup.

Outputs:
  - groups_with_pdb.txt
  - groups_without_pdb.txt
"""

DEFAULT_METHODS = {"X-ray diffraction", "Electron microscopy", "Solution NMR"}


def load_first_col_as_set(path: Path) -> Set[str]:
    if not path.exists():
        return set()
    out: Set[str] = set()
    with path.open("r") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            out.add(s.split("\t")[0])
    return out


def load_oma_groups(path: Path, fingerprints: Set[str]) -> Dict[str, List[str]]:
    fp_to_oma: Dict[str, List[str]] = {}
    with path.open("r") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            _group_num, fp, *oma_entries = parts
            if fp in fingerprints:
                fp_to_oma.setdefault(fp, []).extend(oma_entries)
    return fp_to_oma


def load_oma_uniprot(path: Path) -> Dict[str, str]:
    m: Dict[str, str] = {}
    with path.open("r") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) == 2:
                m[parts[0]] = parts[1]
    return m


def check_pdb_for_uniprot(
    uniprot_id: str,
    methods: Set[str],
    session: requests.Session,
    timeout: int = 20,
    sleep_s: float = 0.0,
) -> Optional[str]:
    """
    Returns a short string like '1abc (X-ray diffraction)' if any experimental structure exists,
    otherwise None.
    """
    url = f"https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/{uniprot_id}"
    try:
        r = session.get(url, timeout=timeout)
        if sleep_s:
            time.sleep(sleep_s)
        if r.status_code != 200:
            return None
        data = r.json()
        if uniprot_id not in data:
            return None
        for entry in data[uniprot_id]:
            if entry.get("experimental_method") in methods:
                return f"{entry.get('pdb_id')} ({entry.get('experimental_method')})"
    except requests.exceptions.RequestException as e:
        print(f"Warning: request failed for {uniprot_id}: {e}")
    except ValueError:
        print(f"Warning: JSON decode failed for {uniprot_id}")
    return None


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Screen orthogroups: mark groups that have at least one experimental PDB structure (stop at first hit)."
    )
    ap.add_argument("--fingerprints", required=True, type=Path, help="Fingerprints list (one per line, or first column used).")
    ap.add_argument("--oma-groups", required=True, type=Path, help="OMA groups table (tab-separated).")
    ap.add_argument("--oma-uniprot", required=True, type=Path, help="OMA entry -> UniProt mapping (tab-separated).")
    ap.add_argument("--with-pdb", default=Path("groups_with_pdb.txt"), type=Path, help="Output file to append hits.")
    ap.add_argument("--without-pdb", default=Path("groups_without_pdb.txt"), type=Path, help="Output file to append misses.")
    ap.add_argument(
        "--methods",
        default=";".join(sorted(DEFAULT_METHODS)),
        help="Semicolon-separated experimental methods to accept.",
    )
    ap.add_argument("--sleep", type=float, default=0.0, help="Optional sleep (seconds) between PDBe API requests.")
    args = ap.parse_args()

    methods = set(m.strip() for m in args.methods.split(";") if m.strip())

    # Resume support: don't redo groups already written
    checked = load_first_col_as_set(args.with_pdb) | load_first_col_as_set(args.without_pdb)

    all_fps = load_first_col_as_set(args.fingerprints)
    remaining = all_fps - checked
    print(f"Total fingerprints: {len(all_fps)} | Already processed: {len(checked)} | Remaining: {len(remaining)}")

    if not remaining:
        print("No new fingerprints to process. Exiting.")
        return 0

    fp_to_oma = load_oma_groups(args.oma_groups, remaining)
    oma_to_uniprot = load_oma_uniprot(args.oma_uniprot)

    # Make sure parent dirs exist
    args.with_pdb.parent.mkdir(parents=True, exist_ok=True)
    args.without_pdb.parent.mkdir(parents=True, exist_ok=True)

    session = requests.Session()

    with args.with_pdb.open("a") as with_pdb, args.without_pdb.open("a") as without_pdb:
        for fp, oma_ids in fp_to_oma.items():
            found = False
            for oma_id in oma_ids:
                uniprot = oma_to_uniprot.get(oma_id)
                if not uniprot:
                    continue
                hit = check_pdb_for_uniprot(uniprot, methods, session=session, sleep_s=args.sleep)
                if hit:
                    with_pdb.write(f"{fp}\t{uniprot}\t{hit}\n")
                    found = True
                    break
            if not found:
                without_pdb.write(f"{fp}\n")

    print(f"DONE. Appended results to:\n  - {args.with_pdb}\n  - {args.without_pdb}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

