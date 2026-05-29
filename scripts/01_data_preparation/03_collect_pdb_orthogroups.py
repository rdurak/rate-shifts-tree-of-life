#!/usr/bin/env python3
from __future__ import annotations

import argparse
import time
from pathlib import Path
from typing import Dict, List, Set

import requests

"""
Collect all experimentally resolved PDB structures associated
with orthogroups previously identified as containing structural
evidence.

For each orthogroup:
  - retrieves UniProt mappings
  - queries the PDBe API
  - collects experimental PDB identifiers

Outputs:
  all_pdb_groups_with_pdb.txt
"""

DEFAULT_METHODS = {"X-ray diffraction", "Electron microscopy", "Solution NMR"}


def load_groups_with_pdb(path: Path) -> Set[str]:
    fps: Set[str] = set()
    with path.open("r") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            fps.add(s.split("\t")[0])
    return fps


def load_oma_groups(path: Path, target_fps: Set[str]) -> Dict[str, List[str]]:
    fp_to_oma: Dict[str, List[str]] = {}
    with path.open("r") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            _group_id, fp, *oma_entries = parts
            if fp in target_fps:
                fp_to_oma[fp] = oma_entries
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


def get_experimental_pdb_ids(
    uniprot_id: str,
    methods: Set[str],
    session: requests.Session,
    timeout: int = 20,
    sleep_s: float = 0.0,
) -> List[str]:
    url = f"https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/{uniprot_id}"
    try:
        r = session.get(url, timeout=timeout)
        if sleep_s:
            time.sleep(sleep_s)
        if r.status_code != 200:
            return []
        data = r.json()
        if uniprot_id not in data:
            return []
        out = []
        for entry in data[uniprot_id]:
            if entry.get("experimental_method") in methods and entry.get("pdb_id"):
                out.append(entry["pdb_id"])
        return out
    except requests.exceptions.RequestException as e:
        print(f"Warning: request failed for {uniprot_id}: {e}")
    except ValueError:
        print(f"Warning: JSON decode failed for {uniprot_id}")
    return []


def main() -> int:
    ap = argparse.ArgumentParser(
        description="For groups known to have experimental structure evidence, list all experimental PDB IDs per orthogroup."
    )
    ap.add_argument("--groups-with-pdb", required=True, type=Path, help="Output from step 02_01 (groups_with_pdb.txt).")
    ap.add_argument("--oma-groups", required=True, type=Path, help="OMA groups table (tab-separated).")
    ap.add_argument("--oma-uniprot", required=True, type=Path, help="OMA entry -> UniProt mapping (tab-separated).")
    ap.add_argument("--out", default=Path("all_pdb_groups_with_pdb.txt"), type=Path, help="Output table: fingerprint + PDB IDs.")
    ap.add_argument(
        "--methods",
        default=";".join(sorted(DEFAULT_METHODS)),
        help="Semicolon-separated experimental methods to accept.",
    )
    ap.add_argument("--sleep", type=float, default=0.0, help="Optional sleep (seconds) between PDBe API requests.")
    args = ap.parse_args()

    methods = set(m.strip() for m in args.methods.split(";") if m.strip())

    fps = load_groups_with_pdb(args.groups_with_pdb)
    fp_to_oma = load_oma_groups(args.oma_groups, fps)
    oma_to_uniprot = load_oma_uniprot(args.oma_uniprot)

    print(f"Processing {len(fp_to_oma)} orthogroups (positives) ...")

    args.out.parent.mkdir(parents=True, exist_ok=True)

    session = requests.Session()

    with args.out.open("w") as out:
        for fp, oma_ids in fp_to_oma.items():
            pdb_ids: Set[str] = set()
            for oma_id in oma_ids:
                uniprot = oma_to_uniprot.get(oma_id)
                if not uniprot:
                    continue
                pdb_ids.update(get_experimental_pdb_ids(uniprot, methods, session=session, sleep_s=args.sleep))
            if pdb_ids:
                out.write(fp + "\t" + "\t".join(sorted(pdb_ids)) + "\n")

    print(f"DONE. Wrote: {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

