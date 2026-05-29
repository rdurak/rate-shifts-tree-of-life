#!/usr/bin/env python3

"""
Perform likelihood-ratio tests (LRTs) comparing
RASER and null evolutionary models.

For each orthogroup:
  - extracts log-likelihood values
  - computes the likelihood-ratio statistic
  - calculates chi-square p-values
  - reports statistical significance

Inputs:
  - RASER result files
  - null model result files

Outputs:
  raser_lrt_results.tsv

Example:
  python 03_run_lrt.py \
    --results-dir raser_output \
    --output raser_lrt_results.tsv
"""

from __future__ import annotations

import argparse
import glob
import re
from pathlib import Path

from scipy.stats import chi2


def extract_log_likelihood(
    file_path: Path,
    pattern: str
) -> float | None:

    with file_path.open("r") as f:

        content = f.read()

    match = re.search(
        pattern,
        content
    )

    if match is None:
        return None

    return float(
        match.group(1)
    )


def main() -> None:

    ap = argparse.ArgumentParser(
        description=(
            "Run likelihood-ratio tests "
            "between RASER and null models."
        )
    )

    ap.add_argument(
        "--results-dir",
        default="raser_output",
        help="Directory containing RASER output files"
    )

    ap.add_argument(
        "--output",
        default="raser_lrt_results.tsv",
        help="Output results table"
    )

    ap.add_argument(
        "--df",
        type=int,
        default=3,
        help="Degrees of freedom for chi-square test"
    )

    args = ap.parse_args()

    results_dir = Path(
        args.results_dir
    )

    output_file = Path(
        args.output
    )

    # --------------------------------------------------------
    # checks
    # --------------------------------------------------------

    if not results_dir.exists():

        raise FileNotFoundError(
            f"Results directory not found: "
            f"{results_dir}"
        )

    # --------------------------------------------------------
    # collect result files
    # --------------------------------------------------------

    null_files = glob.glob(
        str(results_dir / "n_*_results_file")
    )

    raser_files = glob.glob(
        str(results_dir / "r_*_outresfile")
    )

    # --------------------------------------------------------
    # match orthogroups
    # --------------------------------------------------------

    null_dict = {

        Path(f).name
        .replace("n_", "")
        .replace("_results_file", ""): f

        for f in null_files
    }

    raser_dict = {

        Path(f).name
        .replace("r_", "")
        .replace("_outresfile", ""): f

        for f in raser_files
    }

    common_keys = sorted(
        set(null_dict.keys())
        &
        set(raser_dict.keys())
    )

    print(
        f"Found {len(common_keys)} "
        f"matching orthogroups."
    )

    # --------------------------------------------------------
    # regex pattern
    # --------------------------------------------------------

    ll_regex = (
        r"Log-likelihood:\s*"
        r"(-?\d+(?:\.\d+)?)"
    )

    results = []

    skipped = 0

    # --------------------------------------------------------
    # process orthogroups
    # --------------------------------------------------------

    for key in common_keys:

        null_file = Path(
            null_dict[key]
        )

        raser_file = Path(
            raser_dict[key]
        )

        null_ll = extract_log_likelihood(
            null_file,
            ll_regex
        )

        if null_ll is None:

            print(
                f"Warning: missing log-likelihood "
                f"in {null_file}"
            )

            skipped += 1

            continue

        raser_ll = extract_log_likelihood(
            raser_file,
            ll_regex
        )

        if raser_ll is None:

            print(
                f"Warning: missing log-likelihood "
                f"in {raser_file}"
            )

            skipped += 1

            continue

        # ----------------------------------------------------
        # likelihood-ratio test
        # ----------------------------------------------------

        D = 2 * (
            raser_ll - null_ll
        )

        p_value = (
            1 - chi2.cdf(
                D,
                df=args.df
            )
        )

        significant = (
            "Yes"
            if p_value < 0.05
            else "No"
        )

        results.append(

            (
                key,
                null_ll,
                raser_ll,
                D,
                p_value,
                significant
            )
        )

    # --------------------------------------------------------
    # write output
    # --------------------------------------------------------

    with output_file.open("w") as out:

        out.write(
            "Orthogroup\t"
            "Null_LL\t"
            "RASER_LL\t"
            "D_statistic\t"
            "p_value\t"
            "Significant\n"
        )

        for row in results:

            (
                key,
                null_ll,
                raser_ll,
                D,
                p_value,
                significant
            ) = row

            out.write(

                f"{key}\t"
                f"{null_ll}\t"
                f"{raser_ll}\t"
                f"{D:.4f}\t"
                f"{p_value:.2e}\t"
                f"{significant}\n"
            )

    # --------------------------------------------------------
    # summary
    # --------------------------------------------------------

    print("\nFinished.")

    print(
        f"Results written: "
        f"{output_file}"
    )

    print(
        f"Orthogroups analyzed: "
        f"{len(results)}"
    )

    print(
        f"Orthogroups skipped: "
        f"{skipped}"
    )


if __name__ == "__main__":
    main()
