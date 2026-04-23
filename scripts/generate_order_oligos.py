#!/usr/bin/env python3
"""Generate TSO and matched FWD PCR oligos to order."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


TAGS = {
    "original": "ATTGCGCAATG",
    "tag_01": "GAATGCGGCCA",
    "tag_02": "CGCATTATGAC",
    "tag_03": "TCGCAATCTGT",
}

TSO_PREFIX_TO_ORDER = "AGAGACAG"
FWD_PREFIX = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
TSO_SUFFIX = "NNNNNNNNrGrGrG"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate final TSO/FWD oligos to order.")
    parser.add_argument(
        "--tso-output",
        type=Path,
        default=Path("results/final_TSO_oligos_to_order.csv"),
    )
    parser.add_argument(
        "--fwd-output",
        type=Path,
        default=Path("results/final_FWD_PCR_primers_to_order.csv"),
    )
    return parser.parse_args()


def write_csv(rows: list[dict[str, str]], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "name",
                "tag_name",
                "tag_sequence",
                "oligo_sequence",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()

    tso_rows = []
    fwd_rows = []
    for tag_name, tag in TAGS.items():
        tso_rows.append(
            {
                "name": f"{tag_name}_TSO",
                "tag_name": tag_name,
                "tag_sequence": tag,
                "oligo_sequence": f"{TSO_PREFIX_TO_ORDER}{tag}{TSO_SUFFIX}",
            }
        )
        fwd_rows.append(
            {
                "name": f"{tag_name}_FWD_PCR",
                "tag_name": tag_name,
                "tag_sequence": tag,
                "oligo_sequence": f"{FWD_PREFIX}{tag}",
            }
        )

    write_csv(tso_rows, args.tso_output)
    write_csv(fwd_rows, args.fwd_output)

    print(f"Wrote: {args.tso_output}")
    for row in tso_rows:
        print(f"{row['name']}\t{row['oligo_sequence']}")
    print(f"\nWrote: {args.fwd_output}")
    for row in fwd_rows:
        print(f"{row['name']}\t{row['oligo_sequence']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
