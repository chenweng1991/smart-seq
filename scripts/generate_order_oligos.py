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
        "--output",
        type=Path,
        default=Path("results/final_oligos_to_order.csv"),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    for tag_name, tag in TAGS.items():
        rows.append(
            {
                "name": f"{tag_name}_TSO",
                "type": "TSO",
                "tag_name": tag_name,
                "tag_sequence": tag,
                "oligo_sequence": f"{TSO_PREFIX_TO_ORDER}{tag}{TSO_SUFFIX}",
            }
        )
        rows.append(
            {
                "name": f"{tag_name}_FWD_PCR",
                "type": "FWD_PCR",
                "tag_name": tag_name,
                "tag_sequence": tag,
                "oligo_sequence": f"{FWD_PREFIX}{tag}",
            }
        )

    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "name",
                "type",
                "tag_name",
                "tag_sequence",
                "oligo_sequence",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote: {args.output}")
    for row in rows:
        print(f"{row['name']}\t{row['oligo_sequence']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
