#!/usr/bin/env python3
"""Design short fixed tags for SMART-seq TSO oligos.

The script generates deterministic candidate tags and filters them for:
- exact length
- GC fraction range
- no long homopolymers
- minimum Hamming distance from the original tag and all selected tags
- no exact reverse-complement duplication among selected tags
- simple complementarity checks against the reverse PCR primer and fixed adapters
- per-base diversity across the original tag plus the selected new tags
"""

from __future__ import annotations

import argparse
import csv
import random
import sys
from dataclasses import dataclass
from pathlib import Path


BASES = "ACGT"
TSO_PREFIX = "AGAGACAG"
FWD_PREFIX = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
TSO_SUFFIX = "NNNNNNNNWWrGrGrG"
REVERSE_PCR_PRIMER = "ACGAGCATCAGCAGCATACGA"
NEXTERA_READ2_PARTIAL = "CTGTCTCTTATACACATCTCCGAGCCCACGAGAC"


@dataclass(frozen=True)
class Tag:
    name: str
    sequence: str
    gc_fraction: float
    min_distance_to_original: int
    min_distance_to_selected: int | str
    max_homopolymer: int
    max_self_complement: int
    max_primer_adapter_complement: int
    worst_primer_adapter_match: str


def gc_fraction(seq: str) -> float:
    return sum(base in "GC" for base in seq) / len(seq)


def hamming_distance(left: str, right: str) -> int:
    if len(left) != len(right):
        raise ValueError("Hamming distance requires equal-length sequences")
    return sum(a != b for a, b in zip(left, right))


def max_homopolymer(seq: str) -> int:
    longest = 1
    current = 1
    for left, right in zip(seq, seq[1:]):
        current = current + 1 if left == right else 1
        longest = max(longest, current)
    return longest


def reverse_complement(seq: str) -> str:
    table = str.maketrans("ACGT", "TGCA")
    return seq.translate(table)[::-1]


def longest_common_substring(left: str, right: str) -> int:
    best = 0
    lengths = [[0] * (len(right) + 1) for _ in range(len(left) + 1)]
    for i, left_base in enumerate(left, start=1):
        for j, right_base in enumerate(right, start=1):
            if left_base == right_base:
                lengths[i][j] = lengths[i - 1][j - 1] + 1
                best = max(best, lengths[i][j])
    return best


def longest_complement_run(left: str, right: str) -> int:
    return longest_common_substring(left, reverse_complement(right))


def worst_primer_adapter_complement(seq: str) -> tuple[int, str]:
    checks = {
        "reverse_pcr_primer": REVERSE_PCR_PRIMER,
        "forward_pcr_prefix": FWD_PREFIX,
        "tso_prefix": TSO_PREFIX,
        "nextera_read2_partial": NEXTERA_READ2_PARTIAL,
    }
    scored = [
        (longest_complement_run(seq, reference), name)
        for name, reference in checks.items()
    ]
    return max(scored)


def has_forbidden_pattern(seq: str) -> bool:
    # Conservative short-tag filters: avoid enzyme-like or low-complexity motifs.
    forbidden = {
        "GGGG",
        "CCCC",
        "AAAA",
        "TTTT",
        "GATC",
        "GAATTC",
        "AAGCTT",
        "GGATCC",
        "CTGCAG",
        "GCGGCCGC",
    }
    return any(pattern in seq for pattern in forbidden)


def passes_basic_filters(
    seq: str,
    *,
    length: int,
    gc_min: float,
    gc_max: float,
    max_run: int,
    max_self_complement: int,
    max_primer_adapter_complement: int,
) -> bool:
    if len(seq) != length:
        return False
    if any(base not in BASES for base in seq):
        return False
    if not (gc_min <= gc_fraction(seq) <= gc_max):
        return False
    if max_homopolymer(seq) > max_run:
        return False
    if has_forbidden_pattern(seq):
        return False
    if longest_complement_run(seq, seq) > max_self_complement:
        return False
    if worst_primer_adapter_complement(seq)[0] > max_primer_adapter_complement:
        return False
    return True


def diversity_score(seqs: list[str]) -> tuple[float, float]:
    """Return a score that favors balanced base fractions at every position."""
    worst_min_fraction = 1.0
    squared_error = 0.0
    for position in range(len(seqs[0])):
        counts = {base: 0 for base in BASES}
        for seq in seqs:
            counts[seq[position]] += 1
        fractions = [counts[base] / len(seqs) for base in BASES]
        worst_min_fraction = min(worst_min_fraction, min(fractions))
        squared_error += sum((fraction - 0.25) ** 2 for fraction in fractions)
    return worst_min_fraction, -squared_error


def design_tags(
    *,
    original: str,
    count: int,
    length: int,
    gc_min: float,
    gc_max: float,
    min_distance: int,
    max_run: int,
    seed: int,
    max_attempts: int,
    max_self_complement: int,
    max_primer_adapter_complement: int,
    pool_original_for_diversity: bool,
) -> list[Tag]:
    original = original.upper()
    if not passes_basic_filters(
        original,
        length=length,
        gc_min=0.0,
        gc_max=1.0,
        max_run=max_run,
        max_self_complement=length,
        max_primer_adapter_complement=length,
    ):
        raise ValueError(f"Original tag failed basic validation: {original}")

    rng = random.Random(seed)
    if count == 3 and pool_original_for_diversity:
        selected = design_perfect_four_total_tags(
            original=original,
            length=length,
            gc_min=gc_min,
            gc_max=gc_max,
            min_distance=min_distance,
            max_run=max_run,
            max_self_complement=max_self_complement,
            max_primer_adapter_complement=max_primer_adapter_complement,
            rng=rng,
            max_attempts=max_attempts,
        )
        return summarize_tags(original, selected)

    blocked = {original, reverse_complement(original)}
    candidates: list[str] = []

    for _ in range(max_attempts):
        seq = "".join(rng.choice(BASES) for _ in range(length))
        if seq in blocked or reverse_complement(seq) in blocked:
            continue
        if not passes_basic_filters(
            seq,
            length=length,
            gc_min=gc_min,
            gc_max=gc_max,
            max_run=max_run,
            max_self_complement=max_self_complement,
            max_primer_adapter_complement=max_primer_adapter_complement,
        ):
            continue
        if hamming_distance(seq, original) < min_distance:
            continue
        if seq in candidates:
            continue
        candidates.append(seq)
        blocked.add(seq)
        blocked.add(reverse_complement(seq))
        if len(candidates) >= max(count * 500, 5000):
            break

    selected: list[str] = []
    diversity_seed = [original] if pool_original_for_diversity else []
    while len(selected) < count:
        eligible = [
            seq
            for seq in candidates
            if seq not in selected
            and all(hamming_distance(seq, prior) >= min_distance for prior in selected)
        ]
        if not eligible:
            break
        selected.append(
            max(
                eligible,
                key=lambda seq: (
                    diversity_score(diversity_seed + selected + [seq]),
                    min([hamming_distance(seq, original)] + [hamming_distance(seq, p) for p in selected]),
                    seq,
                ),
            )
        )

    if len(selected) < count:
        raise RuntimeError(
            f"Only found {len(selected)} tags from {len(candidates)} candidates. "
            "Try lowering --min-distance or widening --gc-min/--gc-max."
        )

    return summarize_tags(original, selected)


def design_perfect_four_total_tags(
    *,
    original: str,
    length: int,
    gc_min: float,
    gc_max: float,
    min_distance: int,
    max_run: int,
    max_self_complement: int,
    max_primer_adapter_complement: int,
    rng: random.Random,
    max_attempts: int,
) -> list[str]:
    """Design 3 tags that perfectly balance the original into a 4-tag pool."""
    for _ in range(max_attempts):
        selected = ["", "", ""]
        for original_base in original:
            remaining = [base for base in BASES if base != original_base]
            rng.shuffle(remaining)
            for index, base in enumerate(remaining):
                selected[index] += base

        if not all(
            passes_basic_filters(
                seq,
                length=length,
                gc_min=gc_min,
                gc_max=gc_max,
                max_run=max_run,
                max_self_complement=max_self_complement,
                max_primer_adapter_complement=max_primer_adapter_complement,
            )
            for seq in selected
        ):
            continue
        if any(hamming_distance(seq, original) < min_distance for seq in selected):
            continue
        if any(
            hamming_distance(left, right) < min_distance
            for index, left in enumerate(selected)
            for right in selected[:index]
        ):
            continue
        if len(set(selected + [reverse_complement(seq) for seq in selected])) != 6:
            continue
        return selected

    raise RuntimeError(
        "No perfectly balanced 4-total-tag set found. "
        "Try increasing --max-attempts or relaxing complementarity filters."
    )


def summarize_tags(original: str, selected: list[str]) -> list[Tag]:
    tags: list[Tag] = []
    for index, seq in enumerate(selected, start=1):
        prior = selected[: index - 1]
        primer_adapter_score, primer_adapter_name = worst_primer_adapter_complement(seq)
        tags.append(
            Tag(
                name=f"tag_{index:02d}",
                sequence=seq,
                gc_fraction=gc_fraction(seq),
                min_distance_to_original=hamming_distance(seq, original),
                min_distance_to_selected=(
                    "NA" if not prior else min(hamming_distance(seq, p) for p in prior)
                ),
                max_homopolymer=max_homopolymer(seq),
                max_self_complement=longest_complement_run(seq, seq),
                max_primer_adapter_complement=primer_adapter_score,
                worst_primer_adapter_match=primer_adapter_name,
            )
        )
    return tags


def base_diversity_rows(seqs: list[str]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for position in range(len(seqs[0])):
        counts = {base: 0 for base in BASES}
        for seq in seqs:
            counts[seq[position]] += 1
        row = {"position": str(position + 1)}
        for base in BASES:
            row[f"{base}_count"] = str(counts[base])
            row[f"{base}_fraction"] = f"{counts[base] / len(seqs):.3f}"
        row["min_fraction"] = f"{min(counts.values()) / len(seqs):.3f}"
        row["max_fraction"] = f"{max(counts.values()) / len(seqs):.3f}"
        rows.append(row)
    return rows


def write_csv(tags: list[Tag], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "name",
                "sequence",
                "gc_fraction",
                "min_distance_to_original",
                "min_distance_to_selected",
                "max_homopolymer",
                "max_self_complement",
                "max_primer_adapter_complement",
                "worst_primer_adapter_match",
                "tso_oligo",
                "forward_pcr_primer",
            ],
        )
        writer.writeheader()
        for tag in tags:
            writer.writerow(
                {
                    "name": tag.name,
                    "sequence": tag.sequence,
                    "gc_fraction": f"{tag.gc_fraction:.3f}",
                    "min_distance_to_original": tag.min_distance_to_original,
                    "min_distance_to_selected": tag.min_distance_to_selected,
                    "max_homopolymer": tag.max_homopolymer,
                    "max_self_complement": tag.max_self_complement,
                    "max_primer_adapter_complement": tag.max_primer_adapter_complement,
                    "worst_primer_adapter_match": tag.worst_primer_adapter_match,
                    "tso_oligo": f"{TSO_PREFIX}{tag.sequence}{TSO_SUFFIX}",
                    "forward_pcr_primer": f"{FWD_PREFIX}{tag.sequence}",
                }
            )


def write_diversity_csv(rows: list[dict[str, str]], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Design balanced 11 bp SMART-seq TSO fixed tags."
    )
    parser.add_argument("--original", default="ATTGCGCAATG")
    parser.add_argument(
        "--count",
        type=int,
        default=3,
        help=(
            "Number of new tags to design. The default designs 3 new tags "
            "for use with the original tag, giving 4 total tags."
        ),
    )
    parser.add_argument("--length", type=int, default=11)
    parser.add_argument("--gc-min", type=float, default=0.45)
    parser.add_argument("--gc-max", type=float, default=0.64)
    parser.add_argument("--min-distance", type=int, default=5)
    parser.add_argument("--max-homopolymer", type=int, default=2)
    parser.add_argument("--max-self-complement", type=int, default=4)
    parser.add_argument("--max-primer-adapter-complement", type=int, default=5)
    parser.add_argument("--seed", type=int, default=20260423)
    parser.add_argument("--max-attempts", type=int, default=1_000_000)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/designed_4total_11bp_tags.csv"),
    )
    parser.add_argument(
        "--diversity-output",
        type=Path,
        default=Path("results/designed_4total_11bp_tags_base_diversity.csv"),
    )
    parser.add_argument(
        "--exclude-original-from-diversity",
        action="store_true",
        help="Optimize/report diversity using only newly designed tags.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    tags = design_tags(
        original=args.original,
        count=args.count,
        length=args.length,
        gc_min=args.gc_min,
        gc_max=args.gc_max,
        min_distance=args.min_distance,
        max_run=args.max_homopolymer,
        seed=args.seed,
        max_attempts=args.max_attempts,
        max_self_complement=args.max_self_complement,
        max_primer_adapter_complement=args.max_primer_adapter_complement,
        pool_original_for_diversity=not args.exclude_original_from_diversity,
    )
    write_csv(tags, args.output)
    diversity_seqs = [] if args.exclude_original_from_diversity else [args.original.upper()]
    diversity_seqs.extend(tag.sequence for tag in tags)
    diversity_rows = base_diversity_rows(diversity_seqs)
    write_diversity_csv(diversity_rows, args.diversity_output)

    print(
        "name\tsequence\tgc_fraction\tmin_distance_to_original\t"
        "min_distance_to_selected\tmax_homopolymer\tmax_self_complement\t"
        "max_primer_adapter_complement\tworst_primer_adapter_match"
    )
    for tag in tags:
        print(
            f"{tag.name}\t{tag.sequence}\t{tag.gc_fraction:.3f}\t"
            f"{tag.min_distance_to_original}\t{tag.min_distance_to_selected}\t"
            f"{tag.max_homopolymer}\t{tag.max_self_complement}\t"
            f"{tag.max_primer_adapter_complement}\t{tag.worst_primer_adapter_match}"
        )
    print("\nPer-base diversity across reported pool:")
    print("position\tA\tC\tG\tT\tmin_fraction\tmax_fraction")
    for row in diversity_rows:
        print(
            f"{row['position']}\t{row['A_fraction']}\t{row['C_fraction']}\t"
            f"{row['G_fraction']}\t{row['T_fraction']}\t"
            f"{row['min_fraction']}\t{row['max_fraction']}"
        )
    print(f"\nWrote: {args.output}")
    print(f"Wrote: {args.diversity_output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
