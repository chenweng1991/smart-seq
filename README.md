# SMART-seq TSO 11 bp Tag Design

This repository records a small SMART-seq template-switching oligo (TSO) tag
design for an Illumina/Nextera-compatible workflow.

The original 11 bp tag is:

```text
ATTGCGCAATG
```

The current recommended design uses the original tag plus 3 additional 11 bp
tags. Together, the 4-tag pool has exactly 25% A/C/G/T at every position across
the 11 bp tag region.

## Recommended 4-tag Set

| name | 11 bp tag |
| --- | --- |
| original | `ATTGCGCAATG` |
| tag_01 | `GAATGCGGCCA` |
| tag_02 | `CGCATTATGAC` |
| tag_03 | `TCGCAATCTGT` |

Each new tag differs from the original at all 11 positions, and the tags are
also pairwise distance 11 from each other.

## Reproduce the Design

Run:

```bash
python scripts/design_11bp_tags.py
```

This writes:

```text
results/designed_4total_11bp_tags.csv
results/designed_4total_11bp_tags_base_diversity.csv
```

Generate the final TSO and matched FWD PCR oligos to order:

```bash
python scripts/generate_order_oligos.py
```

This writes:

```text
results/final_TSO_oligos_to_order.csv
results/final_FWD_PCR_primers_to_order.csv
```

The script uses only the Python standard library.

## Filters

The script checks:

- 11 bp tag length
- GC fraction between 45% and 64%
- maximum homopolymer length of 2
- minimum Hamming distance from the original and from selected tags
- no exact reverse-complement duplication
- simple self-complementarity screen
- simple complementarity screen against the reverse PCR primer and fixed adapter
  sequences
- per-base diversity across the original plus newly designed tags

## Oligo Context

The original RT TSO is:

```text
TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGATTGCGCAATGNNNNNNNNrGrGrG
```

For a designed tag, the corresponding matched TSO is:

```text
TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG<TAG>NNNNNNNNrGrGrG
```

and the corresponding matched forward PCR primer is:

```text
TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG<TAG>
```

Important caveat: because the 11 bp tag is at the 3' end of the forward PCR
primer, using multiple tags requires the matching forward PCR primer mixture.
The tag pool improves sequencing diversity in the tag region, but it does not
remove the primer-mixture tradeoff.

More protocol-specific notes are in
[`SMARTseq_TSO_11mer_design.md`](SMARTseq_TSO_11mer_design.md).
