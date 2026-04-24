# SMART-seq TSO 11-mer design

## Chosen 11-mer

`ATTGCGCAATG`

This matches the fixed 11 nt sequence present at the 3' end of the forward PCR
primer:

`TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGATTGCGCAATG`

## RT TSO

Order the TSO with the fixed 11-mer followed by the random 8 nt UMI/randomer
then the Smart-seq3xpress-style `WW` spacer and terminal template-switching
riboguanosines:

`/5Biosg/AGAGACAGATTGCGCAATGNNNNNNNNWWrGrGrG`

Structure:

`TSO fixed prefix - fixed 11-mer - N8 - WW - rGrGrG`

## PCR primers after RT

Forward PCR primer:

`TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGATTGCGCAA*T*G`

Forward PCR primer without phosphorothioate notation:

`TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGATTGCGCAATG`

Reverse PCR primer:

`ACGAGCATCAGCAGCATAC*G*A`

Reverse PCR primer without phosphorothioate notation:

`ACGAGCATCAGCAGCATACGA`

## Final Illumina/Nextera context

After PCR and tagmentation, the library should carry the standard Nextera-style
P5/P7 sequencing structure with sample indexes (`JJJJJJJJ`) introduced during
indexing PCR:

P5 strand:

`5'-AATGATACGGCGACCACCGAGATCTACACJJJJJJJJTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG<insert>CTGTCTCTTATACACATCTCCGAGCCCACGAGACJJJJJJJJTAGAGCATACGGCAGAAGACGAAC-3'`

P7 strand:

`3'-TTACTATGCCGCTGGTGGCTCTAGATGTGJJJJJJJJAGCAGCCGTCGCAGTCTACACATATTCTCTGTC<insert>GACAGAGAATATGTGTAGAGGCTCGGGTGCTCTGJJJJJJJJATCTCGTATGCCGTCTTCTGCTTG-5'`

## Notes

- The fixed 11-mer has balanced composition: 5 A/T and 6 G/C bases.
- It avoids long homopolymers and directly matches the forward PCR primer
  annealing sequence.
- The `WW` spacer preserves the Smart-seq3xpress improved junction between the
  UMI and the template-switching `rGrGrG`.
- Keep the `N8` bases random during synthesis if they are intended as a UMI or
  molecular randomer.

## Designing additional 11 bp tags

Use the local script:

```bash
python scripts/design_11bp_tags.py
```

Default settings:

- original tag: `ATTGCGCAATG`
- number of new tags: `10`
- tag length: `11`
- GC fraction: `0.45-0.64`
- minimum Hamming distance from the original and other selected tags: `5`
- maximum homopolymer run: `2`

The default output table is:

`results/designed_11bp_tags.csv`

To generate a different set while keeping the same filters:

```bash
python scripts/design_11bp_tags.py --seed 123
```

## Recommended 4-tag balanced set

For a small pool, use the original tag plus 3 new tags:

| name | 11 bp tag |
| --- | --- |
| original | `ATTGCGCAATG` |
| tag_01 | `GAATGCGGCCA` |
| tag_02 | `CGCATTATGAC` |
| tag_03 | `TCGCAATCTGT` |

This 4-tag pool has exactly 25% A/C/G/T at every position across the 11 bp tag.

Generate this table and its diversity report with:

```bash
python scripts/design_11bp_tags.py --count 3 \
  --output results/designed_4total_11bp_tags.csv \
  --diversity-output results/designed_4total_11bp_tags_base_diversity.csv
```

The `--count 3` is intentional here: it designs 3 new tags to use together with
the original `ATTGCGCAATG`, giving 4 total tags.
