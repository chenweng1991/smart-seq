## Goal

Solve the SMART-seq3 low sequencing-diversity issue by redesigning the tag / TSO region while preserving the Sankaran-lab SMART-seq3 workflow for mitochondrial lineage-tracing experiments.

## Source Protocols

- SMART-seq3 protocol source:[https://benchling.com/s/prt-EuK3tH2lKjp2f0v7LMcd?m=slm-hDryr4LasF8dspuFRDEJ](https://benchling.com/s/prt-EuK3tH2lKjp2f0v7LMcd?m=slm-hDryr4LasF8dspuFRDEJ)
- https://dx.doi.org/10.17504/protocols.io.bcq4ivyw
- Design/code repository: https://github.com/chenweng1991/smart-seq
- Potential stagger design reference: [[70 Resources/Sequencing Systems/Staggered Primer Design]]
- Sankaran lab TruSeq adapter sequences: [[Assets/others/Sankaran_lab_truseq.docx]]

## Key Problem

Current TSO:

```text
Smartseq3_N8_TSO | IDT | RNase-Free HPLC | 100 uM | /5Biosg/AGAGACAGATTGCGCAATGNNNNNNNNrGrGrG
```

The fixed `AGAGACAGATTGCGCAATG` segment creates low sequence diversity early in the read. The design goal is to improve sequencing diversity while preserving the required TSO / UMI function and keeping PCR simple.

This redesign now also adopts the Smart-seq3xpress-style `WW` spacer between the
`N8` UMI and terminal `rGrGrG`, so the tag replacement is combined with the
improved TSO junction rather than the original `N8-rGrGrG` junction.

Relevant sequence relationship:

```text
                         AGAGACAGATTGCGCAATGNNNNNNNNWWrGrGrG        (TSO)
TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGATTGCGCAATG               (FWD PCR primer)
TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG        (Nextera Tn5 adapter Read 1 sequencing primer)
```

This explains that `ATTGCGCAATG` is the tag used in this design, after the shared `AGAGACAG` sequence, 
from which the illumina Read1 starts

## Design Option 1: Balanced Tag Replacement

Original tag and proposed replacement tags:

```text
ATTGCGCAATG  original
GAATGCGGCCA  replacement 1
CGCATTATGAC  replacement 2
TCGCAATCTGT  replacement 3
```

Together, these 4 tags are designed to satisfy:

1. Perfect per-base diversity across the 11 bp tag. At every tag position, the 4 tags contain 25% A, 25% C, 25% G, and 25% T, so every cycle across the tag region is balanced.
2. Maximum Hamming distance from the original. Each replacement tag differs from `ATTGCGCAATG` at all 11 positions.
3. High distance between replacement tags. Each replacement pair is also Hamming distance 11, the maximum possible for 11 bp tags.
4. Basic oligo-quality filters. GC fraction 45-64%, max homopolymer <=2, self-complementarity <=4 bp run, adapter / reverse-primer complementarity <=5 bp run, and forbidden motifs avoided.

Pairwise Hamming distances:

```text
GAATGCGGCCA vs ATTGCGCAATG  distance 11
CGCATTATGAC vs ATTGCGCAATG  distance 11
TCGCAATCTGT vs ATTGCGCAATG  distance 11

GAATGCGGCCA vs CGCATTATGAC  distance 11
GAATGCGGCCA vs TCGCAATCTGT  distance 11
CGCATTATGAC vs TCGCAATCTGT  distance 11
```

Important caveat: this design still requires 4 matched forward PCR primers, 
because the 11 bp tag is at the 3' end of the FWD PCR primer. 

The sequence diversity is excellent, but in practice we DO NOT mix them for PCR, instead we use unique 
pairs for different samples. 

eg, some cells use original_TSO --> original_FWD_PCR
    some  other cells use tag_01_TSO  --> tag_01_FWD_PCR

## Exact Oligo Changes For Tag Replacement

TSO oligos:

```text
original_TSO
AGAGACAGATTGCGCAATGNNNNNNNNWWrGrGrG

tag_01_TSO
AGAGACAGGAATGCGGCCANNNNNNNNWWrGrGrG

tag_02_TSO
AGAGACAGCGCATTATGACNNNNNNNNWWrGrGrG

tag_03_TSO
AGAGACAGTCGCAATCTGTNNNNNNNNWWrGrGrG
```

FWD PCR primers:

```text
original_FWD_PCR
TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGATTGCGCAATG

tag_01_FWD_PCR
TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGGAATGCGGCCA

tag_02_FWD_PCR
TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCGCATTATGAC

tag_03_FWD_PCR
TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGTCGCAATCTGT
```

The matched forward PCR primers stay the same length and composition except for
the 11 bp tag itself, because the added `WW` spacer sits downstream of the PCR
primer annealing site and only changes the RT TSO.

## Exact Filters Added In The Script

- tag length fixed at 11 bp
- GC fraction between 45% and 64%
- max homopolymer length <=2
- minimum Hamming distance >=5 from the original and from selected replacement tags
- exact reverse-complement duplicates excluded
- self-complementarity screen by longest complement run <=4 bp
- complementarity screen against the reverse PCR primer, TSO prefix / adapter sequence, and Nextera Read2 partial sequence by longest complement run <=5 bp
- forbidden motifs excluded: `AAAA`, `TTTT`, `CCCC`, `GGGG`, `GATC`, `GAATTC`, `AAGCTT`, `GGATCC`, `CTGCAG`, `GCGGCCGC`
- for the 4-tag design, the script explicitly searches for a set with perfect per-base balance, meaning 25% A / C / G / T at every tag position across the 4 tags
