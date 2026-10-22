[README.md](https://github.com/user-attachments/files/33834942/README.md)
# Long-read nanopore sequencing analysis

R code for classifying Oxford Nanopore long reads by the presence of template-derived ends, and for measuring the size distribution of the resulting fragments. Fragment size classes are then used to quantify gene insertion.

This repository accompanies:

> Nouri R., et al. *Replacement of full-length CFTR in human airway stem cells using a single AAV.* (Manuscript in preparation.)
>
> Vaidyanathan Laboratory, Abigail Wexner Research Institute, Nationwide Children's Hospital.

## Overview

`NanoporeAnalysis_maptotemplate_reversecomplement.R` takes a basecalled FASTQ file and a reference template FASTA, and asks, for every read, whether the first and last several hundred bases of the template are present near the corresponding ends of the read. Because nanopore reads arrive in either orientation, the script evaluates both the forward configuration and the reverse-complement configuration before making a call.

Each read is assigned to one of four classes:

| Class | Meaning |
| --- | --- |
| `both_ends` | Template left and right ends both detected near the read edges |
| `left_only` | Only the template left end detected |
| `right_only` | Only the template right end detected |
| `neither` | No template end detected |

Reads carrying both ends represent intact template-spanning molecules. Their length distribution separates into species of different sizes, and the ratio between those species is the quantitative readout of insertion.

## What the script does

1. Reads a FASTQ file into a data frame using base R. Plain and gzip-compressed files are both accepted.
2. Reads the template FASTA and extracts the terminal `end_len` bases from each end as anchor sequences.
3. Sanitizes read sequences, removing whitespace and replacing any non-IUPAC character with `N`.
4. Searches each read for the anchor sequences using `Biostrings::matchPattern` with mismatches and indels allowed, requiring the hit to fall within `edge_window` bases of the appropriate read edge. Forward and reverse-complement configurations are both tested.
5. Classifies every read into one of the four classes above.
6. Plots stacked length histograms for the four classes on shared bins.
7. Counts reads falling inside two user-defined size windows and reports the ratio between them.
8. Writes one CSV per class.

## Requirements

- R version 4.x
- [Biostrings](https://bioconductor.org/packages/Biostrings) (Bioconductor), which pulls in BiocGenerics and IRanges

Install with:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("Biostrings")
```

Exact package versions used for the analysis in the manuscript are recorded in [`sessionInfo.txt`](sessionInfo.txt).

## Input files

| File | Description |
| --- | --- |
| `*.fastq` or `*.fastq.gz` | Basecalled nanopore reads. The script expects a well-formed FASTQ whose line count is divisible by four. |
| `*.fa` | Reference template FASTA, single sequence, used to define the left and right anchors. |

Example input files sufficient to run the script end to end are provided in [`data/`](data/).

## Running the analysis

1. Clone or download the repository.
2. Place your FASTQ and template FASTA in the working directory, or point the script at them.
3. Open the script and set the three file paths at the top: the working directory passed to `setwd()`, the FASTQ name passed to `fastq_to_df()`, and the FASTA name passed to `readLines()`.
4. Adjust the parameters below if needed.
5. Run the script from top to bottom.

```r
source("NanoporeAnalysis_maptotemplate_reversecomplement.R")
```

## Parameters

| Parameter | Default | Description |
| --- | --- | --- |
| `end_len` | 500 | Length in bases of the template terminus used as an anchor. Values between 300 and 800 are typical. Longer anchors are more specific but less tolerant of truncated reads. |
| `max_mismatch` | 50 | Mismatches tolerated when matching an anchor, accommodating the nanopore error rate. |
| `edge_window` | 800 | How close to the read edge an anchor hit must fall to count. |
| `peak_6kb` | `c(5500, 7000)` | Size window in bases for the smaller fragment species. |
| `peak_10kb` | `c(9500, 13000)` | Size window in bases for the larger fragment species. |

The two size windows should be set by inspecting the `both_ends` histogram for your own template, since they depend on the expected sizes of the edited and unedited products.

## Outputs

Written to the working directory:

| File | Contents |
| --- | --- |
| `reads_both_ends.csv` | Reads with both template ends, one row per read, with ID, sequence, quality, length, end calls, and orientation |
| `reads_left_only.csv` | Reads with the left template end only |
| `reads_right_only.csv` | Reads with the right template end only |

Printed to the console:

- A contingency table of read counts per class
- Read counts inside each size window, and the ratio between windows

Plotted to the active graphics device:

- A four-panel figure of read-length histograms on shared bins, one panel per class

## Repository contents

```
.
├── NanoporeAnalysis_maptotemplate_reversecomplement.R   analysis script
├── data/                                                example input files
├── sessionInfo.txt                                      R and package versions
└── README.md
```

## Citation

If you use this code, please cite the manuscript above and this repository.

## Contact

Questions about the analysis can be directed to the corresponding author of the manuscript, or raised as an issue in this repository.
