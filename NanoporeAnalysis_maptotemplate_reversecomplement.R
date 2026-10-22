rm(list = ls())

## =========================
## 0) Setup
## =========================
setwd("/Users/sxv041/Desktop/LabStuff/NANOPOREANALYSIS")  # <-- change this

# ## Biostrings only (no ShortRead)
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# if (!requireNamespace("Biostrings", quietly = TRUE)) BiocManager::install("Biostrings")

library(Biostrings)
# 
# install.packages("RCurl")
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("GenomeInfoDb","Biostrings"))

## =========================
## 1) Read FASTQ into a data.frame (base R; supports .gz)
## =========================
fastq_to_df <- function(file) {
  is_gz <- grepl("\\.gz$", file, ignore.case = TRUE)
  con <- if (is_gz) gzfile(file, open = "rt") else file(file, open = "rt")
  on.exit(close(con))
  
  lines <- readLines(con)
  n <- length(lines)
  if (n == 0) stop("FASTQ file appears empty: ", file)
  if (n %% 4 != 0) stop("FASTQ malformed: line count not divisible by 4. Lines: ", n)
  
  ids_idx   <- seq(1, n, 4)
  seqs_idx  <- seq(2, n, 4)
  plus_idx  <- seq(3, n, 4)
  quals_idx <- seq(4, n, 4)
  
  data.frame(
    read_id   = sub("^@", "", lines[ids_idx]),
    sequence  = lines[seqs_idx],
    plus_line = lines[plus_idx],
    quality   = lines[quals_idx],
    length    = nchar(lines[seqs_idx]),
    stringsAsFactors = FALSE
  )
}

reads_df <- fastq_to_df("RN173-5b-SV65.fastq")   # <-- your FASTQ (or .fastq.gz)

## =========================
## 2) Read template FASTA and define ends
## =========================
suppressPackageStartupMessages({
  library(BiocGenerics)
  library(IRanges)
  library(Biostrings)
})

fa_lines <- readLines("CFTRedited_SV65Primers.fa")  # <-- your template FASTA
template_seq <- paste(fa_lines[!grepl("^>", fa_lines)], collapse = "")
template <- DNAString(template_seq)

end_len <- 500  # how much of each end to use as an anchor (300–800 typical)
left_end  <- subseq(template, start = 1, width = end_len)
right_end <- subseq(template, start = length(template) - end_len + 1, width = end_len)

## =========================
## 2.5) Sanitize sequences (fix hidden chars / non-IUPAC)
## =========================
# Keep IUPAC DNA letters only; replace anything else with N
sanitize_seq <- function(x) {
  x <- gsub("\\s+", "", x)  # remove whitespace
  x <- toupper(x)
  x <- gsub("[^ACGTRYSWKMBDHVN]", "N", x)
  x
}

reads_df$sequence <- vapply(reads_df$sequence, sanitize_seq, character(1))
reads_dna <- DNAStringSet(reads_df$sequence)



## =========================
## 3) Orientation-aware end-matching
## =========================


end_len      <- 500
max_mismatch <- 50
edge_window  <- 800

left_end_rc  <- reverseComplement(left_end)
right_end_rc <- reverseComplement(right_end)

end_hit_near_edge <- function(read_seq, end_seq, side = c("left","right"),
                              max_mismatch = 50, edge_window = 800) {
  side <- match.arg(side)
  w <- length(read_seq)
  if (w < 100) return(FALSE)
  
  hits <- Biostrings::matchPattern(
    pattern = end_seq,
    subject = read_seq,
    max.mismatch = max_mismatch,
    with.indels = TRUE,
    fixed = FALSE
  )
  if (length(hits) == 0) return(FALSE)
  
  if (side == "left") any(start(hits) <= edge_window)
  else               any(end(hits) >= (w - edge_window + 1))
}

# For each read, evaluate forward and reverse-complement configurations
call_ends <- function(read_seq, max_mismatch, edge_window) {
  # Forward orientation: left_end at left, right_end at right
  f_left  <- end_hit_near_edge(read_seq, left_end,  "left",  max_mismatch, edge_window)
  f_right <- end_hit_near_edge(read_seq, right_end, "right", max_mismatch, edge_window)
  
  # Reverse orientation: (template right end RC) at left, (template left end RC) at right
  r_left  <- end_hit_near_edge(read_seq, right_end_rc, "left",  max_mismatch, edge_window)
  r_right <- end_hit_near_edge(read_seq, left_end_rc,  "right", max_mismatch, edge_window)
  
  fwd_both <- f_left & f_right
  rev_both <- r_left & r_right
  
  if (fwd_both) {
    list(has_left=TRUE, has_right=TRUE, orientation="forward")
  } else if (rev_both) {
    list(has_left=TRUE, has_right=TRUE, orientation="reverse")
  } else {
    # “Left” / “right” here are in TEMPLATE terms, regardless of read direction:
    # template-left could be found as left_end at left (forward) OR left_end_rc at right (reverse)
    template_left_present  <- f_left | r_right
    template_right_present <- f_right | r_left
    
    list(has_left=template_left_present,
         has_right=template_right_present,
         orientation=ifelse((r_left | r_right) & !(f_left | f_right), "reverse_or_partial", "unknown"))
  }
}

res <- lapply(reads_dna, call_ends, max_mismatch=max_mismatch, edge_window=edge_window)

reads_df$has_left  <- vapply(res, `[[`, logical(1), "has_left")
reads_df$has_right <- vapply(res, `[[`, logical(1), "has_right")
reads_df$orientation <- vapply(res, `[[`, character(1), "orientation")

reads_df$class <- ifelse(
  reads_df$has_left & reads_df$has_right, "both_ends",
  ifelse(reads_df$has_left & !reads_df$has_right, "left_only",
         ifelse(!reads_df$has_left & reads_df$has_right, "right_only", "neither"))
)

## =========================
## 4) Split into 3 new structures
## =========================

reads_both_ends  <- reads_df[reads_df$class == "both_ends", ]
reads_left_only  <- reads_df[reads_df$class == "left_only", ]
reads_right_only <- reads_df[reads_df$class == "right_only", ]

## (optional) If you also want the "neither" group:
reads_neither <- reads_df[reads_df$class == "neither", ]

## quick sanity check
table(reads_df$class)
c(both_ends = nrow(reads_both_ends),
  left_only = nrow(reads_left_only),
  right_only = nrow(reads_right_only),
  neither = nrow(reads_neither))

##

reads_both_ends_dna  <- reads_dna[reads_df$class == "both_ends"]
reads_left_only_dna  <- reads_dna[reads_df$class == "left_only"]
reads_right_only_dna <- reads_dna[reads_df$class == "right_only"]


## =========================
## Histograms with common bins (4 classes)
## =========================

all_lengths <- c(
  reads_both_ends$length,
  reads_left_only$length,
  reads_right_only$length,
  reads_neither$length
)

breaks <- pretty(range(all_lengths), n = 100)

par(mfrow = c(4, 1), mar = c(4, 4, 2, 1))

hist(reads_both_ends$length,
     breaks = breaks,
     col = "steelblue",
     main = "Reads with both ends",
     xlab = "Read length (bp)")

hist(reads_left_only$length,
     breaks = breaks,
     col = "darkorange",
     main = "Reads with left end only",
     xlab = "Read length (bp)")

hist(reads_right_only$length,
     breaks = breaks,
     col = "forestgreen",
     main = "Reads with right end only",
     xlab = "Read length (bp)")

hist(reads_neither$length,
     breaks = breaks,
     col = "grey50",
     main = "Reads with neither end",
     xlab = "Read length (bp)")

par(mfrow = c(1, 1))


## =========================
## Define peak windows (bp)
## =========================

peak_6kb  <- c(5500, 7000)
peak_10kb <- c(9500, 13000)


lengths_both <- reads_both_ends$length

n_6kb  <- sum(lengths_both >= peak_6kb[1]  & lengths_both <= peak_6kb[2])
n_10kb <- sum(lengths_both >= peak_10kb[1] & lengths_both <= peak_10kb[2])

n_total <- length(lengths_both)

ratio_10kb_to_6kb <- n_10kb / n_6kb
ratio_10kb_to_total <- n_10kb / n_total

cat("Both-end reads summary:\n")
cat("  Total reads:", n_total, "\n")
cat("  ~6 kb reads:", n_6kb, "\n")
cat("  ~10 kb reads:", n_10kb, "\n\n")

cat("Ratios:\n")
cat("  10 kb / 6 kb:", round(ratio_10kb_to_6kb, 3), "\n")
cat("  10 kb / total:", round(ratio_10kb_to_total, 3), "\n")


write.csv(subset(reads_df, class == "left_only"),
          "reads_left_only.csv", row.names = FALSE)
write.csv(subset(reads_df, class == "right_only"),
          "reads_right_only.csv", row.names = FALSE)
write.csv(subset(reads_df, class == "both_ends"),
          "reads_both_ends.csv", row.names = FALSE)
