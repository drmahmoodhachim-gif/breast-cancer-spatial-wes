#!/usr/bin/env Rscript
# Genomic scar HRD scores from FACETS allele-specific segments (scarHRD).
# Input columns follow scarHRD's allele-specific file: sample, chromosome, start, end,
# total CN, A allele CN, B allele CN, ploidy. Chromosomes are chr-prefixed (GRCh38).
# Integer copy number is the FACETS estimate; ploidy is the FACETS ploidy.
# Output: 10_tmb_signatures/hrd_scar.tsv
args <- commandArgs(trailingOnly = TRUE)
proj <- if (length(args)) args[[1]] else "."
suppressPackageStartupMessages(library(scarHRD))

pp <- read.delim(file.path(proj, "07_cnv_purity", "purity_ploidy.tsv"), stringsAsFactors = FALSE)
work <- file.path(proj, "10_tmb_signatures", "hrd_scar_work")
dir.create(work, recursive = TRUE, showWarnings = FALSE)
rows <- list()
for (i in seq_len(nrow(pp))) {
  sid <- pp$sample_id[i]
  ploidy <- as.numeric(pp$ploidy[i])
  f <- file.path(proj, "07_cnv_purity", "facets", sid, paste0(sid, ".segments.tsv"))
  seg <- read.delim(f, stringsAsFactors = FALSE)
  seg <- seg[!is.na(seg$tcn.em) & !is.na(seg$lcn.em) & seg$tcn.em >= 0, ]
  chrom <- sub("^chr", "", as.character(seg$chrom))
  chrom <- paste0("chr", chrom)
  minor <- as.numeric(seg$lcn.em)
  total <- as.numeric(seg$tcn.em)
  major <- pmax(0, total - minor)
  tab <- data.frame(SampleID = sid, Chromosome = chrom, Start_position = seg$start,
                    End_position = seg$end, total_cn = total, A_cn = minor, B_cn = major,
                    ploidy = ploidy, check.names = FALSE)
  infile <- file.path(work, paste0(sid, ".scar_input.tsv"))
  write.table(tab, infile, sep = "\t", quote = FALSE, row.names = FALSE)
  res <- scar_score(infile, reference = "grch38", seqz = FALSE, ploidy = ploidy,
                    chr.in.names = TRUE, outputdir = work)
  vals <- as.numeric(res[1, ])
  rows[[sid]] <- data.frame(sample = sid, LOH = vals[1], TAI = vals[2], LST = vals[3],
                            HRD_sum = vals[4], stringsAsFactors = FALSE)
  cat(sprintf("[hrd] %s LOH=%g TAI=%g LST=%g HRD_sum=%g\n", sid, vals[1], vals[2], vals[3], vals[4]))
}
out <- do.call(rbind, rows)
outfile <- file.path(proj, "10_tmb_signatures", "hrd_scar.tsv")
write.table(out, outfile, sep = "\t", quote = FALSE, row.names = FALSE)
cat("[hrd] wrote", outfile, "\n")
