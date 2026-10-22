#!/usr/bin/env Rscript
# T-cell fraction per sample from TCRA-locus coverage (T cell ExTRECT; Bentham et al., Nature 2021).
# Rearranged T cells delete part of the TCRA locus (chr14), so exome coverage there dips in proportion
# to the fraction of T cells in the sample.
#
# Usage: Rscript 11_tcell_extrect.R <project_dir> <bam_root_dir>
# Inputs : recalibrated BAMs (found under bam_root_dir), padded Exome 2.0 Plus target BED,
#          optional FACETS purity/segments for T and J (purity + TCRA copy-number adjustment)
# Outputs: 14_immune/tcell/tcell_fraction.tsv, per-sample QC PDFs, provenance in the log
suppressPackageStartupMessages({ library(TcellExTRECT) })
args <- commandArgs(trailingOnly = TRUE)
proj <- args[1]; bam_root <- args[2]
out_dir <- file.path(proj, "14_immune", "tcell"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
bed <- file.path(proj, "03_resources", "exome2plus_targets_hg38.bed")
if (!file.exists(bed)) stop("Target BED not found: ", bed)

data("tcra_seg_hg38")
# Exome 2.0 Plus is not one of the bundled kits -> build the TCRA exon table from our own capture BED
exons <- createExonDFBed(bed, "hg38")
cat(sprintf("[tcell] TCRA exons from capture BED: %d\n", nrow(exons)))
if (nrow(exons) < 10) warning("Few TCRA exons in the capture design: estimates may be unreliable")

pp_file <- file.path(proj, "07_cnv_purity", "purity_reconciled.tsv")
pp <- if (file.exists(pp_file)) read.delim(pp_file) else NULL
tcra_cn <- function(sid) {   # FACETS total CN at the TCRA locus (chr14 ~22.0 Mb, hg38)
  f <- file.path(proj, "07_cnv_purity", "facets", sid, paste0(sid, ".segments.tsv"))
  if (!file.exists(f)) return(NA)
  s <- read.delim(f); s$chrom <- sub("^chr", "", as.character(s$chrom))
  hit <- s[s$chrom == "14" & s$start <= 22.0e6 & s$end >= 22.0e6, ]
  if (nrow(hit)) hit$tcn.em[1] else NA
}

bams <- list.files(bam_root, pattern = "^DH[0-9]{5}_[NJT]\\.recal\\.bam$", recursive = TRUE, full.names = TRUE)
if (!length(bams)) stop("No *.recal.bam found under ", bam_root)
res <- list()
for (b in sort(bams)) {
  sid <- sub("\\.recal\\.bam$", "", basename(b))
  cov_file <- tryCatch(getCovFromBam(bamPath = b, outPath = paste0(out_dir, "/"), vdj.seg = tcra_seg_hg38),
                       error = function(e) NA)
  if (is.na(cov_file) || !file.exists(cov_file) || file.size(cov_file) == 0) {
    # Fallback: compute depth ourselves (chr-prefixed GRCh38) and give ExTRECT a pos/reads table
    cov_file <- file.path(out_dir, paste0(sid, "_TCRA.txt"))
    system2("samtools", c("depth", "-a", "-Q", "20", "-r", "chr14:21000000-23000000", b), stdout = cov_file)
    cov <- read.delim(cov_file, header = FALSE, col.names = c("chr", "pos", "reads"))[, c("pos", "reads")]
  } else {
    cov <- loadCov(cov_file)
  }
  out <- tryCatch(runTcellExTRECT(cov, exons, tcra_seg_hg38, "hg38", sample_name = sid),
                  error = function(e) runTcellExTRECT(cov, exons, tcra_seg_hg38, "hg38"))
  out$sample <- sid
  specimen <- sub(".*_", "", sid)
  out$purity <- NA; out$TCRA_cn <- NA; out$adjusted <- FALSE
  if (specimen %in% c("T", "J") && !is.null(pp) && sid %in% pp$sample) {
    pur <- pp$used_purity[pp$sample == sid]; cn <- tcra_cn(sid)
    if (!is.na(pur) && !is.na(cn) && cn > 0) {
      qf <- out$qcFit  # adjustTcellExTRECT drops qcFit; keep it
      out <- adjustTcellExTRECT(out, purity = pur, TCRA.cn = cn)
      out$qcFit <- qf
      out$purity <- pur; out$TCRA_cn <- cn; out$adjusted <- TRUE
    }
  }
  pdf(file.path(out_dir, paste0(sid, "_tcell_extrect_qc.pdf")), width = 7, height = 4)
  print(plotTcellExTRECT(cov, exons, tcra_seg_hg38, "hg38", sample_name = sid))
  dev.off()
  res[[sid]] <- out
  cat(sprintf("[tcell] %s fraction=%.4f qcFit=%.2f adjusted=%s\n", sid,
              out$TCRA.tcell.fraction[1], out$qcFit[1], out$adjusted[1]))
}
all_cols <- unique(unlist(lapply(res, names)))
tab <- do.call(rbind, lapply(res, function(x) { for (c in setdiff(all_cols, names(x))) x[[c]] <- NA; as.data.frame(x)[, all_cols] }))
write.table(tab, file.path(out_dir, "tcell_fraction.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("[tcell] wrote", file.path(out_dir, "tcell_fraction.tsv"), "\n")
