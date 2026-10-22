#!/usr/bin/env Rscript
# Usage: Rscript 05_facets.R <snp_pileup.csv.gz> <outdir> <sample_id> <cval>
suppressPackageStartupMessages(library(facets))
args <- commandArgs(trailingOnly = TRUE)
pile <- args[1]; outdir <- args[2]; sid <- args[3]; cval <- as.numeric(args[4])
set.seed(1234)  # FACETS uses random SNP sampling; fixed seed for reproducibility

rcmat <- readSnpMatrix(pile)
xx    <- preProcSample(rcmat, gbuild = "hg38", ndepth = 25, snp.nbhd = 250)
oo    <- procSample(xx, cval = cval)
fit   <- emcncf(oo)

seg <- fit$cncf
seg$chrom <- ifelse(seg$chrom == 23, "X", as.character(seg$chrom))
seg$sample_id <- sid
write.table(seg[, c("sample_id","chrom","start","end","num.mark","cnlr.median","mafR","cf.em","tcn.em","lcn.em")],
            file.path(outdir, paste0(sid, ".segments.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

# Genome-doubling heuristic: >50% of genome with major CN >= 2
seglen <- seg$end - seg$start
major  <- seg$tcn.em - ifelse(is.na(seg$lcn.em), 0, seg$lcn.em)
wgd    <- sum(seglen[major >= 2], na.rm = TRUE) / sum(seglen, na.rm = TRUE) > 0.5
loh    <- sum(seglen[!is.na(seg$lcn.em) & seg$lcn.em == 0 & seg$tcn.em > 0]) / sum(seglen)

write.table(data.frame(sample_id = sid, purity = fit$purity, ploidy = fit$ploidy,
                       dipLogR = oo$dipLogR, frac_genome_LOH = round(loh, 4), WGD_heuristic = wgd, cval = cval),
            file.path(outdir, paste0(sid, ".purity_ploidy.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

pdf(file.path(outdir, paste0(sid, ".facets_qc.pdf")), width = 8, height = 8)
plotSample(x = oo, emfit = fit, sname = sid)
logRlogORspider(oo$out, oo$dipLogR)
dev.off()
cat(sprintf("[facets] %s purity=%.2f ploidy=%.2f\n", sid, fit$purity, fit$ploidy))
