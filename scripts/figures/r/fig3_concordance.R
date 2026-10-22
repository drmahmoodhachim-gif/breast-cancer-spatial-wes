#!/usr/bin/env Rscript
# SNV sharing versus breakpoint sharing, and MEDICC2 trees with WGD branches.
.script_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
.fig_r_dir <- if (length(.script_file)) dirname(normalizePath(.script_file[1])) else file.path(Sys.getenv("PROJECT_DIR", "."), "scripts", "figures", "r")
source(file.path(.fig_r_dir, "fig_style.R"))
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggtree)
  library(ape)
  library(patchwork)
  library(ggrepel)
})

root <- PROJECT_DIR
var <- read.delim(file.path(root, "06_variants/somatic_filtered.tsv"), stringsAsFactors = FALSE)
bp <- read.delim(file.path(root, "07_cnv_purity/cn_phylogeny/breakpoint_sharing.tsv"), stringsAsFactors = FALSE)
snv <- var[var$variant_type == "SNV" | grepl("^SNV$|SNP", var$variant_type) | (nchar(var$ref) == 1 & nchar(var$alt) == 1), ]
if (!nrow(snv)) snv <- var
patients <- c("DH00018","DH00019","DH00022","DH00023","DH00024")
rows <- lapply(patients, function(p) {
  g <- snv[snv$patient == p, ]
  t_present <- g$T_present %in% c(TRUE, "True", "true")
  j_present <- g$J_present %in% c(TRUE, "True", "true")
  n_t <- sum(t_present)
  data.frame(patient = p,
             snv_shared_fraction = if (n_t) sum(t_present & j_present) / n_t else NA_real_,
             stringsAsFactors = FALSE)
})
d <- merge(do.call(rbind, rows), bp, by = "patient")
d$hjust <- ifelse(d$snv_shared_fraction > 0.75, 1.08, -0.08)
d$vjust <- ifelse(d$fraction_T_within_1Mb_of_J > 0.82, 1.3, -0.4)
d$excluded <- d$patient == "DH00024"
p_sc <- ggplot(d, aes(snv_shared_fraction, fraction_T_within_1Mb_of_J, label = patient)) +
  geom_point(aes(shape = excluded), size = 2, colour = COL_SHARED) +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 1), guide = "none") +
  ggrepel::geom_text_repel(size = 2, family = FONT, colour = INK, min.segment.length = 0,
                           segment.size = 0.2, box.padding = 0.6, point.padding = 0.45, force = 4,
                           max.overlaps = 20) +
  coord_cartesian(xlim = c(0, 1.02), ylim = c(0, 1.02)) +
  labs(x = "Fraction of tumour SNVs also in J", y = "Fraction of T breakpoints within 1 Mb of J", tag = "a") +
  theme_paper() +
  theme(plot.tag = element_text(face = "bold", size = 8))

trees <- lapply(seq_along(patients), function(i) {
  p <- patients[i]
  path <- file.path(root, "07_cnv_purity/cn_phylogeny", p, paste0(p, ".medicc2_input_final_tree.new"))
  summary <- file.path(root, "07_cnv_purity/cn_phylogeny", p, paste0(p, ".medicc2_input_summary.tsv"))
  tr <- ape::read.tree(path)
  wgd <- FALSE
  if (file.exists(summary)) {
    sm <- read.delim(summary, header = FALSE, stringsAsFactors = FALSE)
    wgd <- any(grepl("^WGD", sm$V2))
  }
  tip_col <- setNames(ifelse(grepl("_N$", tr$tip.label), COL_N,
                             ifelse(grepl("_J$", tr$tip.label), COL_J, COL_T)), tr$tip.label)
  depth <- max(node.depth.edgelength(tr))
  off <- max(0.08, 0.06 * depth)
  gt <- ggtree(tr) +
    geom_tiplab(aes(colour = label), size = 2.2, offset = off, show.legend = FALSE) +
    geom_treescale(x = 0, y = 0.2, fontsize = 2, linesize = 0.3, offset = 0.15) +
    scale_colour_manual(values = tip_col) +
    xlim(0, depth + off + max(1.5, 0.45 * depth)) +
    labs(title = p, tag = if (i == 1) "b" else NULL) +
    theme_void(base_family = FONT) +
    theme(plot.title = element_text(size = 7, family = FONT, colour = INK, hjust = 0.5),
          plot.tag = element_text(face = "bold", size = 8))
  if (wgd) {
    node <- ape::getMRCA(tr, c(paste0(p, "_J"), paste0(p, "_T")))
    gt <- gt + geom_label2(aes(subset = (node == !!node), label = "WGD"),
                           size = 1.8, colour = "white", fill = COL_T, label.size = 0)
  }
  gt
})

panel <- p_sc + wrap_plots(trees, ncol = 1) + plot_layout(widths = c(0.7, 1.3))
inputs <- c(
  file.path(root, "06_variants/somatic_filtered.tsv"),
  file.path(root, "07_cnv_purity/cn_phylogeny/breakpoint_sharing.tsv")
)
save_fig(function() print(panel), "Fig3_concordance", width_mm = 183, height_mm = 140,
         inputs = inputs, script = "fig3_concordance.R",
         note = "SNV share is the fraction of T-present SNVs also present in J. Breakpoint share is D034. WGD labels are MEDICC2 events on the J+T stem.")
