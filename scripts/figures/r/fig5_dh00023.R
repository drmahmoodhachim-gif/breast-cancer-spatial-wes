#!/usr/bin/env Rscript
# DH00023: VAF groups, clonal spheres, chr7 copy number, KMT2C lollipop.
source("/home/chimeric/WES_TJN_project/scripts/figures/r/fig_style.R")
suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
  library(gridExtra)
  library(ggrepel)
})

root <- PROJECT_DIR
ccf <- read.delim(file.path(root, "09_ccf_clonality/ccf.tsv"), stringsAsFactors = FALSE)
ccf <- ccf[ccf$patient == "DH00023", ]
lin <- read.delim(file.path(root, "09_ccf_clonality/DH00023_lineages.tsv"), stringsAsFactors = FALSE)
vars <- lin[lin$record == "variant", c("var_id", "lineage")]
ccf <- merge(ccf, vars, by = "var_id", all.x = TRUE)
ccf$lineage <- ifelse(ccf$lineage %in% c("trunk", "J-lineage", "T-lineage"), ccf$lineage, "unassigned")
group_col <- c(trunk = COL_SHARED, "J-lineage" = COL_J, "T-lineage" = COL_T, unassigned = "#bdbdbd")
p_a <- ggplot(ccf, aes(T_vaf, J_vaf, colour = lineage)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.3, colour = "#e4e3df") +
  geom_point(size = 1.2, alpha = 0.9) +
  scale_colour_manual(values = group_col, name = NULL) +
  coord_cartesian(xlim = c(0, 0.3), ylim = c(0, 0.3)) +
  labs(x = "Tumour VAF", y = "Junction VAF", tag = "a") +
  theme_paper() +
  theme(legend.position = "bottom", plot.tag = element_text(face = "bold", size = 8))

bars <- read.delim(file.path(root, "09_ccf_clonality/dh00023_vaf_ratio_bars.tsv"), stringsAsFactors = FALSE)
bars$region <- factor(bars$region, levels = c("J", "T"))
bars$lineage <- factor(bars$lineage, levels = c("trunk", "J-lineage", "T-lineage"))
p_b <- ggplot(bars, aes(region, vaf_ratio_capped, fill = lineage)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.64, colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = c(trunk = COL_SHARED, "J-lineage" = COL_J, "T-lineage" = COL_T), name = NULL) +
  scale_y_continuous(limits = c(0, 1.15), expand = c(0, 0)) +
  labs(x = NULL, y = "VAF-ratio estimate", tag = "b") +
  theme_paper() +
  theme(legend.position = "bottom", plot.tag = element_text(face = "bold", size = 8))
nobs <- read.delim(file.path(root, "09_ccf_clonality/dh00023_n_observed_expected.tsv"), stringsAsFactors = FALSE)
nobs$lineage <- factor(nobs$lineage, levels = c("trunk", "J-lineage", "T-lineage"))
nlong <- rbind(
  data.frame(lineage = nobs$lineage, which = "Observed", value = nobs$observed_alt, lo = NA_real_, hi = NA_real_),
  data.frame(lineage = nobs$lineage, which = "Expected", value = nobs$expected_alt, lo = nobs$betabinom_lo, hi = nobs$betabinom_hi)
)
p_in <- ggplot(nlong, aes(lineage, value, fill = which)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.64, colour = "white", linewidth = 0.2) +
  geom_errorbar(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.72), width = 0.15, linewidth = 0.3, na.rm = TRUE) +
  scale_fill_manual(values = c(Observed = "#3d3d3a", Expected = "#c9c8c2"), name = NULL) +
  labs(x = NULL, y = "Normal alt reads") +
  theme_paper() +
  theme(legend.position = "bottom", axis.text.x = element_text(size = 6))

seg_files <- list.files(file.path(root, "07_cnv_purity/facets"), pattern = "DH00023_[JT].*segments.tsv$", recursive = TRUE, full.names = TRUE)
seg <- do.call(rbind, lapply(seg_files, read.delim, stringsAsFactors = FALSE))
seg <- seg[sub("^chr", "", seg$chrom) == "7", ]
seg$sample <- ifelse(grepl("_T$", seg$sample_id), "T", "J")
seg$sample <- factor(seg$sample, levels = c("J", "T"))
steps <- rbind(
  data.frame(sample = seg$sample, x = seg$start / 1e6, tcn = seg$tcn.em, lcn = seg$lcn.em),
  data.frame(sample = seg$sample, x = seg$end / 1e6, tcn = seg$tcn.em, lcn = seg$lcn.em)
)
steps <- steps[order(steps$sample, steps$x), ]
p_c <- ggplot(steps, aes(x)) +
  geom_step(aes(y = tcn), direction = "hv", linewidth = 0.4) +
  geom_step(aes(y = lcn), direction = "hv", linewidth = 0.4, linetype = "dashed") +
  geom_vline(xintercept = 152.16, linewidth = 0.3, colour = INK) +
  annotate("text", x = 150.2, y = 3.55, label = "KMT2C", size = 2.6, family = FONT, hjust = 1) +
  facet_grid(sample ~ .) +
  coord_cartesian(ylim = c(0, 4), xlim = c(0, 159.3)) +
  labs(x = "chr7 (Mb)", y = "Copy number", tag = "c") +
  theme_paper() +
  theme(strip.text = element_text(size = 8, face = "bold"), plot.tag = element_text(face = "bold", size = 8))
g_c <- ggplotGrob(p_c)

# UniProt Q8NEZ4 isoform 1, MANE-Select ENST00000262189.11 / NM_170606.3.
# Sequence check: residue 2197 is Q and residue 4553 is L. 4553 is inside FYRN 4545-4605.
dom <- read.delim(file.path(root, "06_variants/kmt2c_Q8NEZ4_domains.tsv"), stringsAsFactors = FALSE)
dom$kind <- ifelse(grepl("^PHD", dom$name), "PHD",
                   ifelse(grepl("^prePHD", dom$name), "pre-PHD", dom$name))
dom_main <- dom[dom$name != "RING", ]
dom_ring <- dom[dom$name == "RING", ]
dom_col <- c(PHD = "#5C6B7A", "pre-PHD" = "#8E99A4", RING = "#1F4E79",
             DHHC = "#6B8F71", FYRN = "#C47B2B", FYRC = "#7A4E62",
             SET = "#3d3d3a", PostSET = "#9A9A94")
pins <- data.frame(x = c(2197, 4553), label = c("p.Gln2197Ter", "p.Leu4553Arg"),
                   fill = c(COL_J, COL_T), hjust = c(0.5, 1), stringsAsFactors = FALSE)
p_d <- ggplot() +
  geom_rect(aes(xmin = 1, xmax = 4911, ymin = -0.04, ymax = 0.04), fill = "#e4e3df", colour = NA) +
  geom_rect(data = dom_main, aes(xmin = start, xmax = end, ymin = -0.08, ymax = 0.08, fill = kind), colour = NA) +
  geom_rect(data = dom_ring, aes(xmin = start, xmax = end, ymin = 0.10, ymax = 0.18, fill = kind), colour = NA) +
  scale_fill_manual(values = dom_col, name = NULL) +
  geom_segment(data = pins, aes(x = x, xend = x, y = 0.22, yend = 0.48), linewidth = 0.35, colour = INK) +
  geom_point(data = pins, aes(x = x, y = 0.52), size = 2.6, colour = pins$fill) +
  geom_label(data = pins, aes(x = x, y = 0.78, label = label), hjust = pins$hjust, fill = pins$fill, colour = "white",
             size = 3.2, family = FONT, label.size = 0, label.padding = unit(0.15, "lines")) +
  coord_cartesian(xlim = c(1, 4911), ylim = c(-0.55, 1.05), expand = FALSE) +
  labs(x = "KMT2C residue (Q8NEZ4, NM_170606.3)", y = NULL, tag = "d") +
  theme_paper() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "bottom", plot.tag = element_text(face = "bold", size = 8))
g_d <- ggplotGrob(p_d)

g_a <- ggplotGrob(p_a)
g_b <- ggplotGrob(p_b)
g_in <- ggplotGrob(p_in)
g_bstack <- arrangeGrob(g_b, g_in, ncol = 1, heights = c(1.35, 1.15))
draw_all <- function() {
  grid.newpage()
  gridExtra::grid.arrange(g_a, g_bstack, g_c, g_d, ncol = 2,
                          top = textGrob("DH00023", gp = gpar(fontsize = 7, fontfamily = FONT)))
}
inputs <- c(
  file.path(root, "09_ccf_clonality/ccf.tsv"),
  file.path(root, "09_ccf_clonality/DH00023_lineages.tsv"),
  file.path(root, "09_ccf_clonality/dh00023_vaf_ratio_bars.tsv"),
  file.path(root, "09_ccf_clonality/dh00023_n_observed_expected.tsv"),
  file.path(root, "07_cnv_purity/facets/DH00023_T/DH00023_T.segments.tsv"),
  file.path(root, "07_cnv_purity/facets/DH00023_J/DH00023_J.segments.tsv"),
  file.path(root, "06_variants/kmt2c_Q8NEZ4_domains.tsv")
)
save_fig(draw_all, "Fig5_DH00023", width_mm = 183, height_mm = 175, inputs = inputs,
         script = "fig5_dh00023.R",
         note = "Panel b upper bars are VAF-ratio estimates, not cell fractions. Domain bars are every UniProt Q8NEZ4 Domain and zinc-finger feature on NM_170606.3. RING is drawn above PHD2 because the intervals overlap. Residue 4553 is inside FYRN 4545-4605. Residue 2197 is Q and is not inside a named domain.")
