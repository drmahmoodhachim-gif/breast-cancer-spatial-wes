#!/usr/bin/env Rscript
# Figure 2 oncoprint. Rows: drivers, core pathway genes hit in >=2 patients, and
# named copy-number events. Columns: junction and tumour.
.script_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
.fig_r_dir <- if (length(.script_file)) dirname(normalizePath(.script_file[1])) else file.path(Sys.getenv("PROJECT_DIR", "."), "scripts", "figures", "r")
source(file.path(.fig_r_dir, "fig_style.R"))
suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

root <- PROJECT_DIR
var <- read.delim(file.path(root, "06_variants/somatic_filtered.tsv"), stringsAsFactors = FALSE)
pur <- read.delim(file.path(root, "07_cnv_purity/purity_ploidy.tsv"), stringsAsFactors = FALSE)
ccf <- read.delim(file.path(root, "09_ccf_clonality/ccf.tsv"), stringsAsFactors = FALSE)
hrd <- read.delim(file.path(root, "10_tmb_signatures/hrd_scar.tsv"), stringsAsFactors = FALSE)
tmb <- read.delim(file.path(root, "12_tables/main/tmb.tsv"), stringsAsFactors = FALSE)
clin <- read.delim(file.path(root, "config/clinical_metadata.tsv"), stringsAsFactors = FALSE)
seg_files <- list.files(file.path(root, "07_cnv_purity/facets"), pattern = "segments.tsv$", recursive = TRUE, full.names = TRUE)
seg <- do.call(rbind, lapply(seg_files, read.delim, stringsAsFactors = FALSE))
seg$chrom <- sub("^chr", "", as.character(seg$chrom))

core <- c("PIK3CA","PIK3R1","AKT1","PTEN","MTOR","TSC1","TSC2",
          "BRCA2","BRCA1","PALB2","FANCM","RAD54B","DMC1","POLE2","LIG3","UNG","ATM","CHEK2",
          "TP53","RBL1","RB1","CDKN2A","CCND1",
          "KMT2C","KMT2D","NSD2","BPTF","SETDB2","SIN3B","BCORL1","BCL11A","RBBP4","ARID1A",
          "AMER1","WNT9A","FZD3","WLS","CXXC5","MAP3K4","EPHA5","NRG2","PRKCB","RALGAPA2",
          "RHOA","DLC1","WWC2","MGA","MYCL","ELF5","RUNX2")
pa <- c("missense_variant","stop_gained","frameshift_variant","start_lost","inframe_deletion",
        "inframe_insertion","splice_acceptor_variant","splice_donor_variant","stop_lost","protein_altering_variant")
var$csq1 <- sub("&.*", "", var$consequence)
sev <- function(c) {
  c <- as.character(c)
  if (any(grepl("stop_gained|frameshift|splice_acceptor|splice_donor|start_lost|stop_lost", c))) "Truncating"
  else if (any(grepl("missense", c))) "Missense"
  else if (any(grepl("inframe|protein_altering", c))) "In-frame"
  else ""
}

patients <- c("DH00018","DH00019","DH00022","DH00023","DH00024")
samples <- unlist(lapply(patients, function(p) paste0(p, c("_J","_T"))))
drivers <- sort(unique(var$gene[var$is_driver_gene %in% c(TRUE, "True", "true") & var$csq1 %in% pa]))
core_n <- tapply(var$patient[var$gene %in% core & var$csq1 %in% pa], var$gene[var$gene %in% core & var$csq1 %in% pa], function(x) length(unique(x)))
core_keep <- sort(names(core_n)[core_n >= 2 & !names(core_n) %in% drivers])
snv_rows <- c(drivers, core_keep)

mut_mat <- matrix("", nrow = length(snv_rows), ncol = length(samples), dimnames = list(snv_rows, samples))
for (g in snv_rows) {
  sub <- var[var$gene == g & var$csq1 %in% pa, ]
  for (p in unique(sub$patient)) {
    sp <- sub[sub$patient == p, ]
    for (spec in c("J","T")) {
      present <- if (spec == "J") sp$J_present %in% c(TRUE, "True", "true") else sp$T_present %in% c(TRUE, "True", "true")
      if (!any(present)) next
      mut_mat[g, paste0(p, "_", spec)] <- sev(sp$csq1[present])
    }
  }
}

cn_events <- read.delim(file.path(root, "07_cnv_purity/cn_events.tsv"), stringsAsFactors = FALSE)
cn_events <- cn_events[cn_events$called %in% c(TRUE, "True", "true"), ]
cn_label <- c(amp = "Amp", gain = "Gain", loss = "Loss", loh = "LOH")
cn_rows <- c("CCND1 amp", "FGFR1 amp", "8q gain", "1q gain", "16q loss", "17p LOH")
cn_mat <- matrix("", nrow = length(cn_rows), ncol = length(samples), dimnames = list(cn_rows, samples))
for (i in seq_len(nrow(cn_events))) {
  ev <- cn_events$event[i]
  if (!ev %in% cn_rows) next
  kind <- if (grepl("amp", ev)) "Amp" else if (grepl("gain", ev)) "Gain" else if (grepl("loss", ev)) "Loss" else "LOH"
  cn_mat[ev, cn_events$sample[i]] <- kind
}
rownames(mut_mat)[rownames(mut_mat) == "BRCA2"] <- "BRCA2 (VUS)"
mat <- rbind(mut_mat, cn_mat)

purity <- ploidy <- wgd <- hrd_sum <- tmb_v <- er <- pr <- her2 <- grade <- setNames(rep(NA, length(samples)), samples)
for (sid in samples) {
  p <- substr(sid, 1, 7); spec <- substr(sid, 9, 9)
  hit <- ccf[ccf$patient == p, ]
  purity[sid] <- if (spec == "T") hit$T_purity[1] else hit$J_purity[1]
  ploidy[sid] <- pur$ploidy[match(sid, pur$sample_id)]
  wgd[sid] <- ifelse(pur$WGD_heuristic[match(sid, pur$sample_id)] %in% c(TRUE, "True", "true"), "Yes", "No")
  hrd_sum[sid] <- hrd$HRD_sum[match(sid, hrd$sample)]
  tmb_v[sid] <- tmb$TMB_coding_per_mb[match(sid, tmb$sample_id)]
  cl <- clin[clin$patient == p, ]
  er[sid] <- as.character(cl$ER)
  pr[sid] <- as.character(cl$PR)
  her2[sid] <- as.character(cl$HER2)
  grade[sid] <- as.character(cl$grade)
}
keep <- function(x, cols) cols[names(cols) %in% unique(as.character(x))]
er <- as.character(er); pr <- as.character(pr); her2 <- as.character(her2); grade <- as.character(grade)
col_fun_samples <- ifelse(grepl("_T$", samples), COL_T, COL_J)

col <- c(Missense = "#2a78d6", Truncating = "#eb6834", "In-frame" = "#1baf7a",
         Amp = "#B2182B", Gain = "#F4A582", Loss = "#92C5DE", LOH = "#4D4D4D")
col <- col[names(col) %in% unique(mat[nzchar(mat)])]

top <- HeatmapAnnotation(
  Purity = anno_barplot(purity, gp = gpar(fill = col_fun_samples), ylim = c(0, 1), height = unit(8, "mm")),
  Ploidy = anno_barplot(ploidy, gp = gpar(fill = col_fun_samples), ylim = c(0, 4), height = unit(8, "mm")),
  WGD = wgd,
  HRD = anno_barplot(hrd_sum, gp = gpar(fill = col_fun_samples), height = unit(8, "mm")),
  TMB = anno_barplot(tmb_v, gp = gpar(fill = col_fun_samples), height = unit(8, "mm")),
  col = list(WGD = keep(wgd, c(Yes = "#4D4D4D", No = "#e4e3df"))),
  annotation_name_gp = gpar(fontsize = 7, fontfamily = FONT),
  annotation_name_side = "left",
  simple_anno_size = unit(4, "mm"),
  show_legend = c(WGD = TRUE),
  annotation_legend_param = list(WGD = list(
    title_gp = gpar(fontsize = 7, fontfamily = FONT),
    labels_gp = gpar(fontsize = 7, fontfamily = FONT)))
)
bottom <- HeatmapAnnotation(
  ER = er, PR = pr, HER2 = her2, Grade = grade,
  col = list(
    ER = keep(er, c(Pos = "#4D4D4D", Neg = "#d9d9d9")),
    PR = keep(pr, c(Pos = "#4D4D4D", Neg = "#d9d9d9")),
    HER2 = keep(her2, c(Neg = "#d9d9d9", Equivocal = "#878787", Pos = "#4D4D4D")),
    Grade = keep(grade, c("2" = "#d9d9d9", "3" = "#4D4D4D"))
  ),
  annotation_name_gp = gpar(fontsize = 7, fontfamily = FONT),
  annotation_name_side = "left",
  simple_anno_size = unit(4, "mm"),
  annotation_legend_param = list(
    ER = list(title_gp = gpar(fontsize = 7, fontfamily = FONT), labels_gp = gpar(fontsize = 7, fontfamily = FONT)),
    PR = list(title_gp = gpar(fontsize = 7, fontfamily = FONT), labels_gp = gpar(fontsize = 7, fontfamily = FONT)),
    HER2 = list(title_gp = gpar(fontsize = 7, fontfamily = FONT), labels_gp = gpar(fontsize = 7, fontfamily = FONT)),
    Grade = list(title_gp = gpar(fontsize = 7, fontfamily = FONT), labels_gp = gpar(fontsize = 7, fontfamily = FONT))
  )
)

mat[mat == ""] <- "none"
heat_col <- c(none = "#f7f6f4", col)
heat_col <- heat_col[names(heat_col) %in% unique(as.vector(mat))]
ht <- Heatmap(
  mat, col = heat_col, name = "Alteration",
  top_annotation = top, bottom_annotation = bottom,
  rect_gp = gpar(col = "white", lwd = 0.4),
  row_names_gp = gpar(fontsize = 7, fontfamily = FONT, fontface = "italic"),
  column_labels = substr(colnames(mat), 9, 9),
  column_names_gp = gpar(fontsize = 7, fontfamily = FONT),
  cell_fun = function(j, i, x, y, w, h, fill) {
    if (rownames(mat)[i] == "KMT2C" && colnames(mat)[j] == "DH00023_T") {
      grid.rect(x - w * 0.18, y, width = w * 0.58, height = h * 0.82,
                gp = gpar(fill = "#2a78d6", col = "white", lwd = 0.4))
      grid.rect(x + w * 0.28, y, width = w * 0.32, height = h * 0.82,
                gp = gpar(fill = "#eb6834", col = "white", lwd = 0.4))
    }
  },
  column_names_rot = 0,
  column_title_gp = gpar(fontsize = 8, fontfamily = FONT),
  column_split = factor(substr(colnames(mat), 1, 7), levels = patients),
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_names_side = "left",
  heatmap_legend_param = list(
    title = "Alteration", at = names(heat_col)[names(heat_col) != "none"],
    title_gp = gpar(fontsize = 7, fontfamily = FONT),
    labels_gp = gpar(fontsize = 7, fontfamily = FONT)
  )
)
message(paste(capture.output(cn_mat), collapse = "\n"))

inputs <- c(
  file.path(root, "06_variants/somatic_filtered.tsv"),
  file.path(root, "07_cnv_purity/purity_ploidy.tsv"),
  file.path(root, "09_ccf_clonality/ccf.tsv"),
  file.path(root, "10_tmb_signatures/hrd_scar.tsv"),
  file.path(root, "12_tables/main/tmb.tsv"),
  file.path(root, "config/clinical_metadata.tsv"),
  file.path(root, "07_cnv_purity/cn_events.tsv")
)
save_fig(function() draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right"),
         "Fig2_oncoprint", width_mm = 183, height_mm = 150, inputs = inputs, script = "fig2_oncoprint.R",
         note = "Copy-number rows are cn_events.tsv (D044). DH00023 T KMT2C is a split cell. The BRCA2 row is labelled BRCA2 (VUS). WGD Yes is dark grey. Patient titles are 8 pt.")
