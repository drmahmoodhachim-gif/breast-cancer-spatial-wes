# Publication style shared with scripts/figures/fig_style.py.
# Widths are journal columns in millimetres. Type is 7 pt Helvetica (Arial metric).

PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = "/home/chimeric/WES_TJN_project")
WIDTH_MM <- c(single = 89, onehalf = 120, double = 183)
COL_N <- "#0072B2"
COL_J <- "#E69F00"
COL_T <- "#882255"
COL_SHARED <- "#3d3d3a"
INK <- "#1a1a1a"
INK2 <- "#52514e"
FONT <- "Helvetica"

mm_in <- function(mm) mm / 25.4

git_hash <- function() {
  hash <- tryCatch(
    system2("git", c("-C", PROJECT_DIR, "rev-parse", "--short", "HEAD"), stdout = TRUE, stderr = FALSE),
    error = function(e) "no-git"
  )
  if (length(hash) == 0 || !nzchar(hash[[1]])) "no-git" else hash[[1]]
}

file_md5 <- function(path) {
  if (!file.exists(path) || file.info(path)$size > 2e8) return(NA_character_)
  unname(tools::md5sum(path))
}

theme_paper <- function(base_size = 7) {
  ggplot2::theme_bw(base_size = base_size, base_family = FONT) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = 0.2, colour = "#e4e3df"),
      axis.text = ggplot2::element_text(colour = INK2, size = 6),
      axis.title = ggplot2::element_text(colour = INK, size = 7),
      plot.title = ggplot2::element_text(colour = INK, size = 7, face = "bold"),
      legend.text = ggplot2::element_text(size = 6),
      legend.title = ggplot2::element_text(size = 7),
      legend.key.size = grid::unit(3, "mm"),
      strip.text = ggplot2::element_text(size = 7, colour = INK),
      strip.background = ggplot2::element_rect(fill = "white", colour = NA)
    )
}

panel_label <- function(label) {
  ggplot2::labs(tag = label) +
    ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 8, family = FONT),
                   plot.tag.position = c(0.01, 0.98))
}

#' Draw `expr` (a function, or an expression) to PDF, SVG and 600 dpi PNG.
save_fig <- function(draw, name, kind = c("main", "supplementary"),
                     width_mm = 183, height_mm = 120,
                     inputs = character(), note = "", script = "fig_style.R") {
  kind <- match.arg(kind)
  draw <- if (is.function(draw)) draw else function() eval(draw, envir = parent.frame())
  out <- file.path(PROJECT_DIR, "11_figures", kind)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  w <- mm_in(width_mm)
  h <- mm_in(height_mm)
  pdf_path <- file.path(out, paste0(name, ".pdf"))
  svg_path <- file.path(out, paste0(name, ".svg"))
  png_path <- file.path(out, paste0(name, ".png"))
  grDevices::cairo_pdf(pdf_path, width = w, height = h, family = FONT, pointsize = 7)
  draw()
  grDevices::dev.off()
  svglite::svglite(svg_path, width = w, height = h, pointsize = 7)
  draw()
  grDevices::dev.off()
  grDevices::png(png_path, width = width_mm / 25.4 * 600, height = height_mm / 25.4 * 600,
                 res = 600, type = "cairo", family = FONT, pointsize = 7)
  draw()
  grDevices::dev.off()
  rec <- list(
    output = file.path("11_figures", kind, paste0(name, ".pdf")),
    created = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
    script = script,
    git_commit = git_hash(),
    r = paste(R.version$major, R.version$minor, sep = "."),
    params = list(note = note, width_mm = width_mm, height_mm = height_mm, font = FONT, pointsize = 7),
    inputs = lapply(inputs, function(p) list(path = p, md5 = file_md5(p)))
  )
  # `%||%` may be absent
  jsonlite::write_json(rec, paste0(pdf_path, ".provenance.json"), auto_unbox = TRUE, pretty = TRUE, null = "null")
  message("[figure] ", kind, "/", name, " (.pdf .svg .png)")
  invisible(c(pdf_path, svg_path, png_path))
}

# Avoid a hard dependency on rlang for %||%
if (!exists("%||%", mode = "function")) {
  `%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || (is.character(a) && !nzchar(a))) b else a
}
