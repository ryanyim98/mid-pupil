# Resolve figures/ under repo root whether getwd() is repo root, step2_pooled_analyses,
# analyses_by_sections, or this helper/.
mid_pupil_figures_dir <- function() {
  tails <- c(
    "figures",
    file.path("..", "figures"),
    file.path("..", "..", "figures"),
    file.path("..", "..", "..", "figures")
  )
  wd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  for (rel in tails) {
    cand <- normalizePath(file.path(wd, rel), winslash = "/", mustWork = FALSE)
    if (!is.na(cand) && dir.exists(cand)) return(cand)
  }
  repo_figures <- file.path(
    normalizePath(file.path(wd, "..", "..", ".."), winslash = "/", mustWork = FALSE),
    "figures"
  )
  if (!dir.exists(repo_figures)) dir.create(repo_figures, recursive = TRUE, showWarnings = FALSE)
  repo_figures
}

mid_pupil_figure_path <- function(name, ext = "png") {
  file.path(mid_pupil_figures_dir(), paste0(name, ".", ext))
}

# Save a ggplot (or patchwork) as PNG and TIFF with a systematic base name.
save_figure <- function(name, plot = ggplot2::last_plot(), dpi = 300, ...) {
  png_path <- mid_pupil_figure_path(name, "png")
  tiff_path <- mid_pupil_figure_path(name, "tiff")
  ggplot2::ggsave(png_path, plot = plot, dpi = dpi, ...)
  ggplot2::ggsave(tiff_path, plot = plot, dpi = dpi, ...)
  invisible(list(png = png_path, tiff = tiff_path))
}

# Save a magick image object as PNG and TIFF (e.g. after panel labels).
save_figure_image <- function(name, img) {
  png_path <- mid_pupil_figure_path(name, "png")
  tiff_path <- mid_pupil_figure_path(name, "tiff")
  magick::image_write(img, png_path)
  magick::image_write(img, tiff_path)
  invisible(list(png = png_path, tiff = tiff_path))
}
