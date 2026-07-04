# Generate svyder hex logo + favicon assets from the chosen SVG source.
# Source of truth (edit here, then re-run): pkgdown/hex/option-1-profile.svg
# All rendering is local (rsvg + magick); no external favicon service is used.
suppressMessages({library(magick); library(rsvg)})

src <- "pkgdown/hex/option-1-profile.svg"
stopifnot(file.exists(src))

# --- man/figures/logo.* (pkgdown navbar + README) ---
dir.create("man/figures", showWarnings = FALSE, recursive = TRUE)
file.copy(src, "man/figures/logo.svg", overwrite = TRUE)
rsvg::rsvg_png(src, "man/figures/logo.png", width = 480)

# Remove stale logo variants from the previous design
for (f in c("logo-small.png", "logo-large.png", "favicon-16.png", "favicon-32.png"))
  if (file.exists(file.path("man/figures", f))) file.remove(file.path("man/figures", f))

# --- favicons: square, hex centred on a transparent background ---
dir.create("pkgdown/favicon", showWarnings = FALSE, recursive = TRUE)
file.copy(src, "pkgdown/favicon/favicon.svg", overwrite = TRUE)

sq <- function(px) {
  img <- magick::image_read_svg(src, height = px)         # keeps aspect (taller than wide)
  magick::image_extent(img, geometry = sprintf("%dx%d", px, px),
                       gravity = "center", color = "none")
}
# pkgdown's modern favicon set references favicon-96x96.png (not 16/32/48)
magick::image_write(sq(96),  "pkgdown/favicon/favicon-96x96.png",  format = "png")
magick::image_write(sq(180), "pkgdown/favicon/apple-touch-icon.png", format = "png")
# drop the older standalone sizes if present (only the multi-res .ico needs them)
for (f in c("favicon-16x16.png", "favicon-32x32.png", "favicon-48x48.png"))
  if (file.exists(file.path("pkgdown/favicon", f))) file.remove(file.path("pkgdown/favicon", f))
magick::image_write(sq(192), "pkgdown/favicon/web-app-manifest-192x192.png", format = "png")
magick::image_write(sq(512), "pkgdown/favicon/web-app-manifest-512x512.png", format = "png")

ico <- magick::image_join(sq(16), sq(32), sq(48))
magick::image_write(ico, "pkgdown/favicon/favicon.ico", format = "ico")

# --- web app manifest ---
writeLines(c(
  '{',
  '  "name": "svyder",',
  '  "short_name": "svyder",',
  '  "icons": [',
  '    { "src": "web-app-manifest-192x192.png", "sizes": "192x192", "type": "image/png" },',
  '    { "src": "web-app-manifest-512x512.png", "sizes": "512x512", "type": "image/png" }',
  '  ],',
  '  "theme_color": "#4393C3",',
  '  "background_color": "#ffffff",',
  '  "display": "standalone"',
  '}'
), "pkgdown/favicon/site.webmanifest")

cat("hex + favicon assets generated from", src, "\n")
