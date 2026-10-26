# make_logos.R -- round 2: star-themed hex-sticker candidates for stareal.
#
# In the package: Bob chose design C, the tessellated star with its time
# slices (2026-10-25). stareal_hex.png and stareal_hex.svg here are that
# design at full size; man/figures/logo.png is its 240 x 278 version. To
# redraw it, run in this folder:
#   LOGOS_ONLY=C Rscript make_logos.R
# then copy star_C_240x278.png to man/figures/logo.png (and star_C.png and
# star_C.svg over stareal_hex.*). The candidates' other letters are kept in
# the script for reference. The note about ../logos below refers to the
# first round, which is not kept.
#
# Run from anywhere:   Rscript make_logos.R
# Everything is written next to this script (no absolute paths).
#
# Geometry follows the hexb.in / hexSticker standard: a pointy-top hexagon
# 2 in vertex-to-vertex and 1.732 in (2 cos 30 deg) flat-to-flat. Drawing units
# are inches with the origin at the centre of the hexagon, so the hexagon's
# vertices are at (0, +-1) and its flat sides at x = +-0.866.
#
# Outputs, per design X:
#   star_X.png          520 x 600 px (1.733 x 2 in at 300 dpi), transparent corners
#   star_X.svg          same physical size (svglite); text is converted to outlines
#   star_X_240x278.png  pkgdown size (what usethis::use_logo() writes)
# plus contact_sheet.png and contact_sheet_small.png.
#
# Text is drawn as polygons read straight from the macOS system font files
# (TrueType 'glyf' outlines, laid out with textshaping, so kerning is kept).
# That makes the PNG and SVG identical and the SVG independent of the
# viewer's installed fonts.
#
# The helpers (hexagon, text outlines, rendering, edge checks) are the
# round-1 script's (../logos/make_logos.R). Changes: text_poly() can return
# the label in pieces (`parts`) and name_parts() places such pieces; the star
# helpers are new. Each run prints every design's clearance from the edge.
#
# Needs: ggplot2, sf, ragg, svglite, magick, png, colorspace, viridisLite,
# systemfonts, textshaping; macOS fonts Avenir Next, Charter, Futura,
# Gill Sans, Optima and Helvetica Neue. No data, no downloads.
# To redraw only some designs: LOGOS_ONLY=CE Rscript make_logos.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(sf)
})
sf_use_s2(FALSE)

args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grep("^--file=", args)])
out_dir <- if (length(script)) dirname(normalizePath(script)) else getwd()

# ---- geometry helpers -------------------------------------------------------

COS30 <- cos(pi / 6)

# Pointy-top hexagon of circumradius r centred at (cx, cy), as an sfc polygon.
hexagon <- function(r = 1, cx = 0, cy = 0, rot = 90) {
  a <- (rot + 60 * 0:5) * pi / 180
  p <- cbind(cx + r * cos(a), cy + r * sin(a))
  st_sfc(st_polygon(list(rbind(p, p[1, ]))))
}

# Circumradius of the hexagon whose sides are `d` inside the unit hexagon.
inset_r <- function(d) 1 - d / COS30

# Line width in inches -> ggplot linewidth.
lw <- function(inches) inches * 96 / .pt

rot2 <- function(a) matrix(c(cos(a), -sin(a), sin(a), cos(a)), 2)  # row-vector form

# Affine map of every coordinate, (x, y) -> (x, y) %*% M + shift, done
# elementwise. (sf's own `sfc * matrix` goes through BLAS, which on this Mac
# returns last-bit differences between a ring's first and last point, so the
# rings come out unclosed and GEOS refuses them.)
affine <- function(g, M = diag(2), shift = c(0, 0)) {
  f <- function(m) {
    if (is.matrix(m)) {
      x <- m[, 1]; y <- m[, 2]
      m[, 1] <- M[1, 1] * x + M[2, 1] * y + shift[1]
      m[, 2] <- M[1, 2] * x + M[2, 2] * y + shift[2]
    } else {
      x <- m[1]; y <- m[2]
      m[1] <- M[1, 1] * x + M[2, 1] * y + shift[1]
      m[2] <- M[1, 2] * x + M[2, 2] * y + shift[2]
    }
    m
  }
  tr <- function(s) {
    if (!is.list(s)) return(f(s))
    a <- attributes(s); s <- lapply(s, tr); attributes(s) <- a; s
  }
  st_sfc(lapply(st_geometry(g), tr), crs = st_crs(g))
}

gauss <- function(x, y, cx, cy, s) exp(-((x - cx)^2 + (y - cy)^2) / (2 * s^2))

# Index a palette vector by values in [0, 1].
pal_at <- function(pal, v) {
  v <- pmin(1, pmax(0, v))
  pal[round(v * (length(pal) - 1)) + 1]
}

# Lloyd-relaxed Voronoi tessellation of `region` with n cells.
voronoi_cells <- function(region, n, seed, iter = 4) {
  set.seed(seed)
  region <- st_union(region)
  pts <- st_sample(region, n, exact = TRUE)
  env <- st_as_sfc(st_bbox(st_buffer(region, 1)))
  for (k in 0:iter) {
    v <- st_collection_extract(st_voronoi(st_union(pts), env[[1]]), "POLYGON")
    v <- st_intersection(v, region)
    v <- v[as.numeric(st_area(v)) > 0]
    if (any(st_geometry_type(v) == "GEOMETRYCOLLECTION"))
      v <- st_collection_extract(v, "POLYGON")
    if (k < iter) pts <- st_centroid(v)
  }
  v
}

# ---- text as outlines -------------------------------------------------------

.font_cache <- new.env()

ttf_load <- function(path, index = 0L) {
  key <- paste(path, index)
  if (!is.null(.font_cache[[key]])) return(.font_cache[[key]])
  r <- readBin(path, "raw", n = file.size(path))
  u16 <- function(o) as.integer(r[o + 1L]) * 256L + as.integer(r[o + 2L])
  u32 <- function(o) as.numeric(u16(o)) * 65536 + u16(o + 2L)
  off <- if (identical(r[1:4], charToRaw("ttcf"))) u32(12L + 4L * index) else 0
  tab <- list()
  for (i in seq_len(u16(off + 4L)) - 1L) {
    rec <- off + 12L + 16L * i
    tab[[rawToChar(r[rec + 1:4])]] <- u32(rec + 8L)
  }
  if (is.null(tab$glyf)) stop("font has no TrueType outlines: ", path)
  n_glyph <- u16(tab$maxp + 4L)
  loca <- if (u16(tab$head + 50L) == 1L) {
    vapply(0:n_glyph, function(g) u32(tab$loca + 4L * g), 0)
  } else {
    vapply(0:n_glyph, function(g) 2 * u16(tab$loca + 2L * g), 0)
  }
  f <- list(r = r, glyf = tab$glyf, loca = loca, upem = u16(tab$head + 18L))
  assign(key, f, envir = .font_cache)
  f
}

# Contours of glyph g (0-based id): list of data frames (x, y, on-curve).
ttf_glyph <- function(f, g) {
  r <- f$r
  u8  <- function(o) as.integer(r[o + 1L])
  u16 <- function(o) as.integer(r[o + 1L]) * 256L + as.integer(r[o + 2L])
  s16 <- function(o) { v <- u16(o); if (v >= 32768L) v - 65536L else v }
  s8  <- function(o) { v <- u8(o); if (v >= 128L) v - 256L else v }
  if (f$loca[g + 2L] <= f$loca[g + 1L]) return(list())
  start <- f$glyf + f$loca[g + 1L]
  nc <- s16(start)
  if (nc >= 0L) {
    end_pts <- vapply(seq_len(nc) - 1L, function(i) u16(start + 10L + 2L * i), 0L)
    n <- end_pts[nc] + 1L
    p <- start + 10L + 2L * nc
    p <- p + 2L + u16(p)
    flags <- integer(n); i <- 1L
    while (i <= n) {
      fl <- u8(p); p <- p + 1L
      flags[i] <- fl; i <- i + 1L
      if (bitwAnd(fl, 8L)) {
        k <- u8(p); p <- p + 1L
        if (k > 0) { flags[i:(i + k - 1L)] <- fl; i <- i + k }
      }
    }
    coord <- function(short_bit, same_bit) {
      v <- numeric(n); cur <- 0
      for (k in seq_len(n)) {
        fl <- flags[k]
        if (bitwAnd(fl, short_bit)) {
          d <- u8(p); p <<- p + 1L
          cur <- cur + if (bitwAnd(fl, same_bit)) d else -d
        } else if (!bitwAnd(fl, same_bit)) {
          cur <- cur + s16(p); p <<- p + 2L
        }
        v[k] <- cur
      }
      v
    }
    x <- coord(2L, 16L)
    y <- coord(4L, 32L)
    on <- bitwAnd(flags, 1L) > 0
    starts <- c(1L, end_pts[-nc] + 2L)
    lapply(seq_len(nc), function(k) {
      ix <- starts[k]:(end_pts[k] + 1L)
      data.frame(x = x[ix], y = y[ix], on = on[ix])
    })
  } else {
    out <- list(); p <- start + 10L
    repeat {
      fl <- u16(p); gi <- u16(p + 2L); p <- p + 4L
      if (bitwAnd(fl, 1L)) { a1 <- s16(p); a2 <- s16(p + 2L); p <- p + 4L }
      else { a1 <- s8(p); a2 <- s8(p + 1L); p <- p + 2L }
      m <- diag(2)
      if (bitwAnd(fl, 8L)) { m <- diag(s16(p) / 16384, 2); p <- p + 2L }
      else if (bitwAnd(fl, 64L)) { m <- diag(c(s16(p), s16(p + 2L)) / 16384); p <- p + 4L }
      else if (bitwAnd(fl, 128L)) {
        m <- matrix(c(s16(p), s16(p + 2L), s16(p + 4L), s16(p + 6L)) / 16384, 2)
        p <- p + 8L
      }
      dx <- if (bitwAnd(fl, 2L)) a1 else 0
      dy <- if (bitwAnd(fl, 2L)) a2 else 0
      for (cc in ttf_glyph(f, gi)) {
        xy <- cbind(cc$x, cc$y) %*% m
        out[[length(out) + 1L]] <- data.frame(x = xy[, 1] + dx, y = xy[, 2] + dy, on = cc$on)
      }
      if (!bitwAnd(fl, 32L)) break
    }
    out
  }
}

# Flatten a quadratic TrueType contour into a closed ring.
flatten_contour <- function(cc, steps = 12L) {
  x <- cc$x; y <- cc$y; on <- cc$on; n <- length(x)
  X <- Y <- numeric(0); ON <- logical(0)
  for (i in seq_len(n)) {
    j <- if (i == n) 1L else i + 1L
    X <- c(X, x[i]); Y <- c(Y, y[i]); ON <- c(ON, on[i])
    if (!on[i] && !on[j]) {
      X <- c(X, (x[i] + x[j]) / 2); Y <- c(Y, (y[i] + y[j]) / 2); ON <- c(ON, TRUE)
    }
  }
  k <- which(ON)[1]
  idx <- c(k:length(X), seq_len(k - 1L))
  X <- X[idx]; Y <- Y[idx]; ON <- ON[idx]
  m <- length(X); tt <- seq_len(steps) / steps
  px <- X[1]; py <- Y[1]; i <- 2L
  while (i <= m + 1L) {
    ii <- if (i > m) 1L else i
    if (ON[ii]) {
      px <- c(px, X[ii]); py <- c(py, Y[ii]); i <- i + 1L
    } else {
      jj <- if (i + 1L > m) 1L else i + 1L
      x0 <- px[length(px)]; y0 <- py[length(py)]
      px <- c(px, (1 - tt)^2 * x0 + 2 * (1 - tt) * tt * X[ii] + tt^2 * X[jj])
      py <- c(py, (1 - tt)^2 * y0 + 2 * (1 - tt) * tt * Y[ii] + tt^2 * Y[jj])
      i <- i + 2L
    }
  }
  ring <- cbind(px, py)
  ring <- ring[c(TRUE, rowSums(abs(diff(ring))) > 1e-9), , drop = FALSE]
  if (any(ring[1, ] != ring[nrow(ring), ])) ring <- rbind(ring, ring[1, ])
  ring
}

font_face <- function(family, style) {
  fonts <- systemfonts::system_fonts()
  hit <- fonts[fonts$family == family & fonts$style == style, ]
  if (!nrow(hit)) stop("font not found: ", family, " ", style)
  list(path = hit$path[1], index = hit$index[1])
}

# A label as one sfc polygon (inches). The ink box is scaled to `width`
# (or the em to `size`) and placed so that (x, y) is at (hjust, vjust) of it.
# tracking: extra space between letters, in em.
# parts: if given, e.g. c(2, 2, 3), the label is returned as a list of that
# many pieces (2 glyphs, 2 glyphs, 3 glyphs), placed as the whole label is.
text_poly <- function(label, family, style = "Regular", width = NULL,
                      size = 0.3, x = 0, y = 0, hjust = 0.5, vjust = 0.5,
                      tracking = 0, parts = NULL) {
  ff <- font_face(family, style)
  sh <- textshaping::shape_text(label, path = ff$path, index = ff$index, size = 100)
  f <- ttf_load(sh$shape$font_path[1], sh$shape$font_index[1])
  s <- 100 / f$upem
  glyphs <- lapply(seq_len(nrow(sh$shape)), function(i) {
    cs <- Filter(function(cc) nrow(cc) >= 3, ttf_glyph(f, sh$shape$index[i]))
    if (!length(cs)) return(NULL)
    g <- NULL
    for (cc in cs) {
      r <- flatten_contour(cc)
      pg <- st_polygon(list(cbind(r[, 1] * s + sh$shape$x_offset[i] + (i - 1) * tracking * 100,
                                  r[, 2] * s)))
      g <- if (is.null(g)) pg else st_sym_difference(g, pg)
    }
    g
  })
  geom <- st_union(st_sfc(Filter(Negate(is.null), glyphs)))
  bb <- st_bbox(geom)
  sc <- if (is.null(width)) size / 100 else width / unname(bb["xmax"] - bb["xmin"])
  place <- function(g) affine(g, diag(sc, 2), -unname(bb[c("xmin", "ymin")]) * sc)
  bb2 <- st_bbox(place(geom))
  shift <- c(x - hjust * unname(bb2["xmax"]), y - vjust * unname(bb2["ymax"]))
  if (is.null(parts)) return(affine(place(geom), shift = shift))
  grp <- rep(seq_along(parts), parts)
  lapply(seq_along(parts), function(k)
    affine(place(st_union(st_sfc(Filter(Negate(is.null), glyphs[grp == k])))), shift = shift))
}

# ---- layers and rendering ---------------------------------------------------

# Polygons with identity colours. fill / col may be vectors (one per feature).
poly <- function(g, fill = NA, col = NA, lwd = 0, alpha = 1) {
  g <- st_geometry(g)
  d <- st_sf(fill = rep_len(fill, length(g)), col = rep_len(col, length(g)), geometry = g)
  geom_sf(data = d, aes(fill = fill, colour = col), linewidth = lw(lwd),
          alpha = alpha, lineend = "round", linejoin = "round")
}

# Polyline from a data frame (x, y, col), drawn segment by segment.
path <- function(df, lwd) {
  n <- nrow(df)
  seg <- data.frame(x = df$x[-n], y = df$y[-n], xend = df$x[-1], yend = df$y[-1],
                    col = df$col[-n])
  geom_segment(data = seg, aes(x = x, y = y, xend = xend, yend = yend, colour = col),
               linewidth = lw(lwd), lineend = "round")
}

# The standard frame: background hexagon, art, border ring on top.
sticker <- function(bg, border, d, art) {
  ring <- st_difference(hexagon(1), hexagon(inset_r(d)))
  c(list(poly(hexagon(1), fill = bg)), art, list(poly(ring, fill = border)))
}

render <- function(layers, file, w_px, h_px, res, type = c("png", "svg")) {
  type <- match.arg(type)
  w_in <- w_px / res; h_in <- h_px / res
  p <- ggplot() + layers +
    scale_fill_identity() + scale_colour_identity() +
    coord_sf(xlim = c(-w_in, w_in) / 2, ylim = c(-h_in, h_in) / 2,
             expand = FALSE, datum = NA) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0),
          plot.background = element_rect(fill = NA, colour = NA),
          panel.background = element_rect(fill = NA, colour = NA),
          legend.position = "none")
  if (type == "png") {
    ragg::agg_png(file, width = w_px, height = h_px, res = res, background = "transparent")
  } else {
    svglite::svglite(file, width = w_in, height = h_in, bg = "transparent")
  }
  print(p)
  invisible(dev.off())
}

# Distance (in) from the nearest point of g to the hexagon's outer edge;
# negative if g pokes outside. The border itself is ~0.04 in of this.
edge_margin <- function(g) {
  xy <- st_coordinates(st_union(st_geometry(g)))[, 1:2, drop = FALSE]
  a <- seq(0, 300, by = 60) * pi / 180
  COS30 - max(xy %*% rbind(cos(a), sin(a)))
}

margins <- new.env()
check_safe <- function(g, name, what, min_margin = 0.12) {
  m <- edge_margin(g)
  margins[[paste(name, what)]] <- m
  if (m < min_margin) warning(sprintf("design %s: %s is %.3f in from the edge (want >= %.2f)",
                                      name, what, m, min_margin), call. = FALSE)
  invisible(m)
}

# A transform that scales and moves art (keeping its shape) so that the art's
# bottom sits at `bottom`, centred left-right, as large as possible while
# staying `margin` in from the hexagon edge and within max_w x max_h.
# Returns a function to apply to every layer of that art (sfc or sf).
fitter <- function(g, bottom, margin, max_w = Inf, max_h = Inf) {
  g <- st_union(st_geometry(g))
  bb <- st_bbox(g)
  org <- unname(c((bb["xmin"] + bb["xmax"]) / 2, bb["ymin"]))
  make <- function(s) function(x) {
    if (inherits(x, "sf")) { st_geometry(x) <- make(s)(st_geometry(x)); return(x) }
    affine(x, diag(s, 2), c(0, bottom) - org * s)
  }
  lo <- 0
  hi <- min(max_w / unname(bb["xmax"] - bb["xmin"]), max_h / unname(bb["ymax"] - bb["ymin"]), 1e3)
  if (edge_margin(make(hi)(g)) >= margin) return(make(hi))
  for (i in 1:40) {
    mid <- (lo + hi) / 2
    if (edge_margin(make(mid)(g)) >= margin) lo <- mid else hi <- mid
  }
  make(lo)
}

# The package name, `width` wide, moved down (or up) until its nearest point
# is `clearance` in from the hexagon's outer edge.
name_poly <- function(family, style, width = 1.0, clearance = 0.13,
                      side = c("bottom", "top"), ...) {
  side <- match.arg(side)
  w0 <- text_poly("stareal", family, style, width = width, y = 0, ...)
  gap <- function(dy) edge_margin(affine(w0, shift = c(0, dy))) - clearance
  dy <- uniroot(gap, if (side == "bottom") c(-0.95, 0) else c(0, 0.95), tol = 1e-6)$root
  affine(w0, shift = c(0, dy))
}
# The same, split into pieces of `parts` glyphs (see text_poly()).
name_parts <- function(family, style, parts, width = 1.0, clearance = 0.13,
                       side = c("bottom", "top"), ...) {
  side <- match.arg(side)
  w0 <- text_poly("stareal", family, style, width = width, y = 0, ...)
  gap <- function(dy) edge_margin(affine(w0, shift = c(0, dy))) - clearance
  dy <- uniroot(gap, if (side == "bottom") c(-0.95, 0) else c(0, 0.95), tol = 1e-6)$root
  lapply(text_poly("stareal", family, style, width = width, y = 0, parts = parts, ...),
         affine, shift = c(0, dy))
}
ymax <- function(g) unname(st_bbox(g)["ymax"])
ymin <- function(g) unname(st_bbox(g)["ymin"])

# ---- star helpers -----------------------------------------------------------

# Vertices of an n-pointed star: outer radius R, inner radius r (r / R = 0.382
# is the classic pentagram outline), first tip at angle `rot` (degrees).
star_xy <- function(n = 5, R = 1, r = 0.382 * R, cx = 0, cy = 0, rot = 90) {
  a <- (rot + 180 / n * 0:(2 * n - 1)) * pi / 180
  rad <- rep(c(R, r), n)
  cbind(x = cx + rad * cos(a), y = cy + rad * sin(a))
}
star_poly <- function(...) {
  p <- star_xy(...)
  st_sfc(st_polygon(list(rbind(p, p[1, ]))))
}

# Round off the convex corners (a star's tips) of g with radius k.
round_tips <- function(g, k) st_buffer(st_buffer(g, -k, nQuadSegs = 12), k, nQuadSegs = 12)

# Discs of radius `rad` (recycled) centred at the rows of xy.
discs <- function(xy, rad, nq = 16) {
  rad <- rep_len(rad, nrow(xy))
  st_sfc(lapply(seq_len(nrow(xy)), function(i)
    st_buffer(st_point(as.numeric(xy[i, 1:2])), rad[i], nQuadSegs = nq)))
}

# A four-pointed sparkle (concave sides) of half-length s centred at (cx, cy).
sparkle <- function(cx, cy, s, waist = 0.22, rot = 0, n = 4) {
  st <- star_xy(n, R = s, r = waist * s, cx = cx, cy = cy, rot = 90 + rot)
  st_polygon(list(rbind(st, st[1, ])))
}

# Line segments (columns x, y, xend, yend) as polygons `w` inches wide.
segs_poly <- function(seg, w, cap = "ROUND") {
  st_buffer(st_sfc(lapply(seq_len(nrow(seg)), function(i)
    st_linestring(rbind(c(seg$x[i], seg$y[i]), c(seg$xend[i], seg$yend[i]))))),
    w / 2, endCapStyle = cap, nQuadSegs = 6)
}

# Edges of the Delaunay triangulation of the points xy, as a data frame of
# index pairs (i, j).
delaunay_edges <- function(xy) {
  tri <- st_collection_extract(st_triangulate(st_sfc(st_multipoint(as.matrix(xy[, 1:2])))), "POLYGON")
  key <- paste(round(xy[, 1], 9), round(xy[, 2], 9))
  e <- do.call(rbind, lapply(tri, function(t) {
    v <- match(paste(round(t[[1]][1:3, 1], 9), round(t[[1]][1:3, 2], 9)), key)
    rbind(v[1:2], v[2:3], v[c(3, 1)])
  }))
  e <- unique(t(apply(e, 1, sort)))
  data.frame(i = e[, 1], j = e[, 2])
}

# Mix colour(s) a towards b by fraction p (in sRGB).
mix <- function(a, b, p) {
  A <- grDevices::col2rgb(a) / 255; B <- grDevices::col2rgb(b) / 255
  p <- rep_len(p, ncol(A))
  M <- A * rep(1 - p, each = 3) + B * rep(p, each = 3)
  grDevices::rgb(M[1, ], M[2, ], M[3, ])
}

# ---- design A: constellation (dark, plasma) ---------------------------------
# A five-pointed star drawn as a constellation: nodes on and inside the star,
# joined by their neighbourhood graph. Most nodes are faint (unsurveyed
# cells); a few are bright sparkles (the surveyed ones). Node and edge
# colours follow one smooth surface. Faint background stars fill the sky.

design_A <- function() {
  bg <- "#0b1230"; border <- "#f89540"; ink <- "#f6f1ff"
  pal <- viridisLite::plasma(256)
  shape <- star_poly(5, R = 1, r = 0.47)
  # nodes: the 10 vertices, one point midway along each side, and interior
  # points spread by Lloyd relaxation
  v <- star_xy(5, R = 1, r = 0.47)
  mids <- (v + v[c(2:10, 1), ]) / 2
  inner <- st_coordinates(st_centroid(voronoi_cells(st_buffer(shape, -0.16), 14, seed = 4, iter = 6)))
  xy <- rbind(v, mids, inner[, 1:2])
  colnames(xy) <- c("x", "y")
  e <- delaunay_edges(xy)
  seg <- data.frame(x = xy[e$i, 1], y = xy[e$i, 2], xend = xy[e$j, 1], yend = xy[e$j, 2])
  keep <- st_covered_by(st_sfc(lapply(seq_len(nrow(seg)), function(k)
    st_linestring(rbind(c(seg$x[k], seg$y[k]), c(seg$xend[k], seg$yend[k]))))),
    st_buffer(shape, 1e-6), sparse = FALSE)[, 1]
  seg <- seg[keep, ]; e <- e[keep, ]
  field <- plogis(-0.9 + 1.3 * xy[, 1] + 0.7 * xy[, 2] +
                    2.2 * gauss(xy[, 1], xy[, 2], 0.55, 0.35, 0.45))
  # surveyed nodes: the nodes nearest four spread-out targets
  targets <- rbind(c(0.95, 0.31), c(-0.05, 0.62), c(-0.45, -0.45), c(0.12, -0.12))
  obs <- unique(apply(targets, 1, function(p) which.min((xy[, 1] - p[1])^2 + (xy[, 2] - p[2])^2)))

  word <- name_poly("Avenir Next", "Demi Bold Italic", width = 1.02)
  tf <- fitter(shape, bottom = ymax(word) + 0.12, margin = 0.165, max_w = 1.25)
  pts <- st_coordinates(tf(st_sfc(lapply(seq_len(nrow(xy)), function(k) st_point(xy[k, ])))))
  sseg <- data.frame(x = pts[e$i, 1], y = pts[e$i, 2], xend = pts[e$j, 1], yend = pts[e$j, 2])
  ecol <- pal_at(pal, 0.25 + 0.7 * (field[e$i] + field[e$j]) / 2)
  # the outline: vertex k -> midpoint 10 + k -> vertex k + 1
  rim_pairs <- rbind(cbind(1:10, 11:20), cbind(11:20, c(2:10, 1)))
  rim <- paste(e$i, e$j) %in% paste(pmin(rim_pairs[, 1], rim_pairs[, 2]), pmax(rim_pairs[, 1], rim_pairs[, 2]))
  edges <- segs_poly(sseg, ifelse(rim, 0.016, 0.009))
  node_col <- pal_at(pal, 0.25 + 0.7 * field)
  faint <- setdiff(seq_len(nrow(pts)), obs)
  dots <- discs(pts[faint, , drop = FALSE], 0.021)
  halo1 <- discs(pts[obs, , drop = FALSE], 0.085, 24)
  halo2 <- discs(pts[obs, , drop = FALSE], 0.055, 24)
  sparks <- st_sfc(lapply(obs, function(k) sparkle(pts[k, 1], pts[k, 2], 0.075, waist = 0.2)))
  cores <- discs(pts[obs, , drop = FALSE], 0.022)

  # background sky: small faint stars away from the constellation and name
  set.seed(9)
  keep_out <- st_union(c(st_buffer(tf(shape), 0.06), st_buffer(st_as_sfc(st_bbox(word)), 0.04)))
  cand <- st_sample(hexagon(inset_r(0.16)), 120)
  cand <- cand[!st_intersects(cand, keep_out, sparse = FALSE)[, 1]]
  cand <- cand[seq_len(min(46, length(cand)))]
  sky_r <- runif(length(cand), 0.004, 0.011)
  sky <- discs(st_coordinates(cand), sky_r, 8)
  sky_a <- runif(length(cand), 0.25, 0.75)

  check_safe(c(halo1, dots, edges), "A", "art", 0.12)
  check_safe(sky, "A", "sky", 0.12)
  check_safe(word, "A", "name")
  art <- list(
    poly(sky, fill = "#cfd8ff", alpha = sky_a),
    poly(edges, fill = ecol, alpha = 0.55),
    poly(halo1, fill = pal_at(pal, 0.25 + 0.7 * field[obs]), alpha = 0.18),
    poly(halo2, fill = pal_at(pal, 0.25 + 0.7 * field[obs]), alpha = 0.30),
    poly(dots, fill = node_col[faint]),
    poly(sparks, fill = "#fff6dc"),
    poly(cores, fill = "#ffffff"),
    poly(word, fill = ink)
  )
  sticker(bg, border, 0.04, art)
}

# ---- design B: a star of cells (light, magma) -------------------------------
# A five-pointed star tessellated into Voronoi cells, each filled from one
# smooth surface. White dots mark the few surveyed cells.

design_B <- function() {
  bg <- "#f8f4fa"; border <- "#2d1152"; ink <- "#2d1152"
  pal <- viridisLite::magma(256)
  shape <- round_tips(star_poly(5, R = 1, r = 0.48), 0.035)
  cells <- voronoi_cells(shape, 40, seed = 3, iter = 8)
  xy <- st_coordinates(st_point_on_surface(cells))
  field <- plogis(-0.4 + 1.1 * xy[, 1] + 0.9 * xy[, 2] +
                    2.4 * gauss(xy[, 1], xy[, 2], 0.45, 0.25, 0.35) -
                    1.8 * gauss(xy[, 1], xy[, 2], -0.55, -0.35, 0.4))
  targets <- rbind(c(0.02, 0.62), c(0.72, 0.2), c(-0.62, 0.12), c(0.38, -0.5), c(-0.12, -0.2))
  obs <- unique(apply(targets, 1, function(p) which.min((xy[, 1] - p[1])^2 + (xy[, 2] - p[2])^2)))

  word <- name_poly("Charter", "Bold Italic", width = 1.0)
  tf <- fitter(shape, bottom = ymax(word) + 0.10, margin = 0.15, max_w = 1.3)
  cells <- tf(cells)
  dots <- discs(st_coordinates(st_point_on_surface(cells[obs]))[, 1:2], 0.029)
  check_safe(cells, "B", "art", 0.12)
  check_safe(word, "B", "name")
  art <- list(
    poly(cells, fill = pal_at(pal, 0.15 + 0.73 * field), col = bg, lwd = 0.014),
    poly(dots, fill = "#ffffff", col = border, lwd = 0.009),
    poly(word, fill = ink)
  )
  sticker(bg, border, 0.04, art)
}

# ---- design C: a stack of time slices (dark, inferno) -----------------------
# The star as a stack of slices, one per time step, the oldest at the back.
# The front slice (the latest step) is tessellated into cells and coloured by
# the estimated surface; the slices behind it step down the palette.

design_C <- function() {
  bg <- "#140a24"; border <- "#fca50a"; ink <- "#fff4d6"
  pal <- viridisLite::inferno(256)
  shape <- round_tips(star_poly(5, R = 1, r = 0.48), 0.035)
  cells <- voronoi_cells(shape, 30, seed = 8, iter = 8)
  xy <- st_coordinates(st_point_on_surface(cells))
  field <- plogis(-0.6 + 1.6 * xy[, 1] + 1.2 * xy[, 2] +
                    2.2 * gauss(xy[, 1], xy[, 2], 0.45, 0.3, 0.35))
  n_back <- 4
  step <- c(-0.085, 0.085)                      # each older slice sits up and to the left
  back <- lapply(n_back:1, function(k) affine(shape, shift = k * step))
  word <- name_poly("Avenir Next", "Heavy Italic", width = 1.04)
  tf <- fitter(st_union(do.call(c, c(back, list(shape)))), bottom = ymax(word) + 0.10,
               margin = 0.15, max_w = 1.35)
  back <- lapply(back, tf); cells <- tf(cells); front <- tf(shape)
  slice_col <- pal_at(pal, c(0.12, 0.19, 0.26, 0.33))
  art <- list()
  for (k in seq_along(back))
    art <- c(art, list(poly(back[[k]], fill = slice_col[k], col = bg, lwd = 0.014)))
  art <- c(art, list(
    poly(front, fill = bg, col = bg, lwd = 0.014),
    poly(cells, fill = pal_at(pal, 0.50 + 0.48 * field), col = bg, lwd = 0.009)
  ))
  check_safe(st_union(do.call(c, back)), "C", "art", 0.12)
  check_safe(word, "C", "name")
  sticker(bg, border, 0.04, c(art, list(poly(word, fill = ink))))
}

# ---- design D: a star on a lattice of cells (light, cold-hot) ---------------
# A lattice of hexagonal cells, the same shape as the sticker. The cells
# inside a five-pointed star carry the estimated surface, cold to hot; the
# rest of the lattice is empty. White dots mark the surveyed cells.

design_D <- function() {
  bg <- "#fbfaf7"; border <- "#1b2a49"; ink <- "#1b2a49"
  empty_fill <- "#eef0f4"
  # cold to hot in HCL: blue, purple, pink, orange, yellow
  pal <- colorspace::sequential_hcl(256, h = c(265, 420), c = c(70, 90, 80), l = c(28, 88),
                                    power = c(0.9, 1.1))
  rc <- 0.047                                  # cell circumradius (in)
  word <- name_poly("Gill Sans", "Bold", width = 1.0)
  tf <- fitter(star_poly(5, R = 1, r = 0.5), bottom = ymax(word) + 0.11,
               margin = 0.17, max_w = 1.3)
  shape <- tf(star_poly(5, R = 1, r = 0.5))
  ctr <- st_coordinates(tf(st_sfc(st_point(c(0, 0)))))[1, 1:2]
  # pointy-top lattice with a cell centred on the star's centre
  qr <- expand.grid(q = -14:14, r = -14:14)
  cx <- ctr[1] + sqrt(3) * rc * (qr$q + qr$r / 2)
  cy <- ctr[2] + 1.5 * rc * qr$r
  pts <- st_sfc(lapply(seq_along(cx), function(i) st_point(c(cx[i], cy[i]))))
  inside <- st_intersects(pts, shape, sparse = FALSE)[, 1]
  cell <- function(i, f = 0.9) hexagon(rc * f, cx[i], cy[i])
  # the empty lattice: cells wholly inside a hexagon inset from the edge, and
  # clear of the name
  room <- st_intersection(hexagon(inset_r(0.17)),
                          st_as_sfc(st_bbox(c(xmin = -1, ymin = ymax(word) + 0.07, xmax = 1, ymax = 1))))
  all_cells <- do.call(c, lapply(seq_along(cx), cell))
  in_room <- st_covered_by(all_cells, room, sparse = FALSE)[, 1]
  u <- (cx - ctr[1]) / 0.6; w <- (cy - ctr[2]) / 0.6
  field <- plogis(-0.3 + 2.4 * u + 1.1 * w + 2.0 * gauss(u, w, 0.45, 0.35, 0.35) -
                    2.4 * gauss(u, w, -0.55, -0.3, 0.4))
  idx <- which(inside)
  targets <- rbind(c(0, 0.75), c(0.72, 0.2), c(-0.7, 0.2), c(0.42, -0.62), c(-0.42, -0.62), c(0, 0))
  obs <- unique(apply(targets, 1, function(p) idx[which.min((u[idx] - p[1])^2 + (w[idx] - p[2])^2)]))
  dots <- discs(cbind(cx[obs], cy[obs]), 0.025)
  check_safe(all_cells[in_room | inside], "D", "art", 0.12)
  check_safe(word, "D", "name")
  art <- list(
    poly(all_cells[in_room & !inside], fill = empty_fill),
    poly(all_cells[inside], fill = pal_at(pal, 0.02 + 0.96 * field[inside])),
    poly(dots, fill = "#ffffff", col = ink, lwd = 0.008),
    poly(word, fill = ink)
  )
  sticker(bg, border, 0.04, art)
}

# ---- design E: star x areal, overprinted (light, blue/pink) -----------------
# Two-ink "risograph" print. A blue star overlaps a pink patch of four cells;
# where they overlap the inks multiply to deep indigo. The name is printed the same way:
# "star" in blue, "areal" in pink, and the letters they share, "ar", in the
# overlap colour.

design_E <- function() {
  bg <- "#f5f0e6"
  blue <- "#3255a4"; pink <- "#ff48b0"
  over <- grDevices::rgb(t(grDevices::col2rgb(blue) * grDevices::col2rgb(pink) / 255^2))
  border <- over
  star <- star_poly(5, R = 1, r = 0.5, cx = -0.25, cy = 0.12)
  # the "areal" part: an irregular hexagon cut into four cells
  ang <- c(10, 68, 130, 192, 248, 312) * pi / 180
  rad <- c(0.98, 0.92, 1.0, 0.9, 1.0, 0.94) * 0.85
  patchxy <- cbind(0.5 + rad * cos(ang), -0.12 + rad * sin(ang))
  patch <- st_sfc(st_polygon(list(rbind(patchxy, patchxy[1, ]))))
  cells <- voronoi_cells(patch, 4, seed = 7, iter = 10)
  parts <- name_parts("Futura", "Bold", parts = c(2, 2, 3), width = 1.12)
  word <- st_union(do.call(c, parts))
  tf <- fitter(st_union(c(star, patch)), bottom = ymax(word) + 0.11, margin = 0.15, max_w = 1.32)
  star <- tf(star)
  cells <- st_buffer(tf(cells), -0.012, joinStyle = "MITRE")
  both <- st_intersection(star, st_union(cells))
  check_safe(c(star, cells), "E", "art", 0.12)
  check_safe(word, "E", "name")
  art <- list(
    poly(cells, fill = pink), poly(star, fill = blue), poly(both, fill = over),
    poly(parts[[1]], fill = blue), poly(parts[[2]], fill = over), poly(parts[[3]], fill = pink)
  )
  sticker(bg, border, 0.04, art)
}

# ---- design F: ripples from the data (light, plasma) ------------------------
# A five-pointed star built from nested star-shaped bands around one surveyed
# cell at its centre, like ripples. Each band is cut into ten cells (the
# halves of the points). Hot at the data, cooling band by band away from it.

design_F <- function() {
  bg <- "#fdf8f1"; border <- "#1f1a5c"; ink <- "#1f1a5c"
  pal <- rev(viridisLite::plasma(256))          # hot (yellow) to cold (blue)
  v <- star_xy(5, R = 1, r = 0.5)
  outline <- function(s) st_polygon(list(rbind(s * v, s * v[1, ])))
  bands <- c(0.22, 0.47, 0.74, 1)               # the core, then three bands
  # ten wedges, split along the rays from the centre to every vertex
  wa <- pi / 2 + seq(0, 2 * pi, length.out = 11)
  pieces <- list(); band <- c()
  for (j in 2:length(bands)) {
    ring <- st_difference(outline(bands[j]), outline(bands[j - 1]))
    for (k in 1:10) {
      wedge <- st_polygon(list(rbind(c(0, 0), 3 * c(cos(wa[k]), sin(wa[k])),
                                     3 * c(cos(wa[k + 1]), sin(wa[k + 1])), c(0, 0))))
      pieces[[length(pieces) + 1]] <- st_cast(st_intersection(ring, wedge), "POLYGON")
      band <- c(band, j - 1)
    }
  }
  cells <- st_sfc(pieces)
  core <- st_sfc(outline(bands[1]))
  word <- name_poly("Optima", "Bold Italic", width = 1.0)
  tf <- fitter(st_sfc(outline(1)), bottom = ymax(word) + 0.09, margin = 0.15, max_w = 1.4)
  cells <- tf(cells); core <- tf(core)
  ctr <- st_coordinates(tf(st_sfc(st_point(c(0, 0)))))[1, 1:2]
  dot <- discs(rbind(ctr), 0.034)
  check_safe(cells, "F", "art", 0.12)
  check_safe(word, "F", "name")
  art <- list(
    poly(cells, fill = pal_at(pal, 0.25 + 0.67 * (band - 1) / 2), col = bg, lwd = 0.013),
    poly(core, fill = pal_at(pal, 0.1), col = bg, lwd = 0.013),
    poly(dot, fill = "#ffffff", col = ink, lwd = 0.009),
    poly(word, fill = ink)
  )
  sticker(bg, border, 0.04, art)
}

# ---- design G: a faceted star (dark, magma) ---------------------------------
# A bevelled five-pointed star, its ten facets cut into a mesh of small
# triangular cells. Each cell takes its colour from one smooth surface,
# shaded by its facet's angle to the light. White dots mark the
# surveyed cells.

design_G <- function() {
  bg <- "#161c2c"; border <- "#fd9a6a"; ink <- "#f2f5fb"
  pal <- viridisLite::magma(256)
  v <- star_xy(5, R = 1, r = 0.5)               # tips are odd rows, notches even rows
  m <- 4                                         # subdivisions per facet side
  h <- 0.55                                      # height of the centre (for shading)
  light <- c(-0.5, 0.6, 0.62); light <- light / sqrt(sum(light^2))
  tris <- list(); shade <- c(); cen <- NULL
  for (k in 1:10) {
    A <- c(0, 0); B <- v[k, ]; C <- v[k %% 10 + 1, ]
    # facet normal from 3D points (A raised by h)
    P <- c(A, h); Q <- c(B, 0); R <- c(C, 0)
    nv <- c((Q[2] - P[2]) * (R[3] - P[3]) - (Q[3] - P[3]) * (R[2] - P[2]),
            (Q[3] - P[3]) * (R[1] - P[1]) - (Q[1] - P[1]) * (R[3] - P[3]),
            (Q[1] - P[1]) * (R[2] - P[2]) - (Q[2] - P[2]) * (R[1] - P[1]))
    nv <- nv / sqrt(sum(nv^2)); if (nv[3] < 0) nv <- -nv
    lam <- sum(nv * light)
    pt <- function(i, j) A + (B - A) * i / m + (C - A) * j / m
    for (i in 0:(m - 1)) for (j in 0:(m - 1 - i)) {
      p <- rbind(pt(i, j), pt(i + 1, j), pt(i, j + 1))
      tris[[length(tris) + 1]] <- st_polygon(list(rbind(p, p[1, ])))
      shade <- c(shade, lam); cen <- rbind(cen, colMeans(p))
      if (i + j < m - 1) {
        p <- rbind(pt(i + 1, j), pt(i + 1, j + 1), pt(i, j + 1))
        tris[[length(tris) + 1]] <- st_polygon(list(rbind(p, p[1, ])))
        shade <- c(shade, lam); cen <- rbind(cen, colMeans(p))
      }
    }
  }
  tris <- st_sfc(tris)
  field <- plogis(-0.2 + 1.6 * cen[, 1] + 1.0 * cen[, 2] +
                    2.2 * gauss(cen[, 1], cen[, 2], 0.5, 0.3, 0.35) -
                    1.8 * gauss(cen[, 1], cen[, 2], -0.5, -0.4, 0.4))
  base <- pal_at(pal, 0.28 + 0.70 * field)
  s <- (shade - mean(range(shade))) / diff(range(shade))   # -0.5 .. 0.5
  cols <- ifelse(s >= 0, colorspace::lighten(base, 0.35 * s, space = "HLS"),
                 colorspace::darken(base, -0.5 * s, space = "HLS"))
  targets <- rbind(c(0, 0.75), c(0.68, 0.2), c(-0.6, 0.18), c(0.35, -0.5), c(-0.38, -0.5))
  obs <- unique(apply(targets, 1, function(p) which.min((cen[, 1] - p[1])^2 + (cen[, 2] - p[2])^2)))

  word <- name_poly("Helvetica Neue", "Bold Italic", width = 1.0)
  tf <- fitter(st_union(tris), bottom = ymax(word) + 0.10, margin = 0.15, max_w = 1.3)
  tris <- tf(tris)
  dots <- discs(st_coordinates(st_centroid(tris[obs]))[, 1:2], 0.025)
  check_safe(tris, "G", "art", 0.12)
  check_safe(word, "G", "name")
  art <- list(
    poly(tris, fill = cols, col = cols, lwd = 0.002),
    poly(dots, fill = "#ffffff", col = bg, lwd = 0.007),
    poly(word, fill = ink)
  )
  sticker(bg, border, 0.04, art)
}

# ---- build ------------------------------------------------------------------

all_designs <- mget(sprintf("design_%s", LETTERS), ifnotfound = list(NULL))
all_designs <- Filter(Negate(is.null), all_designs)
names(all_designs) <- sub("design_", "", names(all_designs))
designs <- all_designs
only <- Sys.getenv("LOGOS_ONLY")
if (nzchar(only)) designs <- designs[strsplit(only, "")[[1]]]

for (nm in names(designs)) {
  message("design ", nm)
  layers <- designs[[nm]]()
  render(layers, file.path(out_dir, sprintf("star_%s.png", nm)), 520, 600, 300, "png")
  render(layers, file.path(out_dir, sprintf("star_%s.svg", nm)), 520, 600, 300, "svg")
  render(layers, file.path(out_dir, sprintf("star_%s_240x278.png", nm)), 240, 278, 139, "png")
}

# ---- verify geometry --------------------------------------------------------
# Alpha along the edges of each PNG: the hexagon's top and bottom vertices
# should just touch the top and bottom rows, its flat sides the left and right
# columns, and the corners should be fully transparent.

verify <- function(file) {
  a <- png::readPNG(file)[, , 4]
  h <- nrow(a); w <- ncol(a); mid_c <- c(floor((w + 1) / 2), ceiling((w + 1) / 2))
  data.frame(
    file = basename(file), width = w, height = h,
    top = max(a[1, mid_c]), bottom = max(a[h, mid_c]),
    left = max(a[floor(h / 2), 1]), right = max(a[floor(h / 2), w]),
    corners = max(a[1, 1], a[1, w], a[h, 1], a[h, w]),
    outside = max(a[1:3, 1:40], a[1:3, (w - 39):w])
  )
}
L <- names(all_designs)
files <- file.path(out_dir, c(sprintf("star_%s.png", L), sprintf("star_%s_240x278.png", L)))
files <- files[file.exists(files)]
print(do.call(rbind, lapply(files, verify)), digits = 2, row.names = FALSE)
m <- unlist(as.list(margins)); m <- m[order(names(m))]
cat("\nclearance from the hexagon's outer edge (in; border is ~0.04):\n")
print(round(m, 3))

# ---- contact sheets ---------------------------------------------------------

sheet <- function(suffix, out, gap, label_px, ncol) {
  imgs <- lapply(L, function(l)
    magick::image_read(file.path(out_dir, sprintf("star_%s%s.png", l, suffix))))
  n <- length(imgs)
  w <- magick::image_info(imgs[[1]])$width; h <- magick::image_info(imgs[[1]])$height
  nrow <- ceiling(n / ncol)
  W <- ncol * w + (ncol + 1) * gap; H <- nrow * (h + label_px * 2) + (nrow + 1) * gap
  canvas <- magick::image_blank(W, H, "#e4e4e4")
  for (i in seq_len(n)) {
    r <- (i - 1) %/% ncol; cc <- (i - 1) %% ncol
    in_row <- min(ncol, n - r * ncol)              # centre a short last row
    x0 <- (W - (in_row * w + (in_row - 1) * gap)) / 2
    x <- round(x0 + cc * (w + gap)); y <- gap + r * (h + label_px * 2 + gap)
    canvas <- magick::image_composite(canvas, imgs[[i]], offset = sprintf("+%d+%d", x, y))
    canvas <- magick::image_annotate(canvas, L[i], size = label_px, color = "#333333",
                                     font = "Helvetica", weight = 700, gravity = "northwest",
                                     location = sprintf("+%d+%d", round(x + w / 2 - label_px * 0.35),
                                                        round(y + h + label_px * 0.4)))
  }
  magick::image_write(canvas, file.path(out_dir, out), format = "png")
}
if (all(file.exists(file.path(out_dir, sprintf("star_%s.png", L))))) {
  sheet("", "contact_sheet.png", 40, 44, ncol = ceiling(length(L) / 2))
  sheet("_240x278", "contact_sheet_small.png", 24, 24, ncol = length(L))
}
