# igad_adm1.R
#
# Makes data/igad_adm1.rda: first-level administrative units of the eight
# IGAD countries (Djibouti, Eritrea, Ethiopia, Kenya, Somalia, South Sudan,
# Sudan, Uganda), for the examples and tutorials.
#
# Source: geoBoundaries gbOpen, release commit 9469f09 (build of
# 2023-12-12), the simplified geometries. Licences differ by country and are
# recorded in the `licence` column and in ?igad_adm1:
#   DJI ADM1  Public Domain (Natural Earth)
#   ERI ADM1  ODbL 1.0 (OpenStreetMap, Wambacher)
#   ETH ADM1  CC BY 4.0 (Open Africa, Code for Ethiopia)
#   KEN ADM1  Public Domain (RCMRD GeoPortal, Africa GeoPortal)
#   SDN ADM1  CC BY 3.0 IGO (IMWG, OCHA ROSEA)
#   SOM ADM1  ODbL 1.0 (OpenStreetMap, Wambacher)
#   SSD ADM1  CC BY 3.0 IGO (IMWG, NBS, OCHA South Sudan)
#   UGA ADM2  Public Domain (geoBoundaries, Wikimedia Commons)
# Uganda is taken at ADM2 (151 districts, 2006 boundaries): geoBoundaries'
# ADM1 for Uganda is 4 statistical regions, far coarser than the other
# countries' first-level units, while districts are Uganda's working
# administrative units.
#
# The national files come from different sources, so they overlap at some
# borders (by up to ~9,000 km2 at the Ilemi triangle). Each overlap is
# given to one country, in the order KEN, UGA, ETH, ERI, DJI, SOM, SDN, SSD
# (each country loses what the countries before it already cover), so that
# every point lies in at most one unit. They also leave gaps at borders, so
# frames built from them need areal_frame(adjacency_snap = ~2000) to make
# neighbouring units across borders adjacent. The boundaries are for
# illustration and imply no position on disputed borders.
#
# Usage, from the package root: Rscript data-raw/igad_adm1.R
#
# ------------------------------------------------------------------

suppressMessages(library(sf))
sf_use_s2(FALSE)

release <- "https://github.com/wmgeolab/geoBoundaries/raw/9469f09/releaseData/gbOpen"
countries <- data.frame(
  iso3 = c("KEN", "UGA", "ETH", "ERI", "DJI", "SOM", "SDN", "SSD"),       # also the order overlaps are resolved in
  level = c("ADM1", "ADM2", "ADM1", "ADM1", "ADM1", "ADM1", "ADM1", "ADM1"),
  country = c("Kenya", "Uganda", "Ethiopia", "Eritrea", "Djibouti", "Somalia", "Sudan", "South Sudan"),
  licence = c("Public Domain", "Public Domain", "CC BY 4.0", "ODbL 1.0", "Public Domain", "ODbL 1.0",
              "CC BY 3.0 IGO", "CC BY 3.0 IGO"),
  stringsAsFactors = FALSE)

tmp <- tempfile(fileext = ".geojson")
read_country <- function(i) {
  cc <- countries[i, ]
  download.file(sprintf("%s/%s/%s/geoBoundaries-%s-%s_simplified.geojson", release, cc$iso3, cc$level,
                        cc$iso3, cc$level), tmp, quiet = TRUE)
  p <- st_make_valid(st_read(tmp, quiet = TRUE))
  p <- p[order(p$shapeName), ]
  data.frame(id = sprintf("%s_%03d", cc$iso3, seq_len(nrow(p))), name = p$shapeName, iso3 = cc$iso3,
             country = cc$country, level = cc$level, licence = cc$licence, gb_id = p$shapeID,
             stringsAsFactors = FALSE) |>
    st_sf(geometry = st_geometry(p))
}
parts <- lapply(seq_len(nrow(countries)), read_country)

# resolve overlaps: each country loses what the countries before it cover
taken <- NULL
for (i in seq_along(parts)) {
  if (!is.null(taken)) {
    g <- st_difference(st_geometry(parts[[i]]), taken)
    g <- st_collection_extract(st_make_valid(g), "POLYGON")
    st_geometry(parts[[i]]) <- st_cast(g, "MULTIPOLYGON")
  }
  outline <- st_union(st_geometry(parts[[i]]))
  taken <- if (is.null(taken)) outline else st_union(taken, outline)
}
igad_adm1 <- do.call(rbind, parts)
igad_adm1 <- igad_adm1[order(igad_adm1$iso3, igad_adm1$id), ]
st_geometry(igad_adm1) <- st_cast(st_geometry(igad_adm1), "MULTIPOLYGON")
rownames(igad_adm1) <- NULL
stopifnot(all(st_is_valid(igad_adm1)), !anyDuplicated(igad_adm1$id), nrow(igad_adm1) == 268)

# no overlaps left (beyond rounding)
ov <- st_overlaps(igad_adm1)
area_ov <- sum(vapply(seq_along(ov), function(i) {
  j <- ov[[i]][ov[[i]] > i]
  if (!length(j)) return(0)
  sum(as.numeric(st_area(st_transform(st_intersection(igad_adm1[i, ], igad_adm1[j, ]), "+proj=laea +lat_0=5 +lon_0=35"))))
}, 1)) / 1e6
message(sprintf("igad_adm1: %d units; remaining overlap %.2f km2", nrow(igad_adm1), area_ov))

usethis::use_data(igad_adm1, overwrite = TRUE, compress = "xz")
