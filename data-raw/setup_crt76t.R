# setup_crt76t.R
#
# Sourced by the data-raw scripts that work on the crt K76T data of the
# first tutorial: the IGAD frame (GADM 4.1 first-level units, with the
# Natural Earth lakes over 1,000 km2 cut out), monthly steps from 1989 to
# 2024, and the package's igad_crt76t surveys. The same code as the
# tutorial's.
#
# GADM comes from geodata::gadm(). If the environment variable
# STAREAL_GADM_CACHE names a folder of files that geodata downloaded
# earlier (gadm41_<ISO3>_<level>_pk_low.rds), they are used instead of the
# server, which is sometimes down.
#
# ------------------------------------------------------------------

pkgload::load_all(".", quiet = TRUE)

gadm_path <- function() {
  path <- tempdir()
  cache <- Sys.getenv("STAREAL_GADM_CACHE")
  if (nzchar(cache)) {
    dir.create(file.path(path, "gadm"), showWarnings = FALSE)
    file.copy(list.files(cache, pattern = "^gadm41_.*\\.rds$", full.names = TRUE), file.path(path, "gadm"))
  }
  path
}

iso3 <- c("DJI", "ERI", "ETH", "KEN", "SOM", "SSD", "SDN", "UGA")
adm1 <- sf::st_as_sf(geodata::gadm(iso3, level = 1, path = gadm_path(), version = "4.1", resolution = 2))
lakes <- rnaturalearth::ne_download(scale = 10, type = "lakes", category = "physical",
                                    returnclass = "sf", destdir = tempdir())
crs <- "+proj=laea +lat_0=8 +lon_0=37 +units=m"
adm1 <- sf::st_make_valid(sf::st_transform(adm1, crs))
lakes <- sf::st_make_valid(sf::st_transform(lakes, crs))
lakes <- lakes[as.numeric(sf::st_area(lakes)) > 1e9, ]               # lakes over 1,000 km2
water <- sf::st_intersection(sf::st_union(lakes), sf::st_union(adm1))
adm1 <- suppressWarnings(sf::st_difference(adm1, water))
adm1 <- adm1[adm1$ENGTYPE_1 != "Water body", ]           # what is left of Uganda's lake units
fr <- areal_frame(adm1, id = "GID_1", groups = c(country = "COUNTRY"), link_within = "country", verbose = FALSE)
tf <- time_frame(1989, 2024, "month")
d <- areal_data(igad_crt76t, fr, tf, y = "positive", k = "tested",
                coords = c("longitude", "latitude"), date = "collection_day", verbose = FALSE)
