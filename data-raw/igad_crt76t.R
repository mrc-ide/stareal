# igad_crt76t.R
#
# Makes data/igad_crt76t.rda: survey counts of the chloroquine-resistance
# marker crt K76T in the eight IGAD countries (Djibouti, Eritrea, Ethiopia,
# Kenya, Somalia, South Sudan, Sudan, Uganda), from the public STAVE data
# set, for the examples.
#
# Source: STAVE object stave_data_2026.03.17.rds in the scrub repository
# (IDEELResearch/scrub, analysis/data-out/), pinned at commit 2d85d021
# (2026-03-17). scrub is MIT-licensed (copyright O. J. Watson); every survey
# keeps the reference of the study it came from. Counts come from STAVE's
# get_prevalence(), which counts every survey that genotyped the codon,
# zeros included. Needs the STAVE package (github.com/mrc-ide/STAVE).
#
# Surveys are selected by location, not by STAVE's country names, which
# are not reliable (e.g. "ETHIOPIA" and "Ethiopia"; "DRC" for surveys in
# Uganda): a survey is in if it lies within 10 km of a GADM 4.1 first-level
# unit of an IGAD country (large lakes cut out, as in tutorial 1), measured
# in the frame's planar projection, which is the rule
# areal_data() uses to snap points (and the rule of tutorial 1). `country`
# is the country of the nearest unit. GADM is downloaded through
# geodata::gadm() and is not shipped (its licence does not allow
# redistribution).
#
# Usage, from the package root: Rscript data-raw/igad_crt76t.R
#
# ------------------------------------------------------------------

suppressPackageStartupMessages(library(STAVE))

url <- "https://github.com/IDEELResearch/scrub/raw/2d85d021efe1c93dc55300610d944c9fbaf4af1a/analysis/data-out/stave_data_2026.03.17.rds"
tmp <- tempfile(fileext = ".rds")
download.file(url, tmp, mode = "wb", quiet = TRUE)
stopifnot(unname(tools::md5sum(tmp)) == "c1fd108b5665b88041217ced87f6a860")
stave <- suppressWarnings(readRDS(tmp))

source("data-raw/setup_crt76t.R")                         # the frame of tutorial 1, lakes cut out
units <- adm1                                            # already in fr$crs_planar

p <- stave$get_prevalence("crt:76:T")
p <- p[p$denominator > 0, ]
pts <- sf::st_transform(sf::st_as_sf(p, coords = c("longitude", "latitude"), crs = 4326), fr$crs_planar)
inside <- lengths(sf::st_is_within_distance(pts, units, dist = 10000)) > 0
p <- p[inside, ]
nearest <- sf::st_nearest_feature(pts[inside, ], units)         # the unit, and so the country

# typographic punctuation in the text fields (non-breaking spaces, curly
# quotes) as plain ASCII, so the data are ASCII; byte-wise, so it works in
# any locale
plain <- function(v) {
  utf8 <- c("\xc2\xa0", "\xe2\x80\x98", "\xe2\x80\x99", "\xe2\x80\x9c", "\xe2\x80\x9d")
  ascii <- c(" ", "'", "'", '"', '"')
  for (i in seq_along(utf8)) v <- gsub(utf8[i], ascii[i], v, fixed = TRUE, useBytes = TRUE)
  Encoding(v) <- "UTF-8"
  v
}
igad_crt76t <- data.frame(study_id = p$study_id, survey_id = p$survey_id, country = units$COUNTRY[nearest],
                          site = plain(p$site_name), longitude = p$longitude, latitude = p$latitude,
                          collection_start = p$collection_start, collection_end = p$collection_end,
                          collection_day = p$collection_day, positive = as.integer(p$numerator),
                          tested = as.integer(p$denominator), reference = plain(p$reference),
                          reference_year = as.integer(p$reference_year), pmid = as.character(p$PMID),
                          stringsAsFactors = FALSE)
igad_crt76t <- igad_crt76t[order(igad_crt76t$collection_day, igad_crt76t$country, igad_crt76t$site), ]
rownames(igad_crt76t) <- NULL
stopifnot(!anyNA(igad_crt76t[c("longitude", "latitude", "collection_day", "positive", "tested")]),
          all(igad_crt76t$positive <= igad_crt76t$tested))

message(sprintf("igad_crt76t: %d surveys from %d studies, %d positive of %d tested, %s to %s",
                nrow(igad_crt76t), length(unique(igad_crt76t$study_id)), sum(igad_crt76t$positive),
                sum(igad_crt76t$tested), min(igad_crt76t$collection_day), max(igad_crt76t$collection_day)))
print(table(igad_crt76t$country))

usethis::use_data(igad_crt76t, overwrite = TRUE, compress = "xz")
