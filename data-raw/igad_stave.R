# igad_stave.R
#
# Makes data/igad_crt76t.rda and data/igad_k13_622i.rda: survey counts of
# two antimalarial drug-resistance markers in the eight IGAD countries,
# from the public STAVE data set, for the examples and tutorials.
#
# Source: STAVE object stave_data_2026.03.17.rds in the scrub repository
# (IDEELResearch/scrub, analysis/data-out/), pinned at commit 2d85d021
# (2026-03-17). scrub is MIT-licensed (copyright O. J. Watson); every survey
# keeps the reference of the study it came from.
#
# Counts come from STAVE's get_prevalence(), not from get_counts() on the
# exact variant string. get_prevalence() also counts surveys that recorded
# the locus in other variant strings (e.g. multi-locus haplotypes), and,
# for a mutant, the surveys where it was absent; get_counts() on "k13:622:I"
# returns only the 23 IGAD surveys where R622I was found, against 423 that
# genotyped the codon. Needs the STAVE package (github.com/mrc-ide/STAVE).
#
# Usage, from the package root: Rscript data-raw/igad_stave.R
#
# ------------------------------------------------------------------

suppressPackageStartupMessages(library(STAVE))

url <- "https://github.com/IDEELResearch/scrub/raw/2d85d021efe1c93dc55300610d944c9fbaf4af1a/analysis/data-out/stave_data_2026.03.17.rds"
tmp <- tempfile(fileext = ".rds")
download.file(url, tmp, mode = "wb", quiet = TRUE)
stopifnot(unname(tools::md5sum(tmp)) == "c1fd108b5665b88041217ced87f6a860")
stave <- suppressWarnings(readRDS(tmp))

igad <- c("Djibouti", "Eritrea", "Ethiopia", "Kenya", "Somalia", "South Sudan", "Sudan", "Uganda")

subset_variant <- function(variant) {
  p <- stave$get_prevalence(variant)
  p <- p[p$country_name %in% igad & p$denominator > 0, ]
  out <- data.frame(study_id = p$study_id, survey_id = p$survey_id, country = p$country_name,
                    site = p$site_name, longitude = p$longitude, latitude = p$latitude,
                    collection_start = p$collection_start, collection_end = p$collection_end,
                    collection_day = p$collection_day, positive = as.integer(p$numerator),
                    tested = as.integer(p$denominator), reference = p$reference,
                    reference_year = as.integer(p$reference_year), pmid = as.character(p$PMID),
                    stringsAsFactors = FALSE)
  chr <- vapply(out, is.character, logical(1))                 # two site names have a non-breaking space
  out[chr] <- lapply(out[chr], function(v) {                   # byte-wise, so it works in any locale
    v <- gsub("\xc2\xa0", " ", v, fixed = TRUE, useBytes = TRUE)
    Encoding(v) <- "UTF-8"
    v
  })
  out <- out[order(out$collection_day, out$country, out$site), ]
  rownames(out) <- NULL
  stopifnot(!anyNA(out[c("longitude", "latitude", "collection_day", "positive", "tested")]),
            all(out$positive <= out$tested))
  out
}

igad_crt76t <- subset_variant("crt:76:T")
igad_k13_622i <- subset_variant("k13:622:I")
message(sprintf("igad_crt76t: %d surveys, %d positive of %d tested, %s to %s",
                nrow(igad_crt76t), sum(igad_crt76t$positive), sum(igad_crt76t$tested),
                min(igad_crt76t$collection_day), max(igad_crt76t$collection_day)))
message(sprintf("igad_k13_622i: %d surveys (%d with none positive), %d positive of %d tested",
                nrow(igad_k13_622i), sum(igad_k13_622i$positive == 0), sum(igad_k13_622i$positive),
                sum(igad_k13_622i$tested)))

usethis::use_data(igad_crt76t, igad_k13_622i, overwrite = TRUE, compress = "xz")
