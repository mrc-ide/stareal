# datasets.R
#
# Documentation of the data shipped with the package. The script that makes
# them is in data-raw/ (not shipped).
#
# ------------------------------------------------------------------

#' Surveys of the chloroquine-resistance marker crt K76T in IGAD countries
#'
#' Survey counts of the *Plasmodium falciparum* crt K76T mutation (the
#' marker of chloroquine resistance) in the eight member states of the
#' Intergovernmental Authority on Development (IGAD): Djibouti, Eritrea,
#' Ethiopia, Kenya, Somalia, South Sudan, Sudan and Uganda, 1989 to 2024.
#' 328 surveys from 106 studies.
#'
#' Counts are from STAVE's `get_prevalence("crt:76:T")`, which counts every
#' survey that genotyped codon 76, including those where 76T was absent and
#' those that recorded the codon within multi-locus haplotypes. `positive`
#' counts infections carrying 76T, mixed infections included.
#'
#' Surveys were selected by location, not by the country names recorded in
#' STAVE, which are not reliable: a survey is included if it lies within
#' 10 km of the GADM 4.1 outline of an IGAD country, and `country` is that
#' country. Boundaries are not shipped with the package; the tutorials
#' download GADM with `geodata::gadm()`.
#'
#' @format A data frame with 328 rows and columns `study_id`, `survey_id`,
#'   `country` (from the location), `site`, `longitude`, `latitude` (WGS84),
#'   `collection_start` and `collection_end` (the collection window; `NA`
#'   for surveys that give only a date), `collection_day` (the date used for
#'   the survey: the middle of the collection window, or the date the study
#'   gave), `positive`, `tested`, `reference` (a link to the source study),
#'   `reference_year` and `pmid`.
#' @source The public STAVE data set (`stave_data_2026.03.17.rds`, scrub
#'   repository, IDEELResearch/scrub, commit 2d85d02; MIT licence), compiled
#'   from published studies and the WWARN molecular surveyors. Cite the
#'   original studies (`reference`) when using these data.
#' @examples
#' head(igad_crt76t)
#' table(igad_crt76t$country)
#' sum(igad_crt76t$positive) / sum(igad_crt76t$tested)
"igad_crt76t"
