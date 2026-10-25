# datasets.R
#
# Documentation of the data shipped with the package. The scripts that make
# them are in data-raw/ (not shipped).
#
# ------------------------------------------------------------------

#' Administrative units of the eight IGAD countries
#'
#' First-level administrative units of the member states of the
#' Intergovernmental Authority on Development (IGAD): Djibouti, Eritrea,
#' Ethiopia, Kenya, Somalia, South Sudan, Sudan and Uganda. Uganda is given
#' by district (ADM2, 151 units, 2006 boundaries), because its first level
#' in the source is four statistical regions, much coarser than the other
#' countries' first-level units.
#'
#' The national files come from different sources. Where they overlapped,
#' the overlap was given to one country, so every point lies in at most one
#' unit; where they leave gaps, neighbouring units across a border do not
#' quite touch. Build frames with `adjacency_snap` (e.g. 2000 m) so that
#' they count as neighbours. The boundaries are for illustration and imply
#' no position on disputed borders.
#'
#' @format An `sf` data frame (needs the `sf` package) with 268 rows and
#'   columns `id` (e.g. `"KEN_001"`), `name`, `iso3`, `country`, `level`
#'   (`"ADM1"`, or `"ADM2"` for Uganda), `licence`, `gb_id` (the
#'   geoBoundaries shape id) and `geometry` (MULTIPOLYGON, WGS84 longitude
#'   and latitude).
#' @source geoBoundaries (Runfola et al. 2020, PLoS ONE 15(4): e0231866),
#'   gbOpen release of 2023-12-12, simplified geometries. Licences differ by
#'   country (the `licence` column): Public Domain (Djibouti, from Natural
#'   Earth; Kenya, from RCMRD; Uganda); ODbL 1.0 (Eritrea and Somalia, from
#'   OpenStreetMap, (c) OpenStreetMap contributors); CC BY 4.0 (Ethiopia,
#'   from Open Africa / Code for Ethiopia); CC BY 3.0 IGO (Sudan and South
#'   Sudan, from OCHA and the Information Management Working Groups).
#' @seealso [igad_crt76t], [igad_k13_622i], and the tutorials on the
#'   package website.
#' @examples
#' if (requireNamespace("sf", quietly = TRUE)) {
#'   table(igad_adm1$country)
#' }
"igad_adm1"

#' Surveys of the chloroquine-resistance marker crt K76T in IGAD countries
#'
#' Survey counts of the *Plasmodium falciparum* crt K76T mutation (the
#' marker of chloroquine resistance) in the eight IGAD countries, 1989 to
#' 2024: 272 surveys from 90 studies.
#'
#' Counts are from STAVE's `get_prevalence("crt:76:T")`, which also counts
#' surveys that recorded codon 76 within multi-locus haplotypes.
#' `positive` counts infections carrying 76T, mixed infections included.
#'
#' @format A data frame with 272 rows and columns `study_id`, `survey_id`,
#'   `country`, `site`, `longitude`, `latitude` (WGS84), `collection_start`
#'   and `collection_end` (the collection window; `NA` for the 42 surveys
#'   that give only a date), `collection_day` (the date used for the survey:
#'   the middle of the collection window, or the date the study gave),
#'   `positive`, `tested`, `reference` (a link to the source study),
#'   `reference_year` and `pmid`.
#' @source The public STAVE data set (`stave_data_2026.03.17.rds`, scrub
#'   repository, IDEELResearch/scrub, commit 2d85d02; MIT licence), compiled
#'   from published studies and the WWARN molecular surveyors. Cite the
#'   original studies (`reference`) when using these data.
#' @seealso [igad_k13_622i], [igad_adm1]
#' @examples
#' head(igad_crt76t)
#' sum(igad_crt76t$positive) / sum(igad_crt76t$tested)
"igad_crt76t"

#' Surveys of the artemisinin partial-resistance marker k13 R622I in IGAD countries
#'
#' Survey counts of the *Plasmodium falciparum* kelch13 R622I mutation in
#' the eight IGAD countries, 1994 to 2024: 423 surveys that genotyped codon
#' 622, of which 400 found no R622I. The mutation has emerged in Eritrea and
#' Ethiopia.
#'
#' Counts are from STAVE's `get_prevalence("k13:622:I")`, which includes the
#' surveys where the mutation was absent (a count of 0), as a model of
#' prevalence needs.
#'
#' @format A data frame with 423 rows and the same columns as [igad_crt76t]
#'   (39 surveys give only a date, with no collection window).
#' @source As [igad_crt76t].
#' @seealso [igad_crt76t], [igad_adm1]
#' @examples
#' table(igad_k13_622i$country, igad_k13_622i$positive > 0)
"igad_k13_622i"
