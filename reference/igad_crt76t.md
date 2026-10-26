# Surveys of the chloroquine-resistance marker crt K76T in IGAD countries

Survey counts of the *Plasmodium falciparum* crt K76T mutation (the
marker of chloroquine resistance) in the eight member states of the
Intergovernmental Authority on Development (IGAD): Djibouti, Eritrea,
Ethiopia, Kenya, Somalia, South Sudan, Sudan and Uganda, 1989 to 2024.
328 surveys from 106 studies.

## Usage

``` r
igad_crt76t
```

## Format

A data frame with 328 rows and columns `study_id`, `survey_id`,
`country` (from the location), `site`, `longitude`, `latitude` (WGS84),
`collection_start` and `collection_end` (the collection window; `NA` for
surveys that give only a date), `collection_day` (the date used for the
survey: the middle of the collection window, or the date the study
gave), `positive`, `tested`, `reference` (a link to the source study),
`reference_year` and `pmid`.

## Source

The public STAVE data set (`stave_data_2026.03.17.rds`, scrub
repository, IDEELResearch/scrub, commit 2d85d02; MIT licence), compiled
from published studies and the WWARN molecular surveyors. Cite the
original studies (`reference`) when using these data.

## Details

Counts are from STAVE's `get_prevalence("crt:76:T")`, which counts every
survey that genotyped codon 76, including those where 76T was absent and
those that recorded the codon within multi-locus haplotypes. `positive`
counts infections carrying 76T, mixed infections included.

Surveys were selected by location, not by the country names recorded in
STAVE, which are not reliable: a survey is included if it lies within 10
km of the GADM 4.1 outline of an IGAD country, and `country` is that
country. Boundaries are not shipped with the package; the tutorials
download GADM with
[`geodata::gadm()`](https://rdrr.io/pkg/geodata/man/gadm.html).

## Examples

``` r
head(igad_crt76t)
#>                                 study_id
#> 1                         WWARN_11319692
#> 2                         WWARN_15031793
#> 3                         WWARN_19450282
#> 4                         WWARN_15031793
#> 5 MalariaGEN_1132_PF_K1000G_DBS_KE_BEJON
#> 6                         WWARN_19450282
#>                                            survey_id country         site
#> 1                   WWARN_11319692_asar_village_1989   Sudan Asar village
#> 2                           WWARN_15031793_asar_1990   Sudan         Asar
#> 3                         WWARN_19450282_kilifi_1993   Kenya       Kilifi
#> 4                           WWARN_15031793_asar_1993   Sudan         Asar
#> 5 MalariaGEN_1132_PF_K1000G_DBS_KE_BEJON_kilifi_1994   Kenya       Kilifi
#> 6                         WWARN_19450282_kilifi_1994   Kenya       Kilifi
#>   longitude  latitude collection_start collection_end collection_day positive
#> 1  25.34250 13.720400       1989-01-01     1990-12-31     1989-12-31       26
#> 2  30.08740 16.085800       1990-01-01     1990-12-31     1990-07-02       23
#> 3  39.85000 -3.633330       1993-01-01     1993-12-31     1993-07-02       17
#> 4  30.08740 16.085800       1993-01-01     1993-12-31     1993-07-02       36
#> 5  39.68689 -3.174404             <NA>           <NA>     1994-07-01       42
#> 6  39.85000 -3.633330       1994-01-01     1994-12-31     1994-07-02       35
#>   tested                                        reference reference_year
#> 1     39        https://pubmed.ncbi.nlm.nih.gov/11319692/           2001
#> 2     44        https://pubmed.ncbi.nlm.nih.gov/15031793/           2004
#> 3     18        https://pubmed.ncbi.nlm.nih.gov/19450282/           2009
#> 4     50        https://pubmed.ncbi.nlm.nih.gov/15031793/           2004
#> 5     53 https://doi.org/10.12688/wellcomeopenres.24031.1           2025
#> 6     38        https://pubmed.ncbi.nlm.nih.gov/19450282/           2009
#>       pmid
#> 1 11319692
#> 2 15031793
#> 3 19450282
#> 4 15031793
#> 5 41267740
#> 6 19450282
table(igad_crt76t$country)
#> 
#>    Djibouti     Eritrea    Ethiopia       Kenya     Somalia South Sudan 
#>           1           1          46         143           2           6 
#>       Sudan      Uganda 
#>          42          87 
sum(igad_crt76t$positive) / sum(igad_crt76t$tested)
#> [1] 0.712045
```
