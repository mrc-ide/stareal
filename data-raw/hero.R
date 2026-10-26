# hero.R
#
# Makes man/figures/hero-crt76t.png, the home page figure: predicted crt
# K76T prevalence in the eight IGAD countries in July of four years, from
# the fit of tutorial 1 (the default settings, seed 1). From the package
# root:
#   Rscript data-raw/hero.R
#
# ------------------------------------------------------------------

source("data-raw/setup_crt76t.R")
library(ggplot2)
fit <- fit_areal(d, control = list(verbose = FALSE))
stopifnot(diagnose(fit)$pass)
pr <- predict(fit, verbose = FALSE)
years <- c("1995", "2005", "2015", "2024")
m <- do.call(rbind, lapply(years, function(y) {
  x <- sf::st_as_sf(pr, time = paste0(y, "-07"))
  x$year <- y
  x
}))
countries <- sf::st_as_sf(geodata::gadm(iso3, level = 0, path = gadm_path(), version = "4.1", resolution = 2))
countries <- suppressWarnings(sf::st_difference(sf::st_make_valid(sf::st_transform(countries, crs)), water))
g <- ggplot(m) +
  geom_sf(aes(fill = mean), colour = NA) +
  geom_sf(data = sf::st_sf(water), fill = "#bcd6ea", colour = NA) +
  geom_sf(data = countries, fill = NA, colour = "white", linewidth = 0.25) +
  facet_wrap(~year, nrow = 2) +
  scale_fill_viridis_c("K76T\nprevalence", limits = c(0, 1), option = "magma") +
  theme_void(base_size = 13) +
  theme(strip.text = element_text(size = 14, margin = margin(b = 3)), legend.title = element_text(size = 11),
        plot.margin = margin(4, 6, 4, 4))
dir.create("man/figures", showWarnings = FALSE, recursive = TRUE)
ggsave("man/figures/hero-crt76t.png", g, width = 7.5, height = 6.6, dpi = 150, bg = "white")
message("saved man/figures/hero-crt76t.png")
