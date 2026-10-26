# walkthrough.R - a typical run through stareal
#
# Run with demo("walkthrough", package = "stareal"), or source this file.
# Uses only data shipped with packages: a regular lattice, and the North
# Carolina counties that come with sf (plots of predictions also need
# ggplot2). Sections:
#   1. the spatial frame        5. fitting and checking a fit
#   2. the time frame           6. prediction: areas, groups, maps
#   3. simulated data           7. covariates
#   4. real-style survey data   8. cross-validation
#                               9. long fits: checkpoints and resuming
# Everything runs in a couple of minutes.

library(stareal)
set.seed(42)                     # the demo's own random choices (survey points)

# ---- 1. The spatial frame ---------------------------------------------------
# The areal units, their adjacency graph, and any groupings that predictions
# will be aggregated to. Three ways to build one:

# (a) a regular lattice, handy for experiments: 10 x 10 cells, grouped into
#     four 5 x 5 blocks
lat <- areal_lattice(10, 10, blocks = 5)
lat

# (b) a ready-made adjacency matrix, with ids and groupings in row order
W <- matrix(0, 4, 4)
W[cbind(1:3, 2:4)] <- 1
W <- W + t(W)
fr_adj <- areal_frame(adjacency = W, ids = c("a", "b", "c", "d"),
                      groups = data.frame(region = c("north", "north", "south", "south")))
fr_adj

# (c) polygons (an sf object): adjacency is computed (queen by default),
#     islands are linked to the mainland, and survey coordinates can later be
#     assigned to units. Here: the 100 counties of North Carolina, with an
#     invented grouping into three regions by longitude.
nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
lon <- sf::st_coordinates(suppressWarnings(sf::st_centroid(sf::st_geometry(nc))))[, 1]
nc$region <- cut(lon, 3, labels = c("west", "central", "east"))
fr <- areal_frame(nc, id = "FIPS", groups = c(county = "NAME", region = "region"))
fr
fr$links            # island edges added, if any (NULL: the graph was connected)
plot(fr)            # polygons, adjacency edges, island links in red
# A real frame is built once and saved: saveRDS(fr, "frame.rds")

# ---- 2. The time frame ------------------------------------------------------
# A regular grid of months, quarters or years. A year (or month) given as
# the end runs to the end of that period.
tf <- time_frame(2005, 2024, "year")
tf
time_frame("2015-01", "2024-12", "month")
time_frame(n_steps = 12)                            # a plain index, 1..12
time_index(tf, as.Date(c("2004-12-31", "2005-06-15", "2024-12-31")))

# ---- 3. Simulated data ------------------------------------------------------
# One draw of the latent surface (intercept + trend + Leroux x AR1 field),
# then binomial surveys at randomly chosen county-years, with a nugget. The
# result is ordinary survey data, ready to fit, with the truth attached.
sim <- simulate_areal(fr, tf, sigma = 1, rho_s = 0.9, rho_t = 0.8, sigma_n = 0.3,
                      beta = c(-1, 1.5), n_cells = 300, k = 60, seed = 42)
sim                                                 # cells, and the truth
head(sim$surveys)                                   # the simulated surveys
head(sim$cells)                                     # the pooled cells the model will see
dim(sim$truth$prevalence)                           # counties x years, the truth to recover
plot(sim)                                           # where and when the data are

# ---- 4. Real-style survey data ------------------------------------------------
# Real surveys usually come with coordinates and dates. Here each simulated
# survey gets a random point inside its county; then, to show how problems
# are handled, two points are moved offshore east of the state (one ~3 km out,
# snapped back to the nearest county; one ~90 km out, excluded), one count is
# made impossible and one survey is dated before the window.
pts <- sf::st_sample(nc[match(sim$surveys$id, nc$FIPS), ], size = rep(1, nrow(sim$surveys)), by_polygon = TRUE)
xy <- sf::st_coordinates(pts)
surveys <- data.frame(lon = xy[, 1], lat = xy[, 2], date = sim$surveys$date,
                      positive = sim$surveys$y, tested = sim$surveys$k)
east <- sf::st_coordinates(nc)[which.max(sf::st_coordinates(nc)[, 1]), 1:2]
surveys[1:2, c("lon", "lat")] <- rbind(east + c(0.03, 0), east + c(1, 0))
surveys$positive[3] <- surveys$tested[3] + 1                # a data error: y > k
surveys$date[4] <- as.Date("1999-06-01")                    # before the window
d <- areal_data(surveys, fr, tf, y = "positive", k = "tested",
                coords = c("lon", "lat"), date = "date", snap_max_m = 10000)
d                                                   # every exclusion, with its reason
d$surveys[1:4, c("lon", "lat", ".id", ".snap_m", ".excluded")]

# ---- 5. Fitting and checking a fit --------------------------------------------
# The defaults are a pilot run (4 chains x 3,000 iterations); a short run
# keeps the demo quick. Settings can be given as a plain list: anything
# left out takes its default (see ?mcmc_control and ?areal_priors).
# Chains run in parallel, one per worker process, with progress lines.
fit <- fit_areal(sim, control = list(n_iter = 1500, n_burn = 500, n_store = 200, n_chains = 2,
                                     ess_sweeps = 5))
fit
summary(fit)                # hyperparameters and fixed effects, with R-hat and ESS
sim$truth$theta             # ... against the truth
sim$truth$beta
plot(fit)                   # traces, burn-in left of the dashed line

# diagnose() checks EVERY parameter, including each observed cell's latent
# value: hyperparameters can look converged while the cells mix poorly.
# This short demo run is not expected to pass the default thresholds.
dg <- diagnose(fit)
dg
head(dg$table[dg$table$type == "latent", ])

# the draws, for your own tools (the posterior package)
draws <- posterior::as_draws_df(fit, latent = FALSE)
posterior::summarise_draws(draws)

# ---- 6. Prediction: areas, groups, maps ----------------------------------------
# The full county x year surface, one surface per posterior draw, observed
# cells or not; summarised per county-year, and aggregated draw by draw to
# the regions (so their uncertainty accounts for the correlations between
# counties). Weights could be population: here, the 1974 births of each
# county.
pr <- predict(fit, groups = "region", weights = setNames(nc$BIR74, nc$FIPS),
              exceedance = 50, keep_draws = "37119")         # prevalence in percent
pr
head(pr$area)               # county-years: data (y, k) and the posterior of prevalence (%)
head(pr$groups$region)      # region-years
dim(pr$draws)               # all draws for one county, e.g. for its own summaries

# the truth is known here: how often do the 95% intervals cover it?
truth <- sim$truth$prevalence[cbind(match(pr$area$id, fr$ids), pr$area$t)]
mean(pr$area$q2.5 <= truth & truth <= pr$area$q97.5)

plot(pr, time = "2024")                       # a map of one year
plot(pr, time = "2024", stat = "p_gt_50")     # where prevalence is probably above 50%
plot(pr)                                      # time series of the best-surveyed counties
plot(pr, group = "region")                    # and of the regions

# for maps of your own: the table for one year, as polygons
nc_2024 <- sf::st_as_sf(pr, time = "2024")
head(nc_2024[, c("id", "county", "mean", "q2.5", "q97.5")])

# ---- 7. Covariates ------------------------------------------------------------
# The formula can include covariates by area, by time, or by both, each in a
# table keyed by `id` and/or `t` (step number) or `time` (label).
# Covariates must cover every area and time step, since predictions are
# made everywhere. Numeric covariates are standardised, so effects are per
# sd; a character column becomes a factor. Here: the log of each county's
# 1974 births (by area), and a programme that starts in 2015 (by time).
# Data are simulated with known effects, then fitted with the same formula.
covs <- list(data.frame(id = nc$FIPS, log_births = log(nc$BIR74)),
             data.frame(time = tf$labels, programme = ifelse(2005:2024 >= 2015, "on", "off")))
f_cov <- ~ trend() + log_births + programme
sim2 <- simulate_areal(fr, tf, formula = f_cov, covariates = covs,
                       beta = c(intercept = -1, trend = 1, log_births = 0.5, programmeon = -0.7),
                       sigma = 0.8, rho_s = 0.9, rho_t = 0.8, sigma_n = 0.3, n_cells = 300, k = 60, seed = 7)
fit2 <- fit_areal(sim2, formula = f_cov, covariates = covs,
                  priors = list(beta_sd = c(programmeon = 1)),        # priors can be named by design column
                  control = list(n_iter = 1500, n_burn = 500, n_store = 200, n_chains = 2, ess_sweeps = 5))
summary(fit2)               # the covariate effects ...
sim2$truth$beta             # ... against the truth
pr2 <- predict(fit2, groups = "region")       # covariates carry through to every cell
head(pr2$groups$region)

# ---- 8. Cross-validation ------------------------------------------------------
# Hold out each fold of cells in turn, refit, and score the held-out cells:
# the mean log score compares models fitted to the same cells (higher is
# better), and about 5% of observations should fall outside their 95%
# predictive interval. Each fold is a full refit; here only the first of
# five is run.
cv <- crossval(fit, folds = 5, fold = 1)
cv
head(cv$cells)

# ---- 9. Long fits: checkpoints and resuming -------------------------------------
# A long fit can save each chain's state as it goes. If it is interrupted
# (or R is closed), running the same call again resumes every chain where it
# stopped, and gives exactly the draws of an uninterrupted run; once
# finished, running it again just reads the result back. Use a folder
# outside Dropbox/OneDrive. (A temporary folder here.)
ck <- file.path(tempdir(), "nc_chains")
fit_ck <- fit_areal(sim, control = list(n_iter = 600, n_burn = 200, n_store = 100, n_chains = 2,
                                        ess_sweeps = 5, checkpoint_every = 100),
                    checkpoint_dir = ck)
list.files(ck)
fit_again <- fit_areal(sim, control = list(n_iter = 600, n_burn = 200, n_store = 100, n_chains = 2,
                                           ess_sweeps = 5, checkpoint_every = 100),
                       checkpoint_dir = ck)                           # read back, not rerun
identical(fit_ck$draws, fit_again$draws)
