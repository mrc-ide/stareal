# walkthrough.R - a typical run through stareal
#
# Run with demo("walkthrough", package = "stareal"), or source this file.
# Uses only data shipped with packages: a regular lattice, and the North
# Carolina counties that come with sf (plots of predictions also need
# ggplot2). It is extended as each part of the
# package is built; so far it covers frames, time frames, survey data,
# simulation, fitting, diagnostics, prediction and cross-validation.
# Covariates come next.

library(stareal)

# ---- 1. The spatial frame ---------------------------------------------------
# The areal units, their adjacency graph, and any groupings that predictions
# will be aggregated to. Three ways to build one:

# (a) a regular lattice, handy for experiments: 10 x 10 cells, grouped into
#     four 5 x 5 blocks
lat <- areal_lattice(10, 10, blocks = c(5, 5))
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

# ---- 2. The time frame ------------------------------------------------------
# A regular grid of months, quarters or years. A year (or month) given as
# the end runs to the end of that period.
tf <- time_frame(2005, 2024, "year")
tf
time_frame("2015-01", "2024-12", "month")
time_frame(T = 12)                                  # a plain index, 1..12
time_index(tf, as.Date(c("2004-12-31", "2005-06-15", "2024-12-31")))

# ---- 3. Simulated data with a known truth -----------------------------------
# One draw of the latent surface (intercept + trend + Leroux x AR1 field),
# then binomial surveys at randomly chosen county-years, with a nugget.
sim <- simulate_areal(fr, tf, sigma = 1, rho_s = 0.9, rho_t = 0.8, sigma_n = 0.3,
                      beta = c(-1, 1.5), n_cells = 300, k = rpois(300, 60) + 1, seed = 42)
sim
head(sim$data)
dim(sim$truth$prevalence)                           # units x years, the truth to recover

# ---- 4. Survey data into cells ----------------------------------------------
# (a) surveys that already carry a unit id and a date
d <- areal_data(sim$data, fr, tf, y = "positive", k = "tested", area = "id", date = "date")
head(as.data.frame(d))                              # the pooled cells the model will see
plot(d)                                             # where and when the data are

# (b) surveys with coordinates instead (longitude / latitude): each survey
#     gets a random point inside its county. Then, to show how problems are
#     handled, two points are moved offshore, east of the easternmost point
#     of the state: one ~3 km out (snapped back to the nearest county, being
#     within snap_max_m) and one ~90 km out (excluded), plus one data error
#     and one survey from before the window.
pts <- sf::st_sample(nc[match(sim$data$id, nc$FIPS), ], size = rep(1, nrow(sim$data)), by_polygon = TRUE)
xy <- sf::st_coordinates(pts)
surveys <- data.frame(lon = xy[, 1], lat = xy[, 2], date = sim$data$date,
                      positive = sim$data$positive, tested = sim$data$tested)
east <- sf::st_coordinates(nc)[which.max(sf::st_coordinates(nc)[, 1]), 1:2]
surveys[1:2, c("lon", "lat")] <- rbind(east + c(0.03, 0), east + c(1, 0))
surveys$positive[3] <- surveys$tested[3] + 1                # a data error: y > k
surveys$date[4] <- as.Date("1999-06-01")                    # before the window
d2 <- areal_data(surveys, fr, tf, y = "positive", k = "tested",
                 coords = c("lon", "lat"), date = "date", snap_max_m = 10000)
table(d2$surveys$.excluded, useNA = "ifany")        # every survey accounted for
d2$surveys[1:4, c("lon", "lat", ".area", ".snap_m", ".excluded")]

# ---- 5. Fit and diagnose -----------------------------------------------------
# The priors and the sampler settings each have a constructor; the defaults
# are those of the Africa analysis.
areal_priors()
ctl <- mcmc_control(n_iter = 1500, n_burn = 500, n_store = 200, n_chains = 2, ess_sweeps = 5)
ctl
# Chains run in parallel, one per worker process. Progress lines say how
# long is left. A fit this size takes well under a minute; real ones take
# hours, so give checkpoint_dir = "<a folder>" to make them resumable.
fit <- fit_areal(d, control = ctl)
fit
summary(fit)                # hyperparameters and fixed effects, with R-hat and ESS
sim$truth$theta             # ... against the truth
plot(fit)                   # traces, burn-in left of the dashed line

# diagnose() checks EVERY parameter, including each observed cell's latent
# value: hyperparameters can look converged while the cells mix poorly.
# This short demo run is not expected to pass the default thresholds.
dg <- diagnose(fit)
dg
head(dg$table[dg$table$type == "latent", ])

# the draws, for the user's own tools
draws <- posterior::as_draws_array(fit, latent = FALSE)
posterior::summarise_draws(draws)

# ---- 6. Predict ---------------------------------------------------------------
# The full county x year surface, one surface per posterior draw, observed
# cells or not; summarised per county-year, and aggregated draw by draw to
# the regions (so their uncertainty accounts for the correlations between
# counties). Weights could be population: here, the 1974 births of each
# county.
pr <- predict(fit, levels = "region", weights = setNames(nc$BIR74, nc$FIPS),
              exceedance = 0.5, keep_draws = "37119")
pr
head(pr$area)               # county-years: data (y, k) and the posterior of prevalence
head(pr$region)             # region-years
dim(pr$draws)               # all draws for one county, e.g. for its own summaries

# the truth is known here: how often do the 95% intervals cover it?
truth <- sim$truth$prevalence[cbind(match(pr$area$id, fr$ids), pr$area$t)]
mean(pr$area$q2.5 <= truth & truth <= pr$area$q97.5)

plot(pr, time = "2024")                       # a map of one year
plot(pr, time = "2024", stat = "p_gt_0.5")    # where prevalence is probably above 50%
plot(pr)                                      # time series of the best-surveyed counties
plot(pr, level = "region")                    # and of the regions

# ---- 7. Cross-validation ------------------------------------------------------
# Hold out each fold of cells in turn, refit, and score the held-out cells:
# the mean log score compares models fitted to the same cells (higher is
# better), and about 5% of observations should fall outside their 95%
# predictive interval. Each fold is a full refit; here only the first of
# five is run.
cv <- crossval(fit, folds = 5, fold = 1)
cv
head(cv$cells)

# ---- 8. Covariates ------------------------------------------------------------
# The mean can include covariates by area, by time, or by both, each in a
# table keyed by `id` and/or `t` (or `time`, the label). Here: the log of
# each county's 1974 births (by area), and a programme that starts in 2015
# (by time; a character column becomes a factor). Covariates must cover
# every area and time step, since predictions are made everywhere. Numeric
# covariates are standardised, so effects are per sd. Data are simulated
# with known effects, then fitted with the same mean.
covs <- list(data.frame(id = nc$FIPS, log_births = log(nc$BIR74)),
             data.frame(time = tf$labels, programme = ifelse(2005:2024 >= 2015, "on", "off")))
mean_cov <- ~ trend() + log_births + programme
sim2 <- simulate_areal(fr, tf, sigma = 0.8, rho_s = 0.9, rho_t = 0.8, sigma_n = 0.3,
                       beta = c(-1, 1, 0.5, -0.7), n_cells = 300, k = 60, seed = 7,
                       mean = mean_cov, covariates = covs)
sim2$truth$beta
d3 <- areal_data(sim2$data, fr, tf, y = "positive", k = "tested", area = "id", t = "t",
                 verbose = FALSE)
fit2 <- fit_areal(d3, mean = mean_cov, covariates = covs, control = ctl)
summary(fit2)               # the covariate effects, against sim2$truth$beta
pr2 <- predict(fit2, levels = "region")   # covariates carry through to every cell
head(pr2$region)
