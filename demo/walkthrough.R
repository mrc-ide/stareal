# walkthrough.R - a typical run through stareal
#
# Run with demo("walkthrough", package = "stareal"), or source this file.
# Uses only data shipped with packages: a regular lattice, and the North
# Carolina counties that come with sf. It is extended as each part of the
# package is built; so far it covers frames, time frames, survey data,
# simulation, fitting and diagnostics. Prediction comes next.

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
# Not yet built. Coming in the next steps:
#   pr <- predict(fit, levels = c("region"))
#   cv <- crossval(fit, folds = 10)
