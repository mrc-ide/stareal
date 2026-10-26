# Method

This page summarises the method. The derivation and its checks are in
the report of the project the package grew out of (to be published).

## Why this method

Areal models in space and time attach a latent value to every area-time
cell, but surveys usually reach only a small fraction of the cells. The
standard tools treat the empty cells in one of two ways. Each makes one
of two operations on a Gaussian vector cheap and the other expensive:

- **Marginalising**, integrating some cells out, is free when the model
  is specified by its covariance: the observed cells’ covariance is just
  their rows and columns. This is how Gaussian-process (geostatistical)
  models work. They carry only the observed locations, and they predict
  the rest afterwards by kriging. **Conditioning**, the rest given the
  observed, then needs dense matrix solves.
- Areal priors (CAR, Leroux, BYM, AR1 and random walks in time) are
  usually specified by their **sparse precision** matrix, a Gaussian
  Markov random field (GMRF). Here conditioning is cheap, needing only
  sparse blocks of the precision. Marginalising destroys the sparsity:
  the precision of the observed cells alone is a dense Schur complement.
  GMRF-based software (INLA, Rue et al. 2009; CARBayes, Lee 2013; ICAR
  models in Stan, Morris et al. 2019) therefore keeps a latent value for
  every cell, empty or not, because that is what keeps the algebra
  sparse.

With hundreds of areas and a few decades of months, the full table holds
tens of thousands of cells. That is a heavy load for MCMC, and coarser
time steps are often chosen to keep it down. *stareal* takes the
Gaussian-process route for areal models. For a *separable* space-time
prior, the covariance between any two cells is the product of a small
spatial and a small temporal matrix. The observed cells’ covariance is
therefore cheap, as for a Gaussian process, while the areal prior is
kept. The model is fitted on the observed cells alone, and the full
surface is rebuilt afterwards, one posterior draw at a time.

This pays off when the data are sparse relative to the space-time grid:
fine time steps, many areas, few surveys, as is typical of molecular
surveillance. It does not when most cells are observed, when the prior
is not separable (for example, temporal correlation that varies by
area), or when there are tens of thousands of areas.

## The model

Space is a set of m areal units (districts, counties, grid cells) with
an adjacency graph; time is a regular grid of T steps (months, quarters,
years). An *area-time cell* is one area in one time step. Surveys report
y positives out of k tested. For each observed cell,

y\_{at} \sim \text{Binomial}\big(k\_{at},\\
\text{logit}^{-1}(x\_{at})\big), \qquad x\_{at} = z\_{at}^\top \beta +
f\_{at} + \epsilon\_{at}.

- z\_{at}^\top \beta are the **fixed effects**: an intercept, by default
  a linear trend in time
  ([`trend()`](https://mrc-ide.github.io/stareal/reference/trend.md),
  from -1 at the first step to 1 at the last), and any covariates, which
  may vary by area, by time or by both. Their prior is Gaussian, \beta
  \sim N(b_0, B).
- f is a **Gaussian field** with separable covariance \Sigma = \Sigma_t
  \otimes \Sigma_s:
  - in space, a Leroux conditional autoregressive model on the graph,
    with precision Q_s = \tau\\\rho_s R + (1 - \rho_s) I\\, where R =
    D - W is the graph Laplacian (W the adjacency matrix, D the numbers
    of neighbours), \sigma^2 = 1/\tau sets the scale and \rho_s \in
    (0, 1) the strength of spatial dependence;
  - in time, a first-order autoregression with unit variance and lag-one
    correlation \rho_t.
- \epsilon\_{at} \sim N(0, \sigma_n^2) is a **cell-level nugget**:
  variation between surveys that the smooth field does not explain. It
  applies to the observed cells; the predicted surface is the field
  without it.

**Surveys in the same cell are pooled** into one count, y\_{at} of
k\_{at}. This rests on an assumption: given the cell’s latent value
x\_{at} (nugget included), the surveys in it are independent binomial
samples of one prevalence. A sum of independent binomial counts with a
common probability is binomial, so under that assumption pooling loses
nothing. The model has no variation between surveys within a cell. Two
surveys of different sites in the same district and month are treated as
one larger survey. If they differ by more than binomial sampling,
through different sites, populations or methods, pooling overstates how
much the cell’s data say. Longer time steps and larger areas put more
surveys together, and make this more likely.

The hyperparameters \theta = (\sigma, \rho_s, \rho_t, \sigma_n) have
half-normal and Beta priors
([`areal_priors()`](https://mrc-ide.github.io/stareal/reference/areal_priors.md)).

## Fitting on the observed cells only

Most area-time cells have no data: in the [crt K76T
tutorial](https://mrc-ide.github.io/stareal/articles/crt76t.md), 264 of
74,304 unit-months. Because everything except the likelihood is
Gaussian, the unobserved cells can be integrated out *exactly*: the n
observed cells’ latent values x_o have the Gaussian prior

x_o \sim N\big(Z_o b_0,\\ \Sigma\_{oo}\big), \qquad \Sigma\_{oo} =
\Sigma_t\[t_o, t_o\] \circ \Sigma_s\[a_o, a_o\] + Z_o B Z_o^\top +
\sigma_n^2 I,

where a_o and t_o are the observed cells’ areas and time steps and \circ
is the elementwise product. The fixed effects are integrated out too
(the Z_o B Z_o^\top term) and drawn afterwards from their exact
conditional distribution. The spatial block \Sigma_s\[a_o, a_o\] comes
from a single eigendecomposition of R, done once when the frame is
built, so a change of \tau or \rho_s only rescales eigenvalues.

The MCMC therefore works on n latent values, not m \times T, and its
cost per iteration grows as n^3, independently of the size of the full
table. This is what makes continental maps over decades, at monthly
resolution, feasible.

## The sampler

Each iteration updates the latent values, then the hyperparameters.

- **Latent values** by elliptical slice sampling (Murray, Adams &
  MacKay, 2010) against a *posterior-scaled* Gaussian (Nishihara, Murray
  & Adams, 2014): the Gaussian prior combined with a Gaussian
  approximation of each cell’s likelihood (pseudo-data), with the
  approximation divided back out, so the update is exact. Proposing on
  this ellipse, rather than the prior’s, is what lets heavily surveyed
  cells mix.
- **Hyperparameters** by two adaptive Metropolis moves: a centred move
  (\theta given x_o) and the surrogate-data move of Murray and Adams
  (2010), in which x_o moves with \theta. The centred move alone gets
  stuck, because n latent values pin \theta.
- **During burn-in** the proposals adapt, and the pseudo-data are
  refitted twice to a Gaussian approximation of the likelihood at the
  chain’s mean latent values. Both are then fixed, so the draws after
  burn-in come from an exact MCMC. The refit matters most for rare
  outcomes, where many cells have no positives.
- **Chains** start from dispersed points and run in parallel;
  [`diagnose()`](https://mrc-ide.github.io/stareal/reference/diagnose.md)
  checks R-hat and effective sample sizes of every parameter, the latent
  value of every observed cell included.

## Prediction and aggregation

For each posterior draw of (\theta, \beta, x_o), the full m \times T
surface is drawn by conditional simulation (Matheron’s rule): an
unconditional draw of the field, corrected by kriging so that it agrees
with x_o. This is exact, and needs no matrix larger than n \times n
beyond the stored eigenvectors; the work is split into blocks of areas
so that memory stays bounded.

Each draw of the surface is then aggregated to any grouping of the areas
(regions, countries) as a weighted mean, and only then summarised. Group
estimates therefore carry the correlation between their areas, which
averaging area-level summaries would lose.

## Cross-validation

[`crossval()`](https://mrc-ide.github.io/stareal/reference/crossval.md)
holds out folds of observed cells, refits, and predicts each held-out
cell from its exact conditional distribution given the fitted cells,
nugget included. The log predictive density integrates that Gaussian
against the binomial likelihood by adaptive Gauss-Hermite quadrature,
and calibration is checked with a randomised probability integral
transform.

## Scale

Two things set the cost of a fit:

- **The eigendecomposition of the graph.**
  [`areal_frame()`](https://mrc-ide.github.io/stareal/reference/areal_frame.md)
  does it once. It is a dense m \times m problem, costing O(m^3) in the
  number of areas m: quick for hundreds of areas, noticeably slower for
  thousands. Build the frame once and save it with
  [`saveRDS()`](https://rdrr.io/r/base/readRDS.html).
- **Each MCMC iteration.** It costs O(n^3) in the number of observed
  cells n, because of a dense Cholesky factorisation. It does not depend
  on the number of areas or time steps.

The sampler is almost entirely dense linear algebra, so an optimised
BLAS (Intel MKL, OpenBLAS or Apple Accelerate) makes a large difference
over R’s reference BLAS. The method was developed on about 1,300 and
1,800 observed cells. Prediction holds one block of areas at a time, so
memory stays bounded however many draws are used.

## References

- Leroux, B. G., Lei, X. & Breslow, N. (2000). Estimation of disease
  rates in small areas: a new mixed model for spatial dependence. In
  *Statistical Models in Epidemiology, the Environment, and Clinical
  Trials*, 179–191.
- Murray, I., Adams, R. P. & MacKay, D. J. C. (2010). Elliptical slice
  sampling. *AISTATS*.
- Murray, I. & Adams, R. P. (2010). Slice sampling covariance
  hyperparameters of latent Gaussian models. *NeurIPS*.
- Nishihara, R., Murray, I. & Adams, R. P. (2014). Parallel MCMC with
  generalized elliptical slice sampling. *JMLR* 15, 2087–2112.
- Lee, D. (2013). CARBayes: an R package for Bayesian spatial modeling
  with conditional autoregressive priors. *Journal of Statistical
  Software* 55(13).
- Morris, M., Wheeler-Martin, K., Simpson, D., Mooney, S. J., Gelman, A.
  & DiMaggio, C. (2019). Bayesian hierarchical spatial models:
  implementing the Besag York Mollié model in Stan. *Spatial and
  Spatio-temporal Epidemiology* 31, 100301.
- Rue, H., Martino, S. & Chopin, N. (2009). Approximate Bayesian
  inference for latent Gaussian models by using integrated nested
  Laplace approximations. *Journal of the Royal Statistical Society B*
  71, 319–392.
- Andrieu, C. & Thoms, J. (2008). A tutorial on adaptive MCMC.
  *Statistics and Computing* 18, 343–373.
