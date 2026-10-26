# *stareal*: spatial-temporal areal estimation from survey data

The model: binomial counts on area-time cells, logit prevalence = fixed
effects + a separable Gaussian field (Leroux in space, AR1 in time) + a
cell-level nugget. The empty cells are integrated out exactly, so the
MCMC runs on the observed cells only; the full table is reconstructed
afterwards by conditional simulation (Matheron's rule).

## See also

Useful links:

- <https://github.com/mrc-ide/stareal>

- <https://mrc-ide.github.io/stareal/>

- Report bugs at <https://github.com/mrc-ide/stareal/issues>

## Author

**Maintainer**: Robert Verity <r.verity@imperial.ac.uk>
