# dyadfast

Fast dyadic-robust variance estimation for R. Implements the
[Aronow, Samii & Assenova (2015)](https://doi.org/10.1093/pan/mpv018)
cluster-robust variance estimator for dyadic data via a single O(nK)
scatter-add pass, replacing the O(m·n·K) agent loop in
[dyadRobust](https://github.com/wpmarble/dyadRobust).

## Installation

```r
# install.packages("devtools")
devtools::install_github("xiangao/dyadfast")
```

## Usage

```r
library(dyadfast)

res <- dyad_vcov(fit, ego = df$agent_i, alter = df$agent_j)
res$bhat    # coefficients
res$sehat   # dyadic-robust standard errors
res$Vhat    # full covariance matrix
```

Works with any model class compatible with `sandwich::estfun()` and
`sandwich::bread()`: `lm`, `feols` (fixest), `felm` (lfe), and more.

### With fixest

```r
library(fixest)
library(dyadfast)

m <- feols(y ~ x1 + x2 | fe_var,
           data = df, lean = FALSE, demeaned = TRUE)

# Align ego/alter with observations used by feols
obs_rem <- m$obs_selection$obsRemoved
if (is.null(obs_rem) || length(obs_rem) == 0L) {
  ego_used   <- df$agent_i
  alter_used <- df$agent_j
} else {
  ego_used   <- df$agent_i[obs_rem]
  alter_used <- df$agent_j[obs_rem]
}

res <- dyad_vcov(m, ego = ego_used, alter = alter_used)
```

## Formula

```
Vhat = B (H'H - S'S) / n² B
```

where:

- **S** (n × K) — score matrix from `sandwich::estfun(fit)`
- **H** (m × K) — agent-level score sums via `rowsum()` scatter-add:
  `H[i,] = Σ_{j: (i,j) ∈ data} S[ij,]`
- **B** (K × K) — bread `n·(X'X)⁻¹` from `sandwich::bread(fit)`

This is algebraically equivalent to the dyadRobust loop formula but
computed in a single pass over the data.

## Performance

Verified on Fisman et al. (2006) speed-dating data
(`feols(dec ~ amb + attr + intel | iid, weights = wts)`,
n = 3,454 dyads, N = 542 agents):

| | dyadRobust | dyadfast |
|---|---:|---:|
| Elapsed | 1.13 s | 0.012 s |
| Speedup | — | **94×** |
| Max \|SE diff\| | — | 1.6e-15 |

Speedup grows with N (number of agents): dyadRobust loops over all N
agents; dyadfast does not.

## References

Aronow, P.M., Samii, C. & Assenova, V.A. (2015). Cluster-robust
variance estimation for dyadic data. *Political Analysis*, 23(4),
564–577. https://doi.org/10.1093/pan/mpv018
