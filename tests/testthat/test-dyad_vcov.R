library(dyadfast)

# Load bundled sim data (lm on dY ~ dX, dyad1/dyad2 as ego/alter)
load(test_path("sim-dat.RData"))   # -> dyad.sim

fit <- lm(dY ~ dX, data = dyad.sim)
res <- dyad_vcov(fit, ego = dyad.sim$dyad1, alter = dyad.sim$dyad2)

# ── Regression test against pre-computed reference values ─────────────────────
# Reference SEs computed once with dyadRobust and locked in here.
# Max |sehat_dyadfast - sehat_dyadRobust| = 4.4e-16 on this dataset.
ref_sehat <- c("(Intercept)" = 0.21039921, dX = 0.08144004)

test_that("sehat matches dyadRobust reference to 1e-6", {
  expect_equal(res$sehat, ref_sehat, tolerance = 1e-6)
})

test_that("bhat matches lm() coefficients", {
  expect_equal(unname(res$bhat), unname(coef(fit)), tolerance = 1e-15)
})

test_that("Vhat is symmetric", {
  expect_equal(res$Vhat, t(res$Vhat))
})

test_that("Vhat is positive semi-definite", {
  evs <- eigen(res$Vhat, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(evs >= -1e-14))
})

test_that("N and n are reported correctly", {
  expect_equal(res$N, length(unique(c(dyad.sim$dyad1, dyad.sim$dyad2))))
  expect_equal(res$n, nrow(dyad.sim))
})
