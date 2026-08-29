# =============================================================================
# test_mspaf_identity.R
# -----------------------------------------------------------------------------
# Standalone verification of the msPAF_EC50 exact-identity checks described in
# Docs/TaxoTox_Technical_Methods.md Section 5.5 ("Consistency check"):
#
#   msPAF_EC50 = Phi(log10(S) / sigma_mix)
#
# should equal exactly 0.05 whenever S (the summed EC50-based toxic units)
# lands exactly on HC5-PTI = 1, in either of the two cases where the identity
# is exact (not merely approximate):
#   1. A single compound, concentration set exactly at its own HC5.
#   2. Multiple compounds that all share one sigma (e.g. all fall back to the
#      same group-level sigma_group_fallback), scaled so the summed HC5-PTI
#      equals exactly 1.
#
# This mirrors app.R's .calc_hc5() msPAF_EC50 block (see the comments there,
# citing de Zwart & Posthuma, 2005, Environ. Toxicol. Chem. 24(11):2665-2679)
# on a hand-built `denom`/`user_data` pair, without needing the Shiny reactive
# context or the full reference database.
#
# Run: Rscript Code/tests/test_mspaf_identity.R
# =============================================================================

suppressMessages({
  library(dplyr)
  library(tidyr)
})

# Mirrors app.R's .calc_hc5() msPAF_EC50 computation exactly (Docs
# TaxoTox_Technical_Methods.md Section 5.5; de Zwart & Posthuma, 2005, Eq. 4-6).
calc_mspaf <- function(denom, user_data) {
  conc_long <- denom %>%
    select(PREFERRED_NAME, median_conc, sigma_ec50) %>%
    left_join(user_data, by = "PREFERRED_NAME") %>%
    pivot_longer(cols = -c(PREFERRED_NAME, median_conc, sigma_ec50),
                 names_to = "Sample", values_to = "C") %>%
    mutate(TU = if_else(is.na(C) | is.na(median_conc), 0, C / median_conc))

  conc_long %>%
    group_by(Sample) %>%
    summarise(
      S       = sum(TU, na.rm = TRUE),
      S_known = sum(TU[!is.na(sigma_ec50)], na.rm = TRUE),
      sigma_mix = if (S_known > 0) {
        sum((TU / S_known) * sigma_ec50, na.rm = TRUE)
      } else {
        NA_real_
      },
      .groups = "drop"
    ) %>%
    mutate(
      msPAF_EC50 = dplyr::case_when(
        S_known <= 0     ~ 0,
        is.na(sigma_mix) ~ NA_real_,
        TRUE             ~ pnorm(log10(S_known) / sigma_mix)
      )
    )
}

TOL <- 1e-9
n_fail <- 0

check <- function(label, actual, expected, tol = TOL) {
  ok <- !is.na(actual) && abs(actual - expected) < tol
  status <- if (ok) "PASS" else "FAIL"
  cat(sprintf("[%s] %s: actual=%.12f expected=%.12f\n", status, label, actual, expected))
  if (!ok) n_fail <<- n_fail + 1
}

# ---------------------------------------------------------------------------
# Case 1: single compound, concentration set exactly at its own HC5.
# HC5 = median_conc * 10^(qnorm(0.05) * sigma) -- the same log-normal SSD
# relationship used to derive sigma_group_fallback in taxotox_install.R.
# ---------------------------------------------------------------------------
sigma1  <- 0.6
median1 <- 1000
hc5_1   <- median1 * 10^(qnorm(0.05) * sigma1)

denom1 <- tibble(
  PREFERRED_NAME = "CompoundA",
  median_conc    = median1,
  sigma_ec50     = sigma1
)
user_data1 <- tibble(PREFERRED_NAME = "CompoundA", Sample1 = hc5_1)

res1 <- calc_mspaf(denom1, user_data1)
check("Single-compound at HC5", res1$msPAF_EC50[res1$Sample == "Sample1"], 0.05)

# ---------------------------------------------------------------------------
# Case 2: multiple compounds sharing one sigma (uniform-sigma mixture),
# concentrations scaled so the summed HC5-based PTI = 1 exactly, i.e.
# sum(C_i / HC5_i) = 1 with HC5_i = median_i * 10^(qnorm(0.05)*sigma) for all i.
# By construction S = sum(C_i/median_i) = 10^(-qnorm(0.05)*sigma) * sum(C_i/HC5_i)
#                    = 10^(-qnorm(0.05)*sigma) [since sum(C_i/HC5_i) = 1]
# so log10(S)/sigma = -qnorm(0.05), and Phi(-qnorm(0.05)) = 0.05.
# ---------------------------------------------------------------------------
sigma2   <- 0.9
medians2 <- c(500, 2000, 50)
hc5s2    <- medians2 * 10^(qnorm(0.05) * sigma2)
tu_shares <- c(0.5, 0.3, 0.2)  # arbitrary shares of the HC5-PTI=1 total
concs2   <- tu_shares * hc5s2  # C_i such that sum(C_i/HC5_i) = sum(tu_shares) = 1

denom2 <- tibble(
  PREFERRED_NAME = c("CompoundB", "CompoundC", "CompoundD"),
  median_conc    = medians2,
  sigma_ec50     = sigma2
)
user_data2 <- tibble(PREFERRED_NAME = denom2$PREFERRED_NAME, Sample1 = concs2)

res2 <- calc_mspaf(denom2, user_data2)
check("Uniform-sigma mixture at HC5-PTI=1", res2$msPAF_EC50[res2$Sample == "Sample1"], 0.05)

# ---------------------------------------------------------------------------
# Case 3 (sanity, not the exact identity): no detections -> msPAF_EC50 = 0
# ---------------------------------------------------------------------------
user_data3 <- tibble(PREFERRED_NAME = denom2$PREFERRED_NAME, Sample1 = c(0, 0, 0))
res3 <- calc_mspaf(denom2, user_data3)
check("Zero concentrations -> msPAF_EC50 = 0", res3$msPAF_EC50[res3$Sample == "Sample1"], 0)

# ---------------------------------------------------------------------------
# Case 4 (sanity): a compound with unknown sigma_ec50, at TU=0, is excluded
# and the remaining weights renormalise correctly (should reduce to Case 1's
# single-known-compound result).
# ---------------------------------------------------------------------------
denom4 <- tibble(
  PREFERRED_NAME = c("CompoundA", "CompoundUnknownSigma"),
  median_conc    = c(median1, 100),
  sigma_ec50     = c(sigma1, NA_real_)
)
user_data4 <- tibble(PREFERRED_NAME = denom4$PREFERRED_NAME, Sample1 = c(hc5_1, 0))
res4 <- calc_mspaf(denom4, user_data4)
check("Unknown-sigma compound at TU=0 excluded cleanly", res4$msPAF_EC50[res4$Sample == "Sample1"], 0.05)

# ---------------------------------------------------------------------------
# Case 4b: same as Case 4, but the unknown-sigma compound now has a NONZERO
# concentration/TU. Regression test for the S vs. S_known fix: msPAF_EC50
# must use S_known (only known-sigma compounds) as the numerator, not S (all
# compounds) -- otherwise the unknown-sigma compound's TU would inflate the
# numerator while sigma_mix (computed over S_known only) stays blind to it,
# giving a numerator and spread that describe two different compound sets.
# Since the formula now uses S_known throughout, adding TU from an
# unknown-sigma compound must NOT change the result from Case 4's.
# ---------------------------------------------------------------------------
user_data4b <- tibble(PREFERRED_NAME = denom4$PREFERRED_NAME, Sample1 = c(hc5_1, 5000))
res4b <- calc_mspaf(denom4, user_data4b)
check("Unknown-sigma compound with nonzero TU still excluded from numerator",
      res4b$msPAF_EC50[res4b$Sample == "Sample1"], 0.05)

# ---------------------------------------------------------------------------
# Case 5: HC5 derived via the ACTUAL ssdtools/qlnorm natural-log convention,
# not the log10 formula Cases 1-2 assume. This is the path that a real
# ssd_fit_dists(dists="lnorm") + ssd_hc() fit produces (confirmed empirically
# against the installed ssdtools package: fitted `sdlog` matches a
# natural-log SD, i.e. log(X) ~ N(meanlog, sdlog), NOT log10(X) ~
# N(meanlog/ln(10)-ish, sdlog)). Regression test for the bug where
# sigma_ssd_lnorm (natural-log) was coalesced straight into sigma_ec50
# instead of being converted to log10 first (sigma_ssd_lnorm / log(10)) --
# Cases 1-2 cannot catch this because they never construct an HC5 via the
# natural-log qlnorm relationship in the first place.
# ---------------------------------------------------------------------------
sdlog_natural <- 1.15   # a plausible ssdtools-fitted natural-log sdlog
meanlog5      <- log(1000)   # median = 1000 on the natural (linear) scale
hc5_5         <- qlnorm(0.05, meanlog = meanlog5, sdlog = sdlog_natural)

# The value that must end up in sigma_ec50 for the msPAF formula to be
# correct is the LOG10 sigma, i.e. sdlog_natural / log(10) -- exactly the
# fix applied in taxotox_install.R's fit_hc5_ssd_full().
sigma_ec50_correct <- sdlog_natural / log(10)

denom5 <- tibble(
  PREFERRED_NAME = "CompoundE",
  median_conc    = 1000,
  sigma_ec50     = sigma_ec50_correct
)
user_data5 <- tibble(PREFERRED_NAME = "CompoundE", Sample1 = hc5_5)

res5 <- calc_mspaf(denom5, user_data5)
check("Single compound at ssdtools/qlnorm-derived HC5 (natural-log sdlog, converted)",
      res5$msPAF_EC50[res5$Sample == "Sample1"], 0.05)

# Sanity: confirm the UNCONVERTED natural-log sdlog would fail this check --
# guards against silently "fixing" the test instead of the bug.
denom5_bug <- denom5 %>% mutate(sigma_ec50 = sdlog_natural)
res5_bug   <- calc_mspaf(denom5_bug, user_data5)
mspaf_if_bugged <- res5_bug$msPAF_EC50[res5_bug$Sample == "Sample1"]
if (abs(mspaf_if_bugged - 0.05) < TOL) {
  cat("[FAIL] Sanity check: unconverted natural-log sigma should NOT satisfy the identity, but did.\n")
  n_fail <- n_fail + 1
} else {
  cat(sprintf("[PASS] Sanity check: unconverted natural-log sigma correctly gives the WRONG msPAF_EC50 (%.4f, not 0.05) -- confirms Case 5 actually exercises the bug.\n",
              mspaf_if_bugged))
}

cat(sprintf("\n%d check(s) failed.\n", n_fail))
if (n_fail > 0) quit(status = 1)
