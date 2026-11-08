test_that("predictEventTime's CDF matches predictRisk and its quantiles invert it", {
  skip_on_cran()
  data(pbc3, envir = environment())
  ds <- pbc3[!duplicated(pbc3$id), ]
  history <- pbc3[pbc3$id %in% c(2, 7, 13) & pbc3$year <= 3, ]
  history$years <- NA
  for (cr in c(FALSE, TRUE)) {
    sf <- survivalSub(ds, Surv(years, status3) ~ age + sex, if (cr) status4 ~ years + age + sex)
    fx <- if (cr) list(serBilir ~ year + age + years + status4) else list(serBilir ~ year + age + years)
    lf <- longitudinalSub(pbc3[pbc3$status3 == 1, ], fx, list(~ year | id))
    out <- predictEventTime(history, lf, sf, 3, c(1, 3, 6), "year", probs = c(0.25, 0.5))
    expect_s3_class(out, "predictEventTime.BJM")
    expect_equal(rownames(out$cdf), c("2", "7", "13"))
    expect_true(all(apply(out$cdf, 1, diff) >= 0))
    for (h in c(1, 3)) {
      r <- predictRisk(history, lf, sf, 3, h, "year", NULL, NULL, bandcount1 = 200, bandcount2 = 100)
      total <- r$risk_prob_1 + if (cr) r$risk_prob_2 else 0
      expect_equal(unname(out$cdf[, paste0("h", h)]), unname(total[rownames(out$cdf)]), tolerance = 0.01)
    }
    if (cr) expect_equal(out$cif[[1]] + out$cif[[2]], out$cdf)
    expect_true(all(out$quantile[, 1] <= out$quantile[, 2]))
    expect_true(all(out$quantile > 3))
    # the median sits where the CDF crosses 1/2 (the Breslow baseline makes
    # the CDF jump at event times, so check either side of it)
    around <- predictEventTime(history, lf, sf, 3, out$quantile[1, 2] - 3 + c(-0.1, 0.1), "year")$cdf[1, ]
    expect_lte(around[[1]], 0.5)
    expect_gte(around[[2]], 0.5)
  }
  expect_output(print(out), "Quantiles of T")
  expect_error(predictEventTime(history, lf, sf, 3, -1, "year"), "positive")
})
