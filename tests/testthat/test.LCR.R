library(stratEst)

test_that("latent class regression",  {
  strats = rbind(ALLD,ALLC,GRIM,TFT,WSLS,T2)
  # treatment 1
  set.seed(1)
  data <- DF2011all[DF2011all$delta==0.5 & DF2011all$r==32 & DF2011all$supergame <= 20,]
  covar_mat = matrix(data$supergame-1,ncol=1)
  lcr1 <- stratEst(data, strats, covariates=covar_mat, select="strategies", se="yes", r.trembles ="global", inner.runs = 10, outer.runs = 10, lcr.runs = 500)
  estimated_shares <- round(lcr1$shares,3)
  estimated_ses <- round(lcr1$shares.se,3)
  estimated_coefficients <- round(lcr1$coefficients,3)
  estimated_coefficients_ses <- round(lcr1$coefficients.se,3)
  estimated_trembles <- round(lcr1$trembles[1],3)
  estimated_trembles_ses <- round(lcr1$trembles.se[1],3)
  expect_equal(0.936,estimated_shares[1])
  expect_equal(0.064,estimated_shares[2])
  expect_equal(0.012,estimated_ses[1])
  expect_equal(0.012,estimated_ses[2])
  expect_equal(0.106,estimated_trembles[1])
  expect_equal(0.007,estimated_trembles_ses[1])
  expect_equal(-1.732,estimated_coefficients[1])
  expect_equal(-0.125,estimated_coefficients[2])
  expect_equal(0.153,estimated_coefficients_ses[1])
  expect_equal(0.021,estimated_coefficients_ses[2])
})

