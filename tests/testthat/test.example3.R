library(stratEst)

test_that("test example 3 of vignette",  {
  set.seed(1)
  ## Latent class regression with data from Dal Bo and Frechette (2011)
  ## For the two treatments with R = 32, introduce a dummy which is one if delta = 3/4
  dummy <- as.numeric(DF2011$treatment > 3 )
  data <- as.data.frame(cbind(DF2011,dummy))
  strats <- rbind(ALLD,TFT)
  model <- stratEst(data,strats,covariates = c("dummy"),lcr.runs = 100,penalty = F,print.messages = F)
  expect_equal(0.52,round(model$shares[1],2))
  expect_equal(0.48,round(model$shares[2],2))
  expect_equal(-1.12,round(model$coefficients[1],2))
  expect_equal(2.19,round(model$coefficients[2],2))
  model <- stratEst(data,strats,covariates = c("dummy"),lcr.runs = 100,penalty = T,print.messages = F)
  expect_equal(0.52,round(model$shares[1],2))
  expect_equal(0.48,round(model$shares[2],2))
  expect_equal(-1.11,round(model$coefficients[1],2))
  expect_equal(2.17,round(model$coefficients[2],2))
})


