library(stratEst)

test_that("Example: rock-paper-scissors" , {
  set.seed(1)
  strategies.mixture = list( "nash" = strategies.RPS$nash, "imitate" = strategies.RPS$imitate )
  model.mixture <- stratEst.model(data.WXZ2014,strategies.mixture)
  expect_equal(-22358.43,round(model.mixture$loglike,2))
})

test_that("Example: Dal Bo & Frechette 2011",  {
  set.seed(1)
  data <- stratEst.data(DF2011,input =c("choice","other.choice"), input.lag = 1 )
  model.DF2011 <- stratEst.model(data,strategies.DF2011,sample.id="treatment",verbose=F)
  expect_equal(0.920,round(as.numeric(model.DF2011$shares$treatment.D5R32[1]),3))
  expect_equal(0.080,round(as.numeric(model.DF2011$shares$treatment.D5R32[4]),3))
  expect_equal(0.043,round(as.numeric(model.DF2011$shares.se[1]),3))
  expect_equal(0.043,round(as.numeric(model.DF2011$shares.se[4]),3))
  expect_equal(0.362,round(as.numeric(model.DF2011$gammas.par[1]),3))
})


test_that("Example: Fudenberg, Rand, Dreber (2012)" , {
  set.seed(1)
  data <- stratEst.data(FRD2012,input =c("last.choice","last.other") )
  model.FRD2012 <- stratEst.model(data,strategies.FRD2012,sample.id="bc",verbose=F)
  expect_equal(0.193,round(as.numeric(model.FRD2012$shares$bc.1.5[2]),3))
  expect_equal(0.139,round(as.numeric(model.FRD2012$shares$bc.1.5[7]),3))
  expect_equal(0.053,round(as.numeric(model.FRD2012$shares.se[2]),3))
  expect_equal(0.035,round(as.numeric(model.FRD2012$shares.se[3]),3))
  expect_equal(0.46,round(as.numeric(model.FRD2012$gammas.par[1]),2))
})







