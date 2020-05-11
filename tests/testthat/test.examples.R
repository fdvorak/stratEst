library(stratEst)

test_that("Example: rock-paper-scissors" , {
  set.seed(1)
  data <- stratEst.data( data = RPS,input.vars = c("choice"), output = "choice", input.lag = 1 )
  mixed = stratEst.strategy( inputs = c(NA,"rock","paper","scissors") , outputs = c("rock","paper","scissors") , num.states = 1 , responses = NA )
  cyclic = stratEst.strategy( inputs = c(NA,"rock","paper","scissors") , outputs = c("rock","paper","scissors") , responses = c(1/3,1/3,1/3,0,0,1,1,0,0,0,1,0)  )
  anti = stratEst.strategy( inputs = c(NA,"rock","paper","scissors") , outputs = c("rock","paper","scissors") , responses = c(1/3,1/3,1/3,0,1,0,0,0,1,1,0,0)  )
  model.mixed <- stratEst( data , list("mixed" = mixed) , outer.runs = 1, verbose = F )
  model.conditional <- stratEst( data.RPS , list( "anti" = anti , "cyclic" = cyclic ) , outer.runs = 1,verbose=F )
  expect_equal(-10974.59,round(model.mixed$loglike,2))
  expect_equal(-10945.63,round(model.conditional$loglike,2))
})

test_that("Example: Dal Bo & Frechette 2011",  {
  set.seed(1)
  data <- stratEst.data(DF2011,input.vars =c("choice","other_choice"),output="choice", input.lag = 1 )
  model.DF2011 <- stratEst(data,strategies.DF2011,sample.id="treatment",verbose=F)
  expect_equal(0.920,round(as.numeric(model.DF2011$shares$treatment.D5R32[1]),3))
  expect_equal(0.080,round(as.numeric(model.DF2011$shares$treatment.D5R32[4]),3))
  expect_equal(0.043,round(as.numeric(model.DF2011$shares.se[1]),3))
  expect_equal(0.043,round(as.numeric(model.DF2011$shares.se[4]),3))
  expect_equal(0.362,round(as.numeric(model.DF2011$gammas.par[1]),3))
})


test_that("Example: Fudenberg, Rand, Dreber (2012)" , {
  set.seed(1)
  data <- stratEst.data(FRD2012,input.vars =c("last_choice","last_other"),output="choice" )
  model.FRD2012 <- stratEst(data,strategies.FRD2012,sample.id="bc",verbose=F)
  expect_equal(0.193,round(as.numeric(model.FRD2012$shares$bc.1.5[2]),3))
  expect_equal(0.139,round(as.numeric(model.FRD2012$shares$bc.1.5[7]),3))
  expect_equal(0.053,round(as.numeric(model.FRD2012$shares.se[2]),3))
  expect_equal(0.035,round(as.numeric(model.FRD2012$shares.se[3]),3))
  expect_equal(0.46,round(as.numeric(model.FRD2012$gammas.par[1]),2))
})




