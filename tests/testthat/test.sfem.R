library(stratEst)

test_that("SFEM for treatment 1, Dal Bo & Frechette 2011",  {
  strats <- rbind(ALLD,ALLC,GRIM,TFT,WSLS,T2)
  # treatment 1
  set.seed(1)
  sfem1 <- stratEst(DF2011[DF2011[,1]==1,],strats,select = "strategies", crit ="aic", se= "bs",print.messages = F)
  estimated_shares <- round(sfem1$shares,3)
  estimated_trembles <- round(sfem1$trembles,3)
  estimated_ses <- round(sfem1$shares.se,3)
  estimated_trembles_ses <- round(sfem1$trembles.se,3)
  expect_equal(0.920,estimated_shares[1])
  expect_equal(0.080,estimated_shares[2])
  expect_equal(0.043,estimated_ses[1])
  expect_equal(0.043,estimated_ses[2])
  expect_equal(0.060,estimated_trembles[1])
  expect_equal(0.020,estimated_trembles_ses[1])
})


