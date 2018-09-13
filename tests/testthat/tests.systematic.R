library(stratEst)

test_that("systematic tests",  {
  # response
  response_pure <- stratEst(DF2011[DF2011[,1]==1,],3,response= "pure")
  response_mixed <- stratEst(DF2011[DF2011[,1]==1,],3,response= "mixed")
  # restrict trembles
  restrict_trembles_no <- stratEst(DF2011[DF2011[,1]==1,],strategies=rbind(ALLD,TFT,WSLS),r.trembles="no")
  expect_equal(5,length(restrict_trembles_no$trembles))
  restrict_trembles_global <- stratEst(DF2011[DF2011[,1]==1,],strategies=rbind(ALLD,TFT,WSLS),r.trembles="global")
  expect_equal(1,length(restrict_trembles_global$trembles))
  restrict_trembles_strats <- stratEst(DF2011[DF2011[,1]==1,],strategies=rbind(ALLD,TFT,WSLS),r.trembles="strategies")
  expect_equal(3,length(restrict_trembles_strats$trembles))
  restrict_trembles_states <- stratEst(DF2011[DF2011[,1]==1,],strategies=rbind(ALLD,TFT,WSLS),r.trembles="states")
  expect_equal(2,length(restrict_trembles_states$trembles))
  # restrict responses
  restrict_responses_no <- stratEst(DF2011[DF2011[,1]==1,],3,r.responses="no")
  expect_equal(15,length(restrict_responses_no$responses))
  restrict_responses_global <- stratEst(DF2011[DF2011[,1]==1,],3,r.responses="global")
  expect_equal(1,length(restrict_responses_global$responses))
  restrict_responses_strategies <- stratEst(DF2011[DF2011[,1]==1,],3,r.responses="strategies")
  expect_equal(3,length(restrict_responses_strategies$responses))
  restrict_responses_states <- stratEst(DF2011[DF2011[,1]==1,],3,r.responses="states")
  expect_equal(5,length(restrict_responses_states$responses))
  # select
  select_strategies <- stratEst(DF2011[DF2011[,1]==1,],4,select = "strategies", crit ="aic")
  expect_equal(3,length(select_strategies$shares))


  #select_responses <- stratEst(DF2011[DF2011[,1]==1,],strategies=rbind(ALLD,TFT,WSLS),select = "all", crit ="icl")
  #expect_equal(3,lenght(select_strats$shares))
})

