library(stratEst)

test_that("systematic tests",  {
  # shares
  one_fixed_share <- stratEst(DF2011[DF2011[,1]==1,],3,shares=c(0.5,NA,NA))
  two_fixed_shares <- stratEst(DF2011[DF2011[,1]==1,],3,shares=c(0.5,0.5,NA))
  # cluster
  cluster_bs <- stratEst(DF2011LCR[DF2011LCR[,1] == 1 & DF2011LCR[,4]==32 & DF2011LCR[,5]==0.5,],2,cluster = DF2011LCR$group[DF2011LCR[,1] == 1 & DF2011LCR[,4]==32 & DF2011LCR[,5]==0.5], bs.samples = 100 )
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
  select_strategies_aic <- stratEst(DF2011[DF2011[,1]==6,],3,select = "strategies", crit ="aic")
  expect_equal(2,length(select_strategies_aic$shares))
  select_strategies_bic <- stratEst(DF2011[DF2011[,1]==6,],3,select = "strategies", crit ="bic")
  expect_equal(2,length(select_strategies_bic$shares))
  select_strategies_icl <- stratEst(DF2011[DF2011[,1]==6,],3,select = "strategies", crit ="icl")
  expect_equal(2,length(select_strategies_icl$shares))
  select_trembles_global <- stratEst(DF2011[DF2011[,1]==6,],3,response="pure",r.trembles = "global",select = "trembles", crit ="icl")
  expect_equal(3,length(select_trembles_global$trembles))
  select_responses_global <- stratEst(DF2011[DF2011[,1]==6,],3,r.responses = "global",select = "responses", crit ="icl")
  expect_equal(4,length(select_responses_global$responses))
  select_both_global <- stratEst(DF2011[DF2011[,1]==6,],3,r.trembles = "global",r.responses = "global",select = "both", crit ="icl")
  expect_equal(4,length(select_both_global$responses))
  select_all_global <- stratEst(DF2011[DF2011[,1]==6,],3,r.trembles = "global",r.responses = "global",select = "all", crit ="icl")
  expect_equal(4,length(select_all_global$responses))
  expect_equal(3,length(select_all_global$shares))
})

