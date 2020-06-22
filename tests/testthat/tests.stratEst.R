library(stratEst)

test_that("core systematic tests",  {
  # unsorted data
  set.seed(1)
  sorted_data <- data.DF2011[data.DF2011$treatment=="D5R48",]
  unsorted_data <- sorted_data[sample(nrow(sorted_data)),]
  unsorted <- stratEst.model(unsorted_data,strategies=strategies.DF2011,inner.runs = 5, outer.runs = 2,sample.id = "treatment", verbose = F)
  set.seed(1)
  sorted <- stratEst.model(sorted_data,strategies=strategies.DF2011,inner.runs = 5, outer.runs = 2,sample.id = "treatment", verbose = F)
  unsorted_shares <- round(as.numeric(unsorted$shares),3)
  sorted_shares <- round(as.numeric(sorted$shares),3)
  expect_equal(sorted_shares[1],unsorted_shares[1])
  expect_equal(sorted_shares[2],unsorted_shares[2])
  expect_equal(sorted_shares[4],unsorted_shares[4])
  expect_equal(sorted_shares[5],unsorted_shares[5])

  # # fix shares
  data <- data.DF2011[data.DF2011$treatment=="D5R40",]
  not_fixed <- stratEst.model(data,3,inner.runs = 10, outer.runs = 2, verbose = F)
  not_fixed_samples <- stratEst.model(data,3,sample.id="treatment",inner.runs = 10, outer.runs = 2,verbose = F)
  one_fixed_share <- stratEst.model(data,3,shares=c(0.5,NA,NA),inner.runs = 10, outer.runs = 5,verbose = F)
  expect_equal(round(as.numeric(one_fixed_share$shares[1]),1),0.5)
  two_fixed_shares <- stratEst.model(data,3,shares=c(0.5,0.2,NA),inner.runs = 10, outer.runs = 5,verbose = F)
  expect_equal(round(as.numeric(two_fixed_shares$shares[c(1,2)]),1),c(0.5,0.2))
  all_fixed_shares <- stratEst.model(data,3,shares=c(0.5,0.2,0.3),inner.runs = 10, outer.runs = 5,verbose = F)
  expect_equal(round(as.numeric(all_fixed_shares$shares),1),c(0.5,0.2,0.3))
  some_fixed_samples <- stratEst.model(data.DF2011[data.DF2011$treatment %in% c("D5R32","D5R40","D5R48"),],3,shares=list("treatment.D5R32" = c(0.5,NA,NA),"treatment.D5R40" = c(NA,0.2,NA), "treatment.D5R48" = c(NA,NA,0.1)),sample.id="treatment",inner.runs = 10, outer.runs = 2,verbose = F)
  expect_equal(round(c(as.numeric(some_fixed_samples$shares$treatment.D5R32[1]),as.numeric(some_fixed_samples$shares$treatment.D5R40[2]),as.numeric(some_fixed_samples$shares$treatment.D5R48[3])),1),c(0.5,0.2,0.1))
  all_fixed_samples <- stratEst.model(data.DF2011[data.DF2011$treatment %in% c("D5R32","D5R40","D5R48"),],3,shares=list("treatment.D5R32" = c(0.5,0.2,0.3),"treatment.D5R40" = c(0.6,0.2,0.2), "treatment.D5R48" = c(0.3,0.6,0.1)),sample.id="treatment",inner.runs = 10, outer.runs = 2,verbose = F)
  expect_equal(round(as.numeric(unlist(all_fixed_samples$shares)),1),c(0.5,0.2,0.3,0.6,0.2,0.2,0.3,0.6,0.1))

  # fix probs

  # fix trembles


  # restrict trembles
  data <- data.DF2011[data.DF2011$treatment=="D5R40",]
  restrict_trembles_no <- stratEst.model(data,strategies.PD[c("ALLD","TFT","WSLS")],r.trembles="no",verbose = F)
  expect_equal(5,length(restrict_trembles_no$trembles.par))
  restrict_trembles_global <- stratEst.model(data,strategies.PD[c("ALLD","TFT","WSLS")],r.trembles="global",verbose = F)
  expect_equal(1,length(restrict_trembles_global$trembles.par))
  restrict_trembles_strats <- stratEst.model(data,strategies.PD[c("ALLD","TFT","WSLS")],r.trembles="strategies",verbose = F)
  expect_equal(3,length(restrict_trembles_strats$trembles.par))
  restrict_trembles_states <- stratEst.model(data,strategies.PD[c("ALLD","TFT","WSLS")],r.trembles="states",verbose = F)
  expect_equal(2,length(restrict_trembles_states$trembles.par))

  # restrict probs
  data <- data.DF2011[data.DF2011$treatment=="D5R40",]
  restrict_probs_no <- stratEst.model(data,3,r.probs="no",verbose = F)
  expect_equal(30,length(restrict_probs_no$probs.par))
  restrict_probs_global <- stratEst.model(data,3,r.probs="global",verbose = F)
  expect_equal(2,length(restrict_probs_global$probs.par))
  restrict_probs_strategies <- stratEst.model(data,3,r.probs="strategies",verbose = F)
  expect_equal(6,length(restrict_probs_strategies$probs.par))
  restrict_probs_states <- stratEst.model(data,3,r.probs="states",verbose = F)
  expect_equal(10,length(restrict_probs_states$probs.par))

  # cluster
  #cluster_bs <- stratEst.model(DF2011LCR[DF2011LCR[,1] == 1 & DF2011LCR[,4]==32 & DF2011LCR[,5]==0.5,],2,cluster.id = "group", bs.samples = 50, verbose = F )

  # select
  set.seed(1)
  data <- data.DF2011[data.DF2011$treatment=="D75R48",]
  select_strategies_icl <- stratEst.model(data,3,select = "strategies", crit ="icl",verbose = F)
  expect_equal(2,length(select_strategies_icl$shares))
  set.seed(1)
  select_trembles_global <- stratEst.model(data,3,response="pure",r.trembles = "global",select = "trembles", crit ="icl",verbose = F)
  expect_equal(3,length(select_trembles_global$trembles.par))
  set.seed(1)
  select_all_global <- stratEst.model(data,3,r.trembles = "global",r.probs = "global",select = c("strategies","trembles","probs"), crit ="icl",verbose = F)
  expect_equal(12,length(select_all_global$probs.par))
  expect_equal(3,length(select_all_global$shares.par))
})

test_that("additional systematic tests",  {
  skip_on_cran()
  # response
  data <- data.DF2011[data.DF2011$treatment=="D5R32",]
  set.seed(1)
  response_pure <- stratEst.model(data,3,response= "pure",verbose = F)
  set.seed(1)
  response_mixed <- stratEst.model(data,3,response= "mixed",verbose = F)
  expect_equal(1,round(as.numeric(all( response_pure$probs.par == 0 | response_pure$probs.par == 1 ))))


  # #select
  # set.seed(1)
  # data <- data.DF2011[data.DF2011$treatment=="D75R48",]
  # select_strategies_aic <- stratEst.model(data,3,select = "strategies", crit ="aic",verbose = F)
  # expect_equal(2,round(as.numeric(length(select_strategies_aic$shares))))
  # set.seed(1)
  # select_strategies_bic <- stratEst.model(data,3,select = "strategies", crit ="bic",verbose = F)
  # expect_equal(2,round(as.numeric(length(select_strategies_bic$shares))))
  # set.seed(1)
  # select_trembles_global <- stratEst.model(data,3,response="pure",r.trembles = "global",select = "trembles", crit ="icl",verbose = F)
  # expect_equal(4,round(as.numeric(length(select_trembles_global$trembles.par))))
  # set.seed(1)
  # select_probs_global <- stratEst.model(data,3,r.probs = "global",select = "probs", crit ="icl",verbose = F)
  # expect_equal(8,round(as.numeric(length(select_probs_global$probs.par))))
  # set.seed(1)
  # select_both_global <- stratEst.model(data,3,r.trembles = "global",r.probs = "global",select = c("probs","trembles"), crit ="icl",verbose = F)
  # expect_equal(16,round(as.numeric(length(select_both_global$probs.par))))
})

# test_that("multivariate output test",  {
#   skip_on_cran()
#   N = 32                                       # number of subjects
#   num_probs = c(2,3,4,5)                   # number of distinct probs of the machines
#   tremble = 0.2
#   shares = runif(5)                            # generate 5 shares
#   shares = shares/sum(shares)                  # normalize shares
#   for( m in 1:length(num_probs)){
#     check = 0
#     while( check == 0 ){
#       strategy_mat = matrix(runif(5*5*num_probs[m]),5*5,num_probs[m])
#       strategy_mat = matrix(as.numeric(strategy_mat == apply(strategy_mat,1,max)),5*5,num_probs[m])
#       if( sum(as.numeric(apply(strategy_mat,1,sum) == 1)) == 5*5 ){ check = 1 }
#     }
#
#     strats = cbind(rep(c(1:5),5),strategy_mat,matrix(unlist(lapply(c(2:5), function(x) rep(x,5*5) )),5*5,5-1))
#
#     response_mat = (1-strategy_mat)*tremble/(num_probs[m]-1) + strategy_mat*(1-tremble)
#
#     strategy_vec = t(rmultinom(N, size = 1, prob = shares))%*%matrix(c(1:5),5,1)
#     strat_id = rep( strategy_vec , each= (15*5))
#
#     real_shares = rep(NA,5)
#     for( i in 1:5 ){
#       real_shares[i] = sum(as.numeric(strategy_vec==i))/length(strategy_vec)
#     }
#
#     id = rep(NA,15*5*N)
#     game = rep( unlist(lapply(c(1:15), function(x) rep(x,5) )) ,N )
#     period = rep(c(1:5),15*N)
#     input =  rep(c(1:5)-1,15*N)
#     output = rep(NA,15*5*N)
#     rand_vec = runif(15*5*N)
#     pr_vec = rep(NA,15*5*N)
#
#     for( i in 1:N ){
#       id[ ((i-1)*15*5 + 1) : ((i-1)*15*5 + 15*5) ] = rep(i,15*5)
#     }
#
#     for( j in 1:(15*5*N) ){
#       output[ j ] = t(rmultinom(1, size = 1, prob = response_mat[1+input[j] + (strat_id[j]-1)*5 , ] ))%*%matrix(c(1:num_probs[m]),num_probs[m],1)
#     }
#
#     data <- as.data.frame(cbind(id,game,period,input,output))
#     data <- stratEst.data(data,input.vars = "input" , output = "output")
#
#     P = stratEst.model(data,5,response="pure",outer.runs = 2, verbose = F)
#     M = stratEst.model(data,5,response="mixed",outer.runs = 2, verbose = F)
#
#     r_string = paste("r",as.character(num_probs[m]),sep="")
#
#     P_l = nrow( P$probs )
#     M_l = nrow( M$probs )
#
#     expect_equal( 25 , round(as.numeric(P_l)) )
#     expect_equal( 25  , round(as.numeric(M_l)) )
#
#   }
# })

