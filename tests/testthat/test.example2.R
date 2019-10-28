library(stratEst)

test_that("test example 2 of vignette",  {
  set.seed(1)
  id <- c(62,62,62,62,87,87,87,87)
  game <- c(4,4,4,4,4,4,4,4)
  period <- c(1,2,3,4,1,2,3,4)
  input <- c(0,1,2,3,0,1,3,2)
  output <- c(2,2,1,2,2,1,2,1)
  data <- as.data.frame(cbind(id,game,period,input,output))

  state <- c(1,2,3,1,2)
  r1 <- c(0.5,0,1,0.1,NA)
  r2 <- c(0.5,1,0,0.9,NA)
  t1 <- c(2,2,2,2,1)
  t2 <- c(3,3,3,2,1)
  t3 <- c(2,2,2,2,1)

  strategies <- as.data.frame(cbind(state,r1,r2,t1,t2,t3),row.names = c("S1.1","S1.2","S1.3","S2.1","S2.2"))

  model <- stratEst(data,strategies,print.messages = F)
  expect_equal(0.5,model$shares[1])
  expect_equal(0.5,model$shares[2])
  expect_equal(1,model$strategies$r2[2])
  expect_equal(0,model$strategies$r2[3])
})


