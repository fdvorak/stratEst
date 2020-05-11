# This function executes the post-processing after the estimation

stratEst.post <- function( cpp.output , stratEst.return , strategies , covariates , response , unique_ids , num_unique_ids , input , output , unique_inputs , unique_outputs , num_unique_inputs , num_unique_outputs , sample , sample.id , sample_factor , num_samples , specific_shares , specific_responses , specific_trembles , sample_is_factor , integer_strategies , LCR , response_mat_col_index , crit , num_obs , se , quantiles ){

  num_responses_to_est = length(stratEst.return$responses.par)
  num_trembles_to_est = length(stratEst.return$trembles.par)

  if( num_samples == 1 | ( (specific_responses == F | num_responses_to_est == 0) & (specific_trembles == F | num_trembles_to_est == 0) ) ){
    # strategies post-processing
    stratEst.return$strategies[ stratEst.return$strategies == -1 ] = NA
    num_strategies = sum(stratEst.return$strategies[,1] == 1)
    num_strategies_sample = num_strategies
    num_unique_inputs <- length(unique(input[input!=0]))
    post_sid <- cpp.output$sid
    unique_post_sids <- unique(post_sid)
    if( integer_strategies ){
      state <- stratEst.return$strategies[,1]
      response_mat <- matrix(stratEst.return$strategies[,2:(1+num_unique_outputs)],nrow(stratEst.return$strategies), num_unique_outputs)
      r_names <- rep(NA,num_unique_outputs)
      for( outs in 1:num_unique_outputs ){
        r_names[outs] <- paste("output.",as.character(unique_outputs[outs]),sep="")
      }
      colnames(response_mat) <- r_names

      transition_mat <- matrix(stratEst.return$strategies[, ((1+num_unique_outputs+1):(1+num_unique_outputs+num_unique_inputs))],nrow(stratEst.return$strategies), num_unique_inputs)
      t_names <- rep(NA,num_unique_inputs)
      for( ins in 1:num_unique_inputs ){
        t_names[ins] <- paste("input.",as.character(unique_inputs[ins]),sep="")
      }
      colnames(transition_mat) <- t_names
      tremble_vec <- cpp.output$trembles
      colnames(tremble_vec) <- "tremble"
      strategies_mat <- as.data.frame(cbind(response_mat,tremble_vec,transition_mat))

      names_of_strategies <- NULL
      for( s in 1:num_strategies ){
        names_of_strategies <- c(names_of_strategies,paste("strategy.",as.character(s),sep=""))
      }
      stratEst.return$strategies <- strategies_mat
    }
    else if( class(strategies) == "list" ){
      strategies_list <- list(NULL)
      names_of_strategies <- names(strategies)
      if( is.null(names_of_strategies) ){
        for( s in 1:length(strategies) ){
          names_of_strategies <- c(names_of_strategies,paste("strategy.",as.character(s),sep=""))
        }
      }
      for( strs in 1:num_strategies){
        strategy <- strategies[[unique_post_sids[strs]]]
        for( out in 1:num_unique_outputs ){
          strategy[,response_mat_col_index[unique_post_sids[strs],out]] <- stratEst.return$strategies[post_sid == unique_post_sids[strs],(1+out)]
        }
        trembles_strategy <- cpp.output$trembles[post_sid == unique_post_sids[strs]]
        if( sum(trembles_strategy > 0) > 0 ){
          strategy$tremble <- trembles_strategy
        }
        strategies_list[[strs]] <- strategy
      }
      stratEst.return$strategies <- strategies_list
    }

    # create list of strategies
    if( class(strategies) != "list" ){
      strategies_list <- list(NULL)
      for( i in 1:num_strategies ){
        strategy <- stratEst.return$strategies[ post_sid == unique_post_sids[i] , ]
        rownames(strategy) <- c(1:nrow(strategy))
        strategies_list[[i]] <- strategy
      }
    }

  }else{  # num_samples > 1
    # strategies post-processing
    stratEst.return$strategies[ stratEst.return$strategies == -1 ] = NA
    num_strategies = sum(stratEst.return$strategies[,1] == 1)
    num_strategies_sample <- (num_strategies/num_samples)
    num_unique_inputs <- length(unique(input[input!=0]))
    original_sid <- cpp.output$sid
    post_sid <- original_sid
    post_sid <- rep( post_sid , num_samples ) + rep( c(0:(num_samples-1)) , each = length(post_sid))*max(post_sid)
    unique_post_sids <- unique(post_sid)
    unique_original_sids <- unique(original_sid)
    if( integer_strategies ){
      #state <- stratEst.return$strategies[,1]
      response_mat <- cpp.output$responses #matrix(stratEst.return$strategies[,2:(1+num_unique_outputs)],nrow(stratEst.return$strategies), num_unique_outputs)
      r_names <- rep(NA,num_unique_outputs)
      for( outs in 1:num_unique_outputs ){
        r_names[outs] <- paste("output.",as.character(unique_outputs[outs]),sep="")
      }
      colnames(response_mat) <- r_names

      transition_mat <- matrix(stratEst.return$strategies[, ((1+num_unique_outputs+1):(1+num_unique_outputs+num_unique_inputs))],nrow(stratEst.return$strategies), num_unique_inputs)
      t_names <- rep(NA,num_unique_inputs)
      for( ins in 1:num_unique_inputs ){
        t_names[ins] <- paste("input.",as.character(unique_inputs[ins]),sep="")
      }
      colnames(transition_mat) <- t_names
      tremble_vec <- cpp.output$trembles
      colnames(tremble_vec) <- "tremble"
      strategies_mat <- as.data.frame(cbind(response_mat,tremble_vec,transition_mat))

      names_of_strategies <- NULL

      for( s in 1:num_strategies_sample ){
        names_of_strategies <- c(names_of_strategies,paste("strategy.",as.character(s),sep=""))
      }
      stratEst.return$strategies <- strategies_mat

    }
    else if( class(strategies) == "list" ){
      names_of_strategies <- names(strategies)
      if( is.null(names_of_strategies) ){
        for( s in 1:length(strategies) ){
          names_of_strategies <- c(names_of_strategies,paste("strategy.",as.character(unique_original_sids[s]),sep=""))
        }
      }
      strategies_list <- list(NULL)
      strategies_list_sample <- list(NULL)
      for( sam in 1:num_samples ){
        for( str in 1:num_strategies_sample){
          strategy <- strategies[[unique_post_sids[str]]]
          strs <- (sam-1)*num_strategies_sample + str
          for( out in 1:num_unique_outputs ){
            strategy[,response_mat_col_index[unique_post_sids[str],out]] <- stratEst.return$strategies[post_sid == unique_post_sids[strs],(1+out)]
          }
          trembles_strategy <- cpp.output$trembles[post_sid == unique_post_sids[strs]]

          if( sum(trembles_strategy > 0) > 0 ){
            strategy$tremble <- trembles_strategy
          }
          strategies_list_sample[[str]] <- strategy
        }
        names(strategies_list_sample) <- names_of_strategies[unique_original_sids]
        strategies_list[[sam]] <- strategies_list_sample
      }
      stratEst.return$strategies <- strategies_list
    }
  }

  if( class(strategies) == "list"){
    names_of_strategies <- names_of_strategies[ cpp.output$selected.strats ]
  }

  # shares post-processing
  shares_mat <- stratEst.return$shares
  names_of_samples <- rep(NA,num_samples)
  if( num_samples > 1 ){
    if( sample_is_factor ){
      unique_factors <- sort(unique(sample_factor))
      for( smps in 1:num_samples ){
        names_of_samples[smps] <- paste(sample.id,as.character(unique_factors[smps]),sep=".")
      }
    }else{
      unique_samples <- sort(unique(sample))
      for( smps in 1:num_samples ){
        names_of_samples[smps] <-paste(sample.id,as.character(unique_samples[smps]),sep=".")
      }
    }
  }
  if( num_samples == 1 | specific_shares == F ){
    frame_shares_sample <- as.data.frame(t(shares_mat[,1]))
    colnames(frame_shares_sample) <- names_of_strategies
    rownames(frame_shares_sample) <- "share"
    stratEst.return$shares <- frame_shares_sample
    stratEst.return$shares = as.matrix(stratEst.return$shares)
  }else{
    shares <- list(NULL)
    for( i in 1:num_samples ){
      frame_shares_sample <- as.data.frame(t(shares_mat[,i]))
      rownames(frame_shares_sample) <- names_of_samples[i]

      colnames(frame_shares_sample) <- names_of_strategies
      shares[[i]] <- frame_shares_sample
    }
    names(shares) <- names_of_samples
    stratEst.return$shares <- shares
  }

  if( length(stratEst.return$shares.par) > 0 ){
    rownames(stratEst.return$shares.indices) = names_of_strategies
    if( num_samples == 1 | specific_shares == F ){
      colnames(stratEst.return$shares.indices) = "share"
    }
    else{
      colnames(stratEst.return$shares.indices) = names_of_samples
    }
  }

  if( num_samples > 1 & ( (specific_responses & num_responses_to_est > 0) | (specific_trembles & num_trembles_to_est > 0) ) ){
    # create list of strategies
    if( class(strategies) != "list" ){
      strategies_list <- list(NULL)
      for( j in 1:num_samples ){
        strategies_list_sample <- list(NULL)
        counter_names <- 0
        for( i in 1:num_strategies_sample ){
          strategy <- stratEst.return$strategies[ post_sid == unique_post_sids[(num_strategies_sample*(j-1)+i)] , ]
          rownames(strategy) <- c(1:nrow(strategy))
          strats_sample_list <- strategy
          counter_names <- counter_names + nrow(strats_sample_list)
          strategies_list_sample[[i]] <- strats_sample_list
        }
        names(strategies_list_sample) <- names_of_strategies
        strategies_list[[j]] <- strategies_list_sample
      }
    }else{
      for( j in 1:num_samples ){
        strat_list_sample <- stratEst.return$strategies[[j]]
        names(strat_list_sample)  <- names_of_strategies
        stratEst.return$strategies[[j]] <- strat_list_sample
      }
    }
  }

  # assignments post-processing
  assignment_row_names <- NULL
  for( i in 1:num_unique_ids ){
    assignment_row_names <- c(assignment_row_names,paste( "id" , as.character(unique_ids[i]) , sep = "" ) )
  }
  rownames(stratEst.return$posterior.assignment) <- assignment_row_names
  colnames(stratEst.return$posterior.assignment) <- names_of_strategies

  # coefficients and priors post-processing
  stratEst.return$coefficients.par <- cpp.output$coefficients.list$coefficients.par
  if( LCR & num_strategies > 1 ){
    colnames(stratEst.return$coefficients) <- names_of_strategies
    rownames(stratEst.return$coefficients) <- covariates

    colnames(stratEst.return$prior.assignment) <- names_of_strategies
    rownames(stratEst.return$prior.assignment) <- assignment_row_names
  }

  # num.obs num.ids and free.par and residual degrees
  # num_obs = rep( NA , num_samples )
  # num_unique_ids = rep( NA , num_samples )
  # for( i in 1:num_samples ){
  #   num_obs[i] <- length(id[sample.id == i])
  #   num_unique_ids[i]  <- length(unique(id[ sample.id == i ]))
  # }
  stratEst.return$num.obs <- num_obs
  stratEst.return$num.ids <- num_unique_ids
  stratEst.return$num.par = length(stratEst.return$shares.par) + length(stratEst.return$responses.par) + length(stratEst.return$trembles.par) + length(stratEst.return$coefficients.par)
  stratEst.return$free.par = cpp.output$n.par
  stratEst.return$res.degrees = stratEst.return$num.ids - stratEst.return$free.par
  if( stratEst.return$res.degrees < 0 ){
    warning("stratEst warning: Residual degrees of freedom are negative. Parameter inference invalid.")
  }

  if( length(stratEst.return$shares.par) > 0 ){
    if( se == "bootstrap" ){
      stratEst.return$shares.quantiles = cpp.output$shares.list$shares.quantiles
    }
    else{
      shares.quantiles = matrix(NA,length(stratEst.return$shares.par),length(quantiles))
      for( i in 1:length(quantiles)){
        shares.quantiles[,i] =  stratEst.return$shares.par + stats::qt( quantiles[i] , df = stratEst.return$res.degrees )*stratEst.return$shares.se
      }
      stratEst.return$shares.quantiles = shares.quantiles
    }
    rownames(stratEst.return$shares.quantiles) <- paste("par.",as.character(seq(1,nrow(stratEst.return$shares.quantiles),by = 1)),sep="")
    colnames(stratEst.return$shares.quantiles) <- paste(as.character(quantiles*100),"%",sep="")
  }
  if( length(stratEst.return$responses.par) > 0 ){
    if( se == "bootstrap" ){
      stratEst.return$responses.quantiles = cpp.output$responses.list$responses.quantiles
    }
    else{
      responses.quantiles = matrix(NA,length(stratEst.return$responses.par),length(quantiles))
      for( i in 1:length(quantiles)){
        responses.quantiles[,i] =  stratEst.return$responses.par + stats::qt( quantiles[i] , df = stratEst.return$res.degrees )*stratEst.return$responses.se
      }
      stratEst.return$responses.quantiles = responses.quantiles
    }
    rownames(stratEst.return$responses.quantiles) <- paste("par.",as.character(seq(1,nrow(stratEst.return$responses.quantiles),by = 1)),sep="")
    colnames(stratEst.return$responses.quantiles) <- paste(as.character(quantiles*100),"%",sep="")
  }
  if( length(stratEst.return$trembles.par) > 0 ){
    if( se == "bootstrap" ){
      stratEst.return$trembles.quantiles = cpp.output$trembles.list$trembles.quantiles
    }
    else{
      trembles.quantiles = matrix(NA,length(stratEst.return$trembles.par),length(quantiles))
      for( i in 1:length(quantiles)){
        trembles.quantiles[,i] =  stratEst.return$trembles.par + stats::qt( quantiles[i] , df = stratEst.return$res.degrees )*stratEst.return$trembles.se
      }
      stratEst.return$trembles.quantiles = trembles.quantiles
    }
    rownames(stratEst.return$trembles.quantiles) <- paste("par.",as.character(seq(1,nrow(stratEst.return$trembles.quantiles),by = 1)),sep="")
    colnames(stratEst.return$trembles.quantiles) <- paste(as.character(quantiles*100),"%",sep="")
  }
  if( length(stratEst.return$coefficients.par) > 0 ){
    if( se == "bootstrap" ){
      stratEst.return$coefficients.quantiles = cpp.output$coefficients.list$coefficients.quantiles
    }
    else{
      coefficients.quantiles = matrix(NA,length(stratEst.return$coefficients.par),length(quantiles))
      for( i in 1:length(quantiles)){
        coefficients.quantiles[,i] =  stratEst.return$coefficients.par + stats::qt( quantiles[i] , df = stratEst.return$res.degrees )*stratEst.return$coefficients.se
      }
      stratEst.return$coefficients.quantiles = coefficients.quantiles
    }
    rownames(stratEst.return$coefficients.quantiles) <- paste("par.",as.character(seq(1,nrow(stratEst.return$coefficients.quantiles),by = 1)),sep="")
    colnames(stratEst.return$coefficients.quantiles) <- paste(as.character(quantiles*100),"%",sep="")
  }

  # gammas
  stratEst.return$gammas = -1/log(stratEst.return$trembles/(1-stratEst.return$trembles))
  if( length(stratEst.return$trembles.par) > 0 ){
    stratEst.return$gammas.par = -1/log(stratEst.return$trembles.par/(1-stratEst.return$trembles.par))
    stratEst.return$gammas.se = -1/log(stratEst.return$trembles.se/(1-stratEst.return$trembles.se))
  }

  # #convergence post-processing
  stratEst.return$convergence[stratEst.return$convergence == -1] = NA
  if( num_strategies == 1 | LCR ){
    stratEst.return$convergence[1] = NA
  }
  if( response == "pure"){
    stratEst.return$convergence[2] = NA
  }
  if( sum(is.na(stratEst.return$convergence)) < length(stratEst.return$convergence) ){
    shares <- stratEst.return$convergence[1]
    responses <- stratEst.return$convergence[2]
    trembles <- stratEst.return$convergence[3]
    coefficients <- stratEst.return$convergence[4]
    convergence_vec <- cbind(shares,responses,trembles,coefficients)
    stratEst.return$convergence <- matrix(convergence_vec,1,length(convergence_vec))
    rownames(stratEst.return$convergence) <- "max.abs.score"
    colnames(stratEst.return$convergence) <- c("shares","responses","trembles","coefficients")
    #stratEst.return$convergence <- stratEst.return$convergence[, colSums(is.na(stratEst.return$convergence)) != nrow(stratEst.return$convergence)]
  }else{
    stratEst.return$convergence = NULL
  }

  convergence_string <- "no parameters estimated"

  # add names to strategies list
  if( class(strategies) != "list" ){
    stratEst.return$strategies <- strategies_list
    if( num_samples == 1 | ( (specific_responses == F | num_responses_to_est == 0) & (specific_trembles == F | num_trembles_to_est == 0) ) ){
      names(stratEst.return$strategies) <- names_of_strategies
    }else{
      names(stratEst.return$strategies) <- names_of_samples

    }
  }
  else{
    if( num_samples == 1 | ( (specific_responses == F | num_responses_to_est == 0) & (specific_trembles == F | num_trembles_to_est == 0) ) ){
      names(stratEst.return$strategies) <- names_of_strategies
    }
    else{
      names(stratEst.return$strategies) <- names_of_samples

    }
  }

  # strategy matrix row names
  if( "list" %in% class(stratEst.return$strategies[[1]])  ){
    strategies_sample_list <- NULL
    strategies_print <- stratEst.return$strategies
    names_samples <- names(stratEst.return$strategies)
    for( i in 1:length(stratEst.return$strategies) ){
      strategies_sample <- round(do.call(rbind,strategies_print[[i]]),3)
      row_names_strategies_sample <- row.names(strategies_sample)
      row.names(strategies_sample) =  paste( names_samples[i], "." , row_names_strategies_sample, sep="")
      strategies_sample_list <- rbind( strategies_sample_list , strategies_sample )
    }
    strategies_matrix_names = rownames(strategies_sample_list)
  }else{
    strategies_matrix_names = rownames(do.call(rbind,stratEst.return$strategies))
  }

  # state.obs post-processing
  rownames(stratEst.return$state.obs) <- strategies_matrix_names
  colnames(stratEst.return$state.obs) <- "weighted.obs"

  rownames(stratEst.return$responses) = strategies_matrix_names
  r_names <- rep(NA,num_unique_outputs)
  for( outs in 1:num_unique_outputs ){
    r_names[outs] <- paste("output.",as.character(unique_outputs[outs]),sep="")
  }
  colnames(stratEst.return$responses) = r_names

  rownames(stratEst.return$trembles) = strategies_matrix_names
  colnames(stratEst.return$trembles) = "tremble"

  rownames(stratEst.return$gammas) = strategies_matrix_names
  colnames(stratEst.return$gammas) = "gamma"

  if( length( stratEst.return$shares.par ) > 0 ){
    row_names_shares = paste("par.",as.character(seq(1,length(stratEst.return$shares.par),by=1)),sep="")
    rownames(stratEst.return$shares.par) = row_names_shares
    rownames(stratEst.return$shares.se) = row_names_shares
    rownames(stratEst.return$shares.score) = row_names_shares
    rownames(stratEst.return$shares.covar) = row_names_shares
    rownames(stratEst.return$shares.fisher) = row_names_shares

    colnames(stratEst.return$shares.par) = "probability"
    colnames(stratEst.return$shares.se) = "standard error"
    colnames(stratEst.return$shares.score) = "score"
    colnames(stratEst.return$shares.covar) = row_names_shares
    colnames(stratEst.return$shares.fisher) = row_names_shares
  }
  if( length( stratEst.return$responses.par ) > 0 ){
    rownames(stratEst.return$responses.indices) = strategies_matrix_names
    colnames(stratEst.return$responses.indices) = r_names

    row_names_responses = paste("par.",as.character(seq(1,length(stratEst.return$responses.par),by=1)),sep="")
    rownames(stratEst.return$responses.par) = row_names_responses
    rownames(stratEst.return$responses.se) = row_names_responses
    rownames(stratEst.return$responses.score) = row_names_responses
    rownames(stratEst.return$responses.covar) = row_names_responses
    rownames(stratEst.return$responses.fisher) = row_names_responses

    colnames(stratEst.return$responses.par) = "probability"
    colnames(stratEst.return$responses.se) = "standard error"
    colnames(stratEst.return$responses.score) = "score"
    colnames(stratEst.return$responses.covar) = row_names_responses
    colnames(stratEst.return$responses.fisher) = row_names_responses
  }
  if( length( stratEst.return$trembles.par ) > 0 ){
    rownames(stratEst.return$trembles.indices) = strategies_matrix_names
    colnames(stratEst.return$trembles.indices) = "tremble"

    row_names_trembles = paste("par.",as.character(seq(1,length(stratEst.return$trembles.par),by=1)),sep="")
    rownames(stratEst.return$trembles.par) = row_names_trembles
    rownames(stratEst.return$trembles.se) = row_names_trembles
    rownames(stratEst.return$trembles.score) = row_names_trembles
    rownames(stratEst.return$trembles.covar) = row_names_trembles
    rownames(stratEst.return$trembles.fisher) = row_names_trembles

    rownames(stratEst.return$gammas.par) = row_names_trembles
    rownames(stratEst.return$gammas.se) = row_names_trembles

    colnames(stratEst.return$trembles.par) = "probability"
    colnames(stratEst.return$trembles.se) = "standard error"
    colnames(stratEst.return$trembles.score) = "score"
    colnames(stratEst.return$trembles.covar) = row_names_trembles
    colnames(stratEst.return$trembles.fisher) = row_names_trembles

    colnames(stratEst.return$gammas.par) = "gamma"
    colnames(stratEst.return$gammas.se) = "standard error"
  }
  if( length( stratEst.return$coefficients.par ) > 0 ){
    row_names_coefficients = paste("par.",as.character(seq(1,length(stratEst.return$coefficients.par),by=1)),sep="")
    rownames(stratEst.return$coefficients.par) = row_names_coefficients
    rownames(stratEst.return$coefficients.se) = row_names_coefficients
    rownames(stratEst.return$coefficients.score) = row_names_coefficients
    rownames(stratEst.return$coefficients.covar) = row_names_coefficients
    rownames(stratEst.return$coefficients.fisher) = row_names_coefficients

    colnames(stratEst.return$coefficients.par) = "coefficient"
    colnames(stratEst.return$coefficients.se) = "standard error"
    colnames(stratEst.return$coefficients.score) = "score"
    colnames(stratEst.return$coefficients.covar) = row_names_coefficients
    colnames(stratEst.return$coefficients.fisher) = row_names_coefficients
  }

  crit <- matrix(cpp.output$fit[1,4:6],1,3)
  colnames(crit) = c("aic","bic","icl")
  rownames(crit) = "crit value"
  stratEst.return$crit = crit

  # delete empty list entries
  if( length(stratEst.return$shares.par) == 0 ){
    stratEst.return$shares.indices = NULL
  }
  if( length(stratEst.return$responses.par) == 0 ){
    stratEst.return$responses.indices = NULL
  }
  else{
    if( all( stratEst.return$responses.par == 0 | stratEst.return$responses.par == 1 ) ){
      stratEst.return$responses.se = NULL
      stratEst.return$responses.score = NULL
      stratEst.return$responses.covar = NULL
      stratEst.return$responses.fisher = NULL
    }
  }
  if( length(stratEst.return$trembles.par) == 0 ){
    stratEst.return$trembles.indices = NULL
  }

  stratEst.return <- stratEst.return[lapply(stratEst.return,length)>0]

  # return result
  return(stratEst.return)
}

