// #define ARMA_NO_DEBUG
#include <RcppArmadillo.h>
using namespace Rcpp;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// stratEst_EM
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

arma::field<arma::mat> stratEst_EM(arma::cube& output_cube, arma::cube& sum_outputs_cube, arma::vec& strat_id, arma::mat shares, arma::vec responses, arma::vec trembles, arma::mat response_mat, arma::mat tremble_mat, arma::mat indices_shares, arma::mat indices_responses, arma::mat indices_trembles, arma::mat& responses_to_sum, std::string& response, int eval_pre , double tol_eval, int max_eval, arma::vec& sample_of_ids, arma::vec& num_ids_sample ) {

  arma::field<arma::mat> F(23,1);
  int num_ids = output_cube.n_slices;
  int num_rows_response_mat = response_mat.n_rows;
  int num_cols_response_mat = response_mat.n_cols;
  arma::uvec shares_to_est = find( indices_shares != 0 );
  arma::vec unique_indices_responses = unique(indices_responses( find( indices_responses != 0 ) ));
  int num_responses_to_est = unique_indices_responses.n_elem;
  arma::vec unique_indices_trembles = unique(indices_trembles( find( indices_trembles != 0 ) ) );
  int num_trembles_to_est = unique_indices_trembles.n_elem;
  int k = shares.n_rows;
  int num_samples = shares.n_cols;
  arma::mat remaining_shares_mat( k , num_samples , arma::fill::zeros );
  arma::vec num_fixed_shares_col( num_samples , arma::fill::zeros );
  for (int i = 0; i < num_samples; i++) {
    arma::vec shares_col = shares( arma::span::all , i );
    arma::vec indices_shares_col = indices_shares( arma::span::all , i );
    arma::uvec fixed_shares_col = find( indices_shares_col == 0 );
    double remaining_shares = 1 - accu( shares_col( fixed_shares_col ) );
    remaining_shares_mat( arma::span::all , i ).fill( remaining_shares );
    num_fixed_shares_col(i) = fixed_shares_col.n_elem;
  }
  int free_params = num_responses_to_est + num_trembles_to_est + shares_to_est.n_elem;
  double eps = 1;
  double eps_now = 1;
  arma::vec ll_val( 1 , 1 , arma::fill::zeros );
  arma::vec entropy_k( 1 , arma::fill::ones);
  arma::vec aic_val( 1 , arma::fill::ones);
  arma::vec bic_val( 1 , arma::fill::ones);
  arma::vec icl_val( 1 , arma::fill::ones);
  arma::vec new_ll_val( 1 , 1 , arma::fill::zeros );
  arma::vec new_entropy_k( 1 , arma::fill::ones);
  arma::mat i_shares_mat( num_ids , k , arma::fill::zeros );
  arma::cube i_shares_cube( num_rows_response_mat , num_cols_response_mat , num_ids , arma::fill::zeros );
  arma::mat state_obs( num_rows_response_mat , 1 , arma::fill::zeros );
  arma::mat new_shares = shares;
  arma::mat new_responses = responses;
  arma::vec new_trembles = trembles;
  arma::umat indices_responses_to_sum = find( responses_to_sum == 1 );
  arma::umat indices_non_fixed_responses = find( indices_responses == 0 && responses_to_sum == 0 );
  int eval = eval_pre;

  // create indices response cube & indidices trembles ctube
  arma::cube indices_responses_cube( num_rows_response_mat , num_cols_response_mat , num_ids , arma::fill::zeros );
  for (int i = 0; i < num_ids; i++) {
    indices_responses_cube.slice(i) = indices_responses;
  }
  arma::cube indices_trembles_tube( num_rows_response_mat , 1 , num_ids , arma::fill::zeros );
  for (int i = 0; i < num_ids; i++) {
    indices_trembles_tube.slice(i) = indices_trembles.col(0);
  }

  // calculate remaining non-fixed response to distribute
  arma::mat incomplete_response_mat = response_mat;
  incomplete_response_mat( indices_non_fixed_responses ).fill(0);
  arma::vec remaining_response_vec = 1 - sum( incomplete_response_mat , 1 );
  arma::mat remaining_response_mat = repmat( remaining_response_vec , 1 , num_cols_response_mat );

  while (  eval < max_eval+eval_pre && eps != tol_eval && eps != arma::datum::nan ) {
    eval++;
    Rcpp::checkUserInterrupt();
    // parameters are assigned to updated parameter values
    shares = new_shares;
    trembles = new_trembles;
    responses = new_responses;
    ll_val = new_ll_val;
    entropy_k = new_entropy_k;

    // calculate emission probabilities with trembles
    arma::mat pr_mat = response_mat % (1 - tremble_mat) + ( 1 - response_mat ) % ( tremble_mat / (tremble_mat.n_cols - 1) );   // probability for emission with trembles

    // new parameters are zero
    new_shares.elem( shares_to_est ).fill(0);
    new_responses.fill(0);
    new_trembles.fill(0);
    new_ll_val(0) = 0;
    new_entropy_k(0) = 0;

    // loop through ids
    for (int i = 0; i < num_ids; i++) {

      // calculate the probability for each outcome in each state of each strategy
      arma::mat pr_outcomes_states_strategies_entities_mat( num_rows_response_mat , num_cols_response_mat , arma::fill::zeros );
      for ( int l = 0; l < num_cols_response_mat; l++){
        for (int j = 0; j < num_rows_response_mat; j++){
          pr_outcomes_states_strategies_entities_mat(j,l) = pow( pr_mat(j,l) , output_cube(j,l,i) );
        }
      }
      arma::vec pr_states_strategies_entities_mat = prod( pr_outcomes_states_strategies_entities_mat , 1 );
      arma::vec pr_entity_k( k , arma::fill::zeros );
      for ( int j = 0; j < k; j++){
        pr_entity_k(j) = prod( pr_states_strategies_entities_mat( find( strat_id == j+1 ) ) );
      }

      // log likelihood contribution of subject
      pr_entity_k %= shares( arma::span::all , sample_of_ids(i) - 1 );
      new_ll_val += -log( sum( pr_entity_k , 0 ) );

      // share contribution of subject as posterior probability of i to use k / N
      arma::vec i_shares = pr_entity_k.each_row() / sum( pr_entity_k , 0 );
      arma::uvec shares_to_est_col = find( indices_shares( arma::span::all , sample_of_ids(i) - 1 ) != 0 );
      arma::vec new_shares_col = new_shares( arma::span::all , sample_of_ids(i) - 1 );
      new_shares_col( shares_to_est_col ) += ( i_shares( shares_to_est_col )  / num_ids_sample( sample_of_ids(i) - 1 ) );
      new_shares( arma::span::all , sample_of_ids(i) - 1 ) = new_shares_col;

      i_shares_mat.row(i) = i_shares.t();
      arma::mat lines_entity_k = repmat( i_shares , 1 , num_cols_response_mat );
      arma::mat entity_slice( num_rows_response_mat , num_cols_response_mat );
      for ( int l = 0; l < num_rows_response_mat; l++){
        entity_slice.row(l) = lines_entity_k.row( strat_id( l ) - 1 );
      }
      i_shares_cube.slice(i) = entity_slice;

      // entropy k contribution
      arma::rowvec i_entropy_k = i_shares.t() % log( i_shares.t() );
      i_entropy_k.replace(arma::datum::nan, 0);
      new_entropy_k -= sum( i_entropy_k , 1 );

    }

    // correct shares for remainder
    for (int j = 0; j < num_samples; j++) {
      if( num_fixed_shares_col(j) > 0 && num_fixed_shares_col(j) != k ){
        arma::uvec shares_to_est_col = find( indices_shares( arma::span::all , j ) != 0 );
        arma::vec new_shares_col = new_shares( arma::span::all , j );
        arma::vec remaining_shares_col = remaining_shares_mat( arma::span::all , j );
        arma::vec shares_of_remaining_shares = new_shares_col( shares_to_est_col );
        new_shares_col( shares_to_est_col ) = remaining_shares_col( shares_to_est_col ) %  ( shares_of_remaining_shares / accu( shares_of_remaining_shares ) );
        new_shares( arma::span::all , j ) = new_shares_col;
      }
    }

    // update responses
    arma::cube weigthed_output_cube = i_shares_cube % output_cube;
    arma::cube weigthed_sums_cube = i_shares_cube % sum_outputs_cube;
    for ( int j = 0; j < num_responses_to_est; j++){
      new_responses(j) = accu( weigthed_output_cube( find( indices_responses_cube == j+1 ) ) ) / accu( weigthed_sums_cube( find( indices_responses_cube == j+1 ) ) );
    }
    new_responses.replace(arma::datum::nan, -1 );                             // clean responses (-1 indicates no obs)

    // calculate state observations
    state_obs = sum( weigthed_sums_cube.tube( 0 , 0 , num_rows_response_mat-1 , 0 ) , 2);

    // fill response mat with normalized new values
    for (int i = 0; i < num_responses_to_est ; i++) {
      arma::umat responses_to_fill = find( indices_responses == i+1 );
      if ( new_responses(i) == -1 ){
        response_mat( responses_to_fill ).fill(-1);
      }
      else{
        arma::mat response_value_mat = response_mat;
        response_value_mat.fill( new_responses(i) );
        response_mat( responses_to_fill ) = response_value_mat( responses_to_fill );
      }
    }
    response_mat( indices_responses_to_sum ).fill(0);
    arma::vec value_to_sum_vec = 1 - sum( response_mat , 1 );
    arma::mat value_to_sum_mat = repmat( value_to_sum_vec , 1 , num_cols_response_mat );
    response_mat( indices_responses_to_sum ) = value_to_sum_mat( indices_responses_to_sum );

    // transform into discrete responses
    if ( response == "pure" && num_responses_to_est > 0 ){
      arma::umat indices_rows_with_estimated = find( sum( indices_responses , 1 )  > 0 );
      arma::mat target_mat = response_mat.rows( indices_rows_with_estimated );
      arma::mat response_mat_with_estimated = response_mat.rows( indices_rows_with_estimated );
      int num_response_mat_with_estimated = response_mat_with_estimated.n_rows;
      for (int i = 0; i < num_response_mat_with_estimated; i++) {
        arma::rowvec target_row = response_mat_with_estimated.row( i );
        arma::uword index_of_max = index_max( target_row );
        target_row.fill(0);
        target_row( index_of_max ) = 1;
        target_mat.row(i) = target_row;
      }
      response_mat.rows( indices_rows_with_estimated ) = target_mat;
      for (int i = 0; i < num_responses_to_est; i++) {
        arma::mat new_discrete_response_value = unique( response_mat( find( indices_responses_cube.slice(0) == i+1 ) ) );
        new_responses(i) = new_discrete_response_value(0,0);
      }
    }

    //update trembles
    arma::cube response_cube( num_rows_response_mat , num_cols_response_mat , num_ids , arma::fill::zeros );
    response_cube.each_slice() = response_mat;
    arma::mat tremble_correction_factor_mat( num_rows_response_mat , num_cols_response_mat , arma::fill::ones );
    arma::cube tremble_correction_factor_cube( num_rows_response_mat , num_cols_response_mat , num_ids, arma::fill::ones );
    tremble_correction_factor_mat( find( response_mat == 0 ) ).fill( num_cols_response_mat-1 );
    tremble_correction_factor_mat( find( response_mat == 1 ) ).fill( -1 );
    tremble_correction_factor_cube.each_slice() = tremble_correction_factor_mat;
    arma::cube weigthed_output_sum_diffs_cube = tremble_correction_factor_cube % i_shares_cube % ( output_cube - ( sum_outputs_cube % response_cube )  );
    arma::cube weigthed_output_sum_diffs_tube = sum( weigthed_output_sum_diffs_cube , 1 );
    arma::cube weigthed_sums_tube = sum( weigthed_sums_cube , 1 );
    for ( int j = 0; j < num_trembles_to_est; j++){
      new_trembles(j) = accu( weigthed_output_sum_diffs_tube( find( indices_trembles_tube == j+1 ) ) ) / accu( weigthed_sums_tube( find( indices_trembles_tube == j+1 ) ) );
    }
    new_trembles.replace(arma::datum::nan, -1 );

    // fill tremble mat with new values
    for (int i = 0; i < num_trembles_to_est ; i++) {
      arma::mat tremble_value_mat = tremble_mat;
      arma::umat trembles_to_fill = find( indices_trembles == i+1 );
      tremble_value_mat.fill( new_trembles(i) );
      tremble_mat( trembles_to_fill ) = tremble_value_mat( trembles_to_fill );
    }

    // check overshooting and calculate eps for tolerance
    if (eval > eval_pre+1 ) { eps_now = (1 - (new_ll_val(0) / ll_val(0))); }          // current epsilon
    if ( new_ll_val(0) == 0 ){ eps_now = 0; }
    if ( new_ll_val.is_finite() ){  eps = eps_now; }                                  // only continue if no overshoot
    else { eps = arma::datum::nan; }                                                  // if overshooting occured report results from last eval

  } // end while

  // calculate selection criteria
  aic_val = new_ll_val + free_params;                                                 // update AIC value
  bic_val = new_ll_val + ( free_params * log( (double) num_ids ) )/2;                           // update BIC value
  icl_val = bic_val + new_entropy_k;                                                  // update ICL value

  // prepare output
  double LL = new_ll_val(0);
  if( LL == arma::datum::nan ){
    LL = ll_val(0);
  }
  double AIC = aic_val(0);
  double BIC = bic_val(0);
  double ICL = icl_val(0);
  double E = entropy_k(0);

  F(0,0) = new_shares;
  F(1,0) = new_responses;
  F(2,0) = new_trembles;
  F(3,0) = response_mat;
  F(4,0) = tremble_mat;
  F(5,0) = LL;
  F(6,0) = AIC;
  F(7,0)  = BIC;
  F(8,0) = ICL;
  F(9,0) = eval;
  F(10,0) = eps;
  F(11,0) = E;
  F(12,0) = i_shares_mat;
  F(20,0) = state_obs;
  F(21,0) = indices_responses;
  F(22,0) = indices_trembles;
  return(F);
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// stratEst_LCR_EM
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

arma::field<arma::mat> stratEst_LCR_EM(arma::cube& output_cube, arma::cube& sum_outputs_cube, arma::vec& strat_id, arma::mat& covariate_mat, arma::mat shares_mat, arma::vec responses, arma::vec trembles, arma::vec coefficients, arma::mat response_mat, arma::mat tremble_mat, arma::mat coefficient_mat , arma::uvec shares_to_est, arma::mat indices_responses, arma::mat indices_trembles, bool estimate_coefficients, arma::uvec coefficients_to_est, arma::mat& responses_to_sum, std::string& response, int eval_pre , double tol_eval, int max_eval, double newton_stepsize, bool penalty) {

  arma::field<arma::mat> F(23,1);
  arma::vec shares = shares_mat.col(0);
  int num_ids = output_cube.n_slices;
  int num_rows_response_mat = response_mat.n_rows;
  int num_cols_response_mat = response_mat.n_cols;
  arma::vec unique_indices_responses = unique(indices_responses( find( indices_responses != 0 ) ));
  int num_responses_to_est = unique_indices_responses.n_elem;
  arma::vec unique_indices_trembles = unique(indices_trembles( find( indices_trembles != 0 ) ) );
  int num_trembles_to_est = unique_indices_trembles.n_elem;
  int num_coefficients_to_est = coefficients.n_elem;
  int num_coefficients = num_coefficients_to_est + coefficient_mat.n_rows;
  int num_rows_coefficient_mat = coefficient_mat.n_rows;
  int k = shares.n_elem;
  int free_params = num_responses_to_est + num_trembles_to_est + shares_to_est.n_elem;
  double eps = arma::datum::inf;
  double eps_now = arma::datum::inf;
  arma::vec ll_val( 1 , 1 , arma::fill::zeros );
  arma::vec entropy_k( 1 , arma::fill::ones);
  arma::vec aic_val( 1 , arma::fill::ones);
  arma::vec bic_val( 1 , arma::fill::ones);
  arma::vec icl_val( 1 , arma::fill::ones);
  arma::vec new_ll_val( 1 , 1 , arma::fill::zeros );
  arma::vec new_entropy_k( 1 , arma::fill::ones);
  arma::mat i_shares_mat( num_ids , k , arma::fill::zeros );
  arma::cube i_shares_cube( num_rows_response_mat , num_cols_response_mat , num_ids , arma::fill::zeros );
  arma::mat state_obs( num_rows_response_mat , 1 , arma::fill::zeros );
  arma::vec new_shares = shares;
  arma::mat new_responses = responses;
  arma::vec new_trembles = trembles;
  arma::vec new_coefficients = coefficients;
  arma::vec changes_coefficients( num_coefficients_to_est , arma::fill::zeros );
  arma::vec score_vec( num_coefficients , arma::fill::zeros );
  arma::vec short_score_vec( num_coefficients_to_est , arma::fill::zeros );
  arma::mat hessian_mat( num_coefficients , num_coefficients , arma::fill::zeros );
  arma::mat fisher_info( num_coefficients , num_coefficients , arma::fill::zeros );
  arma::mat penalized_fisher_info( num_coefficients_to_est , num_coefficients_to_est , arma::fill::zeros );
  arma::cube d_hessian_cube( num_coefficients , num_coefficients , num_coefficients , arma::fill::zeros );
  arma::cube delta_lk( num_coefficients , num_coefficients , num_coefficients , arma::fill::zeros );
  arma::cube delta_lh( num_coefficients , num_coefficients , num_coefficients , arma::fill::zeros );
  arma::cube delta_kh( num_coefficients , num_coefficients , num_coefficients , arma::fill::zeros );
  arma::mat score_contribution_mat( num_ids , num_coefficients , arma::fill::zeros );
  arma::mat penalized_score_contribution_mat( num_ids , num_coefficients_to_est , arma::fill::zeros );
  arma::mat priors_entities_mat( num_ids , k , arma::fill::zeros );
  arma::umat indices_responses_to_sum = find( responses_to_sum == 1 );
  arma::umat indices_non_fixed_responses = find( indices_responses == 0 && responses_to_sum == 0 );
  arma::vec stepsize_vec( num_coefficients_to_est , arma::fill::ones );
  arma::vec penalty_vec( num_coefficients_to_est , arma::fill::zeros );
  stepsize_vec.fill(newton_stepsize);
  int eval = eval_pre;
  bool coefficients_changed = true;

  // create indices response cube & indidices trembles ctube
  arma::cube indices_responses_cube( num_rows_response_mat , num_cols_response_mat , num_ids , arma::fill::zeros );
  for (int i = 0; i < num_ids; i++) {
    indices_responses_cube.slice(i) = indices_responses;
  }
  arma::cube indices_trembles_tube( num_rows_response_mat , 1 , num_ids , arma::fill::zeros );
  for (int i = 0; i < num_ids; i++) {
    indices_trembles_tube.slice(i) = indices_trembles.col(0);
  }

  //create prior entities mat
  priors_entities_mat.col(0).fill(1);
  priors_entities_mat.cols( 1 , k-1 ) = exp( covariate_mat * coefficient_mat );
  for ( int i = 0; i < num_ids; i++){
    priors_entities_mat.row(i) /= accu( priors_entities_mat.row(i) );
  }

  while (  eval < max_eval+eval_pre && eps != tol_eval && eps != arma::datum::nan && coefficients_changed ) {
    eval++;

    // parameters are assigned to updated parameter values
    shares = new_shares;
    trembles = new_trembles;
    responses = new_responses;
    coefficients = new_coefficients;
    ll_val = new_ll_val;
    entropy_k = new_entropy_k;

    // calculate emission probabilities with trembles
    arma::mat pr_mat = response_mat % (1 - tremble_mat) + ( 1 - response_mat ) % ( tremble_mat / (tremble_mat.n_cols - 1) );   // probability for emission with trembles

    // new parameters are zero
    new_shares( shares_to_est ).fill(0);
    new_responses.fill(0);
    new_trembles.fill(0);
    new_ll_val(0) = 0;
    new_entropy_k(0) = 0;
    score_vec.fill(0);
    hessian_mat.fill(0);
    fisher_info.fill(0);
    d_hessian_cube.fill(0);

      // loop through ids
      for (int i = 0; i < num_ids; i++) {

        // calculate the probability for each outcome in each state of each strategy
        arma::mat pr_outcomes_states_strategies_entities_mat( num_rows_response_mat , num_cols_response_mat , arma::fill::zeros );
        for ( int l = 0; l < num_cols_response_mat; l++){
          for (int j = 0; j < num_rows_response_mat; j++){
            pr_outcomes_states_strategies_entities_mat(j,l) = pow( pr_mat(j,l) , output_cube(j,l,i) );
          }
        }
        arma::vec pr_states_strategies_entities_mat = prod( pr_outcomes_states_strategies_entities_mat , 1 );
        arma::vec pr_entity_k( k , arma::fill::zeros );
        for ( int j = 0; j < k; j++){
          pr_entity_k(j) = prod( pr_states_strategies_entities_mat( find( strat_id == j+1 ) ) );
        }

        // log likelihood contribution of subject
        pr_entity_k %= priors_entities_mat.row(i).t();
        new_ll_val += -log( sum( pr_entity_k , 0 ) );

        // share contribution of subject as posterior probability of i to use k / N
        arma::vec i_shares = pr_entity_k.each_row() / sum( pr_entity_k , 0 );
        new_shares( shares_to_est ) += ( i_shares( shares_to_est )  / num_ids );

        i_shares_mat.row(i) = i_shares.t();
        arma::mat lines_entity_k = repmat( i_shares , 1 , num_cols_response_mat );
        arma::mat entity_slice( num_rows_response_mat , num_cols_response_mat );
        for ( int l = 0; l < num_rows_response_mat; l++){
          entity_slice.row(l) = lines_entity_k.row( strat_id( l ) - 1 );
        }
        i_shares_cube.slice(i) = entity_slice;

        // entropy k contribution
        arma::rowvec i_entropy_k = i_shares.t() % log( i_shares.t() );
        i_entropy_k.replace(arma::datum::nan, 0);
        new_entropy_k -= sum( i_entropy_k , 1 );

        // score and hessian contributions
        arma::cube h_covariate_cube( num_coefficients , num_coefficients , num_coefficients , arma::fill::zeros  );
        arma::vec score_contribution( num_coefficients , arma::fill::zeros );
        arma::vec individual_prior = priors_entities_mat.row(i).t();
        arma::vec s_covariate_vec = repmat( covariate_mat.row(i).t() , k , 1 );
        arma::mat h_covariate_mat = repmat( covariate_mat.row(i).t() * covariate_mat.row(i) , k , k );
        arma::mat h_i_shares( num_coefficients , num_coefficients , arma::fill::zeros );
        arma::vec h_i_shares_vec( num_coefficients , arma::fill::zeros );
        arma::mat h_i_prior( num_coefficients , num_coefficients , arma::fill::zeros );
        arma::vec h_i_prior_vec( num_coefficients , arma::fill::zeros );
        arma::vec h_ones_zeros_vec( num_coefficients , arma::fill::zeros );
        arma::mat h_ones_zeros( num_coefficients , num_coefficients , arma::fill::zeros );
        h_covariate_cube.each_slice() = h_covariate_mat;
        arma::cube third_dim_cube( num_coefficients , num_coefficients , num_coefficients , arma::fill::zeros  );
        arma::mat third_dim_mat( num_coefficients , num_coefficients , arma::fill::zeros  );
        for (int s = 0; s < num_coefficients; s++){
          third_dim_mat.fill(s_covariate_vec(s));
          third_dim_cube.slice(s) = third_dim_mat;
        }
        h_covariate_cube %= third_dim_cube;
        for (int m = 0; m < k; m++){
          for (int j = 0; j < k; j++){
            for ( int l = 0; l < num_rows_coefficient_mat ; l++){
              if( m == 0 ){
                h_i_shares_vec( l + j*num_rows_coefficient_mat ) = i_shares(j);
                h_i_prior_vec( l + j*num_rows_coefficient_mat ) = individual_prior(j);
              }
              if( j == m ){
                h_ones_zeros_vec( l + j*num_rows_coefficient_mat ) = 1;
              }
              else{
                h_ones_zeros_vec( l + j*num_rows_coefficient_mat ) = 0;
              }
            }
            if( m == 0 ){
              score_contribution = s_covariate_vec % ( h_i_shares_vec - h_i_prior_vec );
            }
          }
          h_i_shares( arma::span::all , arma::span(m*num_rows_coefficient_mat , (num_rows_coefficient_mat-1) + m*num_rows_coefficient_mat ) ) = repmat( h_i_shares_vec , 1 , num_rows_coefficient_mat );
          h_i_prior( arma::span::all , arma::span(m*num_rows_coefficient_mat , (num_rows_coefficient_mat-1) + m*num_rows_coefficient_mat ) ) = repmat( h_i_prior_vec , 1 , num_rows_coefficient_mat );
          h_ones_zeros( arma::span::all , arma::span(m*num_rows_coefficient_mat , (num_rows_coefficient_mat-1) + m*num_rows_coefficient_mat ) ) = repmat( h_ones_zeros_vec , 1 , num_rows_coefficient_mat );
        }

        // objects for derivative of hessian
        arma::cube l_i_shares_cube( num_coefficients , num_coefficients , num_coefficients , arma::fill::zeros  );
        arma::cube k_i_shares_cube( num_coefficients , num_coefficients , num_coefficients , arma::fill::zeros  );
        arma::cube h_i_shares_cube( num_coefficients , num_coefficients , num_coefficients , arma::fill::zeros  );
        arma::cube l_i_prior_cube( num_coefficients , num_coefficients , num_coefficients , arma::fill::zeros  );
        arma::cube k_i_prior_cube( num_coefficients , num_coefficients , num_coefficients , arma::fill::zeros  );
        arma::cube h_i_prior_cube( num_coefficients , num_coefficients , num_coefficients , arma::fill::zeros  );

        arma::mat l_i_shares_mat( num_coefficients , num_coefficients , arma::fill::zeros  );
        arma::mat k_i_shares_mat( num_coefficients , num_coefficients , arma::fill::zeros  );
        arma::mat h_i_shares_mat( num_coefficients , num_coefficients , arma::fill::zeros  );
        arma::mat l_i_prior_mat( num_coefficients , num_coefficients , arma::fill::zeros  );
        arma::mat k_i_prior_mat( num_coefficients , num_coefficients , arma::fill::zeros  );
        arma::mat h_i_prior_mat( num_coefficients , num_coefficients , arma::fill::zeros  );

        arma::rowvec h_i_shares_row = h_i_shares.col(0).t();
        arma::rowvec h_i_prior_row = h_i_prior.col(0).t();

        l_i_shares_mat.each_col() = h_i_shares_row.t();
        k_i_shares_mat.each_row() = h_i_shares_row;
        l_i_prior_mat.each_col() = h_i_prior_row.t();
        k_i_prior_mat.each_row() = h_i_prior_row;
        l_i_shares_cube.each_slice() = l_i_shares_mat;
        k_i_shares_cube.each_slice() = k_i_shares_mat;
        l_i_prior_cube.each_slice() = l_i_prior_mat;
        k_i_prior_cube.each_slice() = k_i_prior_mat;

        if( i == 0 ){
          delta_lk.each_slice() = h_ones_zeros;
          for (int l = 0; l < num_coefficients; l++){
            for (int k = 0; k < num_coefficients; k++){
              for (int h = 0; h < num_coefficients; h++){
                double elem = delta_lk( l , k , h );
                delta_lh( l , h , k ) = elem;
                delta_kh( h , l ,  k ) = elem;
              }
            }
          }
        }

        for (int s = 0; s < num_coefficients; s++){
          h_i_shares_mat.fill(h_i_shares_row(s));
          h_i_shares_cube.slice(s) = h_i_shares_mat;
          h_i_prior_mat.fill(h_i_prior_row(s));
          h_i_prior_cube.slice(s) = h_i_prior_mat;
        }

        // individual's hessian, score and fisher cotributions
        arma::mat hessian_contribution = h_covariate_mat % ( h_i_shares % ( h_ones_zeros - h_i_shares.t() )  - h_i_prior % ( h_ones_zeros - h_i_prior.t() ) );
        hessian_mat += hessian_contribution;
        score_vec += score_contribution;
        score_contribution_mat.row(i) = score_contribution.t();
        fisher_info += score_contribution * score_contribution.t();

        // individual's derivative of hessian contribution
        if( penalty ){
          arma::cube d_hessian_i_shares = l_i_shares_cube % ( delta_lk - k_i_shares_cube ) % ( delta_lh - h_i_shares_cube) - l_i_shares_cube % k_i_shares_cube % ( delta_kh - h_i_shares_cube );
          arma::cube d_hessian_i_prior = l_i_prior_cube % ( delta_lk - k_i_prior_cube ) % ( delta_lh - h_i_prior_cube) - l_i_prior_cube % k_i_prior_cube % ( delta_kh - h_i_prior_cube );
          arma::cube d_hessian_contribution = h_covariate_cube % ( d_hessian_i_shares + d_hessian_i_prior );
          d_hessian_cube += d_hessian_contribution;

        }

      }

    //update responses
    arma::cube weigthed_output_cube = i_shares_cube % output_cube;
    arma::cube weigthed_sums_cube = i_shares_cube % sum_outputs_cube;
    for ( int j = 0; j < num_responses_to_est; j++){
      new_responses(j) = accu( weigthed_output_cube( find( indices_responses_cube == j+1 ) ) ) / accu( weigthed_sums_cube( find( indices_responses_cube == j+1 ) ) );
    }
    new_responses.replace(arma::datum::nan, -1 );                             // clean responses (-1 indicates no obs)

    // calculate state observations
    state_obs = sum( weigthed_sums_cube.tube( 0 , 0 , num_rows_response_mat-1 , 0 ) , 2);

    // fill response mat with normalized new values
    for (int i = 0; i < num_responses_to_est ; i++) {
      arma::umat responses_to_fill = find( indices_responses == i+1 );
      if ( new_responses(i) == -1 ){
        response_mat( responses_to_fill ).fill(-1);
      }
      else{
        arma::mat response_value_mat = response_mat;
        response_value_mat.fill( new_responses(i) );
        response_mat( responses_to_fill ) = response_value_mat( responses_to_fill );
      }
    }
    response_mat( indices_responses_to_sum ).fill(0);
    arma::vec value_to_sum_vec = 1 - sum( response_mat , 1 );
    arma::mat value_to_sum_mat = repmat( value_to_sum_vec , 1 , num_cols_response_mat );
    response_mat( indices_responses_to_sum ) = value_to_sum_mat( indices_responses_to_sum );

    // transform into discrete responses
    if ( response == "pure" && num_responses_to_est > 0 ){
      arma::umat indices_rows_with_estimated = find( sum( indices_responses , 1 )  > 0 );
      arma::mat target_mat = response_mat.rows( indices_rows_with_estimated );
      arma::mat response_mat_with_estimated = response_mat.rows( indices_rows_with_estimated );
      int num_response_mat_with_estimated = response_mat_with_estimated.n_rows;
      for (int i = 0; i < num_response_mat_with_estimated; i++) {
        arma::rowvec target_row = response_mat_with_estimated.row( i );
        arma::uword index_of_max = index_max( target_row );
        target_row.fill(0);
        target_row( index_of_max ) = 1;
        target_mat.row(i) = target_row;
      }
      response_mat.rows( indices_rows_with_estimated ) = target_mat;
      for (int i = 0; i < num_responses_to_est; i++) {
        arma::mat new_discrete_response_value = unique( response_mat( find( indices_responses == i+1 ) ) );
        new_responses(i) = new_discrete_response_value(0,0);
      }
    }

    //update trembles
    arma::cube response_cube( num_rows_response_mat , num_cols_response_mat , num_ids , arma::fill::zeros );
    response_cube.each_slice() = response_mat;
    arma::mat tremble_correction_factor_mat( num_rows_response_mat , num_cols_response_mat , arma::fill::ones );
    arma::cube tremble_correction_factor_cube( num_rows_response_mat , num_cols_response_mat , num_ids, arma::fill::ones );
    tremble_correction_factor_mat( find( response_mat == 0 ) ).fill( num_cols_response_mat-1 );
    tremble_correction_factor_mat( find( response_mat == 1 ) ).fill( -1 );
    tremble_correction_factor_cube.each_slice() = tremble_correction_factor_mat;
    arma::cube weigthed_output_sum_diffs_cube = tremble_correction_factor_cube % i_shares_cube % ( output_cube - ( sum_outputs_cube % response_cube )  );
    arma::cube weigthed_output_sum_diffs_tube = sum( weigthed_output_sum_diffs_cube , 1 );
    arma::cube weigthed_sums_tube = sum( weigthed_sums_cube , 1 );
    for ( int j = 0; j < num_trembles_to_est; j++){
      new_trembles(j) = accu( weigthed_output_sum_diffs_tube( find( indices_trembles_tube == j+1 ) ) ) / accu( weigthed_sums_tube( find( indices_trembles_tube == j+1 ) ) );
    }
    new_trembles.replace(arma::datum::nan, -1 );

    // fill tremble mat with new values
    for (int i = 0; i < num_trembles_to_est ; i++) {
      arma::mat tremble_value_mat  = tremble_mat;
      arma::umat trembles_to_fill = find( indices_trembles == i+1 );
      tremble_value_mat.fill( new_trembles(i) );
      tremble_mat( trembles_to_fill ) = tremble_value_mat( trembles_to_fill );
    }

    // update coefficients
    if( estimate_coefficients ){
          short_score_vec = score_vec( arma::span( num_rows_coefficient_mat , num_coefficients-1 ) );
          arma::mat lower_hessian_mat = hessian_mat( arma::span( num_rows_coefficient_mat , num_coefficients-1 ) , arma::span( num_rows_coefficient_mat , num_coefficients-1 ) );
          arma::mat lower_fisher = fisher_info( arma::span( num_rows_coefficient_mat , num_coefficients-1 ) , arma::span( num_rows_coefficient_mat , num_coefficients-1 ) );
          arma::mat inverted_mat( num_coefficients_to_est , num_coefficients_to_est , arma::fill::none );
          arma::mat to_invert_mat = -lower_hessian_mat;                       // observed fisher is negative hesiian since we are minimizing the log likelihood
      if( pinv( inverted_mat , to_invert_mat ) ){
              if( penalty == true ){
                new_ll_val = new_ll_val - log(det(-lower_hessian_mat))/2;     // minus because lll_val is negative log likelihood
                for(int c = 0; c < num_coefficients_to_est; c++){
                  arma::mat d_hessian_slice = d_hessian_cube.slice(num_rows_coefficient_mat+c);
                  arma::mat lower_d_hessian_mat = d_hessian_slice( arma::span( num_rows_coefficient_mat , num_coefficients-1 ) , arma::span( num_rows_coefficient_mat , num_coefficients-1 ) );
                  penalty_vec(c) = trace(inverted_mat*lower_d_hessian_mat)/2;
                }
                short_score_vec += penalty_vec;                               // minus because the negative log likelihood is minimized
              }
        changes_coefficients = stepsize_vec % ( inverted_mat*short_score_vec );
        arma::vec updated_coefficients = coefficients + changes_coefficients;
        new_coefficients( coefficients_to_est ) = updated_coefficients( coefficients_to_est );
        //Rcout<< "ll val: \n"  << new_ll_val << "\n";
        //Rcout<< "score vec: \n"  << short_score_vec << "\n";
      }
      else{
        eval = max_eval+eval_pre;
        new_ll_val.fill(arma::datum::nan);
      }
      //create prior entities mat & fisher info
      coefficient_mat = reshape( new_coefficients , num_rows_coefficient_mat , k-1 );
      priors_entities_mat.col(0).fill(1);
      priors_entities_mat.cols( 1 , k-1 ) = exp( covariate_mat * coefficient_mat );
      for ( int i = 0; i < num_ids; i++){
        priors_entities_mat.row(i) /= accu( priors_entities_mat.row(i) );
      }
    }
    else{
      new_coefficients = coefficients;
    }


    // update shares
    new_shares = mean( priors_entities_mat.t() , 1 );

    // check overshooting and calculate eps for tolerance
    if( estimate_coefficients ){
      if( eval > eval_pre+1 ){
        if( max(abs(changes_coefficients)/abs(coefficients)) <= tol_eval ){
          coefficients_changed = false;
        }
      }
      if (eval > eval_pre+1 ) { eps_now = (1 - (new_ll_val(0) / ll_val(0))); }          // current epsilon
      if ( eps_now < tol_eval ){ eps_now = tol_eval; }
      if ( new_ll_val(0) == 0 ){ eps_now = 0; }
      if ( new_ll_val.is_finite() ){  eps = eps_now; }                                  // only continue if no overshoot
      else { eps = arma::datum::nan; }
    }
    else{
      if (eval > eval_pre+1 ) { eps_now = (1 - (new_ll_val(0) / ll_val(0))); }          // current epsilon
      if ( eps_now < tol_eval ){ eps_now = tol_eval; }
      if ( new_ll_val(0) == 0 ){ eps_now = 0; }
      if ( new_ll_val.is_finite() ){  eps = eps_now; }                                  // only continue if no overshoot
      else { eps = arma::datum::nan; }
    }                                                                                   // if overshooting occured report results from last eval

  } // end while

  // calculate selection criteria
  aic_val = new_ll_val + free_params;                                                 // update AIC value
  bic_val = new_ll_val + ( free_params * log( (double) num_ids ) )/2;                           // update BIC value
  icl_val = bic_val + new_entropy_k;                                                  // update ICL value

  // prepare output
  double LL = new_ll_val(0);
  if( LL == arma::datum::nan ){
    LL = ll_val(0);
  }
  double AIC = aic_val(0);
  double BIC = bic_val(0);
  double ICL = icl_val(0);
  double E = entropy_k(0);

  F(0,0) = new_shares;
  F(1,0) = new_responses;
  F(2,0) = new_trembles;
  F(3,0) = response_mat;
  F(4,0) = tremble_mat;
  F(5,0) = LL;
  F(6,0) = AIC;
  F(7,0)  = BIC;
  F(8,0) = ICL;
  F(9,0) = eval;
  F(10,0) = eps;
  F(11,0) = E;
  F(12,0) = i_shares_mat;
  F(13,0) = new_coefficients;
  F(14,0) = coefficient_mat;
  F(15,0) = priors_entities_mat;
  F(16,0) = short_score_vec;
  F(17,0) = hessian_mat;
  F(18,0) = fisher_info;
  F(19,0) = score_contribution_mat;
  F(20,0) = state_obs;
  F(21,0) = indices_responses;
  F(22,0) = indices_trembles;
  return(F);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// stratEst_SE
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

arma::field<arma::mat> stratEst_SE(arma::cube& output_cube, arma::cube& sum_outputs_cube, arma::vec& strat_id, arma::mat& covariate_mat, arma::mat shares, arma::vec responses, arma::vec trembles, arma::vec coefficients, arma::mat response_mat, arma::mat tremble_mat , arma::mat coefficient_mat , arma::mat i_shares_mat, arma::mat fisher_info_coefficients, arma::vec score_coefficients, arma::mat hessian_mat, arma::mat individual_priors_mat, arma::mat indices_shares, arma::mat indices_responses, arma::mat indices_trembles, std::string response, arma::mat score_contribution_mat, bool CL, arma::vec cluster_id_vec, bool LCR, arma::vec& sample_of_ids, arma::vec& num_ids_sample) {

  arma::field<arma::mat> F(17,1);
  int num_ids = output_cube.n_slices;
  int num_responses_to_est = responses.n_elem;
  int num_trembles_to_est = trembles.n_elem;
  int num_cols_response_mat = response_mat.n_cols;
  int num_rows_response_mat = response_mat.n_rows;
  int num_rows_coefficient_mat = coefficient_mat.n_rows;
  int num_coefficients_to_est = coefficients.n_elem;
  int num_coefficients = num_coefficients_to_est + coefficient_mat.n_rows;
  int k = shares.n_rows;
  int num_samples = shares.n_cols;
  arma::vec shares_vec = vectorise( shares );
  arma::uvec shares_to_est = find( indices_shares > 0 );
  arma::mat convergence( 1 , 4 , arma::fill::ones );
  convergence.fill(-1);
  arma::mat score_responses_slice( num_rows_response_mat , num_cols_response_mat , arma::fill::zeros );
  arma::mat score_trembles_slice( num_rows_response_mat , num_cols_response_mat , arma::fill::zeros );

  arma::mat pr_mat = response_mat % (1 - tremble_mat) + ( 1 - response_mat ) % ( tremble_mat / (tremble_mat.n_cols - 1) );

  // initialize objects
  arma::mat SE_shares( k , num_samples , arma::fill::zeros );
  arma::mat shares_covar( k*num_samples , k*num_samples , arma::fill::zeros );
  arma::vec SE_responses( num_responses_to_est , arma::fill::zeros );
  arma::mat responses_covar( num_responses_to_est , num_responses_to_est , arma::fill::zeros );
  arma::vec SE_trembles( num_trembles_to_est , arma::fill::zeros );
  arma::mat trembles_covar( num_trembles_to_est  , num_trembles_to_est  , arma::fill::zeros );
  arma::vec SE_coefficients( num_coefficients_to_est , arma::fill::zeros );
  arma::mat coefficients_covar( num_coefficients_to_est  , num_coefficients_to_est  , arma::fill::zeros );

  arma::vec score_shares( k*num_samples , arma::fill::zeros );
  arma::vec score_responses( num_responses_to_est , arma::fill::zeros );
  arma::vec score_trembles( num_trembles_to_est , arma::fill::zeros );

  arma::mat fisher_info_shares( k*num_samples , k*num_samples , arma::fill::zeros );
  arma::mat fisher_info_responses( num_responses_to_est , num_responses_to_est , arma::fill::zeros );
  arma::mat fisher_info_trembles( num_trembles_to_est , num_trembles_to_est , arma::fill::zeros );

  // SE of mixed responses & trembles
  arma::mat r_weight_mat( num_rows_response_mat , num_cols_response_mat , arma::fill::zeros );
  arma::mat one_weight_mat( num_rows_response_mat , num_cols_response_mat , arma::fill::ones );
  r_weight_mat.fill( num_cols_response_mat - 1 );
  arma::mat tremble_weight_mat = ( one_weight_mat - response_mat)/r_weight_mat - response_mat;

  for (int i = 0; i < num_ids ; i++) {
    arma::vec score_contribution_responses( num_responses_to_est , arma::fill::zeros );
    arma::vec score_contribution_trembles( num_trembles_to_est , arma::fill::zeros );
    arma::mat lines_entity_k = repmat( i_shares_mat.row(i).t() , 1 , num_cols_response_mat );
    arma::mat entity_slice( num_rows_response_mat , num_cols_response_mat );
    for ( int l = 0; l < num_rows_response_mat; l++){
      entity_slice.row(l) = lines_entity_k.row( strat_id( l ) - 1 );
    }

    // fisher information of responses
    if( response == "mixed" && num_responses_to_est > 0 ){
      score_responses_slice = entity_slice % ( output_cube.slice(i) - ( response_mat % sum_outputs_cube.slice(i) ) );
      for (int j = 0; j < num_responses_to_est ; j++) {
        score_contribution_responses(j) =   accu( score_responses_slice( find( indices_responses == j+1 ) ) );
      }
      score_responses += score_contribution_responses;
      fisher_info_responses += score_contribution_responses * score_contribution_responses.t();
    }

    // fisher information for trembles
    if( num_trembles_to_est > 0 ){
      arma::mat indices_free_trembles = indices_trembles;
      indices_free_trembles.col(0).fill(0);
      score_trembles_slice =  entity_slice % tremble_weight_mat % ( output_cube.slice(i)  /  pr_mat  ) ;
      for (int j = 0; j < num_trembles_to_est ; j++) {
        score_contribution_trembles(j) +=  accu( score_trembles_slice( find( indices_trembles == j+1 ) ) );
      }
      score_trembles += score_contribution_trembles;
      fisher_info_trembles += score_contribution_trembles * score_contribution_trembles.t();
    }
  }

  // SEs of responses
  if( response == "mixed" && num_responses_to_est > 0 ){
    arma::mat inverse_fisher_info_responses = fisher_info_responses;
    if( pinv( inverse_fisher_info_responses , fisher_info_responses ) ){
      arma::mat identity_responses( num_responses_to_est , num_responses_to_est , arma::fill::eye );
      arma::mat rows_mat_responses = repmat( responses , 1 , num_responses_to_est );
      arma::mat cols_mat_responses = repmat( responses.t() , num_responses_to_est , 1 );
      arma::mat no_sum_constraint_responses = cols_mat_responses;
      for (int i = 0; i < num_responses_to_est ; i++) {
        for (int j = 0; j < num_responses_to_est ; j++) {
          arma::rowvec unique_indices_responses_row( num_cols_response_mat , arma::fill::zeros );
          for (int l = 0; l < num_rows_response_mat ; l++) {
            arma::rowvec indices_responses_row = indices_responses.row(l);
            if( any( indices_responses_row == i+1 ) ){
              indices_responses_row( find( indices_responses_row == i+1 ) );
              unique_indices_responses_row = indices_responses_row;
              l = num_rows_response_mat;
            }
          }
          if( any( unique_indices_responses_row == j+1 ) ){
            no_sum_constraint_responses(i,j) = 0;
          }
        }
      }
      identity_responses += no_sum_constraint_responses;
      arma::mat jacobian_responses = rows_mat_responses % ( identity_responses - cols_mat_responses );
      arma::mat covar_responses_mat = jacobian_responses * inverse_fisher_info_responses * jacobian_responses.t() ;
      responses_covar = covar_responses_mat;
      SE_responses = sqrt( covar_responses_mat.diag() );
    }
    else{
      SE_responses.fill(arma::datum::nan);
      responses_covar.fill(arma::datum::nan);
    }

  }

  // SEs of trembles
  if( num_trembles_to_est > 0 ){
    arma::mat inverse_fisher_info_trembles = fisher_info_trembles;
    if( pinv( inverse_fisher_info_trembles , fisher_info_trembles ) ){
      SE_trembles = sqrt( inverse_fisher_info_trembles.diag() );
    }
    else{
      SE_trembles.fill(arma::datum::nan);
    }
  }

  if( LCR ){
    // SEs of coefficients
    arma::mat lower_hessian_mat = hessian_mat( arma::span( num_rows_coefficient_mat , num_coefficients-1 ) , arma::span( num_rows_coefficient_mat , num_coefficients-1 ) );
    arma::mat inverse_neg_lower_hessian_mat = -lower_hessian_mat;
    if( CL ){
      arma::vec unique_clusters = unique( cluster_id_vec );
      int num_clusters = unique_clusters.n_elem;
      arma::cube score_contributions_clusters( num_coefficients , num_coefficients , num_clusters , arma::fill::zeros );
      for (int c = 0; c < num_clusters ; c++) {
        arma::mat score_contribution_cluster = sum( score_contribution_mat.rows( find( cluster_id_vec == unique_clusters(c) ) ) , 0 );
        score_contributions_clusters.slice(c) = score_contribution_cluster.t() * score_contribution_cluster;
      }
      arma::mat meat_mat = sum( score_contributions_clusters , 2 );
      arma::mat lower_meat_mat = meat_mat( arma::span( num_rows_coefficient_mat , num_coefficients-1 ) , arma::span( num_rows_coefficient_mat , num_coefficients-1 ) );
      if( pinv( inverse_neg_lower_hessian_mat , -lower_hessian_mat ) ){
        coefficients_covar = inverse_neg_lower_hessian_mat * lower_meat_mat * inverse_neg_lower_hessian_mat;
        SE_coefficients = sqrt( coefficients_covar.diag() );
      }
      else{
        coefficients_covar.fill(arma::datum::nan);
        SE_coefficients.fill(arma::datum::nan);
      }
    }
    else{
      arma::mat inverse_neg_lower_hessian_mat = lower_hessian_mat;
      if( pinv( inverse_neg_lower_hessian_mat , -lower_hessian_mat ) ){
        arma::mat covar_lower_mat_coefficients = inverse_neg_lower_hessian_mat;
        coefficients_covar = covar_lower_mat_coefficients;
        SE_coefficients = sqrt( covar_lower_mat_coefficients.diag() );
      }
      else{
        coefficients_covar.fill(arma::datum::nan);
        SE_coefficients.fill(arma::datum::nan);
      }
    }

    // SE of shares from covar mat coefficients via delta method
    arma::mat covar_mat_coefficients( num_coefficients , num_coefficients, arma::fill::zeros );
    arma::mat covar_lower_mat_coefficients = covar_mat_coefficients( arma::span( num_rows_coefficient_mat , num_coefficients-1 ) , arma::span( num_rows_coefficient_mat , num_coefficients-1 ) );
    arma::mat inverse_neg_hessian_mat = -hessian_mat;
    if( pinv( inverse_neg_hessian_mat , -hessian_mat ) ){
      covar_mat_coefficients = inverse_neg_hessian_mat;
    }
    else{
      covar_mat_coefficients.fill(arma::datum::nan);
    }
    int length_covariates = coefficient_mat.n_rows;
    arma::mat jacobian_mat_priors( k , num_coefficients , arma::fill::zeros );
    arma::mat zero_mat( k , num_coefficients , arma::fill::zeros );
    arma::rowvec ones_vec( length_covariates , arma::fill::ones );
    for (int j = 0; j < k ; j++) {
      zero_mat( j , arma::span( j*length_covariates , ( j*length_covariates + length_covariates - 1) ) ) = ones_vec;
    }
    for (int i = 0; i < num_ids ; i++) {
      arma::rowvec prior_i = individual_priors_mat.row(i);
      arma::mat prior_i_mat = repmat( prior_i.t() , 1 , num_coefficients );
      arma::mat prior_i_column_mat = repmat( prior_i.t() , 1 , length_covariates );
      arma::mat prior_i_column_vec( 1 , num_coefficients, arma::fill::zeros );
      arma::mat covariate_mat_i = repmat( covariate_mat.row(i) , k , k );
      for (int j = 0; j < k ; j++) {
        prior_i_column_vec( 0 , arma::span( (j*length_covariates) , (j*length_covariates + length_covariates - 1) ) ) = prior_i_column_mat.row(j);
      }
      arma::mat prior_i_column = repmat( prior_i_column_vec , k , 1 );
      jacobian_mat_priors += ( covariate_mat_i % (prior_i_mat % (zero_mat - prior_i_column) ) )/num_ids;
    }
    arma::mat covar_mat_shares = jacobian_mat_priors * covar_mat_coefficients * jacobian_mat_priors.t();
    SE_shares = sqrt( covar_mat_shares.diag() );
    shares_covar = covar_mat_shares;
  }
  else{  // LCR == false
    if( k > 1 ){
      for (int i = 0; i < num_ids ; i++) {
        arma::rowvec individual_shares = shares.col(sample_of_ids(i)-1).t();
        arma::rowvec score_contributions_shares = i_shares_mat.row(i) - individual_shares;
        int start_index_id_fisher = k*(sample_of_ids(i)-1);
        score_shares(arma::span(start_index_id_fisher,(start_index_id_fisher+k-1))) += score_contributions_shares.t();
        fisher_info_shares(arma::span(start_index_id_fisher,(start_index_id_fisher+k-1)),arma::span(start_index_id_fisher,(start_index_id_fisher+k-1))) += score_contributions_shares.t() * score_contributions_shares;
      }
      arma::mat inverse_fisher_info_shares = fisher_info_shares;
      if( pinv( inverse_fisher_info_shares , fisher_info_shares ) ){
        arma::mat identity_shares( k*num_samples , k*num_samples , arma::fill::eye );
        arma::mat rows_mat_shares = repmat( shares_vec , 1 , k*num_samples );
        arma::mat cols_mat_shares = repmat( shares_vec.t() , k*num_samples , 1 );
        arma::mat jacobian_shares = rows_mat_shares % ( identity_shares - cols_mat_shares ) ;
        arma::mat covar_shares_mat = jacobian_shares * inverse_fisher_info_shares * jacobian_shares.t() ;
        SE_shares = sqrt( covar_shares_mat.diag() );
        shares_covar = covar_shares_mat;
      }
      else{
        SE_shares.fill(arma::datum::nan);
        shares_covar.fill(arma::datum::nan);
      }
    }
  }

  // convergence checks based on first order condition
  if( shares_to_est.n_elem > 0 ){
    convergence(0,0) = max(abs( score_shares ));
  }
  if( num_responses_to_est > 0 ){
    convergence(0,1) = max(abs( score_responses ));
  }
  if( num_trembles_to_est > 0 ){
    convergence(0,2) = max(abs( score_trembles ));
  }
  if( LCR ){
    if( num_coefficients_to_est > 0 ){
      convergence(0,3) = max(abs( score_coefficients ));
    }
  }

  F(0,0) = reshape( SE_shares, k , num_samples );
  F(1,0) = SE_responses;
  F(2,0) = SE_trembles;
  F(3,0) = SE_coefficients;
  F(4,0) = convergence;
  F(5,0) = shares_covar;
  F(6,0) = score_shares;
  F(7,0) = fisher_info_shares;
  F(8,0) = responses_covar;
  F(9,0) = score_responses;
  F(10,0) = fisher_info_responses;
  F(11,0) = trembles_covar;
  F(12,0) = score_trembles;
  F(13,0) = fisher_info_trembles;
  F(14,0) = coefficients_covar;
  F(15,0) = score_coefficients;
  F(16,0) = fisher_info_coefficients;

  return(F);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Main Function
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
List stratEst_cpp(arma::mat data, arma::mat strategies, arma::vec sid, arma::mat shares , arma::vec trembles , arma::mat coefficients, arma::mat covariates, bool LCR, arma::vec cluster, std::string response = "mixed", std::string r_responses = "no", std::string r_trembles = "global", std::string select = "no", int min_strategies = 1 , std::string crit = "bic", std::string SE = "yes", int outer_runs = 10, double outer_tol_eval = 0, int outer_max_eval = 1000, int inner_runs = 100, double inner_tol_eval = 0, int inner_max_eval = 100, int LCR_runs = 100, double LCR_tol_eval = 0, int LCR_max_eval = 1000, int BS_samples = 1000, bool print_messages = true, bool integer_strategies = true, double newton_stepsize = 1 , bool penalty = false ) {

  arma::field<arma::mat> R(23,1);
  int rows_data = data.n_rows;
  int cols_data = data.n_cols;
  int rows_strategies = strategies.n_rows;
  int cols_strategies = strategies.n_cols;
  arma::vec state = strategies.col(0);
  int k = accu( state( find( state == 1 ) ) );
  int max_state = max( state );
  arma::vec n_ones( rows_data , 1 , arma::fill::ones );
  int LL_index = 5;
  int crit_index = 7;
  if( crit == "aic" ){
    crit_index = 6;
  }
  else if ( crit == "icl" ){
    crit_index = 8;
  }
  arma::mat complete_share_mat = shares;
  arma::mat failed_inner_runs( 1, 1 , arma::fill::zeros );
  arma::mat failed_outer_runs( 1, 1 , arma::fill::zeros );
  arma::mat failed_inner_LCR_runs( 1, 1 , arma::fill::zeros );
  arma::mat failed_outer_LCR_runs( 1, 1 , arma::fill::zeros );

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // check and transform function inputs
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // generate ids, supergames & rounds 1,2,...,max
  arma::vec id_vec = data.col(0);
  arma::vec match_vec = data.col(1);
  arma::vec round_vec = data.col(2);
  arma::uvec zeros_id = find( id_vec <= 0 );
  int num_zeros_id = zeros_id.n_elem;
  arma::uvec zeros_match = find( match_vec <= 0 );
  int num_zeros_match = zeros_match.n_elem;
  arma::uvec zeros_round = find( round_vec <= 0 );
  int num_zeros_round = zeros_round.n_elem;
  if ( num_zeros_id != 0 || num_zeros_match != 0 || num_zeros_round != 0 ){
    stop("The variables id, game (supergame), and period in of the data frame must contain values greater than zero.");
  }
  arma::vec unique_ids = unique( data.col(0) ) ;
  int num_ids = unique_ids.n_elem;
  for(int i = 0; i < num_ids; i++) {
    data.col(0).replace( unique_ids(i), i+1 );
  }
  id_vec = data.col(0);
  for(int i = 1; i <= num_ids; i++) {
    arma::vec matches_sbj = match_vec( find( data.col(0) == i ) );
    arma::vec unique_matches_sbj = unique( matches_sbj );
    int num_matches_sbj = unique_matches_sbj.n_elem;
    for(int j = 0; j < num_matches_sbj; j++) {
      matches_sbj.replace( unique_matches_sbj(j), j+1 );
      arma::vec rounds_sbj = round_vec( find( id_vec == i && match_vec == unique_matches_sbj(j) ) );
      arma::vec unique_rounds_sbj = unique( rounds_sbj );
      int num_unique_rounds_sbj = unique_rounds_sbj.n_elem;
      int num_rounds_sbj = rounds_sbj.n_elem;
      if( num_unique_rounds_sbj != num_rounds_sbj ){
        stop("The same round number cannot occur several times within the same game for the same subject.");
      }
      int num_rounds_sbj_match = unique_rounds_sbj.n_elem;
      for(int k = 0; k < num_rounds_sbj_match; k++) {
        rounds_sbj.replace( unique_rounds_sbj(k), k+1 );
      }
      round_vec( find( id_vec == i && match_vec == unique_matches_sbj(j)  ) ) = rounds_sbj;
    }
    match_vec( find( data.col(0) == i ) ) = matches_sbj;
  }
  data.col(1) = match_vec;
  data.col(2) = round_vec;

  // generate inputs 1,2,...,N
  arma::vec unique_inputs = unique( data.col(3) );
  arma::vec non_zero_unique_inputs = unique_inputs( find( unique_inputs != 0 ) );
  int num_non_zero_inputs = non_zero_unique_inputs.n_elem;
  for(int i = 0; i < num_non_zero_inputs; i++) {
    data.col(3).replace( non_zero_unique_inputs(i), i+1 );
  }

  // generate outputs 1,2,...,N
  arma::vec unique_outputs = unique( data.col(4) ) ;
  arma::vec non_zero_unique_outputs = unique_outputs( find( unique_outputs != 0 ) );
  int num_non_zero_outputs = non_zero_unique_outputs.n_elem;
  for(int i = 0; i < num_non_zero_outputs; i++) {
    data.col(4).replace( non_zero_unique_outputs(i) , i+1 );
  }

  //check fixed responses
  arma::vec output = data.col(4);
  arma::vec outputs = unique( output );
  int num_outputs = outputs.n_elem;
  arma::mat complete_response_mat( rows_strategies , num_outputs, arma::fill::zeros );
  arma::mat complete_tremble_mat( rows_strategies , num_outputs , arma::fill::zeros );
  arma::vec complete_tremble_vec( rows_strategies , arma::fill::zeros );
  complete_response_mat.fill(arma::datum::nan);
  complete_response_mat.tail_cols( num_non_zero_outputs ) = strategies.cols( 1 , num_non_zero_outputs );
  arma::uvec index_first_free( rows_strategies );
  for(int i = 0; i < rows_strategies; i++) {
    arma::rowvec responses_row = complete_response_mat.row( i );
    arma::uvec indices_free_responses_row = find_nonfinite( responses_row );
    if( indices_free_responses_row.n_elem > 0 ){
      index_first_free(i) = indices_free_responses_row(0);
      arma::vec fixed_responses_row = responses_row( find_finite( responses_row ) );
      arma::vec fixed_mixed_responses_row = fixed_responses_row( find( fixed_responses_row != 0 && fixed_responses_row != 1 ) );
      int num_fixed_responses_row =  fixed_responses_row.n_elem;
      if( num_fixed_responses_row > 0 ){
        if( arma::max( fixed_responses_row ) > 1  ){
          stop("Fixed responses cannot exceed one.");
        }
        if( arma::max( fixed_responses_row ) < 0  ){
          stop("Fixed responses cannot be negative.");
        }
        if( accu( fixed_responses_row ) > 1 ){
          stop("The sum of fixed shares cannot exceed one. stratEst cannot proceed with the current values.");
        }
        else if(  num_fixed_responses_row == num_outputs && accu( fixed_responses_row ) != 1 ){
          stop("The fixed shares must sum to one. stratEst cannot proceed with the current values.");
        }
        if( num_fixed_responses_row < (num_outputs-1)  && fixed_mixed_responses_row.n_elem != 0 && response == "pure" ){
          stop("It is not possible to estimate pure strategy parameters in a row where other parameters are fixed at mixed values. Estimate mixed parameters or change the fixed values to zero or one.");
        }
        if( num_fixed_responses_row < num_non_zero_outputs && ( r_responses != "no" || r_trembles != "no" ) ){
          stop("It is not possible to fix only a subset of response probabilities in the same state if restrictions apply. Fix either all or no probability in each row of strategies.");
        }
      }
    }
  }

  // check samples
  arma::vec sample_id = data.col(5);
  arma::uvec zeros_samples = find( sample_id <= 0 );
  int num_zeros_samples = zeros_samples.n_elem;
  if ( num_zeros_samples != 0  ){
    stop("The variable sample in of the data frame must contain values greater than zero.");
  }
  for(int j = 0; j < num_ids; j++) {
    arma::vec unique_sample = unique( sample_id( find( id_vec == j+1 ) ) );
    int num_unique_sample = unique_sample.n_elem;
      if( num_unique_sample > 1 ){
        stop("The variable sample contains different values for the same subject.");
      }
  }

  // generate samples 1,2,...,N
  arma::vec unique_samples = unique( sample_id ) ;
  int num_samples = unique_samples.n_elem ;
  for(int j = 0; j < num_samples; j++) {
    data.col(5).replace( unique_samples(j) , j+1 );
  }
  sample_id = data.col(5);

  //check fixed shares
  int num_rows_shares = complete_share_mat.n_rows;
  int num_cols_shares = complete_share_mat.n_cols;
  if( k !=  num_rows_shares ){
    stop("Matrix 'shares' must have as many rows as there are strategies.");
  }
  if( num_samples !=  num_cols_shares ){
    stop("Matrix 'shares' must have as many columns as there are samples.");
  }
  arma::vec fixed_shares = complete_share_mat( find_finite( complete_share_mat ) );
  int num_fixed_shares = fixed_shares.n_elem;
  if ( num_fixed_shares > 0 ){
    if( arma::max( fixed_shares ) > 1  ){
      stop("Fixed shares cannot exceed one.");
    }
    if( arma::min( fixed_shares ) < 0  ){
      stop("Fixed shares cannot be negative.");
    }
    arma::mat with_zero_share_mat = complete_share_mat;
    with_zero_share_mat.replace(arma::datum::nan, 0);
    arma::rowvec sums_of_shares = sum( with_zero_share_mat , 0 );
    if( any( sums_of_shares > 1 ) ){
      stop("The sum of all fixed values in a column of 'shares' cannot exceed one. It is not possible to proceed with the current values.");
    }
    for(int j = 0; j < num_samples; j++) {
      arma::vec complete_share_col = complete_share_mat(arma::span::all,j);
      arma::vec fixed_shares_col = complete_share_col( find_finite( complete_share_col ) );
      int num_fixed_shares_col = fixed_shares_col.n_elem;
      if( num_fixed_shares_col == (k-1) ){
        complete_share_col.replace( arma::datum::nan , 1 - accu( fixed_shares_col ) );
      }
      complete_share_mat(arma::span::all,j) = complete_share_col;
    }
  }

  // check covariates
  arma::mat covariate_mat( num_ids , covariates.n_cols + 1 , arma::fill::ones );
  if( LCR ){
    if( num_rows_shares == 1 ){
      stop("Latent class regression requires more than one strategy.");
    }
    if( num_samples > 1 ){
      stop("Latent class regression cannot be run with more than one sample.");
    }
    int cols_covariates = covariates.n_cols;
    int rows_covariates = covariates.n_rows;
    if( rows_covariates != rows_data  ){
      stop("Covariate matrix must have as many rows as data.");
    }
    if( num_fixed_shares > 0 ){
      stop("Shares cannot be fixed with covariates.");
    }
    arma::mat incomplete_covariate_mat( num_ids , covariates.n_cols , arma::fill::ones );
    for(int j = 0; j < num_ids; j++) {
      for(int c = 0; c < cols_covariates; c++) {
        arma::vec covariate = covariates.col(c);
        arma::vec unique_covariate = unique( covariate( find( id_vec == j+1 ) ) );
        int num_unique_covariate = unique_covariate.n_elem;
        if( num_unique_covariate > 1 ){
          stop("Covariate matrix contains different values of the same variable for the same id.");
        }
        else{
          incomplete_covariate_mat( j , c ) = unique_covariate(0);
        }
      }
    }
    arma::mat intercept( num_ids , 1 , arma::fill::ones );
    covariate_mat = join_rows( intercept , incomplete_covariate_mat );
  }

  // check cluster
  int num_clusters = 0;
  arma::vec cluster_id_vec( num_ids , arma::fill::zeros );
  bool CL = cluster.n_elem > 1;
  if( CL ){
    SE = "bs";
    arma::uvec zeros_cluster = find( cluster <= 0 );
    int num_zeros_cluster = zeros_cluster.n_elem;
    if ( num_zeros_cluster != 0 ){
      stop("Cluster must contain values greater than zero.");
    }
    int elem_cluster = cluster.n_elem;
    if( elem_cluster != rows_data  ){
      stop("The cluster vector must have the same number of elements as there are rows in data.");
    }
    for(int j = 0; j < num_ids; j++) {
      arma::vec unique_cluster = unique( cluster( find( id_vec == j+1 ) ) );
      int num_unique_cluster = unique_cluster.n_elem;
      if( num_unique_cluster > 1 ){
        stop("Indivdiuals must be nested within clusters, i.e. the data of one individual appears only in one cluster.");
      }
      else{
        cluster_id_vec(j) = unique_cluster(0);
      }
    }
    // generate cluster 1,2,...,C
    arma::vec unique_cluster = unique( cluster );
    num_clusters = unique_cluster.n_elem;
    for(int j = 0; j < num_clusters; j++) {
      cluster.replace( unique_cluster(j), j+1 );
      cluster_id_vec.replace( unique_cluster(j), j+1 );
    }
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // sort data matrix
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // sort by id, supergame, round
  arma::mat sorted_data( rows_data , cols_data , arma::fill::zeros );
  int line = 0;
  for( int i = 1; i <= num_ids; i++) {
    arma::vec matches_sbj = match_vec( find( data.col(0) == i ) );
    arma::vec unique_matches_sbj = unique( matches_sbj );
    int num_matches_sbj = unique_matches_sbj.n_elem;
    for(int j = 0; j < num_matches_sbj; j++) {
      matches_sbj.replace( unique_matches_sbj(j), j+1 );
      arma::vec rounds_sbj = round_vec( find( id_vec == i && match_vec == unique_matches_sbj(j) ) );
      arma::vec unique_rounds_sbj = unique( rounds_sbj );
      int num_rounds_sbj_match = unique_rounds_sbj.n_elem;
      for(int l = 0; l < num_rounds_sbj_match; l++) {
        sorted_data.row( line ) = data.rows( find( data.col(0) == i && data.col(1) == j+1 && data.col(2) == l+1 ) );
        line = line + 1;
      }
    }
  }
  data = sorted_data;

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // initialize objects
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  arma::vec supergame = data.col(1);
  arma::vec period = data.col(2);
  arma::vec input = data.col(3);
  arma::vec inputs = unique( input );

  // calculate sample of ids
  arma::vec sample_of_ids( num_ids , arma::fill::zeros );
  arma::vec num_ids_sample( num_samples , arma::fill::zeros );
  for (int i = 0; i < num_ids; i++) {
    arma::vec unique_sample_of_ids = unique( sample_id( find( id_vec == i+1 ) ) );
    sample_of_ids(i) = unique_sample_of_ids(0);
    num_ids_sample( unique_sample_of_ids(0) - 1 ) += 1;
  }

  // calculate strategy ids
  arma::colvec complete_strat_id( rows_strategies , arma::fill::zeros );
  complete_strat_id( find( state == 1 ) ).ones();
  double current_id = 0;
  for (int i = 0; i < rows_strategies; i++) {
    current_id += complete_strat_id(i);
    complete_strat_id(i) = current_id;
  }

  // calculate state matrix for strategies
  arma::mat state_mat( rows_data , k , arma::fill::ones );
  arma::cube response_cube( max_state , num_non_zero_outputs , k , arma::fill::zeros );
  response_cube.fill(-1);
  for (int i = 0; i < k; i++) {
    arma::mat strategy = strategies.rows( find( complete_strat_id == (i+1) ) );
    response_cube( arma::span( 0 , ( strategy.n_rows - 1 ) ) , arma::span( 0 , (num_non_zero_outputs-1) ) , arma::span( i , i ) ) = strategy.cols( 1 , num_non_zero_outputs );
    for (int j = 0; j < rows_data; j++) {
      if( period(j) == 1){
        if( input(j) != 0 ){
          state_mat(j,i) = strategy( 0 , ( num_non_zero_outputs + input(j) ) );
        }
      }
      else{
        state_mat(j,i) = strategy( state_mat( j-1 , i ) - 1 , ( num_non_zero_outputs + input(j) ) ); //
      }
    }
  }

  // accumulate number of observed responses conditional on input (rows ids, cols states, slices strats)
  arma::cube complete_output_cube( rows_strategies , num_outputs , num_ids , arma::fill::zeros);
  arma::mat complete_sum_outputs_mat( rows_strategies , num_outputs , arma::fill::zeros);
  arma::cube complete_sum_outputs_cube( rows_strategies , num_outputs , num_ids , arma::fill::zeros);
  arma::mat complete_output_ones( rows_data, output.n_cols , arma::fill::ones );
  for( int l = 0; l < rows_strategies; l++) {
    for( int m = 0; m < num_outputs; m++) {
      for (int i = 0; i < num_ids; i++) {
        complete_output_cube(l,m,i) += accu( complete_output_ones.rows( find( id_vec == (i+1) && state_mat.col( complete_strat_id(l) - 1 ) == state(l) && output == unique_outputs(m) ) ) );
      }
    }
  }
  complete_sum_outputs_mat = sum( complete_output_cube , 1 );
  for (int i = 0; i < num_ids; i++) {
    complete_sum_outputs_cube.slice(i) = repmat( complete_sum_outputs_mat.col(i) , 1 , num_outputs );
  }

  // ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // indices and selection matrices
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // generate complete indices shares
  arma::mat complete_indices_shares( k , num_samples , arma::fill::zeros );
  arma::mat copy_shares_mat = complete_share_mat;
  copy_shares_mat.replace( arma::datum::nan , 1 );
  int num_finite_shares = 1;
  for (int i = 0; i < k; i++) {
    for (int j = 0; j < num_samples; j++) {
      if( copy_shares_mat(i,j) == 1 ){
        complete_indices_shares(i,j) = num_finite_shares;
        num_finite_shares = num_finite_shares + 1;
      }
    }
  }

  arma::mat restriction_base_row( 1 , num_outputs , arma::fill::ones );
  for (int i = 0; i < num_outputs; i++) {
    restriction_base_row(0,i) += i;
  }
  arma::mat complete_indices_responses = repmat( restriction_base_row, rows_strategies , 1 );
  arma::mat complete_indices_trembles( rows_strategies , num_outputs , arma::fill::ones );
  if( r_responses == "no" ){
    for (int i = 0; i < rows_strategies; i++) {
      complete_indices_responses.row(i) += i*num_outputs;
    }
  }
  else if( r_responses == "strategies" ){
    for (int i = 0; i < rows_strategies; i++) {
      complete_indices_responses.row(i) += ( complete_strat_id(i) - 1 )*num_outputs;
    }
  }
  else if( r_responses == "states" ){
    for (int i = 0; i < rows_strategies; i++) {
      complete_indices_responses.row(i) += (state(i)-1)*num_outputs;
    }
  }
  if( r_trembles == "no" ){
    for (int i = 0; i < rows_strategies; i++) {
      complete_indices_trembles.row(i) += i;
    }
  }
  else if( r_trembles == "strategies" ){
    for (int i = 0; i < rows_strategies; i++) {
      complete_indices_trembles.row(i) += ( complete_strat_id(i) - 1 );
    }
  }
  else if( r_trembles == "states" ){
    for (int i = 0; i < rows_strategies; i++) {
      complete_indices_trembles.row(i) += (state(i)-1);
    }
  }
  arma::mat select_responses = complete_indices_responses;
  arma::mat select_trembles = complete_indices_trembles;

  // restriction matrices in case of selection
  if( select == "responses" ){
    complete_indices_responses = repmat( restriction_base_row, rows_strategies , 1 );
    for (int i = 0; i < rows_strategies; i++) {
      complete_indices_responses.row(i) += i*num_outputs;
    }
  }
  else if( select == "trembles"){
    complete_indices_trembles.fill(1);
    for (int i = 0; i < rows_strategies; i++) {
      complete_indices_trembles.row(i) += i;
    }
  }
  else if( select == "both" || select == "all" ){
    complete_indices_trembles.fill(1);
    complete_indices_responses = repmat( restriction_base_row, rows_strategies , 1 );
    for (int i = 0; i < rows_strategies; i++) {
      complete_indices_trembles.row(i) += i;
      complete_indices_responses.row(i) += i*num_outputs;
    }
  }

  // delete fixed responses and responses to sum from response restriction mat
  arma::umat indices_fixed_responses = find_finite( complete_response_mat );
  complete_indices_responses( indices_fixed_responses ).fill(0);
  arma::mat complete_responses_to_sum( complete_indices_responses.n_rows , complete_indices_responses.n_cols , arma::fill::zeros  );
  for (int i = 0; i < rows_strategies; i++) {
    for (int j = 0; j < num_outputs; j++) {
      if( complete_indices_responses(i,j) > 0 ){
        complete_responses_to_sum(i,j) = 1;
        j = complete_indices_responses.n_cols;
      }
    }
  }
  arma::umat complete_indices_responses_to_sum = find( complete_responses_to_sum == 1 );
  complete_indices_responses( complete_indices_responses_to_sum ).fill(0);
  arma::vec unique_non_zero_restrictions = unique( complete_indices_responses( find( complete_indices_responses != 0 ) ) );
  int num_unique_non_zero_restrictions = unique_non_zero_restrictions.n_elem;
  for (int i = 0; i < num_unique_non_zero_restrictions; i++) {
    complete_indices_responses.replace( unique_non_zero_restrictions(i) , i+1 );
  }

  // delete mixed responses from tremble restriction vec
  arma::mat row_with_mixed_responses( rows_strategies , num_outputs , arma::fill::ones );
  if( response == "mixed" ){
    for (int i = 0; i < rows_strategies; i++) {
      for (int j = 0; j < num_outputs; j++) {
        if( complete_response_mat(i,j) == 1 || complete_response_mat(i,j) == 0 ){
          row_with_mixed_responses.row(i).fill(0);
          j = complete_response_mat.n_cols;
        }
      }
    }
    complete_indices_trembles( find( row_with_mixed_responses == 1 ) ).fill(0);
  }
  else{
    arma::mat clean_response_mat = complete_response_mat;
    clean_response_mat.replace( arma::datum::nan, 0 );
    for (int i = 0; i < rows_strategies; i++) {
      for (int j = 0; j < num_outputs; j++) {
        if( clean_response_mat(i,j) == 1 || clean_response_mat(i,j) == 0 ){
          row_with_mixed_responses.row(i).fill(0);
          j = clean_response_mat.n_cols;
        }
      }
    }
    complete_indices_trembles( find( row_with_mixed_responses == 1 ) ).fill(0);
  }
  arma::vec unique_non_zero_trembles = unique( complete_indices_trembles( find( complete_indices_trembles != 0 ) ) );
  int num_unique_non_zero_trembles = unique_non_zero_trembles.n_elem;
  for (int i = 0; i < num_unique_non_zero_trembles; i++) {
    complete_indices_trembles.replace( unique_non_zero_trembles(i) , i+1 );
  }

  // preallocate result field
  int R_num_rows_strategies = 0;
  int R_num_strategies = 0;
  arma::cube R_output_cube( rows_strategies , num_outputs , num_ids , arma::fill::zeros );
  arma::cube R_sum_outputs_cube( rows_strategies , num_outputs , num_ids , arma::fill::zeros );
  arma::vec R_strat_id( rows_strategies , arma::fill::zeros );
  arma::mat R_shares_mat( k , num_samples , arma::fill::zeros );
  arma::mat R_indices_shares( rows_strategies , num_samples , arma::fill::zeros );
  arma::mat R_indices_responses( rows_strategies , num_outputs , arma::fill::zeros );
  arma::mat R_indices_trembles( rows_strategies , num_outputs , arma::fill::zeros );
  arma::mat R_responses_to_sum( rows_strategies , num_outputs , arma::fill::zeros );

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // strategy selection procedure
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  arma::mat complete_strat_id_mat = repmat( complete_strat_id , 1 , num_outputs );
  arma::mat crit_selected_min( 1 , 1 , arma::fill::zeros );
  crit_selected_min.fill( arma::datum::inf );
  arma::vec survivors = unique( complete_strat_id );
  int kill_start = 0;

  // start strategies loop
  for (int K = k; K > 0; K--) {
    if( print_messages == true ){
      Rcout<< "model with " << K << " strategies\n";
    }
    int killed = 0;
    for (int kill = kill_start; kill <= K; kill++) {
      arma::vec survived = survivors( find( survivors != 0 ) );
      arma::vec zero_survivors( survived.n_elem + 1 , arma::fill::zeros );
      zero_survivors( arma::span( 1 , survived.n_elem ) ) = survived;
      arma::vec candidates = survived( find( survived != zero_survivors( kill ) ) );
      int num_candidates = candidates.n_elem;

      // create survivor objects
      arma::mat candidates_mat( rows_strategies , num_outputs , arma::fill::zeros );
      arma::cube candidates_cube( rows_strategies , num_outputs , num_ids , arma::fill::zeros );
      for (int j = 0; j < num_candidates; j++) {
        candidates_mat( find( complete_strat_id_mat == candidates(j) ) ).fill(1);
      }
      int num_rows_response_mat = accu( candidates_mat.col(0) );
      candidates_cube.each_slice() = candidates_mat;

      // create strat_id & share_mat & indices_shares
      arma::vec strat_id( num_rows_response_mat , arma::fill::zeros );
      arma::mat share_mat( num_candidates , num_samples, arma::fill::zeros );
      arma::mat indices_shares( num_candidates , num_samples , arma::fill::zeros );

      if( K == k && kill == 0 ){
        strat_id = complete_strat_id;
        share_mat = complete_share_mat;
        indices_shares = complete_indices_shares;
      }
      else{
        strat_id = complete_strat_id( find ( candidates_mat.col(0) == 1 ) );
        arma::vec unique_strat_id = unique( strat_id );
        int num_unique_strat_id = unique_strat_id.n_elem;
        for (int i = 0; i < num_unique_strat_id; i++) {
          strat_id.replace( unique_strat_id(i) , i+1 );
          share_mat( i , arma::span::all ) = complete_share_mat( unique_strat_id(i) - 1 , arma::span::all );
          indices_shares( i , arma::span::all ) = complete_indices_shares( unique_strat_id(i) - 1 , arma::span::all );
        }
      }

      arma::cube output_cube( num_rows_response_mat , num_outputs , num_ids , arma::fill::zeros );
      if( K == k && kill == 0 ){
        output_cube = complete_output_cube;
      }
      else{
        for (int i = 0; i < num_ids; i++) {
          arma::mat output_slice = complete_output_cube.slice(i);
          arma::vec output_cube_vec = output_slice( find ( candidates_mat == 1 ) );
          arma::mat output_mat = reshape( output_cube_vec , num_rows_response_mat , num_outputs );
          output_cube.slice(i) = output_mat;
        }
      }

      // create sum_outputs_cube
      arma::cube sum_outputs_cube( num_rows_response_mat , num_outputs , num_ids , arma::fill::zeros );
      if( K == k && kill == 0 ){
        sum_outputs_cube = complete_sum_outputs_cube;
      }
      else{
        for (int i = 0; i < num_ids; i++) {
          arma::mat sum_outputs_slice = complete_sum_outputs_cube.slice(i);
          arma::vec sum_outputs_cube_vec = sum_outputs_slice( find ( candidates_mat == 1 ) );
          arma::mat sum_outputs_mat = reshape( sum_outputs_cube_vec , num_rows_response_mat , num_outputs );
          sum_outputs_cube.slice(i) = sum_outputs_mat;
        }
      }

      // update indices_shares
      int num_finite_shares = 1;
      for (int i = 0; i < num_candidates; i++) {
        for (int j = 0; j < num_samples; j++) {
          if( indices_shares(i,j) == 1 ){
            indices_shares(i,j) = num_finite_shares;
            num_finite_shares = num_finite_shares + 1;
          }
        }
      }

      // create indices_responses
      arma::mat indices_responses( num_rows_response_mat , num_outputs , arma::fill::zeros );
      if( K == k && kill == 0 ){
        indices_responses = complete_indices_responses;
      }
      else{
        arma::vec indices_responses_vec = complete_indices_responses( find ( candidates_mat == 1 ) );
        indices_responses = reshape( indices_responses_vec , num_rows_response_mat , num_outputs );
        arma::vec unique_non_zero_indices_responses = unique( indices_responses( find( indices_responses != 0 ) ) );
        int num_unique_non_zero_indices_responses = unique_non_zero_indices_responses.n_elem;
        for (int i = 0; i < num_unique_non_zero_indices_responses; i++) {
          indices_responses.replace( unique_non_zero_indices_responses(i) , i+1 );
        }
      }

      // create indices_trembles
      arma::mat indices_trembles( num_rows_response_mat , num_outputs , arma::fill::zeros );
      if( K == k && kill == 0 ){
        indices_trembles = complete_indices_trembles;
      }
      else{
        arma::vec indices_trembles_vec = complete_indices_trembles( find ( candidates_mat == 1 ) );
        indices_trembles = reshape( indices_trembles_vec , num_rows_response_mat , num_outputs );
        arma::vec unique_non_zero_indices_trembles = unique( indices_trembles( find( indices_trembles != 0 ) ) );
        int num_unique_non_zero_indices_trembles = unique_non_zero_indices_trembles.n_elem;
        for (int i = 0; i < num_unique_non_zero_indices_trembles; i++) {
          indices_trembles.replace( unique_non_zero_indices_trembles(i) , i+1 );
        }
      }

      // create response_mat
      arma::mat response_mat( num_rows_response_mat , num_outputs , arma::fill::zeros );
      if( K == k && kill == 0 ){
        response_mat = complete_response_mat;
      }
      else{
        arma::vec response_mat_vec = complete_response_mat( find ( candidates_mat == 1 ) );
        response_mat = reshape( response_mat_vec , num_rows_response_mat , num_outputs );
      }

      // create tremble_mat
      arma::mat tremble_mat( num_rows_response_mat , num_outputs , arma::fill::zeros );
      if( K == k && kill == 0 ){
        tremble_mat = complete_tremble_mat;
      }
      else{
        arma::vec tremble_mat_vec = complete_tremble_mat( find ( candidates_mat == 1 ) );
        tremble_mat = reshape( tremble_mat_vec , num_rows_response_mat , num_outputs );
      }

      // create responses_to_sum
      arma::mat responses_to_sum( num_rows_response_mat , num_outputs , arma::fill::zeros );
      if( K == k && kill == 0 ){
        responses_to_sum = complete_responses_to_sum;
      }
      else{
        arma::vec responses_to_sum_vec = complete_responses_to_sum( find ( candidates_mat == 1 ) );
        responses_to_sum = reshape( responses_to_sum_vec , num_rows_response_mat , num_outputs );
      }
      arma::umat indices_responses_to_sum = find( responses_to_sum == 1 );

      // identify number of responses, trembles and shares to estimate
      int num_shares_to_est = indices_shares.max();
      int num_responses_to_est = indices_responses.max();
      int num_trembles_to_est = indices_trembles.max();

      // generate indices
      arma::vec num_shares_to_est_col( num_samples , arma::fill::zeros );
      for (int j = 0; j < num_samples; j++) {
        arma::uvec shares_to_est_col = find_nonfinite( share_mat( arma::span::all , j ) );
        num_shares_to_est_col(j) = shares_to_est_col.n_elem;
      }

      // incomplete response mat where values of estimated parameters are added
      response_mat( find_nonfinite( response_mat ) ).fill(0);
      arma::mat incomplete_response_mat = response_mat;
      arma::vec remaining_response_vec = 1 - sum( incomplete_response_mat , 1 );
      arma::mat remaining_response_mat = repmat( remaining_response_vec , 1 , num_outputs );

      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // start EM
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      // best result for k
      arma::field<arma::mat> R_k(23,1);

      // prepare outer runs (index or)
      arma::field<arma::mat> R_outer(23,1);
      double LL_outer_min = 0;

      // outer runs (index ro)
      for (int ro = 0; ro < outer_runs; ro++) {
        arma::field<arma::mat> R_inner(23,1);
        double LL_inner_min = 0;

        // inner runs (index ri)
        for (int ri = 0; ri < inner_runs; ri++) {
          if( print_messages == true ){
            if( kill == 0 ){
              if( select == "strategies" || select == "all" ){
                Rcout<< "estimate complete model (" << ro+1  << "/" << ri+1 << ")     \r";
              }
              else{
                Rcout<< "estimate model (" << ro+1  << "/" << ri+1 << ")     \r";
              }
            }
            else{
              if( integer_strategies == true ){
                Rcout<< "estimate model with " << K - 1 << " strategies " <<  " (" << ro+1  << "/" << ri+1 << ")     \r";
              }
              else{
                Rcout<< "estimate nested model " << kill <<  " (" << ro+1  << "/" << ri+1 << ")     \r";
              }
            }
          }
          //random starting points
          arma::mat start_shares = share_mat;
          arma::uvec shares_to_est = find( indices_shares > 0 );
          if( num_shares_to_est > 0 ){
            start_shares.elem( shares_to_est ).randu();
            for (int j = 0; j < num_samples; j++) {
              arma::vec start_shares_col = start_shares( arma::span::all , j );
              arma::vec indices_shares_col = indices_shares( arma::span::all , j );
              arma::uvec shares_to_est_col = find( indices_shares_col > 0 );
              if(  num_shares_to_est_col(j) == num_candidates ){
                start_shares( arma::span::all , j ) = start_shares_col / accu( start_shares_col );
              }
              else{
                arma::vec fixed_shares_col = start_shares_col( find( indices_shares_col == 0 ) );
                double sum_fixed_shares_col = accu( fixed_shares_col );
                double remaining_shares_col = 1 - sum_fixed_shares_col;
                arma::vec unfixed_shares_col = start_shares_col( shares_to_est_col );
                unfixed_shares_col = unfixed_shares_col / accu( unfixed_shares_col );
                start_shares_col( shares_to_est_col ) = unfixed_shares_col*remaining_shares_col;
                start_shares( arma::span::all, j ) = start_shares_col;
              }
            }
          }

          // // random responses to start (normalized to remaining response)
          arma::vec start_responses = arma::randu( num_responses_to_est );
          for (int i = 0; i < num_responses_to_est; i++) {
            response_mat( find( indices_responses == i+1 ) ).fill( start_responses(i) );
          }
          arma::mat normalized_response_mat = response_mat;
          int num_rows_response_mat = response_mat.n_rows;
          for (int i = 0; i < num_rows_response_mat; i++) {
            arma::rowvec row_to_normalize = normalized_response_mat.row(i);
            arma::rowvec remaining_row = remaining_response_mat.row(i);
            arma::rowvec row_indices_responses = indices_responses.row(i);
            arma::uvec indices_to_normalize = find( row_indices_responses > 0 );
            if ( indices_to_normalize.n_elem > 0 ){
              arma::vec part_one = remaining_row( indices_to_normalize ) % row_to_normalize( indices_to_normalize );
              row_to_normalize( indices_to_normalize ) = part_one.each_row() / ( accu( row_to_normalize( indices_to_normalize ) ) + arma::randu(1) );
            }
            response_mat.row(i) = row_to_normalize;
          }
          for (int i = 0; i < num_responses_to_est; i++) {
            arma::mat normalized_response = unique( response_mat( find( indices_responses == i+1 ) ) );
            start_responses(i) = normalized_response(0,0);
          }
          response_mat( indices_responses_to_sum ).fill(0);
          arma::vec value_to_sum_vec = 1 - sum( response_mat , 1 );
          arma::mat value_to_sum_mat = repmat( value_to_sum_vec , 1 , num_outputs );
          response_mat( indices_responses_to_sum ) = value_to_sum_mat( indices_responses_to_sum );


          // tremble mat to start
          arma::vec start_trembles = arma::randu( num_trembles_to_est );
          start_trembles /= 4;
          for (int i = 0; i < num_trembles_to_est ; i++) {
            arma::mat tremble_value_mat  = tremble_mat;
            arma::umat trembles_to_fill = find( indices_trembles == i+1 );
            tremble_value_mat.fill( start_trembles(i) );
            tremble_mat( trembles_to_fill ) = tremble_value_mat( trembles_to_fill );
          }
          arma::field<arma::mat> R_temp = stratEst_EM( output_cube, sum_outputs_cube, strat_id, start_shares, start_responses, start_trembles, response_mat, tremble_mat, indices_shares, indices_responses, indices_trembles, responses_to_sum, response, 0 , inner_tol_eval, inner_max_eval, sample_of_ids, num_ids_sample );
          arma::mat LL_inner_temp = R_temp(LL_index,0);
          if( LL_inner_temp.has_nan() ){
            failed_inner_runs(0,0) += 1;
            LL_inner_temp.replace(arma::datum::nan, arma::datum::inf);
          }
          if ( ri == 0 || LL_inner_temp(0,0) < LL_inner_min ){
            R_inner = R_temp;
            LL_inner_min = LL_inner_temp(0,0);
          }
        }
        Rcpp::checkUserInterrupt();
        arma::mat pre_eval_vec = R_inner(9,0);
        R_outer = stratEst_EM( output_cube, sum_outputs_cube, strat_id, R_inner(0,0), R_inner(1,0), R_inner(2,0), R_inner(3,0), R_inner(4,0), indices_shares, indices_responses, indices_trembles, responses_to_sum, response, pre_eval_vec(0,0) , outer_tol_eval, outer_max_eval, sample_of_ids, num_ids_sample );
        arma::mat LL_outer_mat = R_outer(LL_index,0);
        if( LL_outer_mat.has_nan() ){
          failed_outer_runs(0,0) += 1;
          LL_outer_mat.replace(arma::datum::nan, arma::datum::inf);
        }
        if ( ro == 0 || LL_outer_mat(0,0) < LL_outer_min ){
          R_k = R_outer;
          LL_outer_min = LL_outer_mat(0,0);
        }
      }

      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // response and tremble fusion procedures
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      arma::mat fused_indices_responses = indices_responses;
      arma::mat fused_indices_trembles = indices_trembles;
      arma::mat new_indices_responses = indices_responses;
      arma::mat new_indices_trembles = indices_trembles;

      if( select != "no" && select != "strategies" && ( num_responses_to_est > 1 || num_trembles_to_est > 1 ) ){
        arma::field<arma::mat> R_temp(23,1);
        arma::field<arma::mat> R_fused(23,1);
        R_fused = R_k;
        arma::mat crit_min = R_k( crit_index , 0 );

        arma::mat pre_eval_vec = R_k(9,0);
        int num_it_responses = 0;
        int num_it_trembles = 0;
        if( select == "responses" ){
          num_it_responses = R_k(1,0).n_elem;
        }
        else if( select == "trembles" ){
          num_it_trembles = R_k(2,0).n_elem;
        }
        else if( select == "both" ){
          num_it_responses = R_k(1,0).n_elem;
          num_it_trembles = R_k(2,0).n_elem;
        }
        else if( select == "all" ){
          num_it_responses = R_k(1,0).n_elem;
          num_it_trembles = R_k(2,0).n_elem;
        }
        int trigger = 1;

        // start response and tremble fusion loop
        while ( trigger == 1 ) {
          trigger = 0;

          // check every combination of responses
          if( select != "strategies" && select != "trembles" ){
            arma::vec unique_non_zero_responses = unique( indices_responses( find( indices_responses != 0 ) ) );
            num_it_responses = unique_non_zero_responses.n_elem;
            for (int r1 = 0; r1 < (num_it_responses-1); r1++) {
              for (int r2 = r1+1; r2 < num_it_responses; r2++) {
                arma::vec c1 = select_responses( find( indices_responses == r1+1 ) );
                arma::vec c2 = select_responses( find( indices_responses == r2+1 ) );
                if(  c1(0) == c2(0) ){
                  arma::vec fused_responses( num_it_responses - 1 , arma::fill::zeros );
                  fused_indices_responses = indices_responses;
                  fused_indices_responses.replace( (r2+1) , (r1+1) );
                  arma::vec unique_non_zero_restrictions = unique( fused_indices_responses( find( fused_indices_responses != 0 ) ) );
                  int num_unique_non_zero_restrictions = unique_non_zero_restrictions.n_elem;

                  for (int i = 0; i < num_unique_non_zero_restrictions; i++) {
                    fused_indices_responses.replace( unique_non_zero_restrictions(i) , i+1 );
                  }
                  arma::mat pre_eval_mat = R_k(11,0);
                  R_temp = stratEst_EM(output_cube, sum_outputs_cube, strat_id, R_k(0,0), fused_responses, R_k(2,0), R_k(3,0), R_k(4,0), indices_shares, fused_indices_responses, indices_trembles, responses_to_sum, response, pre_eval_mat(0,0), outer_tol_eval, outer_max_eval, sample_of_ids, num_ids_sample );
                  arma::mat crit_fused = R_temp( crit_index , 0 );
                  crit_fused.replace( arma::datum::nan , arma::datum::inf );
                  if ( crit_fused(0,0) < crit_min(0,0)  ){
                    R_fused = R_temp;
                    crit_min(0,0) = crit_fused(0,0);
                    new_indices_responses = fused_indices_responses;
                    trigger = 1;
                  }
                }
              }
            }
          }

          // check every combination of trembles
          if( select != "strategies" && select != "responses" ){
            arma::vec unique_non_zero_trembles = unique( indices_trembles( find( indices_trembles != 0 ) ) );
            num_it_trembles = unique_non_zero_trembles.n_elem;
            for (int t1 = 0; t1 < (num_it_trembles-1); t1++) {
              for (int t2 = t1+1; t2 < num_it_trembles; t2++) {
                arma::vec c1 = select_trembles( find( indices_trembles == t1+1 ) );
                arma::vec c2 = select_trembles( find( indices_trembles == t2+1 ) );
                if(  c1(0) == c2(0) ){
                  arma::vec fused_trembles( num_it_trembles - 1 , arma::fill::zeros );
                  fused_indices_trembles = indices_trembles;
                  fused_indices_trembles.replace( (t2+1) , (t1+1) );
                  arma::vec unique_non_zero_restrictions = unique( fused_indices_trembles( find( fused_indices_trembles != 0 ) ) );
                  int num_unique_non_zero_restrictions = unique_non_zero_restrictions.n_elem;
                  for (int i = 0; i < num_unique_non_zero_restrictions; i++) {
                    fused_indices_trembles.replace( unique_non_zero_restrictions(i) , i+1 );
                  }
                  arma::mat pre_eval_mat = R_k(11,0);
                  R_temp = stratEst_EM(output_cube, sum_outputs_cube, strat_id, R_k(0,0), R_k(1,0), fused_trembles, R_k(3,0), R_k(4,0), indices_shares, indices_responses, fused_indices_trembles, responses_to_sum, response, pre_eval_mat(0,0), outer_tol_eval, outer_max_eval, sample_of_ids, num_ids_sample );
                  arma::mat crit_fused = R_temp( crit_index , 0 );
                  crit_fused.replace( arma::datum::nan , arma::datum::inf );
                  if ( crit_fused(0,0) < crit_min(0,0)  ){
                    R_fused = R_temp;
                    crit_min(0,0) = crit_fused(0,0);
                    new_indices_responses = indices_responses;
                    new_indices_trembles = fused_indices_trembles;
                    trigger = 1;
                  }
                }
              }
            }
          }
          R_k = R_fused;
          indices_responses = new_indices_responses;
          indices_trembles = new_indices_trembles;
        }  // end response and tremble fusion loop

      } // end if response and tremble fusion

      //check if R_k is better than current best R
      arma::mat crit_selected = R_k( crit_index , 0 );
      crit_selected.replace( arma::datum::nan , arma::datum::inf );
      if( print_messages == true && ( select == "strategies" || select == "all" ) ){
        Rcout<< "\n";
        Rcout<< crit << ": " << crit_selected(0,0);
        Rcout<< "\n";
      }
      if ( crit_selected(0,0) < crit_selected_min(0,0) || ( K == k && kill == 0 ) ){
        R = R_k;
        killed = kill;
        crit_selected_min(0,0) = crit_selected(0,0);
        R_num_rows_strategies = num_rows_response_mat;
        R_num_strategies = num_candidates;
        R_output_cube( arma::span( 0 , num_rows_response_mat - 1 ) , arma::span::all , arma::span::all ) = output_cube;
        R_sum_outputs_cube( arma::span( 0 , num_rows_response_mat - 1 ) , arma::span::all , arma::span::all ) = sum_outputs_cube;
        R_strat_id( arma::span( 0 , num_rows_response_mat - 1 ) ) = strat_id;
        R_shares_mat( arma::span( 0 , num_candidates - 1 ) , arma::span::all ) = share_mat;
        R_indices_shares( arma::span( 0 , num_candidates - 1 ) , arma::span::all ) = indices_shares;
        R_indices_responses( arma::span( 0 , num_rows_response_mat - 1 ) , arma::span::all ) = indices_responses;
        R_indices_trembles( arma::span( 0 , num_rows_response_mat - 1 ) , arma::span::all ) = indices_trembles;
        R_responses_to_sum( arma::span( 0 , num_rows_response_mat - 1 ) , arma::span::all ) = responses_to_sum;
      }
      if( ( select != "strategies" && select != "all" ) || K == 1 ){
        kill = K+1;
      }
      if( kill == 0 && ( select == "strategies" || select == "all" ) ){
        arma::vec all_strats = survivors( find( survivors != 0 ) );
        arma::mat shares_strats = R(0,0);
        int num_strats = shares_strats.n_rows;
        int num_zero_share = 0;
        int first_zero_share = 0;
        for (int s = 0; s < num_strats; s++){
          if( arma::max( shares_strats(s,arma::span::all) ) < 0.01 ){
            num_zero_share = num_zero_share + 1;
            killed = s;
            survivors( find( survivors == all_strats( killed ) ) ).fill(0);
            K = K - 1;      // reduce number of strategies
            killed = 0;     // s already killed
            kill = -1;      // start all over again with reduced set
            if( print_messages == true && num_zero_share == 1 ){
              first_zero_share = s+1;
              Rcout<< "\r";
              Rcout<< "eliminate strategy " << s+1 << " (estimated share(s) smaller than 1 percent)";
            }
            else if( print_messages == true && num_zero_share == 2 ){
              Rcout<< "\r";
              Rcout<< "eliminate strategies " << first_zero_share << ", " << s+1;
            }
            else if( print_messages == true && num_zero_share > 2 ){
              Rcout<< ", " << s+1;
            }
          }
        }
        if( print_messages == true && kill == -1 ){
          if( num_zero_share > 1 ){
            Rcout<< " (estimated population share(s) smaller than 1 percent)";
          }
          Rcout<< "\nDONE\n";
          Rcout<< "model with " << K << " strategies\n";
        }
      }

      if( integer_strategies == true && kill == 1  ){
        kill = K;
      }

    } // end kill loop


    if( killed == 0 || ( ( select != "strategies" && select != "all" ) )){
      K = 0;
      if( print_messages == true && (select == "strategies" || select == "all") ){
        Rcout<< "\r";
        Rcout<< "end of strategy selection: no strategy can be eliminated";
      }
    }
    else if ( K <= 2 ){
      arma::vec survived = survivors( find( survivors != 0 ) );
      survivors( find( survivors == survived( killed - 1 ) ) ).fill(0);
      K = 0;
    }
    else{
      arma::vec survived = survivors( find( survivors != 0 ) );
      survivors( find( survivors == survived( killed - 1 ) ) ).fill(0);
      kill_start = 1;
      if( print_messages == true ){
        Rcout<< "\r";
        if( integer_strategies == true ){
          Rcout<< "eliminate one strategy (information criterion)";
        }
        else{
          Rcout<< "eliminate strategy " << survived( killed - 1 ) << " (information criterion)";
        }
      }
    }
    if( print_messages == true ){
      Rcout<< "\nDONE\n";
    }
  } // end strategy selection loop

  // dicument which strategy survived the selection
  int num_elem_sid = sid.n_elem;
  int num_survivors = survivors.n_elem;
  arma::vec sid_zeros( num_elem_sid , arma::fill::zeros );
  for (int suv = 0; suv < num_survivors; suv++) {
    sid_zeros( find ( sid == survivors(suv) ) ).fill(1);
  }
  sid = sid( find( sid_zeros == 1 ) );

  //store fused objects
  arma::mat shares_new = R(0,0);
  k = shares_new.n_rows;
  arma::cube output_cube = R_output_cube( arma::span( 0 , R_num_rows_strategies - 1 ) , arma::span::all , arma::span::all );
  arma::cube sum_outputs_cube = R_sum_outputs_cube( arma::span( 0 , R_num_rows_strategies - 1 ) , arma::span::all , arma::span::all );
  arma::vec strat_id = R_strat_id( arma::span( 0 , R_num_rows_strategies - 1 ) );
  arma::umat shares_to_est = find_nonfinite( R_shares_mat( arma::span( 0 , R_num_strategies - 1 ) , arma::span::all ) );
  int num_shares_to_est = shares_to_est.n_elem;
  arma::mat indices_shares = R_indices_shares( arma::span( 0 , k - 1 ) , arma::span::all );
  arma::mat indices_responses = R_indices_responses( arma::span( 0 , R_num_rows_strategies - 1 ) , arma::span::all );
  arma::mat indices_trembles = R_indices_trembles( arma::span( 0 , R_num_rows_strategies - 1 ) , arma::span::all );
  arma::mat responses_to_sum = R_responses_to_sum( arma::span( 0 , R_num_rows_strategies - 1 ) , arma::span::all );

  arma::mat responses_new = R(1,0);
  arma::mat trembles_new = R(2,0);

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // latent class regression
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  if( LCR ){
    if( k == 1 ){
      Rcout<<"warning: only one strategy remains after selection. Cannot proceed with latent class regeression.\n";
      LCR = false;
    }
    else{
      bool replaced = false;
      if( print_messages == true ){
        Rcout<<"estimate latent-class regression model\n";
      }
      //initialize empty coefficients
      arma::mat shares_no_lcr = R(0,0);
      int num_coefficients_to_est = ( covariate_mat.n_cols * (k-1) );
      arma::mat R_shares = R(0,0);
      arma::uvec shares_to_est_lcr = find_nonfinite( shares.col(0) );
      arma::mat pre_eval_vec = R(9,0);
      arma::field<arma::mat> R_LCR_inner(23,1);
      arma::vec LL_inner_min( 1 , 1 , arma::fill::zeros );
      LL_inner_min.fill( arma::datum::inf );
      // inner runs (index ri)
      for (int ri = 0; ri < LCR_runs; ri++) {
        Rcpp::checkUserInterrupt();
        if( print_messages == true ){
          Rcout<< "lcr run " << ri+1  << " (of " << LCR_runs << ")     \r";
        }
        arma::vec start_coefficients( num_coefficients_to_est , arma::fill::zeros );
        arma::mat coefficient_mat( covariate_mat.n_cols , k-1 , arma::fill::zeros );
        //initialize empty coefficients
        int num_coefficients_to_est = ( covariate_mat.n_cols * (k-1) );
        arma::mat R_shares = R(0,0);
        arma::vec start_intercepts( k-1  , arma::fill::zeros );
        if( ri > 0 ){
          start_intercepts = log( R_shares( arma::span( 1 , ( k-1 ) ) , 0 ) / R_shares(0,0) );
          arma::mat random_coefficient_mat( covariate_mat.n_cols , k-1 , arma::fill::randn );
          coefficient_mat = random_coefficient_mat;
        }
        coefficient_mat.row(0) = start_intercepts.t();
        arma::mat start_coefficients_mat  = reshape( coefficient_mat , num_coefficients_to_est , 1 );
        start_coefficients = start_coefficients_mat.col(0);
        arma::mat pre_eval_vec = R(9,0);
        arma::uvec coefficients_to_est = find_finite( start_coefficients );
        arma::field<arma::mat> R_temp = stratEst_LCR_EM( output_cube, sum_outputs_cube, strat_id, covariate_mat, R(0,0), R(1,0), R(2,0), start_coefficients, R(3,0), R(4,0), coefficient_mat, shares_to_est_lcr, indices_responses, indices_trembles, true, coefficients_to_est, responses_to_sum, response, pre_eval_vec(0,0), LCR_tol_eval, LCR_max_eval, newton_stepsize, penalty );
        arma::mat LL_inner_temp = R_temp(LL_index,0);
        arma::mat R_temp_shares = R_temp(0,0);
        arma::mat R_temp_coefficients = R_temp(13,0);
        arma::mat R_temp_responses = R_temp(1,0);
        arma::mat R_temp_trembles = R_temp(2,0);
        if( LL_inner_temp.has_nan() || R_temp_shares.has_nan() || R_temp_coefficients.has_nan() || R_temp_responses.has_nan() || R_temp_trembles.has_nan() ){
          failed_inner_LCR_runs(0,0) += 1;
          LL_inner_temp.fill( arma::datum::inf );
        }
        if (  LL_inner_temp(0,0) > 0 && LL_inner_temp(0,0) < LL_inner_min(0,0) ){
          R = R_temp;
          LL_inner_min(0,0) = LL_inner_temp(0,0);
          replaced = true;
        }
      }
      if( replaced == false ){
        stop("Latent class regression failed. Likelihood values are infinite.");
      }
      if( print_messages == true ){
        Rcout<< "\nDONE\n";
      }
      if( arma::max( vectorise(shares_no_lcr) - vectorise(R(0,0)) ) > 0.05 ){
        Rcout<<"warning: shares of latent class model deviate substantially from the model without covariates.\n";
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // analytical standard erros & convergence
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  arma::field<arma::mat> R_SE(17,1);
  R_SE = stratEst_SE( output_cube, sum_outputs_cube, strat_id, covariate_mat, R(0,0), R(1,0), R(2,0), R(13,0), R(3,0), R(4,0), R(14,0), R(12,0), R(18,0), R(16,0), R(17,0), R(15,0), indices_shares, indices_responses, indices_trembles, response, R(19,0), CL, cluster_id_vec, LCR, sample_of_ids, num_ids_sample );

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // bootstrapped standard errors
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // preallocate reslut matrices
  int num_responses_to_est = R(1,0).n_elem;
  int num_trembles_to_est = R(2,0).n_elem;
  arma::mat BS_shares_SE( num_shares_to_est , 1 , arma::fill::zeros );
  arma::mat BS_coefficients_SE( covariate_mat.n_cols * (k-1) , 1 , arma::fill::zeros  );
  arma::mat BS_responses_SE( num_responses_to_est , 1 , arma::fill::zeros );
  arma::mat BS_trembles_SE( num_trembles_to_est , 1 , arma::fill::zeros  );

  arma::mat indices_shares_BS = indices_shares;

  if( SE == "bs" ){
    arma::mat estimated_shares = R(0,0);
    arma::mat estimated_coefficients = R(13,0);
    arma::mat estimated_responses = R(1,0);
    arma::mat estimated_trembles = R(2,0);
    if( response == "pure" ){
      num_responses_to_est = 0;
    }

    int BS_samples_shares = BS_samples;
    int BS_samples_responses = BS_samples;
    int BS_samples_trembles = BS_samples;
    int BS_samples_coefficients = BS_samples;
    if( print_messages == true ){
      Rcout<<"start bootstrap\n";
    }
    for (int i = 0; i < BS_samples; i++) {
      Rcpp::checkUserInterrupt();
      if( print_messages == true ){
        Rcout<< "sample " << i+1 << " (of " << BS_samples <<")\r";
      }
      arma::field<arma::mat> R_LCR_BS(23,1);
      arma::field<arma::mat> R_NO_LCR_BS(13,1);

      int num_ids_to_sample = num_ids;

      // sample clusters if CL true
      arma::vec sampled_clusters( num_clusters , arma::fill::randu );
      if( CL ){
        sampled_clusters *= num_clusters;
        sampled_clusters = ceil( sampled_clusters );
        int ids_in_sampled_clusters = 0;
        for (int j = 0; j < num_clusters; j++) {
          arma::vec unique_ids_in_sampled_cluster = unique( id_vec( find( cluster == sampled_clusters(j) ) ) );
          ids_in_sampled_clusters = ids_in_sampled_clusters + unique_ids_in_sampled_cluster.n_elem;
        }
        num_ids_to_sample = ids_in_sampled_clusters;
      }
      arma::vec sampled_ids( num_ids_to_sample , arma::fill::zeros );

      // sample ids
      if( CL ){
        arma::mat sampled_ids_mat = unique( id_vec( find( cluster == sampled_clusters(0) ) ) );
        for (int j = 1; j < num_clusters; j++) {
          arma::mat sampled_ids_mat_old = sampled_ids_mat;
          arma::mat sampled_ids_mat_new = unique( id_vec( find( cluster == sampled_clusters(j) ) ) );
          sampled_ids_mat = join_cols(sampled_ids_mat_old,sampled_ids_mat_new);
        }
        sampled_ids = sampled_ids_mat.col(0);
      }
      else{
        sampled_ids.randu();
        sampled_ids *= num_ids_to_sample;
        sampled_ids = ceil( sampled_ids );
      }
      int num_sampled_ids = sampled_ids.n_elem;
      arma::cube output_cube_BS_sample( output_cube.n_rows , output_cube.n_cols , num_sampled_ids , arma::fill::zeros );
      arma::cube sum_outputs_cube_BS_sample( sum_outputs_cube.n_rows , sum_outputs_cube.n_cols , num_sampled_ids , arma::fill::zeros );
      arma::mat covariate_mat_BS_sample( num_sampled_ids , covariate_mat.n_cols , arma::fill::zeros );
      arma::vec sample_of_ids_BS( num_sampled_ids , arma::fill::zeros );
      arma::vec num_ids_sample_BS( num_samples , arma::fill::zeros );
      for (int j = 0; j < num_sampled_ids; j++) {
        output_cube_BS_sample.slice(j) = output_cube.slice( sampled_ids(j) - 1 );
        sum_outputs_cube_BS_sample.slice(j) = sum_outputs_cube.slice( sampled_ids(j) - 1 );
        arma::vec unique_sample_of_ids_BS = unique( sample_id( find( id_vec == sampled_ids(j) ) ) );
        sample_of_ids_BS(j) = unique_sample_of_ids_BS(0);
        num_ids_sample_BS( unique_sample_of_ids_BS(0) - 1 ) += 1;
        if( LCR ){
          covariate_mat_BS_sample.row(j) = covariate_mat.row( sampled_ids(j) - 1 );
        }
      }
      arma::mat pre_eval_mat( 1 , 1 , arma::fill::zeros );
      arma::mat indices_responses_BS( indices_responses.n_rows , indices_responses.n_cols , arma::fill::zeros );
      arma::mat indices_trembles_BS( indices_trembles.n_rows , indices_trembles.n_cols , arma::fill::zeros );
      // bootstrap shares or coefficients
      if( num_shares_to_est > 0 ){
        arma::mat estimated_shares_BS( k , num_samples , arma::fill::zeros );
        if( LCR == false ){
          R_NO_LCR_BS = stratEst_EM( output_cube_BS_sample, sum_outputs_cube_BS_sample, strat_id, R(0,0), R(1,0), R(2,0), R(3,0), R(4,0), indices_shares_BS, indices_responses_BS, indices_trembles_BS, responses_to_sum, response, 0 , outer_tol_eval, outer_max_eval, sample_of_ids_BS, num_ids_sample_BS );
          estimated_shares_BS = R_NO_LCR_BS(0,0);
          if( estimated_shares_BS.is_finite() && estimated_shares_BS.n_elem > 0 ){
            BS_shares_SE += ( ( estimated_shares_BS( shares_to_est ) - estimated_shares( shares_to_est )) % (estimated_shares_BS( shares_to_est ) - estimated_shares( shares_to_est ) ) );
          }
          else{
            BS_samples_shares = BS_samples_shares - 1;
          }
        }
        else{
          arma::mat estimated_coefficients_BS( covariate_mat.n_cols * (k-1) , 1 , arma::fill::zeros );
          arma::uvec coefficients_to_est_bs = find_finite( R(13,0) );
          R_LCR_BS = stratEst_LCR_EM( output_cube_BS_sample, sum_outputs_cube_BS_sample, strat_id, covariate_mat_BS_sample, R(0,0), R(1,0), R(2,0), R(13,0), R(3,0), R(4,0), R(14,0), shares_to_est, indices_responses_BS, indices_trembles_BS, true, coefficients_to_est_bs, responses_to_sum, response, 0 , outer_tol_eval, outer_max_eval, newton_stepsize, true );
          estimated_coefficients_BS = R_LCR_BS(13,0);
          estimated_shares_BS = R_LCR_BS(0,0);
          if( estimated_shares_BS.is_finite() && estimated_coefficients_BS.is_finite() ){
            BS_shares_SE += ( ( estimated_shares_BS( shares_to_est ) - estimated_shares( shares_to_est )) % (estimated_shares_BS( shares_to_est ) - estimated_shares( shares_to_est ) ) );
            BS_coefficients_SE += ( ( estimated_coefficients_BS( coefficients_to_est_bs ) - estimated_coefficients( coefficients_to_est_bs )) % (estimated_coefficients_BS( coefficients_to_est_bs ) - estimated_coefficients( coefficients_to_est_bs ) ) );
          }
          else{
            BS_samples_coefficients = BS_samples_coefficients - 1;
          }
        }

      }
      // bootstrap responses & trembles
      if( response != "pure" ){
        int num_estimated_responses = estimated_responses.n_elem;
        arma::mat estimated_responses_BS( estimated_responses.n_elem , 1 , arma::fill::zeros );
        arma::mat responses_BS( 1 , 1 );
        for (int j = 0; j < num_estimated_responses; j++) {
          responses_BS.fill( estimated_responses(j,0) );
          arma::mat indices_responses_BS_j = indices_responses_BS;
          indices_responses_BS_j( find( indices_responses == j+1 ) ).fill(1);
          if( LCR ){
            arma::uvec coefficients_to_est_bs = find_finite( R(13,0) );
            R_LCR_BS = stratEst_LCR_EM( output_cube_BS_sample, sum_outputs_cube_BS_sample, strat_id, covariate_mat_BS_sample, R(0,0), responses_BS, R(2,0), R(13,0), R(3,0), R(4,0), R(14,0), shares_to_est, indices_responses_BS_j, indices_trembles_BS, false, coefficients_to_est_bs, responses_to_sum, response, 0 , outer_tol_eval, outer_max_eval, newton_stepsize, penalty );
            arma::mat response_estimate_BS = R_LCR_BS(1,0);
            estimated_responses_BS(j,0) = response_estimate_BS(0,0);
          }
          else{
            R_NO_LCR_BS = stratEst_EM( output_cube_BS_sample, sum_outputs_cube_BS_sample, strat_id, R(0,0), responses_BS, R(2,0), R(3,0), R(4,0), indices_shares_BS, indices_responses_BS_j, indices_trembles_BS, responses_to_sum, response, 0 , outer_tol_eval, outer_max_eval, sample_of_ids_BS, num_ids_sample_BS );
            arma::mat response_estimate_BS = R_NO_LCR_BS(1,0);
            estimated_responses_BS(j,0) = response_estimate_BS(0,0);
          }
        }
        if( estimated_responses_BS.is_finite() ){
          BS_responses_SE += ( (estimated_responses_BS - estimated_responses) % (estimated_responses_BS - estimated_responses ) );
        }
        else{
          BS_samples_responses = BS_samples_responses - 1;
        }
      }
      int num_estimated_trembles = estimated_trembles.n_elem;
      arma::mat estimated_trembles_BS( estimated_trembles.n_elem , 1 , arma::fill::zeros );
      arma::mat trembles_BS( 1 , 1 );
      for (int j = 0; j < num_estimated_trembles; j++) {
        trembles_BS.fill( estimated_trembles(j,0) );
        arma::mat indices_trembles_BS_j = indices_trembles_BS;
        indices_trembles_BS_j( find( indices_trembles == j+1 ) ).fill(1);
        if( LCR ){
          arma::uvec coefficients_to_est = find_finite( R(13,0) );
          R_LCR_BS = stratEst_LCR_EM( output_cube_BS_sample, sum_outputs_cube_BS_sample, strat_id, covariate_mat_BS_sample, R(0,0), R(1,0), trembles_BS, R(13,0), R(3,0), R(4,0), R(14,0), shares_to_est, indices_responses_BS, indices_trembles_BS_j, false, coefficients_to_est, responses_to_sum, response, 0, outer_tol_eval, outer_max_eval, newton_stepsize, penalty );
          arma::mat trembles_estimate_BS = R_LCR_BS(2,0);
          estimated_trembles_BS(j,0) = trembles_estimate_BS(0,0);
        }
        else{
          R_NO_LCR_BS = stratEst_EM( output_cube_BS_sample, sum_outputs_cube_BS_sample, strat_id, R(0,0), R(1,0), trembles_BS, R(3,0), R(4,0), indices_shares_BS, indices_responses_BS, indices_trembles_BS_j, responses_to_sum, response, 0, outer_tol_eval, outer_max_eval, sample_of_ids_BS, num_ids_sample_BS );
          arma::mat trembles_estimate_BS = R_NO_LCR_BS(2,0);
          estimated_trembles_BS(j,0) = trembles_estimate_BS(0,0);
        }
      }
      if( estimated_trembles_BS.is_finite() ){
        BS_trembles_SE += ( (estimated_trembles_BS - estimated_trembles) % (estimated_trembles_BS - estimated_trembles ) );
      }
      else{
        BS_samples_trembles = BS_samples_trembles - 1;
      }
    }
    if( LCR ){
      BS_shares_SE = sqrt( BS_shares_SE/BS_samples_coefficients );
      BS_coefficients_SE = sqrt( BS_coefficients_SE/BS_samples_coefficients );
      if( BS_samples_coefficients/BS_samples < 0.9 ){
        BS_shares_SE.fill(-1);
        BS_coefficients_SE.fill(-1);
      }
    }
    else{
      BS_shares_SE = sqrt( BS_shares_SE/BS_samples_shares );
      if( BS_samples_shares/BS_samples < 0.9 ){
        BS_shares_SE.fill(-1);
      }
    }
    BS_responses_SE = sqrt( BS_responses_SE/BS_samples_responses );
    BS_trembles_SE = sqrt( BS_trembles_SE/BS_samples_trembles );
    if( BS_samples_responses/BS_samples < 0.9 ){
      BS_responses_SE.fill(-1);
    }
    if( BS_samples_trembles/BS_samples < 0.9 ){
      BS_trembles_SE.fill(-1);
    }
    if( print_messages == true ){
      Rcout<< "\nDONE\n";
    }
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // output
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // update strategy parameters
  arma::mat strategy_responses = R(3,0);
  arma::vec finally_survived = survivors( find( survivors != 0 ) );
  int num_finally_survived = finally_survived.n_elem;
  arma::vec strat_id_survived( rows_strategies , arma::fill::zeros );
  for (int i = 0; i < num_finally_survived; i++) {
    strat_id_survived( find( complete_strat_id == finally_survived(i) ) ).fill(1);
  }
  arma::mat final_strategies( strategy_responses.n_rows , cols_strategies , arma::fill::zeros  );
  arma::mat strat_id_survived_mat = repmat( strat_id_survived , 1 , cols_strategies );
  arma::vec strategies_alive = strategies( find( strat_id_survived_mat == 1 ) );
  final_strategies = reshape( strategies_alive , strategy_responses.n_rows , cols_strategies );
  final_strategies.cols( 1 , num_non_zero_outputs ) = strategy_responses.cols( num_outputs - num_non_zero_outputs , num_outputs - 1 );

  //prepare output
  arma::mat final_shares = R(0,0);
  arma::mat final_responses = R(1,0);
  arma::mat final_trembles = R(2,0);
  arma::mat final_indices_responses = R(21,0);
  arma::mat final_indices_trembles = R(22,0).col(0);
  arma::mat final_tremble_vec = R(4,0).col(0);
  arma::mat final_LL =  R(5,0);
  arma::mat final_crit = R(crit_index,0);
  arma::mat final_eval = R(9,0);
  arma::mat final_eps = R(10,0);
  arma::mat final_entropy = R(11,0);
  arma::mat final_i_class = R(12,0);
  arma::mat final_coefficients = R(13,0);
  arma::mat final_coefficient_mat = R(14,0);
  arma::mat final_individual_priors = R(15,0);
  arma::mat final_state_obs = R(20,0);

  arma::mat final_SE_shares = R_SE(0,0);
  arma::mat final_SE_responses = R_SE(1,0);
  arma::mat final_SE_trembles = R_SE(2,0);
  arma::mat final_SE_coefficients = R_SE(3,0);
  arma::mat final_convergence = R_SE(4,0);
  arma::mat final_covar_shares = R_SE(5,0);
  arma::mat final_score_shares = R_SE(6,0);
  arma::mat final_fisher_shares = R_SE(7,0);
  arma::mat final_covar_responses = R_SE(8,0);
  arma::mat final_score_responses = R_SE(9,0);
  arma::mat final_fisher_responses = R_SE(10,0);
  arma::mat final_covar_trembles = R_SE(11,0);
  arma::mat final_score_trembles = R_SE(12,0);
  arma::mat final_fisher_trembles = R_SE(13,0);
  arma::mat final_covar_coefficients = R_SE(14,0);
  arma::mat final_score_coefficients = R_SE(15,0);
  arma::mat final_fisher_coefficients = R_SE(16,0);

  arma::mat final_fit( 1 , 3 , arma::fill::zeros );
  final_fit(0,0) = -final_LL(0,0);
  final_fit(0,1) = final_crit(0,0);
  final_fit(0,2) = final_entropy(0,0);
  arma::mat final_solver( 1 , 3 , arma::fill::zeros );
  final_solver(0,0) = final_eval(0,0);
  final_solver(0,1) = final_eps(0,0);
  if( SE == "bs" ){
    if( LCR ){
      final_SE_shares = reshape( BS_shares_SE , k , num_samples );
      final_SE_coefficients = BS_coefficients_SE;
    }
    else{
      if( num_shares_to_est > 0 ){
        final_SE_shares = reshape( BS_shares_SE , k , num_samples );
      }
    }
    if( num_responses_to_est > 0 && response == "mixed"){
      final_SE_responses = BS_responses_SE;
    }
    if( num_trembles_to_est > 0 ){
      final_SE_trembles = BS_trembles_SE;
    }
  }

  List stats_list = List::create( Named("shares.covar") = final_covar_shares, Named("shares.score") =  final_score_shares, Named("shares.fisher") = final_fisher_shares, Named("responses.covar") = final_covar_responses, Named("responses.score") = final_score_responses, Named("responses.fisher") = final_fisher_responses, Named("trembles.covar") = final_covar_trembles, Named("trembles.score") = final_score_trembles, Named("trembles.fisher") = final_fisher_trembles, Named("coefficients.covar") = final_covar_coefficients, Named("coefficients.score") = final_score_coefficients, Named("coefficients.fisher") = final_fisher_coefficients, Named("tremble_vec") = final_tremble_vec );

  List S = List::create( Named("strategies") = final_strategies, Named("shares") = final_shares, Named("shares.se") = final_SE_shares, Named("responses") = final_responses, Named("response.indices") = final_indices_responses, Named("responses.se") = final_SE_responses, Named("trembles") = final_trembles, Named("tremble.indices") = final_indices_trembles, Named("trembles.se") = final_SE_trembles,  Named("coefficients") = final_coefficients, Named("coefficient.mat") = final_coefficient_mat, Named("coefficients.se") = final_SE_coefficients, Named("fit") = final_fit, Named("solver") = final_solver, Named("state.obs") = final_state_obs, Named("assignments") = final_i_class, Named("priors") = final_individual_priors, Named("stats.list") = stats_list, Named("convergence") = final_convergence, Named("sid") = sid );


  return(S);
}






