# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

rcpparma_hello_world <- function() {
    .Call(`_stratEst_rcpparma_hello_world`)
}

rcpparma_outerproduct <- function(x) {
    .Call(`_stratEst_rcpparma_outerproduct`, x)
}

rcpparma_innerproduct <- function(x) {
    .Call(`_stratEst_rcpparma_innerproduct`, x)
}

rcpparma_bothproducts <- function(x) {
    .Call(`_stratEst_rcpparma_bothproducts`, x)
}

stratEst_cpp <- function(data, strategies, sid, shares, trembles, coefficient_mat, covariates, LCR, cluster, quantile_vec, response = "mixed", specific_shares = TRUE, specific_responses = TRUE, specific_trembles = TRUE, specific_coefficients = TRUE, r_responses = "no", r_trembles = "global", select_strategies = FALSE, select_responses = FALSE, select_trembles = FALSE, min_strategies = 1L, crit = "bic", SE = "analytic", outer_runs = 10L, outer_tol_eval = 1e-10, outer_max_eval = 1000L, inner_runs = 100L, inner_tol_eval = 1e-5, inner_max_eval = 100L, LCR_runs = 100L, LCR_tol_eval = 1e-10, LCR_max_eval = 1000L, BS_samples = 1000L, print_messages = TRUE, integer_strategies = TRUE, newton_stepsize = 1, penalty = FALSE) {
    .Call(`_stratEst_stratEst_cpp`, data, strategies, sid, shares, trembles, coefficient_mat, covariates, LCR, cluster, quantile_vec, response, specific_shares, specific_responses, specific_trembles, specific_coefficients, r_responses, r_trembles, select_strategies, select_responses, select_trembles, min_strategies, crit, SE, outer_runs, outer_tol_eval, outer_max_eval, inner_runs, inner_tol_eval, inner_max_eval, LCR_runs, LCR_tol_eval, LCR_max_eval, BS_samples, print_messages, integer_strategies, newton_stepsize, penalty)
}

stratEst_data_cpp <- function(id, game, period, input, lagged_input, lag, num_ids) {
    .Call(`_stratEst_stratEst_data_cpp`, id, game, period, input, lagged_input, lag, num_ids)
}

