#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec transform_pd_cpp(arma::vec id, arma::vec supergame, arma::vec period, arma::vec output, arma::vec p_output, arma::vec input, arma::vec unique_ids) {
  int num_unique_ids = unique_ids.n_elem;
    for(int i = 0; i < num_unique_ids; i++) {
      arma::vec unique_sg = unique( supergame( find( id == unique_ids(i) ) ) );
      int num_unique_sg = unique_sg.n_elem;
      for(int j = 0; j < num_unique_sg; j++) {
        arma::vec hist = output( find( supergame == unique_sg(j) && id == unique_ids(i) ) );
        arma::vec p_hist = p_output( find( supergame == unique_sg(j) && id == unique_ids(i) ) );
        if ( arma::max( period( find( supergame == unique_sg(j) && id == unique_ids(i) ) ) ) >= 2 ){
          int max_period = max( period( find( supergame == unique_sg(j) && id == unique_ids(i) ) ) );
          for( int k = 2; k <= max_period; k++){
            double last = hist(k-2);
            double p_last = p_hist(k-2);
            arma::vec replace( 1 , arma::fill::zeros);
            replace(0) = 1*last*p_last + 2*last*(1 - p_last) + 3*(1 - last)*p_last + 4*(1 - last)*(1 - p_last);
            input( find( id == unique_ids(i) && supergame == unique_sg(j) && period == k ) ) = replace;
          }
        }
      }
    }
  return(input);
}

