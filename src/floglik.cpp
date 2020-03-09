#include <Rcpp.h>
#include <cmath>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;
using namespace std;
typedef Eigen::MappedSparseMatrix<double> MSM;

// [[Rcpp::export]]
double floglik(SEXP & TM, Eigen::Map<Eigen::VectorXd> values, const Eigen::Map<Eigen::VectorXd> rparams,
                         const Eigen::Map<Eigen::MatrixXd> transmatrix, const Eigen::Map<Eigen::VectorXd> init,
                         SEXP & EM, IntegerVector & obs) {
  MSM tmat(as<MSM> (TM));
  MSM emat(as<MSM> (EM));
  SparseVector<double> svector;
  double isum, output = 0;

  values = transmatrix.rightCols(1) + transmatrix.leftCols(rparams.size()) * rparams ;

  svector = emat.col(obs(0)).cwiseProduct(init);
  isum = svector.sum();
  output += log(isum);
  svector = svector / isum;

  for (int i = 1; i < obs.length(); ++i) {
    svector = tmat * svector;
    if (!IntegerVector::is_na(obs(i))) {
      svector = svector.cwiseProduct(emat.col(obs(i)));
      isum = svector.sum();
      output += log(isum);
      svector = svector / isum;
    }
  }

  return(-output);
}
