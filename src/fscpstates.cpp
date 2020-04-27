#include <Rcpp.h>
#include <cmath>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;
using namespace std;
typedef Eigen::MappedSparseMatrix<double, Eigen::RowMajor> MSM;
typedef Eigen::MappedSparseMatrix<double> MSM2;

// [[Rcpp::export]]
Eigen::SparseMatrix<double> fscpstates(const SEXP & TM, const SEXP & alpha, const SEXP & beta,
                                       const NumericVector & sfactors,
                                       const SEXP & EM, const IntegerVector & obs) {
  const MSM tmat(as<MSM> (TM));
  const MSM2 amat(as<MSM2> (alpha));
  const MSM2 bmat(as<MSM2> (beta));
  const MSM2 emat(as<MSM2> (EM));

  SparseMatrix<double> output(amat.rows(), (obs.length() - 1) * amat.rows());

  for (int i = 1; i < obs.length(); ++i)
    if (IntegerVector::is_na(obs(i)))
      output.middleCols(amat.rows() * (i - 1), amat.rows()) =
        tmat.transpose().cwiseProduct(bmat.col(i) * amat.col(i - 1).transpose()).transpose() / sfactors(i);
    else
      output.middleCols(amat.rows() * (i - 1), amat.rows()) =
        tmat.transpose().cwiseProduct((bmat.col(i).cwiseProduct(emat.col(obs(i)))) *
          amat.col(i - 1).transpose()).transpose() / sfactors(i);
  return(output);
}
