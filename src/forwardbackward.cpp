#include <Rcpp.h>
#include <cmath>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;
using namespace std;
typedef Eigen::MappedSparseMatrix<double, Eigen::RowMajor> MSM;
typedef Eigen::MappedSparseMatrix<double> MSM2;

// [[Rcpp::export]]
List fforward(const SEXP & TM, const Eigen::Map<Eigen::VectorXd> init,
               const SEXP & EM, const IntegerVector & obs) {
  const MSM tmat(as<MSM> (TM));
  const MSM2 emat(as<MSM2> (EM));
  SparseVector<double> svector(init.size());
  SparseMatrix<double> alpha(init.size(), obs.length());
  NumericVector sfactors(obs.length());
  VectorXd init2;
  int i = 0;

  if (IntegerVector::is_na(obs(0))) {
    for (init2 = init; (i < obs.length()) && (IntegerVector::is_na(obs(i))); ++i) {
      init2 = tmat.transpose() * init2;
      sfactors(i) = 1;
      alpha.col(i) = init2.sparseView();
    }
    if (i == obs.length())
      return 0;
    svector = emat.col(obs(i)).cwiseProduct(init2);
  }
  else
    svector = emat.col(obs(i)).cwiseProduct(init);

  sfactors(i) = svector.sum();
  svector = svector / sfactors(i);
  alpha.col(i) = svector;

  for (++i; i < obs.length(); ++i) {
    svector = tmat.transpose() * svector;
    if (!IntegerVector::is_na(obs(i))) {
      svector = svector.cwiseProduct(emat.col(obs(i)));
      sfactors(i) = svector.sum();
      svector = svector / sfactors(i);
    }
    else sfactors(i) = 1;
    alpha.col(i) = svector;
  }

  return(List::create(_["alpha"] = wrap(alpha), _["scalefactors"] = sfactors));
}

// [[Rcpp::export]]
Eigen::SparseMatrix<double> fbackward(const SEXP & TM, const SEXP & EM,
                                      const IntegerVector & obs, NumericVector sfactors) {
  const MSM tmat(as<MSM> (TM));
  const MSM2 emat(as<MSM2> (EM));
  SparseVector<double> svector(tmat.cols());
  SparseMatrix<double> beta(tmat.cols(), obs.length());

  svector.setZero();
  svector = svector + VectorXd::Ones(tmat.cols());

  for (int i = obs.length() - 1;i >= 0; --i) {
    beta.col(i) = svector;
    if (!IntegerVector::is_na(obs(i)))
      svector = svector.cwiseProduct(emat.col(obs(i)));
    svector = tmat * svector;
    svector = svector / sfactors(i);
  }

  return(beta);
}
