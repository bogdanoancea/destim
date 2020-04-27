#include <Rcpp.h>
#include <cmath>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;
using namespace std;
typedef Eigen::MappedSparseMatrix<double, Eigen::RowMajor> MSM;


// This is a very simple function to remove duplicates in a matrix
// of float numbers. It only scales well when all but a few rows are
// duplicates

// [[Rcpp::export]]
NumericVector createsteady(const SEXP & TM, bool solverexists = false) {
  const MSM tmat(as<MSM> (TM));
  SparseMatrix<double> A(tmat.rows(), tmat.cols());
  VectorXd x(tmat.rows());
  static SparseLU<SparseMatrix<double, RowMajor>, COLAMDOrdering<int> >   solver;
  A.setIdentity();
  A = tmat.transpose() - A;
  x(x.size() - 1) = 1.0 / A.cols();
  x.head(x.size() - 1) = A.topRightCorner(x.size() - 1, 1) / (-A.cols());
  if (!solverexists)
    solver.analyzePattern(A.topLeftCorner(A.rows() - 1, A.cols() - 1));
  solver.factorize(A.topLeftCorner(A.rows() - 1, A.cols() - 1));
  x.head(x.size() - 1) = solver.solve(x.head(x.size() - 1));
  x = x / x.sum();
  return(wrap(x));
}
