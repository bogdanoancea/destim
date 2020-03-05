#include <Rcpp.h>
#include <cmath>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;
using namespace std;



// This is a very simple function to remove duplicates in a matrix
// of float numbers. It only scales well when all but a few rows are
// duplicates

// [[Rcpp::export]]
Eigen::SparseMatrix<double> createTM(const IntegerMatrix & TL, const NumericVector & transitions, const int & states) {
  int i;

  typedef Triplet<double> T;
  vector<T> tripletList;
  tripletList.reserve(transitions.length());

  for (i = 0; i < transitions.length(); ++i)
    tripletList.push_back(T(TL(0, i) - 1, TL(1, i) - 1, transitions(i)));

  SparseMatrix<double> mat(states, states);
  mat.setFromTriplets(tripletList.begin(), tripletList.end());

  return(mat);
}
