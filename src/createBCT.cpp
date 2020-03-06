#include <Rcpp.h>
#include <cmath>
#include <RcppEigen.h>
#include <algorithm>

using namespace Rcpp;
using namespace Eigen;
using namespace std;



// This is a very simple function to remove duplicates in a matrix
// of float numbers. It only scales well when all but a few rows are
// duplicates
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
Eigen::SparseMatrix<double, Eigen::RowMajor> createBCT(const IntegerMatrix & TL, int S) {
  SparseMatrix<double, RowMajor> mat(S, TL.ncol() + 1);
  int i;

  mat.reserve(VectorXi::Constant(S, TL.ncol() / S + 3));
  for (i = 0; i < S; ++i) {
    mat.insert(i, mat.cols() - 1) = 1;
  }
  for (i = 0; i < TL.ncol(); ++i) {
    mat.insert(TL(0, i) - 1, i) = 1;
  }
  mat.makeCompressed();
  return(mat);
}

// [[Rcpp::export]]
Eigen::SparseMatrix<double, Eigen::RowMajor> createrectangleCT(const IntegerMatrix & TL, int x, int y) {
  SparseMatrix<double, RowMajor> mat(TL.ncol() - 2, TL.ncol() + 1);
  int i, k = 0, move, bmoves[5] = {-1, -1, -1, -1, -1};

  mat.reserve(VectorXi::Constant(TL.ncol() - 1, 4));
  for (i = 0; i < TL.ncol(); ++i) {
    move = TL(0,i) - TL(1,i);
    if (move == 0) {
      if ((TL(0,i) % x > 1) && (TL(0,i) > x) && (TL(0,i) <= x * (y - 1)))
        if (bmoves[0] == -1)
          bmoves[0] = i;
        else {
          mat.insert(k, bmoves[0]) = 1;
          mat.insert(k, i) = -1;
          ++k;
        }
      else if ((TL(0,i) == 1) || (TL(0,i) == x) || (TL(0,i) == x * (y - 1) + 1) || (TL(0,i) == x * y))
        if (bmoves[1] == -1)
          bmoves[1] = i;
        else {
          mat.insert(k, bmoves[1]) = 1;
          mat.insert(k, i) = -1;
          ++k;
        }
      else
        if (bmoves[2] == -1)
          bmoves[2] = i;
        else {
          mat.insert(k, bmoves[2]) = 1;
          mat.insert(k, i) = -1;
          ++k;
        }

    }
    else {
      if ((move == 1) || (move == -1) || (move == x) || (move == -x)) {
        if (bmoves[3] == -1)
          bmoves[3] = i;
        else {
          mat.insert(k, bmoves[3]) = 1;
          mat.insert(k, i) = -1;
          ++k;
        }
      }
      else {
        if (bmoves[4] == -1)
          bmoves[4] = i;
        else {
          mat.insert(k, bmoves[4]) = 1;
          mat.insert(k, i) = -1;
          ++k;
        }
      }
    }
  }
  mat.insert(k, bmoves[0]) = 1;
  mat.insert(k, bmoves[3]) = 4;
  mat.insert(k, bmoves[4]) = 4;
  mat.insert(k++, mat.cols() - 1) = 1;
  mat.insert(k, bmoves[1]) = 1;
  mat.insert(k, bmoves[3]) = 2;
  mat.insert(k, bmoves[4]) = 1;
  mat.insert(k++, mat.cols() -1) = 1;
  mat.insert(k, bmoves[2]) = 1;
  mat.insert(k, bmoves[3]) = 3;
  mat.insert(k, bmoves[4]) = 2;
  mat.insert(k, mat.cols() - 1) = 1;
  mat.makeCompressed();
  return(mat);
}

// [[Rcpp::export]]
bool is_sortedTL(const IntegerMatrix & TL) {
  IntegerMatrix output(2, TL.cols());
  vector<int> idx(TL.cols());

  iota(idx.begin(), idx.end(), 0);

  return is_sorted(idx.begin(), idx.end(), [&TL] (int i, int j) {
    if (TL(1, i) == TL(1,j))
      return (TL(0, i) < TL(0, j));
    else
      return (TL(1, i) < TL(1, j));
  });
}

// [[Rcpp::export]]
IntegerVector orderTL(const IntegerMatrix & TL) {
  IntegerVector output(TL.cols());
  vector<int> idx(TL.cols());

  iota(idx.begin(), idx.end(), 0);

  sort(idx.begin(), idx.end(), [&TL] (int i, int j) {
    if (TL(1, i) == TL(1,j))
      return (TL(0, i) < TL(0, j));
    else
      return (TL(1, i) < TL(1, j));
  });
  for(int i = 0; i < TL.cols(); ++i)
    output(i) = idx.at(i) + 1;

  return output;
}
