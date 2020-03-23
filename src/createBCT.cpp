#include <Rcpp.h>
#include <cmath>
#include <RcppEigen.h>
#include <algorithm>

using namespace Rcpp;
using namespace Eigen;
using namespace std;
typedef Eigen::MappedSparseMatrix<double, Eigen::RowMajor> MSM;



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
Eigen::SparseMatrix<double, Eigen::RowMajor> updateCTaddtran(const SEXP & CT, int newt, int stillt,
                                                             int sizeCT, IntegerVector tCT) {
  int i, wi, j1, j2, eqstillt;

  const MSM cmat(as<MSM> (CT));
  SparseMatrix<double, RowMajor> newcmat(cmat.rows(), cmat.cols() + 1);
  typedef Triplet<double> T;
  vector<T> tripletList;

  tripletList.reserve(sizeCT);
  wi = 0;
  eqstillt = -1;

  for (i = 0; (i < cmat.rows()) && (cmat.row(i).nonZeros() == 2); ++i) {
    MSM::InnerIterator it(cmat,i);
    j1 = it.col();
    ++it;
    j2 = it.col();

    if ((j1 == stillt) || (j2 == stillt)) {
      if (eqstillt == -1)
        eqstillt = j1 + j2 - stillt;
      else {
        tripletList.push_back(T(wi,(eqstillt < newt)?eqstillt:eqstillt+1,1));
        tripletList.push_back(T(wi++,(j1 + j2 - stillt < newt)?j1 + j2 - stillt:
                                  j1 + j2 - stillt + 1,-1));
      }
      continue;
    }

    tripletList.push_back(T(wi,(j1 < newt)?j1:j1+1,1));
    tripletList.push_back(T(wi++,(j2 < newt)?j2:j2+1,-1));
  }

  for (j1 = 0; j1 < tCT.length(); ++j1)
    tripletList.push_back(T(wi, (tCT(j1) < newt)?tCT(j1):tCT(j1)+1,1));
  tripletList.push_back(T(wi, newt, 1));
  tripletList.push_back(T(wi++, newcmat.cols() - 1,1));

  while (i < cmat.rows() && cmat.coeff(i, cmat.cols() - 1) == 1) {
    for (MSM::InnerIterator it(cmat, i);it;++it)
      if (it.col() == stillt)
        tripletList.push_back(T(wi, (eqstillt < newt)?eqstillt:eqstillt + 1, it.value()));
      else
        tripletList.push_back(T(wi, (it.col() < newt)?it.col():it.col() + 1, it.value()));
      ++wi;
      ++i;
  }

  while (i < cmat.rows()) {
    for (MSM::InnerIterator it(cmat, i); it; ++it)
      tripletList.push_back(T(wi, (it.col() < newt)?it.col():it.col() + 1, it.value()));
    ++wi;
    ++i;
  }

  newcmat.setFromTriplets(tripletList.begin(), tripletList.end());
  return(newcmat);
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
    if (TL(0, i) == TL(0,j))
      return (TL(1, i) < TL(1, j));
    else
      return (TL(0, i) < TL(0, j));
  });
}

// [[Rcpp::export]]
IntegerVector orderTL(const IntegerMatrix & TL) {
  IntegerVector output(TL.cols());
  vector<int> idx(TL.cols());

  iota(idx.begin(), idx.end(), 0);

  sort(idx.begin(), idx.end(), [&TL] (int i, int j) {
    if (TL(0, i) == TL(0,j))
      return (TL(1, i) < TL(1, j));
    else
      return (TL(0, i) < TL(0, j));
  });
  for(int i = 0; i < TL.cols(); ++i)
    output(i) = idx.at(i) + 1;

  return output;
}

// [[Rcpp::export]]
int findTorder(const IntegerMatrix & TL, const IntegerVector &T) {
  int i, j, k;

  i = 0;
  j = TL.cols() - 1;
  k = (i + j) / 2;

  while(j - i > 1000) {
    if (TL(0,j) < T(0))
      return(j + 2);
    else if ((TL(0,j) == T(0)) && (TL(1,j) < T(1)))
      return(j + 2);
    else if ((TL(0,j) == T(0)) && (TL(1,j) == T(1)))
      return(0);
    else if (TL(0,i) > T(0))
      return(1);
    else if ((TL(0,i) == T(0)) && (TL(1,i) > T(1)))
      return(1);
    else if ((TL(0,i) == T(0)) && (TL(1,i) == T(1)))
      return(0);
    else if ((TL(0,k) == T(0)) && (TL(1,k) == T(1)))
      return(0);
    else if (TL(0,k) < T(0))
      i = k;
    else if ((TL(0,k) == T(0)) && (TL(1,k) < T(1)))
      i = k;
    else
      j = k;
    k = (i + j) / 2;
  }

  for (int l = i;l <= j; ++l) {
    if ((TL(0,l) == T(0)) && (TL(1,l) == T(1)))
      return(0);
    else if (TL(0,l) > T(0))
      return(l + 1);
    else if ((TL(0,l) == T(0)) && (TL(1,l) > T(1)))
      return(l + 1);
  }
  return(j + 2);
}
