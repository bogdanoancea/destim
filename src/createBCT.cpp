#include <Rcpp.h>
#include <cmath>
#include <RcppEigen.h>
#include <algorithm>

using namespace Rcpp;
using namespace Eigen;
using namespace std;
typedef Eigen::MappedSparseMatrix<double, Eigen::RowMajor> MSM;
typedef Eigen::MappedSparseMatrix<double> MSMC;
typedef MSM::InnerIterator MSMIt;



// This is a very simple function to remove duplicates in a matrix
// of float numbers. It only scales well when all but a few rows are
// duplicates
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

Eigen::SparseMatrix<double, Eigen::RowMajor> extractEQ(const MSM & cmat,
                                                       const IntegerVector & stillt) {
  int i, j1, j2, wi;
  typedef Triplet<double> T;
  vector<T> tripletList;

  tripletList.reserve(cmat.rows() * 2);
  for (i = 0, wi = 0; (i < cmat.rows()) && (cmat.row(i).nonZeros() == 2); ++i) {
    MSMIt it(cmat,i);
    bool is_still = false;
    j1 = 0;
    j2 = stillt.length() - 1;
    while (j2 - j1 > 10) {
      if ((stillt(j1) == it.col()) || (stillt(j2) == it.col())) {
        is_still = true;
        break;
      }
      else if (stillt((j1 + j2) / 2) > it.col())
        j2 = (j1 + j2) / 2;
      else
        j1 = (j1 + j2) / 2;
    }
    if (is_still)
      continue;
    for (;j1 <= j2; ++j1)
      if (stillt(j1) == it.col()) {
        is_still = true;
        break;
      }
    if (is_still)
      continue;
    tripletList.push_back(T(wi,it.col(), it.value()));
    ++it;
    tripletList.push_back(T(wi,it.col(), it.value()));
    ++wi;
  }
  SparseMatrix<double, RowMajor> omat(wi, cmat.cols());
  omat.setFromTriplets(tripletList.begin(), tripletList.end());

  return(omat);
}

Eigen::SparseMatrix<double, Eigen::RowMajor>
  extractRMCT(const SparseMatrix<double, Eigen::RowMajor> & cmat) {
  int i;

  i = cmat.rows() - 1;
  if (cmat.coeff(i, cmat.cols() - 1) == 1) {
    SparseMatrix<double, RowMajor> omat(0, cmat.cols());
    return(omat);
  }
  do --i;
  while (cmat.coeff(i, cmat.cols() - 1) != 1);

  return(cmat.bottomRows(cmat.rows() - i - 1));
}

Eigen::SparseMatrix<double, Eigen::RowMajor>
  createEQBCT(const SparseMatrix<double, RowMajor>  & eqmat,
              const SparseMatrix<double, RowMajor> & bcmat,
              IntegerVector stillt) {
  int i, j, k, l,eqvars[2], bottom;
  IntegerVector eqtran(eqmat.cols() - stillt.length() - 1);
  SparseMatrix<double, RowMajor> trmatrix(eqmat.cols() - 1, eqmat.cols() - eqmat.rows() - 1);

  trmatrix.reserve(VectorXi::Constant(trmatrix.rows(), 1));
  l = 0;
  for (j = 0; j < eqmat.rows(); ++j) {
    k = 0;
    for (SparseMatrix<double, RowMajor>::InnerIterator it(eqmat,j); it; ++it)
      eqvars[k++]=it.index();
    for (k = 0; k <= l; ++k)
      if (k == l) {
        trmatrix.insert(eqvars[0], k) = 1;
        trmatrix.insert(eqvars[1], k) = 1;
        eqtran(l++) = eqvars[0];
        break; // Otherwise the loop would never end
      }
      else {
        if (trmatrix.coeff(eqvars[0], k) != 0) {
          trmatrix.insert(eqvars[1], k) = 1;
          break;
        }
        else if (trmatrix.coeff(eqvars[1], k) != 0) {
          trmatrix.insert(eqvars[0], k) = 1;
          break;
        }
        else
          continue;
      }
  }
  // Fill up the remaining columns
  for (i = 0, j = trmatrix.cols() - 1, k = 0; k < trmatrix.rows(); ++k)
    if (trmatrix.row(k).nonZeros() == 0) {
      if (stillt(i) == k) {
        trmatrix.insert(k, j--) = 1;
        ++i;
      }
      else
        trmatrix.insert(k,l++) = 1;
    }
  trmatrix.makeCompressed();

  SparseMatrix<double, RowMajor> tbcmat(bcmat.rows(), trmatrix.cols());

  tbcmat = bcmat.leftCols(bcmat.cols() - 1) * trmatrix;

  vector<int> idx(tbcmat.rows());
  iota(idx.begin(), idx.end(), 0);

  sort(idx.begin(), idx.end(), [&tbcmat] (int i, int j) {
    SparseMatrix<double, RowMajor>::InnerIterator iti(tbcmat, i);
    SparseMatrix<double, RowMajor>::InnerIterator itj(tbcmat, j);
    for (; iti || itj; ++iti, ++itj) {
      if (!iti)
        return(false);
      else if (!itj)
        return(true);
      else if (iti.col() < itj.col())
        return(true);
      else if (itj.col() < iti.col())
        return(false);
      else if (iti.value() < itj.value())
        return(true);
      else if (itj.value() < iti.value())
        return(false);
    }
    return(true); // already ordered
  });

  SparseMatrix<double, RowMajor> omat(bcmat.rows(), bcmat.cols());
  omat.reserve(VectorXi::Constant(bcmat.rows(), l + 1));

  bottom = tbcmat.rows() - 1;
  for (SparseMatrix<double, RowMajor>::InnerIterator it(tbcmat, idx.at(0));
       it; ++it)
    if (it.col() < l)
      omat.insert(bottom, eqtran(it.col())) = it.value();
    else
      omat.insert(bottom, stillt(it.col() - l)) = it.value();
    omat.insert(bottom--, omat.cols() - 1) = 1;
  for (i = 1, j = 0, k = 0; i < tbcmat.rows(); ++i) {
    bool equal;
    SparseMatrix<double, RowMajor>::InnerIterator iti(tbcmat, idx.at(i));
    SparseMatrix<double, RowMajor>::InnerIterator itj(tbcmat, idx.at(j));
    for(equal = true; (iti.col() < l) || (itj.col() < l); ++iti, ++itj)
      if ((iti.col() != itj.col()) || (iti.value() != itj.value())) {
        equal = false;
        j = i;
        break;
      }
    if (equal) {
      omat.insert(k, stillt(itj.col() - l)) = 1;
      omat.insert(k++, stillt(iti.col() - l)) = -1;
    }
    else {
      SparseMatrix<double, RowMajor>::InnerIterator it(tbcmat, idx.at(i));
      for(;it;++it) {
        if (it.col() < l)
          omat.insert(bottom, eqtran(it.col())) = it.value();
        else
          omat.insert(bottom, stillt(it.col() - l)) = it.value();
      }
      omat.insert(bottom--, omat.cols() - 1) = 1;
    }
  }

  omat.makeCompressed();

  return(omat);

}

// [[Rcpp::export]]
Eigen::SparseMatrix<double, Eigen::RowMajor> frbind(Eigen::SparseMatrix<double, Eigen::RowMajor> mat1,
                                                    Eigen::SparseMatrix<double, Eigen::RowMajor> mat2) {
  SparseMatrix<double, RowMajor> omat(mat1.rows() + mat2.rows(), mat1.cols());

  omat.reserve(VectorXi::Constant(omat.rows(),3));
  for (int i = 0; i < mat1.rows(); ++i)
    for (SparseMatrix<double, RowMajor>::InnerIterator it(mat1, i); it; ++it)
      omat.insert(i, it.col()) = it.value();
  for (int i = 0; i < mat2.rows(); ++i)
    for (SparseMatrix<double, RowMajor>::InnerIterator it(mat2,i); it;++it)
      omat.insert(i + mat1.rows(), it.col()) = it.value();

  omat.makeCompressed();
  return(omat);
}

// [[Rcpp::export]]
Eigen::SparseMatrix<double, Eigen::RowMajor> createEQ(const IntegerVector & tran, const int & ncol) {

  SparseMatrix<double, RowMajor> omat(tran.length() - 1, ncol);
  omat.reserve(VectorXi::Constant(omat.rows(), 2));
  int ftran = tran(0);

  for (int i = 0; i < omat.rows(); ++i) {
    omat.insert(i, ftran) = 1;
    omat.insert(i, tran(i + 1)) = -1;
  }
  omat.makeCompressed();

  return(omat);
}

Eigen::SparseMatrix<double, Eigen::RowMajor>
  canonEQ(const SparseMatrix<double, RowMajor> & eqmat) {
  MatrixXi equiv(eqmat.cols() - 1, 10);
  SparseMatrix<double, RowMajor> omat(eqmat.rows(), eqmat.cols());
  int i, j, k, l, eqvars[2], eqcol[2], maxc;

  equiv.setZero();
  k = -1;
  for (i = 0; i < eqmat.rows(); ++i) {
    l = 0;
    for (SparseMatrix<double, RowMajor>::InnerIterator it(eqmat,i); it; ++it, ++l)
      eqvars[l]=it.index();
    if (k >= 0) {
      maxc = equiv.row(eqvars[0]).head(k + 1).maxCoeff(&eqcol[0]);
      if (maxc == 0)
        eqcol[0] = -1;
      maxc = equiv.row(eqvars[1]).head(k + 1).maxCoeff(&eqcol[1]);
      if (maxc == 0)
        eqcol[1] = -1;
    }
    else {
      eqcol[0] = -1;
      eqcol[1] = -1;
    }
    if ((eqcol[0] == -1) && (eqcol[1] == -1)) {
      if (++k > (equiv.cols() - 1))
        equiv.conservativeResize(equiv.rows(), equiv.cols() + 10);
      equiv(eqvars[0], k) = 1;
      equiv(eqvars[1], k) = 1;
    }
    else if (eqcol[0] == -1)
      equiv(eqvars[0], eqcol[1]) = 1;
    else if (eqcol[1] == -1)
      equiv(eqvars[1], eqcol[0]) = 1;
    else if (eqcol[0] != eqcol[1]) {
      equiv.col(eqcol[0]) = equiv.col(eqcol[0]) + equiv.col(eqcol[1]);
      equiv.col(eqcol[1]) = equiv.col(k);
      equiv.col(k--).setZero();
    }
  }

  omat.reserve(VectorXi::Constant(omat.rows(), 2));
  for (i = 0, l = 0; i <= k; ++i) {
    for(int row = 0, j = -1; row < equiv.rows(); ++row) {
      if (equiv(row, i) == 0)
        continue;
      if (j == -1) {
        j = row;
        continue;
      }
      omat.insert(l,j) = 1;
      omat.insert(l++, row) = -1;
    }
  }

  omat.makeCompressed();

  return(omat);
}

// [[Rcpp::export]]
Eigen::SparseMatrix<double, Eigen::RowMajor>
  faddconstraint(const IntegerVector & tran, const IntegerMatrix & TL,
                 const SEXP & CT, const IntegerVector & stillt, int S) {
  const MSM cmat(as<MSM> (CT));
  SparseMatrix<double, RowMajor> oldeqmat, omat, rmctmat;

  SparseMatrix<double, RowMajor> bctmat(S, TL.ncol() + 1);
  bctmat = createBCT(TL, S);
  SparseMatrix<double, RowMajor> newct(tran.length() - 1, TL.ncol() + 1);
  newct = createEQ(tran, cmat.cols());
  oldeqmat = extractEQ(cmat, stillt);
  if (oldeqmat.rows() != 0)
    omat = frbind(oldeqmat, newct);
  else
    omat = newct;
  omat = canonEQ(omat);
  omat = frbind(omat,createEQBCT(omat, bctmat, stillt));
  rmctmat = extractRMCT(cmat);
  if (rmctmat.rows() != 0)
    omat = frbind(omat, rmctmat);

  return(omat);
}

// [[Rcpp::export]]
Eigen::SparseMatrix<double, Eigen::RowMajor>
  faddconstraintm(const SEXP & CT, const SEXP & CTS) {
    const MSM cmat(as<MSM> (CT));
    const MSMC cmat2(as<MSMC> (CTS));
    SparseMatrix<double, RowMajor> cmat2R(cmat2.rows(), cmat2.cols());

    cmat2R = cmat2;

    return(frbind(cmat,cmat2R));
  }
