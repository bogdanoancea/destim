#include <Rcpp.h>
#include <cmath>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;
using namespace std;
typedef Eigen::MappedSparseMatrix<double, Eigen::RowMajor> MSM;
typedef MSM::InnerIterator MSMIt;


// This is a very simple function to remove duplicates in a matrix
// of float numbers. It only scales well when all but a few rows are
// duplicates

// [[Rcpp::export]]
List createtransmatrix(const SEXP & CT) {
  const MSM mat(as<MSM> (CT));
  int i, j, k, l,eqvars[2];

  // Find where is the first constraint that is not an equality
  for (i = 0; i < mat.rows(); ++i)
    if (mat.row(i).nonZeros() != 2)
      break;

    // Remove equality constraints
  SparseMatrix<double, RowMajor> trmatrix(mat.cols() - 1, mat.cols() - i - 1);
  trmatrix.reserve(VectorXi::Constant(trmatrix.rows(), 1));
  l = 0;
  for (j = 0; j < i; ++j) {
    k = 0;
    for (MSMIt it(mat,j); it; ++it)
      eqvars[k++]=it.index();
    for (k = 0; k <= l; ++k)
      if (k == l) {
        trmatrix.insert(eqvars[0], k) = 1;
        trmatrix.insert(eqvars[1], k) = 1;
        ++l;
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
  for (k = 0; k < trmatrix.rows(); ++k)
    if (trmatrix.row(k).nonZeros() == 0) {
      trmatrix.insert(k, l++) = 1;
    }
  trmatrix.makeCompressed();

  MatrixXd mat2(mat.rows() - i, trmatrix.cols() + 1);
  mat2.leftCols(mat2.cols() - 1) = mat.bottomRows(mat.rows() - i).leftCols(mat.cols() - 1) * trmatrix;
  mat2.rightCols(1) = mat.bottomRows(mat.rows() - i).rightCols(1);


  ColPivHouseholderQR<MatrixXd> decomp(mat2.leftCols(mat2.cols() - 1));
  trmatrix = trmatrix * decomp.colsPermutation();
  IntegerVector parnumber(trmatrix.cols() - decomp.rank());
  for (k = decomp.rank(); k < trmatrix.cols(); ++k)
    for (i = 0; i < trmatrix.rows(); ++i)
      if (trmatrix.innerIndexPtr()[i] == k) {
        parnumber(k - decomp.rank()) = i + 1;
        break;
      }


  mat2.rightCols(1) = - mat2.rightCols(1);
  mat2.leftCols(mat2.cols() - 1) = mat2.leftCols(mat2.cols() - 1) * decomp.colsPermutation();
  mat2.rightCols(mat2.cols() - decomp.rank()) = decomp.householderQ().transpose() * mat2.rightCols(mat2.cols() - decomp.rank());
  MatrixXd transmatrix(mat2.cols() - 1, mat2.cols() - decomp.rank());
  transmatrix.topRows(decomp.rank()) = - decomp.matrixR().topLeftCorner(decomp.rank(), decomp.rank()).template
                  triangularView<Upper>().solve(mat2.topRightCorner(decomp.rank(), mat2.cols() - decomp.rank()));
  transmatrix.bottomLeftCorner(mat2.cols() - decomp.rank() - 1, mat2.cols() - decomp.rank() - 1) =
    MatrixXd::Identity(mat2.cols() - decomp.rank() - 1, mat2.cols() - decomp.rank() - 1);
  transmatrix.bottomRightCorner(mat2.cols() - decomp.rank() - 1, 1) =
    MatrixXd::Zero(mat2.cols() - decomp.rank() - 1, 1);
  return(List::create(wrap(trmatrix * transmatrix), parnumber));
}
