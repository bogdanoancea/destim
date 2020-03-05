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
NumericVector inittransitions(const SEXP & CT) {
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
  MatrixXd eqsys(mat2.rows() + mat2.cols() - 1, mat2.rows() + mat2.cols() - 1);
  eqsys.topLeftCorner(mat2.rows(), mat2.rows()) = MatrixXd::Zero(mat2.rows(), mat2.rows());
  eqsys.topRightCorner(mat2.rows(), mat2.cols() - 1) = mat2.leftCols(mat2.cols() - 1);
  eqsys.bottomLeftCorner(mat2.cols() - 1, mat2.rows()) = mat2.leftCols(mat2.cols() - 1).transpose();
  eqsys.bottomRightCorner(mat2.cols() - 1, mat2.cols() - 1) = 2 * MatrixXd::Identity(mat2.cols() - 1, mat2.cols() - 1);
  ColPivHouseholderQR<MatrixXd> decomp(eqsys);
  VectorXd cterms(eqsys.rows());

  while(true) {
    cterms.head(mat2.rows()) = mat2.rightCols(1);
    cterms.tail(mat2.cols() - 1) = (VectorXd::Random(mat2.cols() - 1) + VectorXd::Constant(mat2.cols() - 1, 1)) / 2;
    cterms = decomp.solve(cterms);
    if ((cterms.tail(mat2.cols() - 1).array() > 0).all())
      break;
    else
      for (i = mat2.rows(); i < cterms.size(); ++i)
        if (cterms(i) < 0)
          cterms.segment(i,1) = (VectorXd::Random(1) + VectorXd::Constant(1,1)) / 2;
  }
  return(wrap(trmatrix * cterms.tail(mat2.cols() - 1)));
}
