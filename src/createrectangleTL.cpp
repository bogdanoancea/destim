#include <Rcpp.h>

using namespace Rcpp;

// This is a very simple function to remove duplicates in a matrix
// of float numbers. It only scales well when all but a few rows are
// duplicates

// [[Rcpp::export]]
IntegerMatrix createrectangleTL(int x, int y) {
  IntegerMatrix mat(2, 9*x*y - 6*x - 6*y + 4);
  int i, j, inci, incj, k = 0;

  for (j = 0; j < y; ++j)
    for (i = 0; i < x; ++i)
      for (inci = -1; inci <= 1; ++inci)
        for (incj = -1; incj <= 1; ++incj)
          if ((i + inci >= 0) && (i + inci < x) && (j + incj >= 0) && (j + incj < y)) {
            mat(0,k) = j * x + i + 1;
            mat(1,k) = (j + incj) * x + i + inci + 1;
            ++k;
          }
  return(mat);
}
