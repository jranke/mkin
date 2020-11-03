#include <R.h>
#include <Rinternals.h>

SEXP SFO_solution(SEXP t, SEXP parent_0_, SEXP k_) {

  int n = length(t);

  double parent_0 = asReal(parent_0_);
  double k = asReal(k_);

  double *pt, *pout;

  SEXP out = PROTECT(allocVector(REALSXP, n));

  pt = REAL(t);
  pout = REAL(out);

  for (int i = 0; i < n; i++) {
    pout[i] = parent_0 * exp(- k * pt[i]);
  }

  UNPROTECT(1);

  return out;
}
