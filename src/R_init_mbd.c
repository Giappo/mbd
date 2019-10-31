#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Fortran calls */
extern void F77_NAME(mbd_fill1d)(double *vec, int *DIMP, double *parms, int *II);
extern void F77_NAME(mbd_initmod)(void (*steadyparms)(int *, double *));
extern void F77_NAME(mbd_initmodpc)(void (*steadyparms)(int *, double *));
extern void F77_NAME(mbd_runmod)(int *neq, double *t, double *Conc, double *dConc, double *yout, int *ip);
extern void F77_NAME(mbd_runmodpc)(int *neq, double *t, double *Conc, double *dConc, double *yout, int *ip);

static const R_FortranMethodDef FortranEntries[] = {
  {"mbd_fill1d", (DL_FUNC) &F77_NAME(mbd_fill1d),  4},
  {"mbd_initmod", (DL_FUNC) &F77_NAME(mbd_initmod),  1},
  {"mbd_initmodpc", (DL_FUNC) &F77_NAME(mbd_initmodpc),  1},
  {"mbd_runmod", (DL_FUNC) &F77_NAME(mbd_runmod),  6},
  {"mbd_runmodpc", (DL_FUNC) &F77_NAME(mbd_runmodpc),  6},
  {NULL, NULL, 0}
};

void R_init_mbd(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
