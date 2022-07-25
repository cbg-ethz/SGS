#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP SGS_all_fam_log_marg_lik(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP SGS_compute_counts(SEXP, SEXP);
extern SEXP SGS_compute_counts_nas(SEXP, SEXP);
extern SEXP SGS_fbp(SEXP);
extern SEXP SGS_fbs(SEXP, SEXP);
extern SEXP SGS_fumt_mask(SEXP, SEXP);
extern SEXP SGS_g2_stat(SEXP, SEXP);
extern SEXP SGS_heom_dist(SEXP, SEXP, SEXP, SEXP);
extern SEXP SGS_in_tabu(SEXP, SEXP);
extern SEXP SGS_is_acyclic(SEXP);
extern SEXP SGS_next_comb(SEXP, SEXP);
extern SEXP SGS_score_node(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"SGS_all_fam_log_marg_lik", (DL_FUNC) &SGS_all_fam_log_marg_lik, 5},
    {"SGS_compute_counts",       (DL_FUNC) &SGS_compute_counts,       2},
    {"SGS_compute_counts_nas",   (DL_FUNC) &SGS_compute_counts_nas,   2},
    {"SGS_fbp",                  (DL_FUNC) &SGS_fbp,                  1},
    {"SGS_fbs",                  (DL_FUNC) &SGS_fbs,                  2},
    {"SGS_fumt_mask",            (DL_FUNC) &SGS_fumt_mask,            2},
    {"SGS_g2_stat",              (DL_FUNC) &SGS_g2_stat,              2},
    {"SGS_heom_dist",            (DL_FUNC) &SGS_heom_dist,            4},
    {"SGS_in_tabu",              (DL_FUNC) &SGS_in_tabu,              2},
    {"SGS_is_acyclic",           (DL_FUNC) &SGS_is_acyclic,           1},
    {"SGS_next_comb",            (DL_FUNC) &SGS_next_comb,            2},
    {"SGS_score_node",           (DL_FUNC) &SGS_score_node,           6},
    {NULL, NULL, 0}
};

void R_init_SGS(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

