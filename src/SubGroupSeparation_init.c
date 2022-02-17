#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP SubGroupSeparation_all_fam_log_marg_lik(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP SubGroupSeparation_compute_counts(SEXP, SEXP);
extern SEXP SubGroupSeparation_compute_counts_nas(SEXP, SEXP);
extern SEXP SubGroupSeparation_fbp(SEXP);
extern SEXP SubGroupSeparation_fbs(SEXP, SEXP);
extern SEXP SubGroupSeparation_fumt_mask(SEXP, SEXP);
extern SEXP SubGroupSeparation_g2_stat(SEXP, SEXP);
extern SEXP SubGroupSeparation_heom_dist(SEXP, SEXP, SEXP, SEXP);
extern SEXP SubGroupSeparation_in_tabu(SEXP, SEXP);
extern SEXP SubGroupSeparation_is_acyclic(SEXP);
extern SEXP SubGroupSeparation_next_comb(SEXP, SEXP);
extern SEXP SubGroupSeparation_score_node(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"SubGroupSeparation_all_fam_log_marg_lik", (DL_FUNC) &SubGroupSeparation_all_fam_log_marg_lik, 5},
    {"SubGroupSeparation_compute_counts",       (DL_FUNC) &SubGroupSeparation_compute_counts,       2},
    {"SubGroupSeparation_compute_counts_nas",   (DL_FUNC) &SubGroupSeparation_compute_counts_nas,   2},
    {"SubGroupSeparation_fbp",                  (DL_FUNC) &SubGroupSeparation_fbp,                  1},
    {"SubGroupSeparation_fbs",                  (DL_FUNC) &SubGroupSeparation_fbs,                  2},
    {"SubGroupSeparation_fumt_mask",            (DL_FUNC) &SubGroupSeparation_fumt_mask,            2},
    {"SubGroupSeparation_g2_stat",              (DL_FUNC) &SubGroupSeparation_g2_stat,              2},
    {"SubGroupSeparation_heom_dist",            (DL_FUNC) &SubGroupSeparation_heom_dist,            4},
    {"SubGroupSeparation_in_tabu",              (DL_FUNC) &SubGroupSeparation_in_tabu,              2},
    {"SubGroupSeparation_is_acyclic",           (DL_FUNC) &SubGroupSeparation_is_acyclic,           1},
    {"SubGroupSeparation_next_comb",            (DL_FUNC) &SubGroupSeparation_next_comb,            2},
    {"SubGroupSeparation_score_node",           (DL_FUNC) &SubGroupSeparation_score_node,           6},
    {NULL, NULL, 0}
};

void R_init_SubGroupSeparation(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

