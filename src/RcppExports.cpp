// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;


// modification_loc
vec modification_loc(vec inx_name, int K, int mx_k);
RcppExport SEXP _GENESIS_modification_loc(SEXP inx_nameSEXP, SEXP KSEXP, SEXP mx_kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vec >::type inx_name(inx_nameSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type mx_k(mx_kSEXP);
    rcpp_result_gen = Rcpp::wrap(modification_loc(inx_name, K, mx_k));
    return rcpp_result_gen;
END_RCPP
}
// cpnorm
long double cpnorm(long double x);
RcppExport SEXP _GENESIS_cpnorm(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< long double >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(cpnorm(x));
    return rcpp_result_gen;
END_RCPP
}
// sumfactorial
long double sumfactorial(int n);
RcppExport SEXP _GENESIS_sumfactorial(SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(sumfactorial(n));
    return rcpp_result_gen;
END_RCPP
}
// sumfactorial_rev
long double sumfactorial_rev(int n, const int& k0);
RcppExport SEXP _GENESIS_sumfactorial_rev(SEXP nSEXP, SEXP k0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const int& >::type k0(k0SEXP);
    rcpp_result_gen = Rcpp::wrap(sumfactorial_rev(n, k0));
    return rcpp_result_gen;
END_RCPP
}
// cdnorm
long double cdnorm(long double x, long double mean, long double sd, bool loglog);
RcppExport SEXP _GENESIS_cdnorm(SEXP xSEXP, SEXP meanSEXP, SEXP sdSEXP, SEXP loglogSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< long double >::type x(xSEXP);
    Rcpp::traits::input_parameter< long double >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< long double >::type sd(sdSEXP);
    Rcpp::traits::input_parameter< bool >::type loglog(loglogSEXP);
    rcpp_result_gen = Rcpp::wrap(cdnorm(x, mean, sd, loglog));
    return rcpp_result_gen;
END_RCPP
}
// cdbinom
long double cdbinom(const int& k, const int& size, long double prob, bool loglog);
RcppExport SEXP _GENESIS_cdbinom(SEXP kSEXP, SEXP sizeSEXP, SEXP probSEXP, SEXP loglogSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type k(kSEXP);
    Rcpp::traits::input_parameter< const int& >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< long double >::type prob(probSEXP);
    Rcpp::traits::input_parameter< bool >::type loglog(loglogSEXP);
    rcpp_result_gen = Rcpp::wrap(cdbinom(k, size, prob, loglog));
    return rcpp_result_gen;
END_RCPP
}
// cdmultinom3
long double cdmultinom3(const int k0, const int k1, const int k2, vec prob, bool loglog);
RcppExport SEXP _GENESIS_cdmultinom3(SEXP k0SEXP, SEXP k1SEXP, SEXP k2SEXP, SEXP probSEXP, SEXP loglogSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type k0(k0SEXP);
    Rcpp::traits::input_parameter< const int >::type k1(k1SEXP);
    Rcpp::traits::input_parameter< const int >::type k2(k2SEXP);
    Rcpp::traits::input_parameter< vec >::type prob(probSEXP);
    Rcpp::traits::input_parameter< bool >::type loglog(loglogSEXP);
    rcpp_result_gen = Rcpp::wrap(cdmultinom3(k0, k1, k2, prob, loglog));
    return rcpp_result_gen;
END_RCPP
}
// loglikelihood
long double loglikelihood(const vec& par, const vec& betahat, const vec& varbetahat, const vec& ldscore, const int& c0, const vec& Nstar, const int& num_threads);
RcppExport SEXP _GENESIS_loglikelihood(SEXP parSEXP, SEXP betahatSEXP, SEXP varbetahatSEXP, SEXP ldscoreSEXP, SEXP c0SEXP, SEXP NstarSEXP, SEXP num_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type par(parSEXP);
    Rcpp::traits::input_parameter< const vec& >::type betahat(betahatSEXP);
    Rcpp::traits::input_parameter< const vec& >::type varbetahat(varbetahatSEXP);
    Rcpp::traits::input_parameter< const vec& >::type ldscore(ldscoreSEXP);
    Rcpp::traits::input_parameter< const int& >::type c0(c0SEXP);
    Rcpp::traits::input_parameter< const vec& >::type Nstar(NstarSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_threads(num_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(loglikelihood(par, betahat, varbetahat, ldscore, c0, Nstar, num_threads));
    return rcpp_result_gen;
END_RCPP
}
// weight
mat weight(const vec& par, const vec& betahat, const vec& varbetahat, const vec& ldscore, const int& c0, const vec& Nstar, const int& num_threads);
RcppExport SEXP _GENESIS_weight(SEXP parSEXP, SEXP betahatSEXP, SEXP varbetahatSEXP, SEXP ldscoreSEXP, SEXP c0SEXP, SEXP NstarSEXP, SEXP num_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type par(parSEXP);
    Rcpp::traits::input_parameter< const vec& >::type betahat(betahatSEXP);
    Rcpp::traits::input_parameter< const vec& >::type varbetahat(varbetahatSEXP);
    Rcpp::traits::input_parameter< const vec& >::type ldscore(ldscoreSEXP);
    Rcpp::traits::input_parameter< const int& >::type c0(c0SEXP);
    Rcpp::traits::input_parameter< const vec& >::type Nstar(NstarSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_threads(num_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(weight(par, betahat, varbetahat, ldscore, c0, Nstar, num_threads));
    return rcpp_result_gen;
END_RCPP
}
// weight_loglikelihood
List weight_loglikelihood(const vec& par, const vec& betahat, const vec& varbetahat, const vec& ldscore, const int& c0, const vec& Nstar, const int& num_threads);
RcppExport SEXP _GENESIS_weight_loglikelihood(SEXP parSEXP, SEXP betahatSEXP, SEXP varbetahatSEXP, SEXP ldscoreSEXP, SEXP c0SEXP, SEXP NstarSEXP, SEXP num_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type par(parSEXP);
    Rcpp::traits::input_parameter< const vec& >::type betahat(betahatSEXP);
    Rcpp::traits::input_parameter< const vec& >::type varbetahat(varbetahatSEXP);
    Rcpp::traits::input_parameter< const vec& >::type ldscore(ldscoreSEXP);
    Rcpp::traits::input_parameter< const int& >::type c0(c0SEXP);
    Rcpp::traits::input_parameter< const vec& >::type Nstar(NstarSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_threads(num_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(weight_loglikelihood(par, betahat, varbetahat, ldscore, c0, Nstar, num_threads));
    return rcpp_result_gen;
END_RCPP
}
// update_pi1
long double update_pi1(const mat& w, const vec& Nstar, const int& num_threads);
RcppExport SEXP _GENESIS_update_pi1(SEXP wSEXP, SEXP NstarSEXP, SEXP num_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const mat& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const vec& >::type Nstar(NstarSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_threads(num_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(update_pi1(w, Nstar, num_threads));
    return rcpp_result_gen;
END_RCPP
}
// onestep_varcomponent
vec onestep_varcomponent(const vec varcomponent, const mat& w, const vec& betahat, const vec& varbetahat, const vec& ldscore, const vec& Nstar, const int& num_threads);
RcppExport SEXP _GENESIS_onestep_varcomponent(SEXP varcomponentSEXP, SEXP wSEXP, SEXP betahatSEXP, SEXP varbetahatSEXP, SEXP ldscoreSEXP, SEXP NstarSEXP, SEXP num_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec >::type varcomponent(varcomponentSEXP);
    Rcpp::traits::input_parameter< const mat& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const vec& >::type betahat(betahatSEXP);
    Rcpp::traits::input_parameter< const vec& >::type varbetahat(varbetahatSEXP);
    Rcpp::traits::input_parameter< const vec& >::type ldscore(ldscoreSEXP);
    Rcpp::traits::input_parameter< const vec& >::type Nstar(NstarSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_threads(num_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(onestep_varcomponent(varcomponent, w, betahat, varbetahat, ldscore, Nstar, num_threads));
    return rcpp_result_gen;
END_RCPP
}
// EM_func
vec EM_func(const vec& par_start, const vec& betahat, const vec& varbetahat, const vec& ldscore, const vec& Nstar, const int& M, int c0, const long double& eps1, const long double& eps2, const long double& eps3, const long double& eps, const int& Meps, const int& steps, const int& num_threads, const bool& print, const int& printfreq, const bool& stratification);
RcppExport SEXP _GENESIS_EM_func(SEXP par_startSEXP, SEXP betahatSEXP, SEXP varbetahatSEXP, SEXP ldscoreSEXP, SEXP NstarSEXP, SEXP MSEXP, SEXP c0SEXP, SEXP eps1SEXP, SEXP eps2SEXP, SEXP eps3SEXP, SEXP epsSEXP, SEXP MepsSEXP, SEXP stepsSEXP, SEXP num_threadsSEXP, SEXP printSEXP, SEXP printfreqSEXP, SEXP stratificationSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type par_start(par_startSEXP);
    Rcpp::traits::input_parameter< const vec& >::type betahat(betahatSEXP);
    Rcpp::traits::input_parameter< const vec& >::type varbetahat(varbetahatSEXP);
    Rcpp::traits::input_parameter< const vec& >::type ldscore(ldscoreSEXP);
    Rcpp::traits::input_parameter< const vec& >::type Nstar(NstarSEXP);
    Rcpp::traits::input_parameter< const int& >::type M(MSEXP);
    Rcpp::traits::input_parameter< int >::type c0(c0SEXP);
    Rcpp::traits::input_parameter< const long double& >::type eps1(eps1SEXP);
    Rcpp::traits::input_parameter< const long double& >::type eps2(eps2SEXP);
    Rcpp::traits::input_parameter< const long double& >::type eps3(eps3SEXP);
    Rcpp::traits::input_parameter< const long double& >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< const int& >::type Meps(MepsSEXP);
    Rcpp::traits::input_parameter< const int& >::type steps(stepsSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_threads(num_threadsSEXP);
    Rcpp::traits::input_parameter< const bool& >::type print(printSEXP);
    Rcpp::traits::input_parameter< const int& >::type printfreq(printfreqSEXP);
    Rcpp::traits::input_parameter< const bool& >::type stratification(stratificationSEXP);
    rcpp_result_gen = Rcpp::wrap(EM_func(par_start, betahat, varbetahat, ldscore, Nstar, M, c0, eps1, eps2, eps3, eps, Meps, steps, num_threads, print, printfreq, stratification));
    return rcpp_result_gen;
END_RCPP
}
// Sk
vec Sk(const vec& par, const long double& betahatk, const long double& varbetahatk, const long double& ldscorek, const int& c0, const int& Nstark);
RcppExport SEXP _GENESIS_Sk(SEXP parSEXP, SEXP betahatkSEXP, SEXP varbetahatkSEXP, SEXP ldscorekSEXP, SEXP c0SEXP, SEXP NstarkSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type par(parSEXP);
    Rcpp::traits::input_parameter< const long double& >::type betahatk(betahatkSEXP);
    Rcpp::traits::input_parameter< const long double& >::type varbetahatk(varbetahatkSEXP);
    Rcpp::traits::input_parameter< const long double& >::type ldscorek(ldscorekSEXP);
    Rcpp::traits::input_parameter< const int& >::type c0(c0SEXP);
    Rcpp::traits::input_parameter< const int& >::type Nstark(NstarkSEXP);
    rcpp_result_gen = Rcpp::wrap(Sk(par, betahatk, varbetahatk, ldscorek, c0, Nstark));
    return rcpp_result_gen;
END_RCPP
}
// Ik
mat Ik(const vec& par, const long double& betahatk, const long double& varbetahatk, const long double& ldscorek, const int& c0, const int& Nstark);
RcppExport SEXP _GENESIS_Ik(SEXP parSEXP, SEXP betahatkSEXP, SEXP varbetahatkSEXP, SEXP ldscorekSEXP, SEXP c0SEXP, SEXP NstarkSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type par(parSEXP);
    Rcpp::traits::input_parameter< const long double& >::type betahatk(betahatkSEXP);
    Rcpp::traits::input_parameter< const long double& >::type varbetahatk(varbetahatkSEXP);
    Rcpp::traits::input_parameter< const long double& >::type ldscorek(ldscorekSEXP);
    Rcpp::traits::input_parameter< const int& >::type c0(c0SEXP);
    Rcpp::traits::input_parameter< const int& >::type Nstark(NstarkSEXP);
    rcpp_result_gen = Rcpp::wrap(Ik(par, betahatk, varbetahatk, ldscorek, c0, Nstark));
    return rcpp_result_gen;
END_RCPP
}
// S
vec S(const vec& par, const vec& betahat, const vec& varbetahat, const vec& ldscore, const int& c0, const vec& Nstar, const int& num_threads);
RcppExport SEXP _GENESIS_S(SEXP parSEXP, SEXP betahatSEXP, SEXP varbetahatSEXP, SEXP ldscoreSEXP, SEXP c0SEXP, SEXP NstarSEXP, SEXP num_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type par(parSEXP);
    Rcpp::traits::input_parameter< const vec& >::type betahat(betahatSEXP);
    Rcpp::traits::input_parameter< const vec& >::type varbetahat(varbetahatSEXP);
    Rcpp::traits::input_parameter< const vec& >::type ldscore(ldscoreSEXP);
    Rcpp::traits::input_parameter< const int& >::type c0(c0SEXP);
    Rcpp::traits::input_parameter< const vec& >::type Nstar(NstarSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_threads(num_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(S(par, betahat, varbetahat, ldscore, c0, Nstar, num_threads));
    return rcpp_result_gen;
END_RCPP
}
// SS
mat SS(const vec& par, const vec& betahat, const vec& varbetahat, const vec& ldscore, const int& c0, const vec& Nstar, const int& num_threads);
RcppExport SEXP _GENESIS_SS(SEXP parSEXP, SEXP betahatSEXP, SEXP varbetahatSEXP, SEXP ldscoreSEXP, SEXP c0SEXP, SEXP NstarSEXP, SEXP num_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type par(parSEXP);
    Rcpp::traits::input_parameter< const vec& >::type betahat(betahatSEXP);
    Rcpp::traits::input_parameter< const vec& >::type varbetahat(varbetahatSEXP);
    Rcpp::traits::input_parameter< const vec& >::type ldscore(ldscoreSEXP);
    Rcpp::traits::input_parameter< const int& >::type c0(c0SEXP);
    Rcpp::traits::input_parameter< const vec& >::type Nstar(NstarSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_threads(num_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(SS(par, betahat, varbetahat, ldscore, c0, Nstar, num_threads));
    return rcpp_result_gen;
END_RCPP
}
// I
mat I(const vec& par, const vec& betahat, const vec& varbetahat, const vec& ldscore, const int& c0, const vec& Nstar, const int& num_threads);
RcppExport SEXP _GENESIS_I(SEXP parSEXP, SEXP betahatSEXP, SEXP varbetahatSEXP, SEXP ldscoreSEXP, SEXP c0SEXP, SEXP NstarSEXP, SEXP num_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type par(parSEXP);
    Rcpp::traits::input_parameter< const vec& >::type betahat(betahatSEXP);
    Rcpp::traits::input_parameter< const vec& >::type varbetahat(varbetahatSEXP);
    Rcpp::traits::input_parameter< const vec& >::type ldscore(ldscoreSEXP);
    Rcpp::traits::input_parameter< const int& >::type c0(c0SEXP);
    Rcpp::traits::input_parameter< const vec& >::type Nstar(NstarSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_threads(num_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(I(par, betahat, varbetahat, ldscore, c0, Nstar, num_threads));
    return rcpp_result_gen;
END_RCPP
}
// mixture_components_marginal
List mixture_components_marginal(const vec& par, const vec& ldscore, const int& c0, const vec& Nstar, const int& num_threads);
RcppExport SEXP _GENESIS_mixture_components_marginal(SEXP parSEXP, SEXP ldscoreSEXP, SEXP c0SEXP, SEXP NstarSEXP, SEXP num_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type par(parSEXP);
    Rcpp::traits::input_parameter< const vec& >::type ldscore(ldscoreSEXP);
    Rcpp::traits::input_parameter< const int& >::type c0(c0SEXP);
    Rcpp::traits::input_parameter< const vec& >::type Nstar(NstarSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_threads(num_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(mixture_components_marginal(par, ldscore, c0, Nstar, num_threads));
    return rcpp_result_gen;
END_RCPP
}
// loglikelihood3
long double loglikelihood3(const vec& par, const vec& betahat, const vec& varbetahat, const vec& ldscore, const int& c0, const vec& Nstar, const int& num_threads);
RcppExport SEXP _GENESIS_loglikelihood3(SEXP parSEXP, SEXP betahatSEXP, SEXP varbetahatSEXP, SEXP ldscoreSEXP, SEXP c0SEXP, SEXP NstarSEXP, SEXP num_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type par(parSEXP);
    Rcpp::traits::input_parameter< const vec& >::type betahat(betahatSEXP);
    Rcpp::traits::input_parameter< const vec& >::type varbetahat(varbetahatSEXP);
    Rcpp::traits::input_parameter< const vec& >::type ldscore(ldscoreSEXP);
    Rcpp::traits::input_parameter< const int& >::type c0(c0SEXP);
    Rcpp::traits::input_parameter< const vec& >::type Nstar(NstarSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_threads(num_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(loglikelihood3(par, betahat, varbetahat, ldscore, c0, Nstar, num_threads));
    return rcpp_result_gen;
END_RCPP
}
// weight3
mat weight3(const vec& par, const vec& betahat, const vec& varbetahat, const vec& ldscore, const int& c0, const vec& Nstar, const int& num_threads);
RcppExport SEXP _GENESIS_weight3(SEXP parSEXP, SEXP betahatSEXP, SEXP varbetahatSEXP, SEXP ldscoreSEXP, SEXP c0SEXP, SEXP NstarSEXP, SEXP num_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type par(parSEXP);
    Rcpp::traits::input_parameter< const vec& >::type betahat(betahatSEXP);
    Rcpp::traits::input_parameter< const vec& >::type varbetahat(varbetahatSEXP);
    Rcpp::traits::input_parameter< const vec& >::type ldscore(ldscoreSEXP);
    Rcpp::traits::input_parameter< const int& >::type c0(c0SEXP);
    Rcpp::traits::input_parameter< const vec& >::type Nstar(NstarSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_threads(num_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(weight3(par, betahat, varbetahat, ldscore, c0, Nstar, num_threads));
    return rcpp_result_gen;
END_RCPP
}
// weight_loglikelihood3
List weight_loglikelihood3(const vec& par, const vec& betahat, const vec& varbetahat, const vec& ldscore, const int& c0, const vec& Nstar, const int& num_threads);
RcppExport SEXP _GENESIS_weight_loglikelihood3(SEXP parSEXP, SEXP betahatSEXP, SEXP varbetahatSEXP, SEXP ldscoreSEXP, SEXP c0SEXP, SEXP NstarSEXP, SEXP num_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type par(parSEXP);
    Rcpp::traits::input_parameter< const vec& >::type betahat(betahatSEXP);
    Rcpp::traits::input_parameter< const vec& >::type varbetahat(varbetahatSEXP);
    Rcpp::traits::input_parameter< const vec& >::type ldscore(ldscoreSEXP);
    Rcpp::traits::input_parameter< const int& >::type c0(c0SEXP);
    Rcpp::traits::input_parameter< const vec& >::type Nstar(NstarSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_threads(num_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(weight_loglikelihood3(par, betahat, varbetahat, ldscore, c0, Nstar, num_threads));
    return rcpp_result_gen;
END_RCPP
}
// update_p3
vec update_p3(const mat& w, const vec& Nstar, const int& num_threads);
RcppExport SEXP _GENESIS_update_p3(SEXP wSEXP, SEXP NstarSEXP, SEXP num_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const mat& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const vec& >::type Nstar(NstarSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_threads(num_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(update_p3(w, Nstar, num_threads));
    return rcpp_result_gen;
END_RCPP
}
// onestep_varcomponent3
vec onestep_varcomponent3(const vec& varcomponent, const mat& w, const vec& betahat, const vec& varbetahat, const vec& ldscore, const vec& Nstar, const int& num_threads);
RcppExport SEXP _GENESIS_onestep_varcomponent3(SEXP varcomponentSEXP, SEXP wSEXP, SEXP betahatSEXP, SEXP varbetahatSEXP, SEXP ldscoreSEXP, SEXP NstarSEXP, SEXP num_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type varcomponent(varcomponentSEXP);
    Rcpp::traits::input_parameter< const mat& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const vec& >::type betahat(betahatSEXP);
    Rcpp::traits::input_parameter< const vec& >::type varbetahat(varbetahatSEXP);
    Rcpp::traits::input_parameter< const vec& >::type ldscore(ldscoreSEXP);
    Rcpp::traits::input_parameter< const vec& >::type Nstar(NstarSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_threads(num_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(onestep_varcomponent3(varcomponent, w, betahat, varbetahat, ldscore, Nstar, num_threads));
    return rcpp_result_gen;
END_RCPP
}
// EM_func3
vec EM_func3(const vec& par_start, const vec& lower_pi, const vec& upper_pi, const vec& betahat, const vec& varbetahat, const vec& ldscore, const vec& Nstar, const int& M, int c0, const long double& eps1, const long double& eps2, const long double& eps3, const long double& eps4, const long double& eps5, const long double& eps, const int& Meps, const int& steps, const int& num_threads, const bool& print, const int& printfreq, const bool& stratification);
RcppExport SEXP _GENESIS_EM_func3(SEXP par_startSEXP, SEXP lower_piSEXP, SEXP upper_piSEXP, SEXP betahatSEXP, SEXP varbetahatSEXP, SEXP ldscoreSEXP, SEXP NstarSEXP, SEXP MSEXP, SEXP c0SEXP, SEXP eps1SEXP, SEXP eps2SEXP, SEXP eps3SEXP, SEXP eps4SEXP, SEXP eps5SEXP, SEXP epsSEXP, SEXP MepsSEXP, SEXP stepsSEXP, SEXP num_threadsSEXP, SEXP printSEXP, SEXP printfreqSEXP, SEXP stratificationSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type par_start(par_startSEXP);
    Rcpp::traits::input_parameter< const vec& >::type lower_pi(lower_piSEXP);
    Rcpp::traits::input_parameter< const vec& >::type upper_pi(upper_piSEXP);
    Rcpp::traits::input_parameter< const vec& >::type betahat(betahatSEXP);
    Rcpp::traits::input_parameter< const vec& >::type varbetahat(varbetahatSEXP);
    Rcpp::traits::input_parameter< const vec& >::type ldscore(ldscoreSEXP);
    Rcpp::traits::input_parameter< const vec& >::type Nstar(NstarSEXP);
    Rcpp::traits::input_parameter< const int& >::type M(MSEXP);
    Rcpp::traits::input_parameter< int >::type c0(c0SEXP);
    Rcpp::traits::input_parameter< const long double& >::type eps1(eps1SEXP);
    Rcpp::traits::input_parameter< const long double& >::type eps2(eps2SEXP);
    Rcpp::traits::input_parameter< const long double& >::type eps3(eps3SEXP);
    Rcpp::traits::input_parameter< const long double& >::type eps4(eps4SEXP);
    Rcpp::traits::input_parameter< const long double& >::type eps5(eps5SEXP);
    Rcpp::traits::input_parameter< const long double& >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< const int& >::type Meps(MepsSEXP);
    Rcpp::traits::input_parameter< const int& >::type steps(stepsSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_threads(num_threadsSEXP);
    Rcpp::traits::input_parameter< const bool& >::type print(printSEXP);
    Rcpp::traits::input_parameter< const int& >::type printfreq(printfreqSEXP);
    Rcpp::traits::input_parameter< const bool& >::type stratification(stratificationSEXP);
    rcpp_result_gen = Rcpp::wrap(EM_func3(par_start, lower_pi, upper_pi, betahat, varbetahat, ldscore, Nstar, M, c0, eps1, eps2, eps3, eps4, eps5, eps, Meps, steps, num_threads, print, printfreq, stratification));
    return rcpp_result_gen;
END_RCPP
}
// Sk3
vec Sk3(const vec& par, const long double& betahatk, const long double& varbetahatk, const long double& ldscorek, const int& c0, const int& Nstark);
RcppExport SEXP _GENESIS_Sk3(SEXP parSEXP, SEXP betahatkSEXP, SEXP varbetahatkSEXP, SEXP ldscorekSEXP, SEXP c0SEXP, SEXP NstarkSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type par(parSEXP);
    Rcpp::traits::input_parameter< const long double& >::type betahatk(betahatkSEXP);
    Rcpp::traits::input_parameter< const long double& >::type varbetahatk(varbetahatkSEXP);
    Rcpp::traits::input_parameter< const long double& >::type ldscorek(ldscorekSEXP);
    Rcpp::traits::input_parameter< const int& >::type c0(c0SEXP);
    Rcpp::traits::input_parameter< const int& >::type Nstark(NstarkSEXP);
    rcpp_result_gen = Rcpp::wrap(Sk3(par, betahatk, varbetahatk, ldscorek, c0, Nstark));
    return rcpp_result_gen;
END_RCPP
}
// Ik3
mat Ik3(const vec& par, const long double& betahatk, const long double& varbetahatk, const long double& ldscorek, const int& c0, const int& Nstark);
RcppExport SEXP _GENESIS_Ik3(SEXP parSEXP, SEXP betahatkSEXP, SEXP varbetahatkSEXP, SEXP ldscorekSEXP, SEXP c0SEXP, SEXP NstarkSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type par(parSEXP);
    Rcpp::traits::input_parameter< const long double& >::type betahatk(betahatkSEXP);
    Rcpp::traits::input_parameter< const long double& >::type varbetahatk(varbetahatkSEXP);
    Rcpp::traits::input_parameter< const long double& >::type ldscorek(ldscorekSEXP);
    Rcpp::traits::input_parameter< const int& >::type c0(c0SEXP);
    Rcpp::traits::input_parameter< const int& >::type Nstark(NstarkSEXP);
    rcpp_result_gen = Rcpp::wrap(Ik3(par, betahatk, varbetahatk, ldscorek, c0, Nstark));
    return rcpp_result_gen;
END_RCPP
}
// S3
vec S3(const vec& par, const vec& betahat, const vec& varbetahat, const vec& ldscore, const int& c0, const vec& Nstar);
RcppExport SEXP _GENESIS_S3(SEXP parSEXP, SEXP betahatSEXP, SEXP varbetahatSEXP, SEXP ldscoreSEXP, SEXP c0SEXP, SEXP NstarSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type par(parSEXP);
    Rcpp::traits::input_parameter< const vec& >::type betahat(betahatSEXP);
    Rcpp::traits::input_parameter< const vec& >::type varbetahat(varbetahatSEXP);
    Rcpp::traits::input_parameter< const vec& >::type ldscore(ldscoreSEXP);
    Rcpp::traits::input_parameter< const int& >::type c0(c0SEXP);
    Rcpp::traits::input_parameter< const vec& >::type Nstar(NstarSEXP);
    rcpp_result_gen = Rcpp::wrap(S3(par, betahat, varbetahat, ldscore, c0, Nstar));
    return rcpp_result_gen;
END_RCPP
}
// SS3
mat SS3(const vec& par, const vec& betahat, const vec& varbetahat, const vec& ldscore, const int& c0, const vec& Nstar, const int& num_threads);
RcppExport SEXP _GENESIS_SS3(SEXP parSEXP, SEXP betahatSEXP, SEXP varbetahatSEXP, SEXP ldscoreSEXP, SEXP c0SEXP, SEXP NstarSEXP, SEXP num_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type par(parSEXP);
    Rcpp::traits::input_parameter< const vec& >::type betahat(betahatSEXP);
    Rcpp::traits::input_parameter< const vec& >::type varbetahat(varbetahatSEXP);
    Rcpp::traits::input_parameter< const vec& >::type ldscore(ldscoreSEXP);
    Rcpp::traits::input_parameter< const int& >::type c0(c0SEXP);
    Rcpp::traits::input_parameter< const vec& >::type Nstar(NstarSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_threads(num_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(SS3(par, betahat, varbetahat, ldscore, c0, Nstar, num_threads));
    return rcpp_result_gen;
END_RCPP
}
// I3
mat I3(const vec& par, const vec& betahat, const vec& varbetahat, const vec& ldscore, const int& c0, const vec& Nstar, const int& num_threads);
RcppExport SEXP _GENESIS_I3(SEXP parSEXP, SEXP betahatSEXP, SEXP varbetahatSEXP, SEXP ldscoreSEXP, SEXP c0SEXP, SEXP NstarSEXP, SEXP num_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type par(parSEXP);
    Rcpp::traits::input_parameter< const vec& >::type betahat(betahatSEXP);
    Rcpp::traits::input_parameter< const vec& >::type varbetahat(varbetahatSEXP);
    Rcpp::traits::input_parameter< const vec& >::type ldscore(ldscoreSEXP);
    Rcpp::traits::input_parameter< const int& >::type c0(c0SEXP);
    Rcpp::traits::input_parameter< const vec& >::type Nstar(NstarSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_threads(num_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(I3(par, betahat, varbetahat, ldscore, c0, Nstar, num_threads));
    return rcpp_result_gen;
END_RCPP
}
// mixture_3components_marginal
List mixture_3components_marginal(const vec& par, const vec& ldscore, const int& c0, const vec& Nstar, const int& num_threads);
RcppExport SEXP _GENESIS_mixture_3components_marginal(SEXP parSEXP, SEXP ldscoreSEXP, SEXP c0SEXP, SEXP NstarSEXP, SEXP num_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type par(parSEXP);
    Rcpp::traits::input_parameter< const vec& >::type ldscore(ldscoreSEXP);
    Rcpp::traits::input_parameter< const int& >::type c0(c0SEXP);
    Rcpp::traits::input_parameter< const vec& >::type Nstar(NstarSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_threads(num_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(mixture_3components_marginal(par, ldscore, c0, Nstar, num_threads));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_GENESIS_modification_loc", (DL_FUNC) &_GENESIS_modification_loc, 3},
    {"_GENESIS_cpnorm", (DL_FUNC) &_GENESIS_cpnorm, 1},
    {"_GENESIS_sumfactorial", (DL_FUNC) &_GENESIS_sumfactorial, 1},
    {"_GENESIS_sumfactorial_rev", (DL_FUNC) &_GENESIS_sumfactorial_rev, 2},
    {"_GENESIS_cdnorm", (DL_FUNC) &_GENESIS_cdnorm, 4},
    {"_GENESIS_cdbinom", (DL_FUNC) &_GENESIS_cdbinom, 4},
    {"_GENESIS_cdmultinom3", (DL_FUNC) &_GENESIS_cdmultinom3, 5},
    {"_GENESIS_loglikelihood", (DL_FUNC) &_GENESIS_loglikelihood, 7},
    {"_GENESIS_weight", (DL_FUNC) &_GENESIS_weight, 7},
    {"_GENESIS_weight_loglikelihood", (DL_FUNC) &_GENESIS_weight_loglikelihood, 7},
    {"_GENESIS_update_pi1", (DL_FUNC) &_GENESIS_update_pi1, 3},
    {"_GENESIS_onestep_varcomponent", (DL_FUNC) &_GENESIS_onestep_varcomponent, 7},
    {"_GENESIS_EM_func", (DL_FUNC) &_GENESIS_EM_func, 17},
    {"_GENESIS_Sk", (DL_FUNC) &_GENESIS_Sk, 6},
    {"_GENESIS_Ik", (DL_FUNC) &_GENESIS_Ik, 6},
    {"_GENESIS_S", (DL_FUNC) &_GENESIS_S, 7},
    {"_GENESIS_SS", (DL_FUNC) &_GENESIS_SS, 7},
    {"_GENESIS_I", (DL_FUNC) &_GENESIS_I, 7},
    {"_GENESIS_mixture_components_marginal", (DL_FUNC) &_GENESIS_mixture_components_marginal, 5},
    {"_GENESIS_loglikelihood3", (DL_FUNC) &_GENESIS_loglikelihood3, 7},
    {"_GENESIS_weight3", (DL_FUNC) &_GENESIS_weight3, 7},
    {"_GENESIS_weight_loglikelihood3", (DL_FUNC) &_GENESIS_weight_loglikelihood3, 7},
    {"_GENESIS_update_p3", (DL_FUNC) &_GENESIS_update_p3, 3},
    {"_GENESIS_onestep_varcomponent3", (DL_FUNC) &_GENESIS_onestep_varcomponent3, 7},
    {"_GENESIS_EM_func3", (DL_FUNC) &_GENESIS_EM_func3, 21},
    {"_GENESIS_Sk3", (DL_FUNC) &_GENESIS_Sk3, 6},
    {"_GENESIS_Ik3", (DL_FUNC) &_GENESIS_Ik3, 6},
    {"_GENESIS_S3", (DL_FUNC) &_GENESIS_S3, 6},
    {"_GENESIS_SS3", (DL_FUNC) &_GENESIS_SS3, 7},
    {"_GENESIS_I3", (DL_FUNC) &_GENESIS_I3, 7},
    {"_GENESIS_mixture_3components_marginal", (DL_FUNC) &_GENESIS_mixture_3components_marginal, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_GENESIS(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
