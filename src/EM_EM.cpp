//------------------------------------------------
//  EM_EM.cpp
//  Goal: Estimate effect size distribution.
//
//  Author: Yan (Dora) Zhang
//  Email: yandorazhang@gmail.com
//------------------------------------------------

#include <cstdlib>
#include <time.h>
#include <iostream>     // std::cout
#include <algorithm>    // std::min
#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]

// --------------------------------//--------------------------------
#include <omp.h>
//[[Rcpp::plugins(openmp)]]
// --------------------------------//--------------------------------

using namespace arma;
using namespace Rcpp;
using namespace std;
using namespace stats;


//--------------------------------
//--------------------------------
// [[Rcpp::export]]
vec modification_loc(vec inx_name, int K, int mx_k){
  vec inx_loc(mx_k);
  inx_loc.fill(0);
  
  for(int i = 0; i < K; i++){
    inx_loc(inx_name(i) - 1) = i+1;
  }
  return(inx_loc);
}



//--------------------------------
//--------------------------------
// [[Rcpp::export]]
long double cpnorm(long double x) // R: pnorm()
{
  // constants
  long double a1 =  0.254829592;
  long double a2 = -0.284496736;
  long double a3 =  1.421413741;
  long double a4 = -1.453152027;
  long double a5 =  1.061405429;
  long double p  =  0.3275911;
  
  // Save the sign of x
  int sign = 1;
  if (x < 0)
    sign = -1;
  x = fabs(x)/sqrt(2.0);
  
  // A&S formula 7.1.26
  long double t = 1.0/(1.0 + p*x);
  long double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);
  
  return 0.5*(1.0 + sign*y);
}

//--------------------------------
//--------------------------------
// [[Rcpp::export]]
long double sumfactorial(int n) // log(factorial(n)) = log(n!)
{
  if(n > 1)
    return log(n) + sumfactorial(n - 1);
  else
    return 0;
}

//--------------------------------
//--------------------------------
// [[Rcpp::export]]
long double sumfactorial_rev(int n, const int &k0) // log(n!/k0!)
{
  if(n > k0)
    return log(n) + sumfactorial_rev(n-1,k0);
  else
    return 0;
}

//--------------------------------
//--------------------------------
// [[Rcpp::export]]
long double cdnorm(long double x, long double mean, long double sd, bool loglog) // R: dnorm()
{
  long double res;
  res = -0.5*log((2.0*PI))  - log(sd)- pow(x-mean, 2)/(2.0*pow(sd,2));
  if(loglog){return res;}
  else{return exp(res);}
}

//--------------------------------
//--------------------------------
// [[Rcpp::export]]
long double cdbinom(const int &k, const int &size, long double prob, bool loglog) // R: dbinom() 
{
  long double res;
  res = sumfactorial_rev(size,size-k) - sumfactorial(k) + k*log(prob) + (size-k)*log(1-prob);
  if(loglog){return res;}
  else{return exp(res);}
}

//--------------------------------
//--------------------------------
// [[Rcpp::export]]
long double cdmultinom3(const int k0, const int k1, const int k2,  vec prob, bool loglog) // R: dbmultinom()
{
    long double res;
    
    res = sumfactorial_rev(k0+k1+k2,k2) - sumfactorial(k0) - sumfactorial(k1)
    + k0*log(prob(0)) + k1*log(prob(1)) + k2*log(prob(2));
    
    if(loglog){return res;}
    else{return exp(res);}
}


//--------------------------------//--------------------------------
// 2-component model
//--------------------------------//--------------------------------

//--------------------------------
//--------------------------------
// [[Rcpp::export]]
long double loglikelihood(const vec & par, const vec &betahat, const vec &varbetahat, const vec &ldscore,const int & c0, const vec &Nstar, const int & num_threads) // loglikelihood function
{
  int K = betahat.n_elem;
  long double marginal_likelihood;
  long double pi1 = par(0);
  long double sigsq = par(1);
  long double a = par(2);
  long double res = 0;
  long double y=0;
  long double tem;
  int k,j;
  // --------------------------------//--------------------------------
  omp_set_num_threads(num_threads);
  #pragma omp parallel for shared(betahat,pi1, sigsq, a, varbetahat,ldscore,c0,Nstar) private(marginal_likelihood,y,k,j,tem) reduction(+:res)
  // --------------------------------//--------------------------------
  for(k=0; k<K; k++){
    marginal_likelihood = 0;
    y = 0;
    for(j=0; j<std::min(int(Nstar(k))+1, c0+1); j++){
      tem = cdbinom(j,Nstar(k),pi1,false);
      y += tem;
      marginal_likelihood += tem*cdnorm(betahat(k),0, sqrt( (j*sigsq)*ldscore(k)/Nstar(k) + varbetahat(k) + a), false);
    }
    res += log(marginal_likelihood/y);
  }
  
  return res;
}

//--------------------------------
// weights function
//--------------------------------
// [[Rcpp::export]]
mat weight(const vec & par, const vec &betahat, const vec &varbetahat, const vec &ldscore, const int & c0, const vec & Nstar, const int &num_threads) // E-step: weights function
{
  int K = betahat.n_elem;
  mat w(K,c0+1);
  w.fill(0.0);
  long double pi1 = par(0);
  long double sigsq = par(1);
  long double a  = par(2);
  vec te;
  long double tem;
  int k,j;
  
  // weights formula.
  // --------------------------------//--------------------------------
  omp_set_num_threads(num_threads);
  #pragma omp parallel for shared(sigsq,a, pi1,betahat,varbetahat,ldscore,Nstar,w,c0) private(k,j,tem)
  // --------------------------------//--------------------------------
  for( k=0; k<K; k++){
    for(j=0; j<std::min(int(Nstar(k))+1, c0+1); j++){
      tem = cdbinom(j,Nstar(k),pi1,false);
      w(k,j) = tem * cdnorm(betahat(k),0, sqrt( (j*sigsq)*ldscore(k)/Nstar(k) + varbetahat(k) + a), false);
    }
  }
  te = sum(w,1);
  te = pow(te, -1.0);
  w = (diagmat(te)) * w;
  
  return w;
}


//--------------------------------
// weights and loglikelihood function
//--------------------------------
// [[Rcpp::export]]
List weight_loglikelihood(const vec & par, const vec &betahat, const vec &varbetahat, const vec &ldscore, const int & c0, const vec & Nstar, const int &num_threads) // E-step: weight-loglikelihood function
{
  int K = betahat.n_elem;
  mat w(K,c0+1);
  w.fill(0.0);
  long double pi1 = par(0);
  long double sigsq = par(1);
  long double a  = par(2);
  vec te;
  long double tem, tem1;
  int k,j;
  long double marginal_likelihood;
  long double res = 0;
  long double y=0;

  // --------------------------------//--------------------------------
  omp_set_num_threads(num_threads);
  #pragma omp parallel for shared(sigsq,a, pi1,betahat,varbetahat,ldscore,Nstar,w,c0) private(marginal_likelihood,y,k,j,tem,tem1) reduction(+:res)
  // --------------------------------//--------------------------------
  for(k=0; k<K; k++){
    marginal_likelihood = 0;
    y = 0;
    for(j=0; j<std::min(int(Nstar(k))+1, c0+1); j++){
      tem = cdbinom(j,Nstar(k),pi1,false);
      tem1 = tem * cdnorm(betahat(k),0, sqrt( (j*sigsq)*ldscore(k)/Nstar(k) + varbetahat(k) + a), false);
      w(k,j) = tem1;
      y += tem;
      marginal_likelihood += tem1;
    }
    res += log(marginal_likelihood/y);
  }
  te = sum(w,1);
  te = pow(te, -1.0);
  w = (diagmat(te)) * w;

  return List::create(
    _["w"] = w,
    _["llk"] = res
  );
}



//--------------------------------
//--------------------------------
// [[Rcpp::export]]
long double update_pi1(const mat & w, const vec & Nstar, const int &num_threads)
{
  long double res1=0;
  long double res2=0;
  long double res;
  int c0 = w.n_cols-1;
  int K = Nstar.n_elem;
  int k,j;
  // --------------------------------//--------------------------------
  omp_set_num_threads(num_threads);
  #pragma omp parallel for shared(w,Nstar,c0) private(k,j) reduction(+:res1,res2)
  // --------------------------------//--------------------------------
  for(k=0; k<K; k++){
    for(j=0; j<std::min(int(Nstar(k))+1, c0+1); j++){
      res1 += w(k,j)*j;
      res2 += w(k,j)*Nstar(k);
    }
  }
  
  res = res1/res2;
  
  return res;
}

//--------------------------------
//--------------------------------
// [[Rcpp::export]]
vec onestep_varcomponent(const vec varcomponent,  const mat & w, const vec &betahat, const vec &varbetahat, const vec &ldscore, const vec & Nstar, const int &num_threads)
{
  
  int K = w.n_rows;
  int c0 = w.n_cols-1;
  long double sigsq = varcomponent(0);
  long double a = varcomponent(1);
  long double tem;
  long double s1=0, dd1=0, s2=0, dd2=0, dd12=0, det=0;
  vec result(2);
  
  int k,j;
  // --------------------------------//--------------------------------
  omp_set_num_threads(num_threads);
  #pragma omp parallel for shared(w, sigsq,a,varbetahat,ldscore,Nstar,c0) private(k,tem,j) reduction(+:s1,dd1,s2,dd2,dd12)
  // --------------------------------//--------------------------------
  for(k=0; k<K; k++)
  {
    for(j=0; j<std::min(int(Nstar(k))+1, c0+1); j++)
    {
      tem = (j*sigsq) * ldscore(k) /double(Nstar(k)) + varbetahat(k) + a;
      s1 += w(k,j)*j*ldscore(k)/(2.0*Nstar(k))*(pow(betahat(k),2)/pow(tem,2) - 1.0/(tem)  );
      s2 += w(k,j)/(2.0)*(pow(betahat(k),2)/pow(tem,2) - 1.0/(tem));
      
      dd1 += w(k, j)*j*ldscore(k)*j*ldscore(k)/(2.0*Nstar(k)*Nstar(k))*(-2.0*pow(betahat(k),2)/pow(tem,3) + 1.0/(pow(tem,2)));
      dd2 += w(k,j)/(2.0)*(-2.0*pow(betahat(k),2)/pow(tem,3) + 1.0/(pow(tem,2)));
      dd12 += w(k,j)*j*ldscore(k)/(2.0*Nstar(k))*(-2.0*pow(betahat(k),2)/pow(tem,3) + 1.0/(pow(tem,2)));
    }
  }
  
  det = (dd1*dd2 - pow(dd12,2));
  
  result(0) = sigsq - (dd2*s1 - dd12*s2)/det;
  result(1) = a - (dd1*s2 - dd12*s1)/det;
  
  return result;
}

//--------------------------------
//--------------------------------
// [[Rcpp::export]]
vec EM_func(const vec &par_start,
            const vec &betahat, const vec & varbetahat, const vec & ldscore,  const vec & Nstar, const int & M,
            int c0, const long double &eps1,const long double &eps2,const long double &eps3,const long double &eps, const int &Meps, const int &steps, const int &num_threads, const bool &print,
            const int &printfreq, const bool &stratification)
{
  vec par(3), prev_par(3);
  par(0) = par_start(0);
  par(1) = par_start(1);
  par(2) = par_start(2);
  long double llk = 0, error_pi, error_sigsq, error_a, increase_ll;
  long double prev_llk=0, pi1;
  vec old_sig(2);
  vec new_sig(2);
  mat w;
  vec result(10);
  List wllk;
  clock_t start_w, finish_w, finish_pi1, finish_sigsq;

  Rcout << "Iteration, prev_loglikelihood, pic, sigmasq, a, Heritability, c0, Seconds(weight_llk), Seconds(pi), Seconds(variance_components)";
  Rcout << endl;

  for(int i=0; i<Meps; i++){
    prev_llk = llk;
    prev_par(0) = par(0); prev_par(1) = par(1); prev_par(2) = par(2);

    // E-step: calculate the weights and log-likelihood
    start_w = clock();
    wllk = weight_loglikelihood(par, betahat, varbetahat, ldscore, c0, Nstar, num_threads);
    w = as<mat>(wllk["w"]);
    llk = as<long double>(wllk["llk"]);
    finish_w = clock();

    // M-step: update the proportion parameters
    pi1 = update_pi1(w,Nstar,num_threads);
    finish_pi1 = clock();

    // M-step: update the variance parameters
    new_sig(0) = par(1); new_sig(1) = par(2);
    if(stratification==false){new_sig(1)=0;}
    for(int j=0; j<steps; j++){
      old_sig = onestep_varcomponent(new_sig, w, betahat,varbetahat, ldscore, Nstar, num_threads);
      if(stratification==false){old_sig(1)=0;}
      if((abs(old_sig(0)-new_sig(0))<1e-20) & (abs(old_sig(1)-new_sig(1))<1e-20)) break;
      new_sig(0) = old_sig(0); new_sig(1) = old_sig(1);
      if(new_sig(0) > 1) new_sig(0) = 1e-5;
      if(new_sig(0) < 0) new_sig(0) = 1e-12;
      if(new_sig(1) > 1) new_sig(1) = 1e-6;
      if(new_sig(1) < -min(varbetahat)) new_sig(1) = -min(varbetahat)/2;
    }
    finish_sigsq = clock();

    // update par vector with the new parameter values
    par(0) = pi1; par(1) = new_sig(0); par(2) = new_sig(1);

    if(isinf(-llk) | isnan(-llk) | isnan(llk)) {par(2) = abs(par(2));}
    if((llk<prev_llk) & (c0<20)) {c0 = c0+1;}
    if((llk<prev_llk) & (c0>=20)) {break;}

    // output results into a vector result
    result(0) = i; result(1) = llk;
    result(2) = par(0); result(3) = par(1); result(4) = par(2);
    result(5) = par(0)*par(1)*M; result(6) = c0;
    result(7) = double(finish_w-start_w)/CLOCKS_PER_SEC;
    result(8) = double(finish_pi1-finish_w)/CLOCKS_PER_SEC;
    result(9) = double(finish_sigsq-finish_pi1)/CLOCKS_PER_SEC;

    if(print==true){
      if(i%printfreq==0){
        for(int r=0; r<10; r++){
          Rcout << result(r)<< ", " ;
        }
        Rcout << endl;
      }
    }

    error_pi = abs(prev_par(0) - par(0)) ;
    error_sigsq = abs(prev_par(1) - par(1)) ;
    error_a = abs(prev_par(2) - par(2)) ;
    increase_ll = (llk - prev_llk)/prev_llk ;

    if(((error_pi< eps1) & (error_sigsq <eps2) & (error_a <eps3)) & (abs(increase_ll) <eps) ) {break;}
  }
  // calcualte the logliklihood under the new par
  llk = loglikelihood(par,betahat,varbetahat,ldscore,c0,Nstar, num_threads);
  result(1) = llk;
  return result;
}

//--------------------------------
//--------------------------------
// [[Rcpp::export]]
vec Sk(const vec & par, const long double &betahatk, const long double &varbetahatk, const long double &ldscorek,const int & c0, const int &Nstark) // score vector for kth SNP 
{
  
  long double pic = par(0);
  long double sigsq = par(1);
  long double a = par(2);
  long double  tem_var, ww, numerator_pic=0, denominator=0, numerator_sigsq=0, numerator_a=0;
  vec res(3);
  int j;
  
  for(j=0; j<std::min(int(Nstark)+1, c0+1); j++){
    tem_var = (j*sigsq)*ldscorek/Nstark + varbetahatk + a;
    ww = cdbinom(j,Nstark,pic,false) * cdnorm(betahatk, 0, sqrt(tem_var), false);
    denominator += ww;
    
    //pic
    numerator_pic += ww*(j/pic - (Nstark-j)/(1-pic));
    
    //sigsq
    numerator_sigsq += ww*(-0.5/tem_var + 0.5*pow(betahatk,2.0)/pow(tem_var,2.0))*j*ldscorek/Nstark ;
    
    //a
    numerator_a += ww*(-0.5/tem_var + 0.5*pow(betahatk,2.0)/pow(tem_var,2.0));
    
  }
  
  res(0) = numerator_pic/denominator;
  res(1) = numerator_sigsq/denominator;
  res(2) = numerator_a/denominator;
  
  return res;
}

//--------------------------------
//--------------------------------
// [[Rcpp::export]]
mat Ik(const vec & par, const long double &betahatk, const long double &varbetahatk, const long double &ldscorek,const int & c0, const int &Nstark) // Information for kth SNP, i.e., -partial(Sk)/partial(theta)
{
  
  long double pic = par(0);
  long double sigsq = par(1);
  long double a = par(2);
  
  long double  tem_var, g, Lk=0;
  vec r(3),  dLk(3);
  dLk.fill(0.0);
  mat dr(3,3), rr(3,3), ddLk(3,3), dLk2(3,3), res(3,3);
  ddLk.fill(0.0);
  dLk2.fill(0.0);
  int j;
  
  for(j=0; j<std::min(int(Nstark)+1, c0+1); j++){
    
    tem_var = (j*sigsq)*ldscorek/Nstark + varbetahatk + a;
    
    g = cdbinom(j,Nstark,pic,false) * cdnorm(betahatk, 0, sqrt(tem_var), false);
    r(0) =  (j/pic - (Nstark-j)/(1-pic));
    r(1)= (-0.5/tem_var + 0.5*pow(betahatk,2.0)/pow(tem_var,2.0))*j*ldscorek/Nstark ;
    r(2) = (-0.5/tem_var + 0.5*pow(betahatk,2.0)/pow(tem_var,2.0));
    
    Lk += g;
    
    dr(0,0) = -j/pow(pic,2.0) - (Nstark-j)/(pow(1-pic, 2.0));
    dr(0,1) = 0;
    dr(0,2) = 0;
    dr(1,0) = 0;
    dr(1,1) = (0.5/pow(tem_var,2.0) - pow(betahatk,2.0)/pow(tem_var,3.0))*pow(j*ldscorek/Nstark,2.0) ;
    dr(1,2) = (0.5/pow(tem_var,2.0) - pow(betahatk,2.0)/pow(tem_var,3.0))*(j*ldscorek/Nstark) ;
    dr(2,0) = 0;
    dr(2,1) = dr(1,2);
    dr(2,2) = (0.5/pow(tem_var,2.0) - pow(betahatk,2.0)/pow(tem_var,3.0));
    
    rr(0,0) = r(0)*r(0);
    rr(0,1) = r(0)*r(1);
    rr(0,2) = r(0)*r(2);
    rr(1,0) = r(1)*r(0);
    rr(1,1) = r(1)*r(1);
    rr(1,2) = r(1)*r(2);
    rr(2,0) = r(2)*r(0);
    rr(2,1) = r(2)*r(1);
    rr(2,2) = r(2)*r(2);
    
    ddLk += g*dr + g*rr;
    dLk += g*r;
  }
  
  dLk2(0,0) = dLk(0)*dLk(0);
  dLk2(0,1) = dLk(0)*dLk(1);
  dLk2(0,2) = dLk(0)*dLk(2);
  dLk2(1,0) = dLk(1)*dLk(0);
  dLk2(1,1) = dLk(1)*dLk(1);
  dLk2(1,2) = dLk(1)*dLk(2);
  dLk2(2,0) = dLk(2)*dLk(0);
  dLk2(2,1) = dLk(2)*dLk(1);
  dLk2(2,2) = dLk(2)*dLk(2);
  
  res = dLk2/pow(Lk,2.0) - ddLk/Lk;
  
  return res;
}

//--------------------------------
//--------------------------------
// [[Rcpp::export]]
vec S(const vec & par, const vec & betahat, const vec &varbetahat, const vec &ldscore,const int & c0, const vec &Nstar, const int &num_threads) // total score for all K SNPs for each parameter
{
  int K = betahat.n_elem;
  vec res(3);
  res.fill(0.0);
  vec tem(3);
  long double res0=0, res1=0, res2=0;
  int k;
  
  // --------------------------------//--------------------------------
  omp_set_num_threads(num_threads);
  #pragma omp parallel for shared(par, betahat, varbetahat,ldscore,Nstar,c0) private(tem,k) reduction(+:res0,res1,res2)
  // --------------------------------//--------------------------------
  for(k=0; k<K; k++){
    tem = Sk(par, (betahat(k)), (varbetahat(k)), (ldscore(k)), c0, (Nstar(k)));
    res0 += tem(0);
    res1 += tem(1);
    res2 += tem(2);
  }
  res(0) = res0;
  res(1) = res1;
  res(2) = res2;
  
  return res;
}


//--------------------------------
//--------------------------------
// [[Rcpp::export]]
mat SS(const vec & par, const vec & betahat, const vec &varbetahat, const vec &ldscore,const int & c0, const vec &Nstar, const int &num_threads) // store the score function for each SNP k, summarize it as a K*3 matrix
{
  int K = betahat.n_elem;
  mat res(K,3); res.fill(0.0);
  vec tem(3);
  int k;
  
  // --------------------------------//--------------------------------
  omp_set_num_threads(num_threads);
  #pragma omp parallel for shared(par, betahat, varbetahat,ldscore,Nstar,c0,res) private(k,tem) 
  // --------------------------------//--------------------------------
  for(k=0; k<K; k++){
    tem = Sk(par,  (betahat(k)), (varbetahat(k)), (ldscore(k)), c0, (Nstar(k)));
    res(k,0) = tem(0);
    res(k,1) = tem(1);
    res(k,2) = tem(2);
  }
  
  return res;
}

//--------------------------------
//--------------------------------
// [[Rcpp::export]]
mat I(const vec & par, const vec & betahat, const vec &varbetahat, const vec &ldscore,const int & c0, const vec &Nstar, const int &num_threads) // 3*3 Information matrix
{
  int K = betahat.n_elem;
  mat res(3,3), tem(3,3);
  res.fill(0.0);
  long double res00=0, res01=0, res02=0, res10=0,res11=0,res12=0,res20=0,res21=0,res22=0;
  int k;
  
  // --------------------------------//--------------------------------
  omp_set_num_threads(num_threads);
  #pragma omp parallel for shared(par, betahat, varbetahat,ldscore,Nstar,c0) private(k,tem) reduction(+:res00,res01,res02,res10,res11,res12,res20,res21,res22)
  // --------------------------------//--------------------------------
  for(k=0; k<K; k++){
    tem = Ik(par,  (betahat(k)), (varbetahat(k)), (ldscore(k)), c0, (Nstar(k)));
    res00 += tem(0,0); res01 += tem(0,1); res02 += tem(0,2);
    res10 += tem(1,0); res11 += tem(1,1); res12 += tem(1,2);
    res20 += tem(2,0); res21 += tem(2,1); res22 += tem(2,2);
  }
  
  res(0,0) = res00; res(0,1) = res01; res(0,2) = res02;
  res(1,0) = res10; res(1,1) = res11; res(1,2) = res12;
  res(2,0) = res20; res(2,1) = res21; res(2,2) = res22;
  
  return res;
}

//--------------------------------
//--------------------------------
// [[Rcpp::export]]
List mixture_components_marginal(const vec & par, const vec &ldscore, const int & c0, const vec &Nstar, const int &num_threads)
{
  int K = ldscore.n_elem;
  
  long double pic = par(0);
  long double sigsq = par(1);
  
  int k,j;
  mat ww(K, (c0+1)); mat vv(K, (c0+1));
  ww.fill(0.0);  vv.fill(0.0);
  
  // weights formula.
  omp_set_num_threads(num_threads);
  #pragma omp parallel for shared(sigsq,pic,ldscore,Nstar,K,c0,ww,vv) private(k,j)
  for(k=0; k<K; k++){
    for(j=0; j<std::min(int(Nstar(k))+1,c0+1); j++){
      ww(k, j) = cdbinom(j,Nstar(k),pic,false);
      vv(k, j) = (j*sigsq)*ldscore(k)/Nstar(k);
    }
  }
  
  return List::create(
    _["proportions"] = ww,
    _["varcomponents"] = vv
  );
  
}



//--------------------------------//--------------------------------
// 3-component model
//--------------------------------//--------------------------------

//--------------------------------
//--------------------------------
// [[Rcpp::export]]
long double loglikelihood3(const vec & par, const vec &betahat, const vec &varbetahat, const vec &ldscore,const int & c0, const vec &Nstar, const int &num_threads) // loglikelihood function
{
    int K = betahat.n_elem;
    long double loginside;
    
    long double pic = par(0);
    long double p0 = par(1);
    long double sig1 = par(2);
    long double sig2 = par(3);
    long double a = par(4);
    
    long double res = 0;
    long double y=0;
    long double tem=0;
    
    long double pi1 = pic*p0;
    long double pi2 = pic*(1-p0);
    vec tem_prob(3);
    tem_prob(0) = pi1;
    tem_prob(1) = pi2;
    tem_prob(2) = 1-pi1-pi2;
    int k,j1,j2;
    
    // -------*-------*-------*-------*-------*-------*-------*
    omp_set_num_threads(num_threads);
    #pragma omp parallel for shared(a,betahat,varbetahat,ldscore,Nstar,tem_prob) private(loginside,y,k,tem,j1,j2) reduction(+:res)
    // -------*-------*-------*-------*-------*-------*-------*
    for(k=0; k<K; k++){
        loginside = 0;
        y = 0;
        for(j1=0; j1<std::min(int(Nstar(k))+1,c0+1); j1++){
            for(j2=0; j2<std::min(int(Nstar(k))+1,c0+1); j2++){
                if(Nstar(k) - j1 - j2 <0) break;
                
                tem = cdmultinom3(j1, j2, Nstar(k)- j1 - j2,tem_prob,false);
                loginside += tem*cdnorm(betahat(k),0, sqrt( (j1*sig1+ j2*sig2 )*ldscore(k)/Nstar(k) + varbetahat(k) + a), false);
                y+=tem;
            }
        }
        res += log(loginside/y);
    }
    
    return res;
}

//--------------------------------
// weights function
//--------------------------------
// [[Rcpp::export]]
mat weight3(const vec & par, const vec &betahat, const vec &varbetahat, const vec &ldscore, const int & c0, const vec & Nstar, const int &num_threads) // E-step: weights function
{
    int K = betahat.n_elem;
    mat w(K, (c0+1)*(c0+1));
    w.fill(0.0);
    
    long double pic = par(0);
    long double p0 = par(1);
    long double sig1 = par(2);
    long double sig2 = par(3);
    long double a = par(4) ;
    
    long double pi1 = pic*p0;
    long double pi2 = pic*(1-p0);
    vec tem_prob(3);
    tem_prob(0) = pi1;
    tem_prob(1) = pi2;
    tem_prob(2) = 1-pi1-pi2;
    
    int k,j1,j2;
    vec te;
    
    // weights formula.
    // -------*-------*-------*-------*-------*-------*-------*
    omp_set_num_threads(num_threads);
    #pragma omp parallel for  shared(a,sig1,sig2,betahat,varbetahat,ldscore,Nstar,w,tem_prob,K,c0) private(k,j1,j2)
    // -------*-------*-------*-------*-------*-------*-------*
    for(k=0; k<K; k++){
        for(j1=0; j1<std::min(int(Nstar(k))+1,c0+1); j1++){
            for(j2=0; j2<std::min(int(Nstar(k))+1,c0+1); j2++){
                if(Nstar(k) - j1 - j2 < 0) break;
                w(k, j1*(c0+1)+j2) =cdmultinom3(j1, j2, Nstar(k)- j1 - j2,tem_prob,false)*cdnorm(betahat(k),0, sqrt( (j1*sig1+ j2*sig2 )*ldscore(k)/Nstar(k) + varbetahat(k) + a), false);
            }
        }
    }
    
    te = sum(w,1);
    te = pow(te, -1.0);
    w = (diagmat(te)) * w;
    return w;
    
}

//--------------------------------
//--------------------------------
// [[Rcpp::export]]
List weight_loglikelihood3(const vec & par, const vec &betahat, const vec &varbetahat, const vec &ldscore,const int & c0, const vec &Nstar, const int &num_threads) // loglikelihood function
{
  int K = betahat.n_elem;
  long double loginside;
  
  long double pic = par(0);
  long double p0 = par(1);
  long double sig1 = par(2);
  long double sig2 = par(3);
  long double a = par(4);
  
  long double res = 0;
  long double y=0;
  long double tem=0,tem1=0;
  
  mat w(K, (c0+1)*(c0+1));
  w.fill(0.0);
  vec te;
  
  long double pi1 = pic*p0;
  long double pi2 = pic*(1-p0);
  vec tem_prob(3);
  tem_prob(0) = pi1;
  tem_prob(1) = pi2;
  tem_prob(2) = 1-pi1-pi2;
  int k,j1,j2;
  
  // -------*-------*-------*-------*-------*-------*-------*
  omp_set_num_threads(num_threads);
  #pragma omp parallel for shared(a,betahat,varbetahat,ldscore,Nstar,tem_prob,w,K,c0) private(tem1,loginside,y,k,tem,j1,j2) reduction(+:res)
  // -------*-------*-------*-------*-------*-------*-------*
  for(k=0; k<K; k++){
    loginside = 0;
    y = 0;
    for(j1=0; j1<std::min(int(Nstar(k))+1,c0+1); j1++){
      for(j2=0; j2<std::min(int(Nstar(k))+1,c0+1); j2++){
        if(Nstar(k) - j1 - j2 <0) break;
        tem = cdmultinom3(j1, j2, Nstar(k)- j1 - j2,tem_prob,false);
        tem1 = tem*cdnorm(betahat(k),0, sqrt( (j1*sig1+ j2*sig2 )*ldscore(k)/Nstar(k) + varbetahat(k) + a), false);
        
        w(k, j1*(c0+1)+j2) =tem1; 
        loginside += tem1;
        y+=tem;
      }
    }
    res += log(loginside/y);
  }
  
  te = sum(w,1);
  te = pow(te, -1.0);
  w = (diagmat(te)) * w;
  
  return List::create(
    _["w"] = w,
    _["llk"] = res
  );
}


//--------------------------------
//--------------------------------
// [[Rcpp::export]]
vec update_p3(const mat & w, const vec & Nstar, const int &num_threads)
{
    
    long double res0=0;
    long double res1=0;
    long double res2=0;
    long double tem;
    vec res(2);
    int c0 = sqrt(w.n_cols) - 1;
    int K = Nstar.n_elem;
    int k,j1,j2;
    
    // -------*-------*-------*-------*-------*-------*-------*
    omp_set_num_threads(num_threads);
    #pragma omp parallel for shared(Nstar,c0,w) private(k,tem,j1,j2) reduction(+:res0,res1,res2)
    // -------*-------*-------*-------*-------*-------*-------*
    for(k=0; k<K; k++){
        for(j1=0; j1<std::min(int(Nstar(k))+1,c0+1); j1++){
            for(j2=0; j2<std::min(int(Nstar(k))+1,c0+1); j2++){
                if(Nstar(k) - j1 - j2 <0) break;
                tem = w(k, j1*(c0+1)+j2);
                res0 += tem*(j1+j2);
                res1 += tem*Nstar(k);
                res2 += tem*j1;
            }
        }
    }
    
    res(0) = res0/res1;
    res(1) = res2/res0;
    
    return res;
}

//--------------------------------
//--------------------------------
// [[Rcpp::export]]
vec onestep_varcomponent3(const vec &varcomponent, const mat & w, const vec &betahat, const vec &varbetahat, const vec &ldscore, const vec & Nstar, const int &num_threads)
{
    int K = betahat.n_elem;
    int c0 = sqrt(w.n_cols) - 1;
    
    long double sig1 = varcomponent(0);
    long double sig2 = varcomponent(1);
    long double a = varcomponent(2);
    
    long double tem, det;
    long double s1 = 0, s2=0, s3 =0, dd11 = 0, dd22=0,dd33=0, dd12=0, dd13=0, dd23=0, dd123=0;
    int k,j1,j2;
    vec result(3);
    
    // -------*-------*-------*-------*-------*-------*-------*
    omp_set_num_threads(num_threads);
    #pragma omp parallel for shared(Nstar,a,sig1,sig2,ldscore,varbetahat,betahat,w) private(tem,k,j1,j2) reduction(+:s1,s2,s3,dd11,dd22,dd33,dd12,dd13,dd23,dd123)
    // -------*-------*-------*-------*-------*-------*-------*
    for(k=0; k<K; k++)
    {
        for(j1=0; j1<std::min(int(Nstar(k))+1,c0+1); j1++){
            for(j2=0; j2<std::min(int(Nstar(k))+1,c0+1); j2++){
                if(Nstar(k) - j1 - j2 <0) break;
                tem = (j1*sig1+j2*sig2) * ldscore(k) /double(Nstar(k)) + varbetahat(k) + a;
                s1 += w(k, j1*(c0+1)+j2)*j1*ldscore(k)/(2.0*Nstar(k))*(pow(betahat(k),2)/pow(tem,2) - 1.0/(tem)  );
                s2 += w(k, j1*(c0+1)+j2)*j2*ldscore(k)/(2.0*Nstar(k))*(pow(betahat(k),2)/pow(tem,2) - 1.0/(tem)  );
                s3 += w(k, j1*(c0+1)+j2)/(2.0)*(pow(betahat(k),2)/pow(tem,2) - 1.0/(tem)  );
                
                dd11 += w(k, j1*(c0+1)+j2)*j1*ldscore(k)*j1*ldscore(k)/(2.0*Nstar(k)*Nstar(k))*(-2.0*pow(betahat(k),2)/pow(tem,3) + 1.0/(pow(tem,2)));
                dd22 += w(k, j1*(c0+1)+j2)*j2*ldscore(k)*j2*ldscore(k)/(2.0*Nstar(k)*Nstar(k))*(-2.0*pow(betahat(k),2)/pow(tem,3) + 1.0/(pow(tem,2)));
                dd33 += w(k, j1*(c0+1)+j2)/(2.0)*(-2.0*pow(betahat(k),2)/pow(tem,3) + 1.0/(pow(tem,2)));
                dd12 += w(k, j1*(c0+1)+j2)*j1*ldscore(k)*j2*ldscore(k)/(2.0*Nstar(k)*Nstar(k))*(-2.0*pow(betahat(k),2)/pow(tem,3) + 1.0/(pow(tem,2)));
                dd13 += w(k, j1*(c0+1)+j2)*j1*ldscore(k)/(2.0*Nstar(k))*(-2.0*pow(betahat(k),2)/pow(tem,3) + 1.0/(pow(tem,2)));
                dd23 += w(k, j1*(c0+1)+j2)*j2*ldscore(k)/(2.0*Nstar(k))*(-2.0*pow(betahat(k),2)/pow(tem,3) + 1.0/(pow(tem,2)));
                dd123 += w(k, j1*(c0+1)+j2)*j1*ldscore(k)*j2*ldscore(k)/(2.0*Nstar(k)*Nstar(k))*(6.0*pow(betahat(k),2)/pow(tem,4) - 2.0/(pow(tem,3)));
            }
        }
    }
    
    det = dd11*(dd33*dd22 - dd23*dd23) - dd12*(dd33*dd12-dd23*dd13) + dd13*(dd23*dd12 - dd22*dd13);
    
    result(0) = sig1 - (1.0/det)*( (dd33*dd22-dd23*dd23)*s1-(dd33*dd12-dd23*dd13)*s2 + (dd23*dd12-dd22*dd13)*s3  );
    result(1) = sig2 -  (1.0/det)*( -(dd33*dd12-dd13*dd23)*s1 + (dd33*dd11-dd13*dd13)*s2 - (dd23*dd11-dd12*dd13)*s3  );
    result(2) = a - (1.0/det)*( (dd23*dd12-dd13*dd22)*s1 - (dd23*dd11-dd13*dd12)*s2 + (dd22*dd11-dd12*dd12)*s3  );
    
    return result;
}

//--------------------------------
//--------------------------------
// [[Rcpp::export]]
vec EM_func3(const vec &par_start, const vec & lower_pi, const vec & upper_pi,
             const vec &betahat, const vec & varbetahat, const vec & ldscore,  const vec & Nstar, const int & M,
             int c0,const long double &eps1,const long double &eps2,const long double &eps3,const long double &eps4,const long double &eps5,
             const long double &eps, const int &Meps, const int &steps, const int &num_threads, const bool &print, const int &printfreq, const bool &stratification)
{
    vec par(5), prev_par(5);
    par(0) = par_start(0);
    par(1) = par_start(1);
    par(2) = par_start(2);
    par(3) = par_start(3);
    par(4) = par_start(4);
    
    long double llk = 0,error_pi, error_p0, error_sig1, error_sig2, error_a, increase_ll;
    long double prev_llk, pic, p0, sig1, sig2, a;
    clock_t start_w, finish_w, finish_p, finish_sigsq; 
    
    List wllk;
    mat w;
    vec tem_p(2);
    vec tem_sig(3);
    vec old(3);
    vec result(12);
    
    Rcout << "Iteration, prev_loglikelihood, pic, p1, sigmasq1, sigmasq2, a, Heritability, c0, Seconds(weight_llk), Seconds(pi), Seconds(variance_components)";
    Rcout << endl;
    
    for(int i=0; i<Meps; i++){
        prev_llk = llk;
        prev_par(0) = par(0); prev_par(1) = par(1); prev_par(2) = par(2); prev_par(3) = par(3);prev_par(4) = par(4);
        
        pic = par(0);
        p0  = par(1);
        sig1 = par(2);
        sig2 = par(3);
        a = par(4);
        
        // update weight and log-likelihood
        start_w = clock();
        wllk = weight_loglikelihood3(par, betahat, varbetahat, ldscore,c0,Nstar,num_threads);
        w = as<mat>(wllk["w"]);
        llk = as<long double>(wllk["llk"]);
        finish_w = clock();
        
        // update proportion
        if ((pic>=lower_pi(0)) & (pic<=upper_pi(0)) & (p0>=lower_pi(1)) & (p0 <= upper_pi(1))) tem_p = update_p3(w,Nstar,num_threads);
        finish_p = clock(); 
        
        if (pic < lower_pi(0)) tem_p(0) = lower_pi(0);
        if (p0  < lower_pi(1)) tem_p(1) = lower_pi(1);
        if (pic > upper_pi(0)) tem_p(0) = upper_pi(0);
        if (p0  > upper_pi(1)) tem_p(1) = upper_pi(1);
        
        // update variance components
        tem_sig(0) = sig1; tem_sig(1) = sig2; tem_sig(2) = a;
        if(stratification==false){tem_sig(2)=0;}
        for(int j=0; j<steps; j++){
            old = onestep_varcomponent3(tem_sig, w, betahat,varbetahat, ldscore, Nstar,num_threads);
            if(stratification==false){old(2)=0;}
            if((abs(old(0)-tem_sig(0))<1e-20) & (abs(old(1)-tem_sig(1))<1e-20) & (abs(old(2)-tem_sig(2))<1e-20) ) break;
            tem_sig(0) = old(0); tem_sig(1) = old(1); tem_sig(2) = old(2);
        }
        if(tem_sig(0) > 1) tem_sig(0) = 1e-5;
        if(tem_sig(0) < 0) tem_sig(0) = 1e-12;
        if(tem_sig(1) > 1) tem_sig(1) = 1e-5;
        if(tem_sig(1) < 0) tem_sig(1) = 1e-12;
        if(tem_sig(2) > 1) tem_sig(2) = 1e-5;
        if(tem_sig(2) < -min(varbetahat)) tem_sig(2) = -min(varbetahat)/2;
        finish_sigsq = clock();
        
        par(0) = tem_p(0); par(1) =tem_p(1);
        par(2) = tem_sig(0); par(3) = tem_sig(1); par(4) = tem_sig(2);
        
        if(isinf(-llk) | isnan(-llk) | isnan(llk)) {par(4) = abs(par(4));}
        if((llk<prev_llk) & (c0<20)) {c0 = c0+1;}
        if((llk<prev_llk) & (c0>=20)) {break;}
        
        result(0) = i; result(1) = llk;
        result(2) = par(0); result(3) = par(1); result(4) = par(2); result(5) = par(3); result(6)= par(4);
        result(7) = M*par(0)*( par(1)*par(2) + (1-par(1))*par(3)); result(8) = c0;
        
        result(9) = double(finish_w-start_w)/CLOCKS_PER_SEC;
        result(10) = double(finish_p-finish_w)/CLOCKS_PER_SEC;
        result(11) = double(finish_sigsq-finish_p)/CLOCKS_PER_SEC;

        if(print==true){
          if(i%printfreq==0){
            for(int r=0; r<12; r++){
              Rcout << result(r)<< ", " ;
            }
            Rcout << endl;
          }
        }
        
        error_pi = abs(prev_par(0) - par(0)) ;
        error_p0 = abs(prev_par(1) - par(1)) ;
        error_sig1 = abs(prev_par(2) - par(2)) ;
        error_sig2 = abs(prev_par(3) - par(3)) ;
        error_a = abs(prev_par(4) - par(4)) ;
        increase_ll = (llk - prev_llk)/prev_llk;
        
        if(((error_pi< eps1) & (error_p0<eps2) & (error_sig1 <eps3) & (error_sig2<eps4) & (error_a<eps5)) & (abs(increase_ll) <eps)){break;}
    }
    
    // update log-likelihood
    llk = loglikelihood3(par,betahat,varbetahat,ldscore,c0,Nstar,num_threads);
    result(1) = llk;
    
    return result;
    
}


//--------------------------------
//--------------------------------
// [[Rcpp::export]]
vec Sk3(const vec & par, const long double &betahatk, const long double &varbetahatk, const long double &ldscorek,const int & c0, const int &Nstark) // score vector for kth SNP
{
    
    long double pic = par(0);
    long double p1 = par(1);
    long double sig1sq = par(2);
    long double sig2sq = par(3);
    long double a = par(4);
    
    vec tem_prob(3);
    tem_prob(0) = pic*p1;
    tem_prob(1) = pic*(1-p1);
    tem_prob(2) = 1-pic;
    int j1,j2;
    
    long double ww, denominator=0,
    numerator_pic=0,numerator_p1=0,numerator_sig1sq=0,
    numerator_sig2sq=0, numerator_a=0;
    vec res(5);
    long double tem_var;
    
    for(j1=0; j1<std::min(Nstark+1,c0+1); j1++){
        for(j2=0; j2<std::min(Nstark+1,c0+1); j2++){
            if(Nstark - j1 - j2 <0) break;
            tem_var = (j1*sig1sq+j2*sig2sq)*ldscorek/Nstark + varbetahatk + a;
            ww = cdmultinom3(j1, j2, Nstark- j1 - j2,tem_prob,false)*cdnorm(betahatk,0, sqrt( tem_var), false);
            denominator += ww;
            
            //pic
            numerator_pic += ww*((j1+j2)/pic - (Nstark-j1-j2)/(1-pic));
            numerator_p1 += ww*(j1/p1 - j2/(1-p1));
            
            //sigma1^2
            numerator_sig1sq += ww*(-0.5/tem_var + 0.5*pow(betahatk,2.0)/pow(tem_var,2.0))*j1*ldscorek/Nstark ;
            numerator_sig2sq += ww*(-0.5/tem_var + 0.5*pow(betahatk,2.0)/pow(tem_var,2.0))*j2*ldscorek/Nstark ;
            
            //sigma1^2
            numerator_a += ww*(-0.5/tem_var + 0.5*pow(betahatk,2.0)/pow(tem_var,2.0));
            
        }
    }
    res(0) = numerator_pic/denominator;
    res(1) = numerator_p1/denominator;
    res(2) = numerator_sig1sq/denominator;
    res(3) = numerator_sig2sq/denominator;
    res(4) = numerator_a/denominator;
    
    return res;
}





//--------------------------------
//--------------------------------
// [[Rcpp::export]]
mat Ik3(const vec & par, const long double &betahatk, const long double &varbetahatk, const long double &ldscorek,const int & c0, const int &Nstark) // Information for kth SNP, i.e., -partial(Sk)/partial(theta)
{
    long double pic = par(0);
    long double p1 = par(1);
    long double sig1sq = par(2);
    long double sig2sq = par(3);
    long double a = par(4);
    
    vec tem_prob(3);
    tem_prob(0) = pic*p1;
    tem_prob(1) = pic*(1-p1);
    tem_prob(2) = 1-pic;
    int j1,j2;
    
    long double  tem_var, g, Lk=0;
    vec r(5),  dLk(5);
    dLk.fill(0.0); r.fill(0.0);
    mat dr(5,5), rr(5,5), ddLk(5,5), dLk2(5,5), res(5,5);
    ddLk.fill(0.0);
    dLk2.fill(0.0);
    dr.fill(0.0);
    rr.fill(0.0);
    res.fill(0.0);
    
    for(j1=0; j1<std::min(Nstark+1,c0+1); j1++){
        for(j2=0; j2<std::min(Nstark+1,c0+1); j2++){
            if(Nstark - j1 - j2 <0) break;
            
            tem_var = (j1*sig1sq+j2*sig2sq)*ldscorek/Nstark + varbetahatk + a;
            g = cdmultinom3(j1, j2, Nstark- j1 - j2,tem_prob,false)*cdnorm(betahatk,0, sqrt(tem_var), false);
            
            r(0) = (j1+j2)/pic - (Nstark-j1-j2)/(1-pic);
            r(1) = j1/p1 - j2/(1-p1);
            r(2) = (-0.5/tem_var + 0.5*pow(betahatk,2.0)/pow(tem_var,2.0))*j1*ldscorek/Nstark ;
            r(3) = (-0.5/tem_var + 0.5*pow(betahatk,2.0)/pow(tem_var,2.0))*j2*ldscorek/Nstark ;
            r(4) = (-0.5/tem_var + 0.5*pow(betahatk,2.0)/pow(tem_var,2.0));
            
            Lk += g;
            
            dr(0,0) = -(j1+j2)/pow(pic,2.0) - (Nstark-j1-j2)/(pow(1-pic, 2.0));
            dr(1,1) = -(j1)/pow(p1,2.0) - (j2)/(pow(1-p1, 2.0));
            dr(2,2) = (0.5/pow(tem_var,2.0) - pow(betahatk,2.0)/pow(tem_var,3.0))*(j1*ldscorek/Nstark)*(j1*ldscorek/Nstark);
            dr(2,3) = (0.5/pow(tem_var,2.0) - pow(betahatk,2.0)/pow(tem_var,3.0))*(j1*ldscorek/Nstark)*(j2*ldscorek/Nstark);
            dr(2,4) = (0.5/pow(tem_var,2.0) - pow(betahatk,2.0)/pow(tem_var,3.0))*(j1*ldscorek/Nstark);
            dr(3,2) = dr(2,3);
            dr(3,3) = (0.5/pow(tem_var,2.0) - pow(betahatk,2.0)/pow(tem_var,3.0))*(j2*ldscorek/Nstark)*(j2*ldscorek/Nstark);
            dr(3,4) = (0.5/pow(tem_var,2.0) - pow(betahatk,2.0)/pow(tem_var,3.0))*(j2*ldscorek/Nstark);
            dr(4,2) = dr(2,4);
            dr(4,3) = dr(3,4);
            dr(4,4) = (0.5/pow(tem_var,2.0) - pow(betahatk,2.0)/pow(tem_var,3.0));
            
            rr(0,0) = r(0)*r(0); rr(0,1) = r(0)*r(1); rr(0,2) = r(0)*r(2); rr(0,3) = r(0)*r(3); rr(0,4) = r(0)*r(4);
            rr(1,0) = rr(0,1);   rr(1,1) = r(1)*r(1); rr(1,2) = r(1)*r(2); rr(1,3) = r(1)*r(3); rr(1,4) = r(1)*r(4);
            rr(2,0) = rr(2,0);   rr(2,1) = rr(1,2);   rr(2,2) = r(2)*r(2); rr(2,3) = r(2)*r(3); rr(2,4) = r(2)*r(4);
            rr(3,0) = rr(0,3);   rr(3,1) = rr(1,3);   rr(3,2) = rr(2,3);   rr(3,3) = r(3)*r(3); rr(3,4) = r(3)*r(4);
            rr(4,0) = rr(0,4);   rr(4,1) = rr(1,4);   rr(4,2) = rr(2,4);   rr(4,3) = rr(3,4);   rr(4,4) = r(4)*r(4);
            
            ddLk += g*dr + g*rr;
            dLk += g*r;
        }
    }
    
    dLk2(0,0) = dLk(0)*dLk(0);dLk2(0,1) = dLk(0)*dLk(1);dLk2(0,2) = dLk(0)*dLk(2);dLk2(0,3) = dLk(0)*dLk(3);dLk2(0,4) = dLk(0)*dLk(4);
    dLk2(1,0) = dLk2(0,1);    dLk2(1,1) = dLk(1)*dLk(1);dLk2(1,2) = dLk(1)*dLk(2);dLk2(1,3) = dLk(1)*dLk(3);dLk2(1,4) = dLk(1)*dLk(4);
    dLk2(2,0) = dLk2(0,2);    dLk2(2,1) = dLk2(1,2);    dLk2(2,2) = dLk(2)*dLk(2);dLk2(2,3) = dLk(2)*dLk(3);dLk2(2,4) = dLk(2)*dLk(4);
    dLk2(3,0) = dLk2(0,3);    dLk2(3,1) = dLk2(1,3);    dLk2(3,2) = dLk2(2,3);    dLk2(3,3) = dLk(3)*dLk(3);dLk2(3,4) = dLk(3)*dLk(4);
    dLk2(4,0) = dLk2(0,4);    dLk2(4,1) = dLk2(1,4);    dLk2(4,2) = dLk2(2,4);    dLk2(4,3) = dLk2(3,4);    dLk2(4,4) = dLk(4)*dLk(4);
    
    res = dLk2/pow(Lk,2.0) - ddLk/Lk;
    
    return res;
}

//--------------------------------
//--------------------------------
// [[Rcpp::export]]
vec S3(const vec & par, const vec & betahat, const vec &varbetahat, const vec &ldscore,const int & c0, const vec &Nstar) // total score for all K SNPs for each parameter
{
    int K = betahat.n_elem;
    vec res(5);
    vec tem(5);
    res.fill(0.0);
    
    for(int k=0; k<K; k++){
        tem = Sk3(par,  (betahat(k)), (varbetahat(k)), (ldscore(k)), c0, (Nstar(k)))/K;
        res(0) += tem(0);
        res(1) += tem(1);
        res(2) += tem(2);
        res(3) += tem(3);
        res(4) += tem(4);
    }
    return res;
}

//--------------------------------
//--------------------------------
// [[Rcpp::export]]
mat SS3(const vec & par, const vec & betahat, const vec &varbetahat, const vec &ldscore,const int & c0, const vec &Nstar, const int &num_threads) // store the score function for each SNP k, summarize it as a K*3 matrix
{
    int K = betahat.n_elem;
    mat res(K,5); res.fill(0.0);
    vec tem(5);
    int k;
    
    // --------------------------------//--------------------------------
    omp_set_num_threads(num_threads);
    #pragma omp parallel for shared(par, betahat, varbetahat,ldscore,Nstar,c0,res) private(k,tem)
    // --------------------------------//--------------------------------
    for(k=0; k<K; k++){
        tem = Sk3(par,  (betahat(k)), (varbetahat(k)), (ldscore(k)), c0, (Nstar(k)));
        res(k,0) = tem(0);
        res(k,1) = tem(1);
        res(k,2) = tem(2);
        res(k,3) = tem(3);
        res(k,4) = tem(4);
    }
    
    return res;
}


//--------------------------------
//--------------------------------
// [[Rcpp::export]]
mat I3(const vec & par, const vec & betahat, const vec &varbetahat, const vec &ldscore,const int & c0, const vec &Nstar, const int &num_threads) // 3*3 Information matrix
{
    int K = betahat.n_elem;
    mat res(5,5), tem(5,5);
    res.fill(0.0), tem.fill(0.0);
    long double res00=0,res01=0,res02=0,res03=0,res04=0,res10=0,res11=0,res12=0,res13=0,res14=0,res20=0,res21=0,res22=0,res23=0,res24=0;
    long double res30=0,res31=0,res32=0,res33=0,res34=0,res40=0,res41=0,res42=0,res43=0,res44=0;
    int k;
    
    // --------------------------------//--------------------------------
    omp_set_num_threads(num_threads);
    #pragma omp parallel for shared(par, betahat, varbetahat,ldscore,Nstar,c0) private(k,tem) reduction(+:res00,res01,res02,res03,res04,res10,res11,res12,res13,res14,res20,res21,res22,res23,res24,res30,res31,res32,res33,res34,res40,res41,res42,res43,res44)
    // --------------------------------//--------------------------------
    for(k=0; k<K; k++){
      tem = Ik3(par,  (betahat(k)), (varbetahat(k)), (ldscore(k)), c0, (Nstar(k)));
      res00 += tem(0,0); res01 += tem(0,1); res02 += tem(0,2); res03 += tem(0,3); res04 += tem(0,4); 
      res10 += tem(1,0); res11 += tem(1,1); res12 += tem(1,2); res13 += tem(1,3); res14 += tem(1,4); 
      res20 += tem(2,0); res21 += tem(2,1); res22 += tem(2,2); res23 += tem(2,3); res24 += tem(2,4); 
      res30 += tem(3,0); res31 += tem(3,1); res32 += tem(3,2); res33 += tem(3,3); res34 += tem(3,4); 
      res40 += tem(4,0); res41 += tem(4,1); res42 += tem(4,2); res43 += tem(4,3); res44 += tem(4,4); 
    }
    
    res(0,0) = res00; res(0,1) = res01; res(0,2) = res02; res(0,3) = res03; res(0,4) = res04;
    res(1,0) = res10; res(1,1) = res11; res(1,2) = res12; res(1,3) = res13; res(1,4) = res14;
    res(2,0) = res20; res(2,1) = res21; res(2,2) = res22; res(2,3) = res23; res(2,4) = res24;
    res(3,0) = res30; res(3,1) = res31; res(3,2) = res32; res(3,3) = res33; res(3,4) = res34;
    res(4,0) = res40; res(4,1) = res41; res(4,2) = res42; res(4,3) = res43; res(4,4) = res44;
    
    return res;
}


//--------------------------------
//--------------------------------
// [[Rcpp::export]]
List mixture_3components_marginal(const vec & par, const vec &ldscore,const int & c0, const vec &Nstar, const int &num_threads)
{
  int K = ldscore.n_elem;
  
  long double pic = par(0);
  long double p0 = par(1);
  long double sig1 = par(2);
  long double sig2 = par(3);

  long double pi1 = pic*p0;
  long double pi2 = pic*(1-p0);
  vec tem_prob(3);
  tem_prob(0) = pi1;
  tem_prob(1) = pi2;
  tem_prob(2) = 1-pi1-pi2;
  int k,j1,j2;
  mat ww(K, (c0+1)*(c0+1)); mat vv(K, (c0+1)*(c0+1));
  ww.fill(0.0);  vv.fill(0.0);
  
  // weights formula.
  omp_set_num_threads(num_threads);
  #pragma omp parallel for  shared(sig1,sig2,ldscore,Nstar,tem_prob,K,c0,ww,vv) private(j1,j2,k)
  for(k=0; k<K; k++){
    for(j1=0; j1<std::min(int(Nstar(k))+1,c0+1); j1++){
      for(j2=0; j2<std::min(int(Nstar(k))+1,c0+1); j2++){
        if(Nstar(k) - j1 - j2 < 0) break;
        ww(k, j1*(c0+1)+j2) = cdmultinom3(j1, j2, Nstar(k)- j1 - j2,tem_prob,false);
        vv(k, j1*(c0+1)+j2) =  (j1*sig1+ j2*sig2 )*ldscore(k)/Nstar(k); 
      }
    }
  }
  
  return List::create(
    _["proportions"] = ww,
    _["varcomponents"] = vv
  );
  
}
