#include <iostream>
#include <string>
#include <RcppArmadillo.h>
#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
using namespace Rcpp;
using namespace std;


//Function to calculate the Whittle likelihood for a given set of parameters
// [[Rcpp::export]]
double Whittle(const vec& omegas, const vec& Pgram, const vec& P, const vec& Z, int Lambda, 
               int L, double MaxFreq){
  //omegas is a vector of fourier frequencies
  //Pgram is a vector of corresponding periodogram values
  //P is an L length vector of BDP weights
  //Z is an L length vector of BDP atoms
  //Lambda is the number of Bernstein polynomial components to use
  //L is the truncation level of the BDP
  //MaxFreq is the maximum frequency to consider for the spectral density
  
  int N_omega=omegas.n_elem;
  vec f(N_omega, fill::zeros);
  vec weights(Lambda, fill::zeros);
  //calculate Bernstein polynomial components weights from BDP weights
  for(int l=0; l<L; l++){
    weights(floor(Z(l)*Lambda)) += P(l);
  }
  //calculate estimated spectral density given parameters and add tiny number to avoid errors from rounding to 0
  for(int n_omega=0; n_omega<N_omega; n_omega++){
    for(int lambda=0; lambda<Lambda; lambda++){
      if(weights(lambda) > 0){f(n_omega) += weights(lambda)*R::dbeta(omegas(n_omega)/MaxFreq,lambda+1,Lambda-lambda,false) + pow(10,-323);}
    }
  }
  //return the log Whittle likelihood value
  return sum(-log(f) - Pgram/f);
}


//function to sample number of Bernstein polynomials, Lambda, using MH step
// [[Rcpp::export]]
int sample_Lambda(int N, int L, const vec& omegas, const mat& Pgrams, const vec& P_cur, 
                  const mat& Z_cur, int Lambda_cur, string Lambda_prior, int Lambda_init, 
                  int Lambda_max, double MaxFreq, string Lambda_prop_dist, bool DBDP){
  //N is the number of subjects
  //L is BDP truncation level
  //omegas is a vector of fourier frequencies
  //Pgrams is a matrix of corresponding periodogram values with a column for each subject
  //P_cur is a vector containing the latest iteration BDP weights
  //Z_cur is a matrix containing the latest iteration BDP atoms with a column for each subject
  //Lambda_cur is the current number of Bernstein polynomial components
  //Lambda_prior is the prior for the number of Bernstein polynomial components
  //Lambda_init is the initial number of Bernstein polynomial components in the MCM chain
  //Lambda_max is the maximum number of Bernstein polynomial components to consider
  //MaxFreq is the maximum frequency to consider for the spectral density
  //Lambda_prop_dist is the proposal distribution to be used for Lambda, up1down1 or poisson. If incorrectly/not specified Lambda will be constant
  //DBDP is a boolean to determine if DP atoms are dependent n covariates
  
  int Lambda_new;
  int Lambda_prop = Lambda_cur;
  double accept_Lambda = 0;
  vec Z_subj(L);
  if(DBDP == false){Z_subj = Z_cur.col(0);}
  //propose new value for Lambda using poisson proposal
  if(Lambda_prop_dist=="up1down1"){
    if(Lambda_cur==1){Lambda_prop=2;}else{
      if(Lambda_cur==Lambda_max){Lambda_prop=Lambda_max-1;}else{
        if(randu()<0.5){Lambda_prop=Lambda_cur-1;}else{Lambda_prop=Lambda_cur+1;}
      }
    }
  }
  if(Lambda_prop_dist=="poisson"){
    Lambda_prop=R::rpois(Lambda_cur);
    if(Lambda_prop==0){Lambda_prop=1;}else{if(Lambda_prop>Lambda_max){Lambda_prop=Lambda_max;}}
    accept_Lambda = R::dpois(Lambda_cur,Lambda_prop,true) - 
      R::dpois(Lambda_prop,Lambda_cur,true);
  }
  
  //calculate MH acceptance probabilities with specifiec prior
  if(Lambda_prop==Lambda_cur){Lambda_new=Lambda_cur;}else{
    if(Lambda_prior=="flat"){accept_Lambda += 0;}
    if(Lambda_prior=="poisson"){accept_Lambda += R::dpois(Lambda_prop,Lambda_init,true) - 
      R::dpois(Lambda_cur,Lambda_init,true);}
    if(Lambda_prior=="expon"){accept_Lambda += -0.05*pow(Lambda_prop,2) + 
      0.05*pow(Lambda_cur,2);}
    if(N >= 1){
      for(int i=0; i<N; i++){
        if(DBDP == true){Z_subj = Z_cur.col(i);}
        accept_Lambda = accept_Lambda +
          Whittle(omegas, Pgrams.col(i), P_cur, Z_subj, Lambda_prop, L, MaxFreq) -
          Whittle(omegas, Pgrams.col(i), P_cur, Z_subj, Lambda_cur, L, MaxFreq);
      }
    }
    //accept or reject proposed Lambda
    if(accept_Lambda > log(randu())){Lambda_new=Lambda_prop;}else{Lambda_new=Lambda_cur;}
  }
  return Lambda_new;
}      


//Function to calculate BDP weights from stick break proportions V for truncation level L
// [[Rcpp::export]]
vec GetWeights(const vec& V, int L){
  vec CumProds=cumprod(1-V);
  vec out(L);
  out(0)=V(0);
  out.tail(L-1)=CumProds.head(L-1)%V.tail(L-1);
  return out;
}


//Function to propose new element of V (weights) or Z (atoms) using uniform proposal
// [[Rcpp::export]]
vec NewVZ(const vec& VZ,const vec& epsilon,double dim){
  //VZ is vector of current weights or atoms
  //epsilon is a vector of uniform halfwidths for proposals
  //dim is the element of the weight or atom vector which is being proposed
  
  vec VZnew=VZ;
  //propose new value
  double prop=VZ(dim)+epsilon(dim)*(2*randu()-1);
  //take modulus 1 to ensure it falls between 0 and 1
  if(prop < 0){VZnew(dim)=1+prop;}else{if(prop < 1){VZnew(dim)=prop;}else{VZnew(dim)=prop-1;}}
  return VZnew;
}


//function to sample stick breaking weights, V, using MH step and uniform proposal
// [[Rcpp::export]]
List sample_V(int N, int L, const vec& omegas, const mat& Pgrams, const vec& V_cur, 
              const vec& P_cur, const mat& Z_cur, int Lambda_cur, const vec& epsilon, 
              double alpha_L_cur, double MaxFreq, bool DBDP){
  //N is the number of subjects
  //L is BDP truncation level
  //omegas is a vector of fourier frequencies
  //Pgrams is a matrix of corresponding periodogram values with a column for each subject
  //V_cur is a vector containing the latest iteration BDP stick breaking weights
  //P_cur is a vector containing the latest iteration BDP weights
  //Z_cur is a matrix containing the latest iteration BDP atoms with a column for each subject
  //Lambda_cur is the current number of Bernstein polynomial components
  //epsilon is a vector of uniform halfwidths for proposals
  //alpha_L_cur is the latest value of the alpha_L concentration parameter
  //MaxFreq is the maximum frequency to consider for the spectral density
  //DBDP is a boolean to determine if DP atoms are dependent n covariates
  
  vec Z_subj(L);
  if(DBDP == false){Z_subj = Z_cur.col(0);}
  vec V_prop(L);
  vec P_prop(L);
  vec V_new=V_cur;
  vec P_new=P_cur;
  double V_new_lik = (alpha_L_cur-1)*sum(log(1-V_new.head(L-1)));
  if(N >= 1){
    for(int i=0; i<N; i++){
      if(DBDP == true){Z_subj = Z_cur.col(i);}
      V_new_lik = V_new_lik +
        Whittle(omegas, Pgrams.col(i), P_new, Z_subj, Lambda_cur, L, MaxFreq);
    }
  }
  double V_prop_lik;
  double accept_V;  
  //loop through updating each of the L-1 random values of stick breaking weights V
  for(int l=0; l<(L-1); l++){
    //propose new V for l-th dimension
    V_prop=NewVZ(V_new, epsilon, l);
    P_prop=GetWeights(V_prop, L);
    //calculate acceptance probability
    V_prop_lik = (alpha_L_cur-1)*sum(log(1-V_prop.head(L-1)));
    if(N >= 1){
      for(int i=0; i<N; i++){
        if(DBDP == true){Z_subj = Z_cur.col(i);}
        V_prop_lik = V_prop_lik +
          Whittle(omegas, Pgrams.col(i), P_prop, Z_subj, Lambda_cur, L, MaxFreq);
      }
    }
    accept_V = V_prop_lik- V_new_lik;
    //accept or reject new V value using MH step
    if(accept_V > log(randu())){
      V_new=V_prop;
      P_new=P_prop;
      V_new_lik = V_prop_lik;
    }
  }
  return List::create(
    _["V_new"] = V_new,
    _["P_new"] = P_new
  ); 
}


//function to sample atoms, Z, using MH step and uniform proposal
// [[Rcpp::export]]
vec sample_Z(int N, int L, const vec& omegas, const mat& Pgrams, const vec& P_cur, 
             const vec& Z_cur, int Lambda_cur, const vec& epsilon, double MaxFreq){
  //N is the number of subjects
  //L is BDP truncation level
  //omegas is a vector of fourier frequencies
  //Pgrams is a matrix of corresponding periodogram values with a column for each subject
  //P_cur is a vector containing the latest iteration BDP weights
  //Z_cur is a vector containing the latest iteration BDP atoms
  //Lambda_cur is the current number of Bernstein polynomial components
  //epsilon is a vector of uniform halfwidths for proposals
  //MaxFreq is the maximum frequency to consider for the spectral density
  
  vec Z_prop(L);
  vec Z_new=Z_cur;
  double Z_new_lik = 0;
  if(N >= 1){
    for(int i=0; i<N; i++){
      Z_new_lik = Z_new_lik +
        Whittle(omegas, Pgrams.col(i), P_cur, Z_new, Lambda_cur, L, MaxFreq);
    }
  }
  double Z_prop_lik;
  double accept_Z;  
  //loop through updating each of the L BDP atoms, Z
  for(int l=0; l<L; l++){
    Z_prop_lik = 0;
    //propose new l-th atom
    Z_prop=NewVZ(Z_new, epsilon, l);
    //calculate acceptance probability
    if(N >= 1){
      for(int i=0; i<N; i++){
        Z_prop_lik = Z_prop_lik +
          Whittle(omegas, Pgrams.col(i), P_cur, Z_prop, Lambda_cur, L, MaxFreq);
      }
    }
    accept_Z = Z_prop_lik - Z_new_lik;
    //acceptance or reject using MH step
    if(accept_Z > log(randu())){Z_new = Z_prop; Z_new_lik = Z_prop_lik;}
  }
  return Z_new;
}


//function to calculate posteior curve samples from posteior variables
// [[Rcpp::export]]
mat GetCurvesBDP(int Nsamp, const vec& omegas, int L, const vec& Lambda_samps, 
                 const mat& P_samps, const mat& Z_samps, double MaxFreq){
  //Nsamp is the number of samples
  //omegas is a vector of fourier frequencies
  //L is BDP truncation level
  //Lambda_samps is a vector of samples of the number of Bernstein polynomial components
  //P_samps is a matrix where each row is an L length sample of BDP weights
  //Z_samps is a matrix where each row is an L length sample of BDP atoms
  //MaxFreq is the maximum frequency to consider for the spectral density
  
  int Lambda_now;
  vec P_now(L);
  vec Z_now(L);
  int nOm=omegas.n_elem;
  mat CurveEsts(Nsamp,nOm,fill::zeros);
  for(int s=0; s<Nsamp; s++){
    Lambda_now=Lambda_samps(s);
    P_now=P_samps.row(s).t();
    Z_now=Z_samps.row(s).t();
    vec w(Lambda_now,fill::zeros);
    for(int l=0; l<L; l++){
      w(floor(Z_now(l)*Lambda_now)) += P_now(l);
    }
    for(int t=0; t<nOm; t++){
      for(int lambda=0; lambda<Lambda_now; lambda++){
        CurveEsts(s, t) += w(lambda)*R::dbeta(omegas(t)/MaxFreq,lambda+1,Lambda_now-lambda,false)/MaxFreq;
      }
    }
  }
  return CurveEsts;
}


//MCMC Sampling Algorithm for single subject spectral BDP
// [[Rcpp::export]]
List SpectralBDP(const vec& Pgram,const vec& omegas , int Nsamp, int L, const vec& epsilon, 
                 double SamplingRate=1, double MaxFreq=0.5, int Sup=500, 
                 int Lambda_init=30, string Lambda_prior="flat", int Lambda_max=300, 
                 string Lambda_prop_dist="poisson", bool alpha_const=true, 
                 double alpha_init=1, double a_alpha_L=1, double b_alpha_L=1){
  //X is a vector containing the time series to be analyzed
  //Nsamp is the number of posterior samples
  //L is the truncation level for the BDP
  //epsilon is a vector of uniform halfwidths for V and Z proposals
  //SamplingRate is the sampling rate of the time series
  //MaxFreq is the maximum frequnecy to consider for the analysis, at most 1/2 the sampling rate
  //Sup is how often status updates are printed to the console, every Sup iterations
  //Lambda_init is the initial value for Lambda and the value for all samples if Lambda is set to constant
  //Lambda_prior is the prior on the number of Bernstein polynomial components, Lambda.
  //either "expon", "poisson", or "flat"
  //Lambda_max is the maximum number of Bernstein polynomial components to consider
  //Lambda_prop_dist is the proposal distribution to be used for Lambda MH steps
  //either "poisson" or "up1down1". If mispecifiec or "const" then lambda will be constant 
  //alpha_const is a boolean to determine if the concetration parameter is set to constant
  //alpha_init is the initial value for concentration parameter and the value for all samples if constant
  //a_alpha_L and b_alpha_L are the gamma hyperparameters to be used if alpha is to be sampled
  
  //determine Fourier frequencies to be considered
  //double T = X.n_elem;
 // int N_freq =omegas.n_elem; //   floor((T-1)/2);
//  int N_freqKeep = 0;
//  vec omegas(N_freq);
//  for(int n_freq=0; n_freq < N_freq; n_freq++){
//    omegas(n_freq) = (n_freq+1)/T*SamplingRate;
//    if((n_freq+1)/T*SamplingRate < MaxFreq){N_freqKeep += 1;}
//  }
//  omegas=omegas.head(N_freqKeep);
  
  //calculate periodograms
//  vec Pgram=square(abs(fft(X)))/T;
 // Pgram=Pgram.head(N_freqKeep);
  
  
  
  //initialize concentration parameter alpha
  vec alpha_L(Nsamp);
  alpha_L(0)=alpha_init;
  
  //initialize BDP (from Choudhuri)
  vec Lambda(Nsamp);
  Lambda(0)=Lambda_init;
  
  mat V(Nsamp, L);
  V.submat(0,0,0,L-2).fill(0.5);
  V(0,L-1)=1;
  
  mat P(Nsamp, L);
  P.row(0)=GetWeights(V.row(0).t(), L).t();
  List NewVP(2);
  
  mat Z(Nsamp, L);
  Z.row(0).randu();
  
  for(int iter=1; iter<Nsamp; iter++){
    //Update number of Bernstein polynomials, Lambda
    Lambda(iter) = sample_Lambda(1, L, omegas, Pgram, P.row(iter-1).t(), Z.row(iter-1).t(), 
                                 Lambda(iter-1), Lambda_prior, Lambda_init, Lambda_max, MaxFreq, Lambda_prop_dist, false);
    
    //Update weights V and P
    NewVP = sample_V(1, L, omegas, Pgram, V.row(iter-1).t(), P.row(iter-1).t(), 
                     Z.row(iter-1).t(), Lambda(iter), epsilon, alpha_L(iter-1), MaxFreq, false);
    V.row(iter) = as<vec>(NewVP("V_new")).t();
    P.row(iter) = as<vec>(NewVP("P_new")).t();
    
    //Update Z
    Z.row(iter) = sample_Z(1, L, omegas, Pgram, P.row(iter).t(), Z.row(iter-1).t(), 
                           Lambda(iter), epsilon, MaxFreq).t();
    
    //Update alpha_L via Gibbs step
    if(alpha_const==false){
      alpha_L(iter)=R::rgamma(a_alpha_L+L-1, b_alpha_L - sum(log(V.row(iter).head(L-1))));
    }
    if(alpha_const==true){
      alpha_L(iter)=alpha_L(iter-1);
    }
    
  //  if((iter-1)/Sup!=iter/Sup){cout << iter << " posterior samples complete" << endl;}
  }
  
  mat Curves=GetCurvesBDP(Nsamp, omegas, L, Lambda, P, Z, MaxFreq);
  
  return List::create(
    _["Lambda"] = Lambda,
    _["P"] = P,
    _["V"] = V,
    _["Z"] = Z,
    _["CurveEsts"] = Curves,
    _["omegas"] = omegas,
    _["Pgram"] = Pgram
  ); 
}  