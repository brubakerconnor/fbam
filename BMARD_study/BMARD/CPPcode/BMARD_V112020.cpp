#include <iostream>
#include <string>
#include <cmath>
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <math.h> 
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
using namespace Rcpp;
using namespace std;

// modification based on spectral splines prior from the SDF estimation method of Choudhury (2014)

//spectral spline function
// compute the mixture of AR(2) spectrums over a set of frequencies 
// [[Rcpp::export]]
vec SSpline(const vec& phi,const vec& psi,const vec& L, const vec& freq){
  //phi is a vector of component weights
  //psi is a vector of components AR2 phases
  //L is a vector of component log modulus each based on a AR2 process  
  double pi=3.14159265358979323846;
  vec fval(freq.size(), fill::zeros);
  for(int j =0; j< freq.size() ; j++ ){
    for(int i =0; i< phi.size() ; i++ ){
      fval[j]+= phi[i]*((((1-exp(-2*L[i]))/(1+exp(-2*L[i])))*(  pow(1+exp(-2*L[i]) ,2)  -4*(  pow( cos(2*pi*psi[i]),2 )  )*exp(-2*L[i]) ))/(  pow( 1-2*exp(-L[i])*cos(2*pi*psi[i])*cos(2*pi*freq[j])+exp(-2*L[i])*cos(4*pi*freq[j]), 2) +  pow (-2*exp(-L[i])*cos(2*pi*psi[i])*sin(2*pi*freq[j])+exp(-2*L[i])*sin(4*pi*freq[j]),2 )   ));  
    }  // end  for i
  } // end for j
  return 2*fval;
}



//Function to calculate the Whittle likelihood for a given set of parameters
// [[Rcpp::export]]
double Whittle(const vec& omegas,const vec& psi_cur,const vec& BW_cur, const vec& Pgram, const vec& P, const vec& Z, int Lambda, 
               int L, double MaxFreq){
  //omegas is a vector of fourier frequencies
  //Pgram is a vector of corresponding periodogram values
  //P is an L length vector of BDP weights
  //Z is an L length vector of BDP atoms
  //Lambda is the number of spectral components to use
  //L is the truncation level of the BDP
  //MaxFreq is the maximum frequency to consider for the spectral density
  //psi_cur is a vector containing the latest iteration of the location parameters
  //Bwidth is a vector containing the latest iteration of the bandwidth (logmodulus) parameters
  
  int N_omega=omegas.size();
  int C= psi_cur.size() ;
  vec f(N_omega, fill::zeros);
  vec weights(Lambda, fill::zeros);
  //calculate  components weights from BDP weights
  for(int l=0; l<L; l++){
    weights[  floor(Z[l]*Lambda)  ] += P[l];
  }
  
  //calculate estimated spectral density given parameters and add tiny number to avoid errors from rounding to 0
  
  f=SSpline(weights, psi_cur, BW_cur,omegas);
  
  //return the log Whittle likelihood value
  return sum(-log(f) - Pgram/f) -sum( log(BW_cur)  ) -.5*pow( C,2) ;
  
  // note: using exponential cuadratic prior for C  -pow( C,2) leads to to a simpler model this prior is found in choudhuri (2014)
  // other possible prior is negative exponential  -C in log scale is a lower penalization could depend in a paraqmter lambda 
  // note using a poisson prior leads to a more complex model and sometimes misleading the components  + lgamma(C+1)   +C*log(.1)
}

//Function to calculate BDP weights from stick break proportions V for truncation level L
// [[Rcpp::export]]
vec GetWeights( vec V, int L){
  
  vec aux=  V.subvec(0,L-2) ;
  vec np(1);
  np(0)=1;
  vec CumProds= join_cols(np, cumprod(1-aux) );
  vec out=CumProds%V;
  return out;
  
}


//Function to propose new location parameter (phase of AR(2)) for a subinterval j with k total components 
// [[Rcpp::export]]
vec Newpsij(const vec& psi ,   const double& epsilonpsi, int j, vec partition ){
  //current phi vector 
  //current psi vector 
  //current L vector
  // j is the index for the subinteval to sample the location parameter psi[j]
  // k total of components considerer
  //epsilon is a vector of uniform halfwidths for proposals chose considering the size of the subinreval that is 1/k
  
  double prop=psi(j);
  vec propvec=psi ;
  double psinew;
  double v= (partition(j+1)-partition(j))/10.0    ; //1.0/static_cast<double>(2*k); 
  //double d=static_cast<double>(j)/static_cast<double>(2*k);
  //double z= static_cast<double>(j-1)/static_cast<double>(2*k);
  //propose new value
  prop = prop +v*(2*randu()-1);
  //take modulus inside the subinterval to ensure it falls inside the subinterval specified
  if(prop < partition(j) ){psinew=v+prop;}else{if(prop < partition(j+1)){psinew=prop;}else{psinew=prop-v;}}
  
  
  propvec(j)=psinew ;
  //propvec[j]= ;
  return propvec;
}


//Function to propose new bandwidth (log-modulus of AR(2) roots) for a subinterval j with k total components 
// [[Rcpp::export]]
vec NewL(const vec& L,const double& epsilon,double dim, double q){
  //VZ is vector of current weights or atoms
  //epsilon is a numerical value of uniform halfwidths for proposals
  //dim is the element of the weight or atom vector which is being proposed
  
  vec Lnew=L;
  //propose new value
  double prop=L(dim)+epsilon*(2*randu()-1);
  //take modulus 1 to ensure it falls between 0 and 1
  if(prop < 0){Lnew(dim)=q+prop;}else{if(prop < q){Lnew(dim)=prop;}else{Lnew(dim)=prop-q;}}
  return Lnew;
}

//Function to calculate BDP weights from stick break proportions V for truncation level L
// [[Rcpp::export]]
vec GetPhi(const vec& Z,const vec& P, int L, int Lambda_now){
  
  vec Phi(Lambda_now,fill::zeros);
  for(int l=0; l<L; l++){
    Phi[ floor(Z[l]*Lambda_now)  ] += P[l];
  }
  
  return Phi;
}



//function to create a new partition points depending if the size grow by 1 or decrease by 1 
// [[Rcpp::export]]
List newpartition(vec partition, int growth, int i, vec BW, vec psi, vec V,const vec& Z,const double& epsBW, int L, double q ){
  vec newpart ;
  vec psiinit;
  vec BWinit;
  vec P(L);
  P=GetWeights(V, L);
  int N =partition.n_elem;
  if(N==2){i=0;}
  if(i==N-2&&growth!=1&&N!=2){i=N-3;}
  int Lambda_now=N-1;
  double b= partition(i+1);
  double a= partition(i);
  vec phi= GetPhi( Z, P,  L, Lambda_now);
  double newpoint=  randu()*(b-a )+a   ;
  vec np(1);
  vec inter(1);
  np(0)=newpoint;
  if (growth==1){//birth
    newpart=join_cols(partition.subvec(0,i), np);
    newpart=join_cols(newpart,  partition.subvec(i+1,N-1));
    if(N!=2){
      if(psi(i)<newpoint ){
        if(i==N-2){
          inter(0)=(newpoint +partition(i+1))/2.0 ;
          psiinit=   join_cols( psi, inter);
          
          
          psi=Newpsij( psiinit ,  (newpart(i+1)-newpart(i))/10.0 , i, newpart );
          
          
          
          inter(0)=BW(i) ;
          BWinit=join_cols(BW ,inter);
          BW=NewL(BWinit,   epsBW  ,i ,q);
        }else{
          inter(0)=(newpoint +partition(i+1))/2.0 ;
          psiinit=   join_cols(psi.subvec(0,i),inter);
          psiinit=   join_cols(psiinit, psi.subvec(i+1,N-2));
          psi=Newpsij( psiinit ,  (newpart(i+1)-newpart(i))/10.0 , i, newpart );
          inter(0)=BW(i) ;
          BWinit=join_cols(BW.subvec(0,i),inter);
          BWinit=join_cols(BWinit ,BW.subvec(i+1,N-2));
          BW=NewL(BWinit,   epsBW  ,i, q);
        }  
        
        
        
      }else{// case psi(i)>=newpoint 
        if(i==0){
          inter(0)=(newpoint +partition(i))/2.0 ;
          psiinit=   inter;
          psiinit=   join_cols(psiinit, psi.subvec(i,N-2));
          psi=Newpsij( psiinit ,  (newpoint-partition(i))/10.0 , i, newpart );
          inter(0)=BW(i) ;
          BWinit= inter;
          BWinit=join_cols(BWinit ,BW.subvec(i,N-2)  );
          BW=NewL(BWinit,   epsBW  ,i , q);
        }else{
          inter(0)=(newpoint +partition(i))/2.0 ;
          psiinit=   join_cols(psi.subvec(0,i-1),inter);
          psiinit=   join_cols(psiinit, psi.subvec(i,N-2));
          psi=Newpsij( psiinit ,  (newpoint-partition(i) )/10.0 , i, newpart );
          inter(0)=BW(i) ;
          BWinit=join_cols(BW.subvec(0,i-1),inter);
          BWinit=join_cols(BWinit ,BW.subvec(i,N-2));
          BW=NewL(BWinit,   epsBW  ,i , q);
        }
      }
    }else{// case N=2
      if(psi(i)<newpoint ){
        
        inter(0)=(newpoint +partition(i+1))/2.0 ;
        psiinit=   join_cols(psi,inter);
        psi=Newpsij( psiinit ,  (newpart(i+2)-newpart(i+1))/10.0 , i+1, newpart );
        inter(0)=BW(i) ;
        BWinit=join_cols(BW,inter);
        BW=NewL(BWinit,   epsBW  ,i+1 , q);        
      }else{
        inter(0)=(newpoint +partition(i))/2.0 ;
        psiinit=   join_cols(inter,psi);
        psi=Newpsij( psiinit ,  (newpart(i+1)-newpart(i))/10.0 , i, newpart );
        inter(0)=BW(i) ;
        BWinit=join_cols(inter, BW);
        BW=NewL(BWinit,   epsBW  ,i , q);
        
      }
      
      
      
    }
  }else{if(N!=2){//death
    if (i==0){
      np(0)=partition(0);
      newpart =join_cols<mat>( np,partition.subvec( 2,N-1) );
      if(phi(i) ==0 &&phi(i+1) ==0 ){ inter(0)= (psi(i)+psi(i+1))/2.0;  }else{
        inter(0)= (phi(i)/(phi(i) + phi(i+1)))*psi(i) +  (phi(i+1)/(phi(i) + phi(i+1)))*psi(i+1)   ;
      }
      if(N==3){
        //   psi = inter;
        psi=  Newpsij( inter ,  (newpart(i+1)-newpart(i))/10.0 , i, newpart );
        
      }else{
        
        psi =join_cols<mat>( inter,psi.subvec( 2,N-2) );
        psi=  Newpsij( psi ,  (newpart(i+1)-newpart(i))/10.0 , i, newpart );
        
        
      }
      inter(0)= sqrt(BW(i)*BW(i+1) )   ;
      
      
      
      
      if(N==3){
        BW=inter;
        BW=NewL(inter, epsBW, i , q);
        
        
      }else{
        BW =join_cols<mat>( inter,BW.subvec( 2,N-2) );
        BW=NewL(BW, epsBW, i , q);
        
      }
    }else{if(i==N-3){
      np(0)=partition(N-1);
      newpart=join_cols<mat>(  partition.subvec(0,N-3), np  );
      if(phi(i) ==0 &&phi(i+1) ==0 ){ inter(0)= (psi(i)+psi(i+1))/2.0;  }else{
        inter(0)= (phi(i)/(phi(i) + phi(i+1)) )*psi(i) +  (phi(i+1)/(phi(i) + phi(i+1)))*psi(i+1);}
      
      psi =join_cols<mat>( psi.subvec( 0,N-4),inter  );
      psi=Newpsij(psi ,  (newpart(i+1)-newpart(i))/10.0 , i,  newpart);
      inter(0)= sqrt(BW(i)*BW(i+1) );
      
      
      BW =join_cols<mat>( BW.subvec( 0,N-4), inter );
      BW=NewL(BW, epsBW, i , q);
      
    }else{
      newpart =join_cols<mat>( partition.subvec(0,i),partition.subvec( i+2,N-1) );
      if(phi(i) ==0 &&phi(i+1) ==0 ){ inter(0)= (psi(i)+psi(i+1))/2.0;  }else{
        inter(0)= (phi(i)/(phi(i) + phi(i+1)))*psi(i) +  (phi(i+1)/(phi(i) + phi(i+1)))*psi(i+1)   ;}
      psiinit =join_cols<mat>( psi.subvec( 0,i-1),inter  );
      psi =join_cols<mat>(psiinit, psi.subvec( i+2,N-2)  );
      psi=Newpsij(psi ,  (newpart(i+1)-newpart(i))/10.0 , i,  newpart);
      
      inter(0)= sqrt(BW(i)*BW(i+1) )   ;
      BWinit =join_cols<mat>(BW.subvec( 0,i-1), inter );
      BW =join_cols<mat>(BWinit, BW.subvec( i+2,N-2) );
      BW=NewL(BW, epsBW, i , q);
      
    }
    }
  }else{newpart=partition;}
  }
  return List::create(
    _["partition"] = newpart,
    _["psi"] = psi,
    _["BW"] = BW);
}


//function to get the width  of a subinterval 
//used to compute uniform densities on the M-H steps
// [[Rcpp::export]]
double Getwidth(vec partition , int i){
  double width= partition[i+1]- partition[i];
  return  width;
}


//function to sample number of components, Lambda, using MH step
// in a death-birth process
// [[Rcpp::export]]
List sample_Lambda(int N, int L, const vec& omegas,const vec& psi_cur,const vec&  BW_cur, const vec& Pgrams, const vec& P_cur, const vec& V_cur,
                   const vec& Z_cur, const double& epsBW, int Lambda_cur, string Lambda_prior, int Lambda_init, 
                   int Lambda_max, double MaxFreq,  vec partition, double q){
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
  // i is the subinterval chosen to make the partition or the death 
  
  List  newpart;
  int M=partition.n_elem-2;
  int Lambda_new;
  int Lambda_prop = Lambda_cur;
  double accept_Lambda = 0;
  int growth;
  int i;
  vec   newp=partition;
  vec   newpsi=psi_cur;
  vec   newBW=BW_cur;
  double width;
  double subwidth;
  double sign;
  
  
  
  if(Lambda_cur==1){Lambda_prop=2;
    growth=1;
    i=0;
    newpart= newpartition( partition,  growth,  i  ,  BW_cur,  psi_cur,  V_cur, Z_cur, epsBW ,  L ,q);
    
    // 0-.5 is the whole interval and set a new point 
    //  int K=randi<int>( distr_param(0, N-2));
    vec psinew=newpart[1];
    vec np=newpart[0];
    width=Getwidth(partition,i);
    if(psi_cur(i)== psinew(i)){
      subwidth= Getwidth(np, i+1 );}else{
        subwidth= Getwidth(np, i );
      }
      
      if(np(1)<psi_cur(i) ){
        sign= -log(np(1))- log(2);
        
      }else{
        sign=-log(.5-np(1)) - log(2) ;
      }
      
      
  }
  else{
    if(Lambda_cur==Lambda_max){Lambda_prop=Lambda_max-1;
      growth=0;
      i=randi<int>( distr_param(0, M));
      newpart= newpartition( partition,  growth,  i  ,  BW_cur,  psi_cur, V_cur ,Z_cur, epsBW  ,  L, q );
      
      vec psinew=newpart[1];
      vec np=newpart[0];
      width=Getwidth(partition,i) +Getwidth(partition,i+1) ;
      if(psi_cur(i)< partition(i+1)){
        subwidth= Getwidth(partition, i+1 );}else{
          subwidth= Getwidth(partition, i );
        }
        if(partition(i+1)<psinew(i+1) ){
          sign= -log( Getwidth(partition, i )  )- log(2) ;
          
        }else{
          sign=-log(Getwidth(partition, i+1 )   ) - log(2) ;
        }
        
    }
    
    else{
      if(randu()<0.5){Lambda_prop=Lambda_cur-1;
        growth=0;
        sign=0;
        i=randi<int>( distr_param(0, M));
        newpart= newpartition( partition, growth,  i  ,  BW_cur,  psi_cur, V_cur, Z_cur ,epsBW,  L ,q );
        
        vec psinew=newpart[1];
        vec np=newpart[0];
        width=Getwidth(partition,i) +Getwidth(partition,i+1) ;
        if(psi_cur(i)< partition(i+1)){
          subwidth= Getwidth(partition, i+1 );}else{
            subwidth= Getwidth(partition, i );
          }
          // cases depends on the case of i
          if(i==0){
            if(partition(i+1)<psinew(i) ){
              sign= -log( Getwidth(partition, i )  ) ;
              
            }else{
              sign=-log(Getwidth(partition, i+1 )   )  ;
            }
            
          }else{      
            if(partition(i)<psinew(i-1) ){
              sign= -log( Getwidth(partition, i-1 )  ) ;
              
            }else{
              sign=-log(Getwidth(partition, i )   )  ;
            }
          }// end else of cases of i
          
          
          
          
      }else{Lambda_prop=Lambda_cur+1;
        growth=1;
        sign=0;
        i=randi<int>( distr_param(0, M));
        newpart= newpartition( partition, growth,  i  ,  BW_cur,  psi_cur, V_cur, Z_cur , epsBW ,  L , q);
        vec psinew=newpart[1];
        vec np=newpart[0];
        width=Getwidth(partition,i);
        if(psi_cur(i)== psinew(i)){
          subwidth= Getwidth(np, i+1 );}else{
            subwidth= Getwidth(np, i );
          }
          
          if(np(i+1)<psi_cur(i) ){
            sign= log( Getwidth(np, i )  ) ;
            
          }else{
            sign=log(Getwidth(np, i+1 )   )  ;
          }
          
      }
      
    }
  }
  
  
  
  accept_Lambda = accept_Lambda +
    Whittle(omegas, newpart[1], newpart[2], Pgrams, P_cur, Z_cur, Lambda_prop, L, MaxFreq) -
    Whittle(omegas,psi_cur, BW_cur, Pgrams, P_cur, Z_cur, Lambda_cur, L, MaxFreq)+  sign  ;  // 
  
  
  //accept or reject proposed Lambda
  if(accept_Lambda > log(randu())){
    Lambda_new=Lambda_prop;
    newp=as<vec>(newpart[0]);
    newpsi=as<vec>(newpart[1]);
    newBW=as<vec>(newpart[2]);
  }else{Lambda_new=Lambda_cur; }
  //  }
  
  return List::create(
    _["Lambda_new"] = Lambda_new,
    _["partition"] = newp,
    _["psi"] = newpsi,
    _["BW"] = newBW ) ;
  
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
List sample_V(int N, int L, const vec& omegas, const vec& psi_cur, const vec& BW_cur,  const mat& Pgrams, const vec& V_cur, 
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
        Whittle(omegas,psi_cur, BW_cur, Pgrams.col(i), P_new, Z_subj, Lambda_cur, L, MaxFreq);
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
          Whittle(omegas, psi_cur, BW_cur, Pgrams.col(i), P_prop, Z_subj, Lambda_cur, L, MaxFreq);
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
vec sample_Z(int N, int L, const vec& omegas, const vec&  psi_cur, const vec&  BW_cur, const mat& Pgrams, const vec& P_cur, 
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
        Whittle(omegas, psi_cur, BW_cur, Pgrams.col(i), P_cur, Z_new, Lambda_cur, L, MaxFreq);
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
          Whittle(omegas, psi_cur, BW_cur, Pgrams.col(i), P_cur, Z_prop, Lambda_cur, L, MaxFreq);
      }
    }
    accept_Z = Z_prop_lik - Z_new_lik;
    //acceptance or reject using MH step
    if(accept_Z > log(randu())){Z_new = Z_prop; Z_new_lik = Z_prop_lik;}
  }
  return Z_new;
}



//function to sample location parameters psi,  given the number of components, using MH step and uniform proposal
// [[Rcpp::export]]
vec sample_psi(int N, int L, const vec& psi_cur, const vec& BW_cur,  const vec& omegas, const vec& Pgrams, const vec& P_cur, const vec& Z_cur,
               int Lambda_cur,  double MaxFreq, vec partition){
  //N is the component to generate a new location parameter psi 
  //L is BDP truncation level
  //omegas is a vector of fourier frequencies
  //Pgrams is a matrix of corresponding periodogram values with a column for each subject
  //P_cur is a vector containing the latest iteration BDP weights
  //psi_cur is a vector containing the latest iteration of the location parameters
  //Bwidth is a vector containing the latest iteration of the bandwidth (logmodulus) parameters
  //Lambda_cur is the current number of spectral components
  //epsilon is a vector of uniform halfwidths for proposals
  //MaxFreq is the maximum frequency to consider for the spectral density
  
  //compute weights from the P weights of the BDP sampling
  
  vec Phi=GetPhi( Z_cur, P_cur,  L,  Lambda_cur);
  double epsilon;
  vec psi_prop(Lambda_cur);
  vec psi_new=psi_cur;
  psi_prop=psi_cur;
  double psi_new_lik = 0;
  
  psi_new_lik = psi_new_lik +
    Whittle(omegas, psi_cur, BW_cur, Pgrams, P_cur, Z_cur,  Lambda_cur, L, MaxFreq);
  
  double psi_prop_lik;
  double accept_psi;  
  //loop through updating each of the Lambda spectral locations
  
    //propose new location for the l-th subinterval
    epsilon=(partition(N+1)-partition(N))/4.0;
    psi_prop=Newpsij( psi_prop , epsilon, N, partition);
    
  
  psi_prop_lik = 0;  
  //calculate acceptance probability
  
  psi_prop_lik = psi_prop_lik +
    Whittle(omegas, psi_prop ,BW_cur, Pgrams, P_cur, Z_cur, Lambda_cur, L, MaxFreq);
  
  accept_psi = psi_prop_lik - psi_new_lik;
  //acceptance or reject using MH step
  if(accept_psi > log(randu())){psi_new = psi_prop; psi_new_lik = psi_prop_lik;}
  
  return psi_new;
}


//function to sample Bandwidth (log modulus) parameters L, using MH step and uniform proposal
// [[Rcpp::export]]
vec sample_L(int N, int L, const vec& psi_cur, const vec& BW_cur,  const vec& omegas, const vec& Pgrams, const vec& P_cur, const vec& Z_cur,
             int Lambda_cur, const double& epsilon, double MaxFreq, double q){
  //N is the component to sample a new bandwidth L
  //L is BDP truncation level
  //omegas is a vector of fourier frequencies
  //Pgrams is a matrix of corresponding periodogram values with a column for each subject
  //P_cur is a vector containing the latest iteration BDP weights
  //psi_cur is a vector containing the latest iteration of the location parameters
  //BW_cur is a vector containing the latest iteration bandwidth (logmodulus) parameters
  //Lambda_cur is the current number of spectral components
  //epsilon is a vector of uniform halfwidths for proposals
  //MaxFreq is the maximum frequency to consider for the spectral density
  
  vec BW_prop(Lambda_cur);
  vec BW_new=BW_cur;
  BW_prop=BW_cur;
  double BW_new_lik = 0;
  
  BW_new_lik = BW_new_lik +
    Whittle(omegas, psi_cur, BW_new,  Pgrams, P_cur, Z_cur,  Lambda_cur, L, MaxFreq);
  
  double BW_prop_lik;
  double accept_BW;  
  //loop through updating each of the L BDP atoms, Z
 
    //propose new l-th atom
    BW_prop=NewL(BW_prop, epsilon, N, q);
  BW_prop_lik = 0;
  //calculate acceptance probability
  
  BW_prop_lik = BW_prop_lik +
    Whittle(omegas, psi_cur, BW_prop, Pgrams, P_cur, Z_cur, Lambda_cur, L, MaxFreq);
  
  accept_BW = BW_prop_lik - BW_new_lik;
  //acceptance or reject using MH step
  if(accept_BW > log(randu())){BW_new = BW_prop; BW_new_lik = BW_prop_lik;}
  
  return BW_new;
}



//function to calculate posteior curve samples from posteior variables
// [[Rcpp::export]]
mat GetCurvesBDP(int Nsamp, const vec& omegas,  const mat& psi_cur, const mat& BW_cur, const vec& Lambda_samps, 
                 const mat& w){
  //Nsamp is the number of samples
  //omegas is a vector of fourier frequencies
  //L is BDP truncation level
  //Lambda_samps is a vector of samples of the number of Bernstein polynomial components
  //P_samps is a matrix where each row is an L length sample of BDP weights
  //Z_samps is a matrix where each row is an L length sample of BDP atoms
  //MaxFreq is the maximum frequency to consider for the spectral density
  
  int Lambda_now;
  
  int nOm=omegas.n_elem;
  mat CurveEsts(Nsamp,nOm,fill::zeros);
  
  for(int s=0; s<Nsamp; s++){
    Lambda_now=Lambda_samps(s);
    
    CurveEsts.row(s)=SSpline(w.row(s).subvec(0,Lambda_now-1).t(), psi_cur.row(s).subvec(0,Lambda_now-1).t(),BW_cur.row(s).subvec(0,Lambda_now-1).t(),  omegas).t();
    
  }
  
  return CurveEsts;
}

//%******************* PickAlpha ******************************

// this function computes the posterior distribution  fot he parameter alpha from the Dirichlet Process
// [[Rcpp::export]]
double  getLogPostAlpha(double x, int iStar, double a0, double b0,int n){
  // n is the sample size , for now will be the original sample size of the time series 
  // until new revision for detemrine how it works the histogram maping VS the DFT maping 
  double pi=3.14159265358979323846;
  
  double alpha = exp(x);
  
  double h= R::beta(alpha +1, n);
  
  double logLike = (iStar-1)*log(alpha) +  log(h) + log(alpha+n) ; 
  
  double logPrior = -.5*log(2*pi)-log(b0)-  .5*(   pow( x - a0,2 )/ pow( b0,2)       )  ;
  
  double logPost = logLike + logPrior;
  
  return  logPost ;
}





//  % This function updates "alpha: based on the old alpha and number of unique
//  % theta's
// [[Rcpp::export]]
double  pickAlpha(double oldAlpha, int iStar, double a0, double b0, int n){
  //iStar is the current number of components 
  int  m = 40;
  int  w = 5;
  //   % slice sampling
  double     x = log(oldAlpha+0.001); 
  double    z = getLogPostAlpha(x, iStar, a0, b0,n) +  log(randu()) ;
  //   the exponential ischanged by the log of a a uniform
  double u = randu();
  double   L = x - w*u;
  double   R = L + w;
  double   v = randu();
  int   J = floor(m*v);
  int   K = (m-1) - J;
  
  while(  J>0 && z < getLogPostAlpha(L, iStar, a0, b0, n)){
    L = L - w;
    J = J - 1;}
  
  
  while( K>0 && z < getLogPostAlpha(R, iStar, a0, b0,n) ){
    R = R+w;
    K = K-1;}
  
  
  u = randu() ;
  double  newX = L + u*(R-L);
  
  while (z > getLogPostAlpha(newX, iStar, a0, b0,n) ){ 
    if (newX < x){
      L = newX;}
    else{
      R = newX;}
    
    u = randu();
    newX = L + u*(R-L);
  }
  
  
  
  
  return  exp(newX);
  
  
  
}  




// [[Rcpp::export]]
double  alphamixgampost( double oldeta, int oldlambda , double a0, double b0, int n){

double newalpha;  
  
double x= (a0 +oldlambda-1)/(static_cast<double>(n)*(b0-log(oldeta)  )  );
  
double pieta=x/(1+x) ;  
  
double u =randu();

if( pieta<=u ){
  newalpha= R::rgamma( a0+oldlambda,  pow(b0-log(oldeta),-1)    )  ;
}else{
  newalpha= R::rgamma(a0+oldlambda-1, pow(b0-log(oldeta),-1)   );
     }
    
  
return(newalpha)   ; 
}





//MCMC Sampling Algorithm for single subject spectral BDP
// [[Rcpp::export]]
List SpectralBDP(double  Tsize ,const vec& omegas,const vec& Pgram, int Nsamp, int L, const vec& epsilon,  const double& epsilon_BW,  
                 double SamplingRate=1, double MaxFreq=0.5, int Sup=500, 
                 int Lambda_init=30, string Lambda_prior="flat", int Lambda_max=300, 
                 string Lambda_prop_dist="poisson", bool alpha_const=true, 
                 double alpha_init=1, double a_alpha_L=1, double b_alpha_L=1  , double q=1  ){
  //X is a vector containing the time series to be analyzed
  //Nsamp is the number of posterior samples
  //L is the truncation level for the BDP
  //epsilon is a vector of uniform halfwidths for V and Z proposals
  //epsilon_BW is a vector of uniform halfwidths for Bandwidth proposals length same as lambda_init
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
  //q is the limit for the uniform distribution of the bandwidth parameters
  //determine Fourier frequencies to be considered
//  double T = X.n_elem;
//  int N_freq = floor((T-1)/2);
//  int N_freqKeep = 0;
//  vec omegas(N_freq);
//  for(int n_freq=0; n_freq < N_freq; n_freq++){
//    omegas(n_freq) = (n_freq+1)/T*SamplingRate;
//    if((n_freq+1)/T*SamplingRate < MaxFreq){N_freqKeep += 1;}
//  }
//  omegas=omegas.head(N_freqKeep);
  
  //calculate periodograms
//  vec Pgram=square(abs(fft(X)))/T;
//  Pgram=Pgram.head(N_freqKeep);
  
  
  
  
  //initialize concentration parameter alpha
  vec alpha_L(Nsamp);
  alpha_L(0)=alpha_init;
  
  //initialize BDP (from Choudhuri)
  vec Lambda(Nsamp);
  Lambda(0)=Lambda_init;
  
  mat w(Nsamp, Lambda_max, fill::zeros);
  
  mat V(Nsamp, L);
  V.submat(0,0,0,L-2).fill(0.5);
  V(0,L-1)=1;
  
  mat P(Nsamp, L);
  P.row(0)=GetWeights(V.row(0).t(), L).t();
  List NewVP(2);
  
  mat Z(Nsamp, L);
  Z.row(0).randu();
  //random start for each parameter in order to evaluate the G-R psrf 
  double eta;
  
  mat psi(Nsamp, Lambda_max, fill::zeros);
  
  mat BW(Nsamp, Lambda_max);
  BW.fill(   randu()*q   );
  //MaxFreq/ static_cast<double>(Lambda_init)
  
  mat partition(Nsamp, Lambda_max+1, fill::zeros);
  
  vec psi_init = zeros<vec>(Lambda_init)  ; 
  //  initialization of the center  psi  
  psi_init(0)= MaxFreq/ static_cast<double>(2*Lambda_init)  ;
  for(int i=1 ; i<Lambda_init;i++){
    psi_init(i)=psi_init(i-1) +  MaxFreq/ static_cast<double>(Lambda_init)   ;
  }
  
  vec partition_init= zeros<vec>(Lambda_init+1)  ; 
  
  // inititalization of the parition 
  
  for(int i=1 ; i<=Lambda_init ;i++){
    partition_init(i)=partition_init(i-1) +  (1.0/static_cast<double>(2*Lambda_init) )  ;
  }  
  
  
  // real initial value after first random movement  
  
  for(int i=0 ; i<Lambda_init;i++){
    psi(0,i) =  Newpsij(psi_init ,  static_cast<double>(1) / static_cast<double>(2*Lambda_init), i, partition_init )(i);
  }
  //keep all the parition vectors?   
  partition.row(0).subvec(0,Lambda_init)=partition_init.t();
  
  vec loglik(Nsamp);
  loglik(0)=Whittle( omegas,   psi.row(0).subvec(0,Lambda_init-1).t() , BW.row(0).subvec(0,Lambda_init-1).t() ,  Pgram,  P.row(0).t(),  Z.row(0).t(), Lambda_init, L, MaxFreq);
  
  
  for(int iter=1; iter<Nsamp; iter++){
    
    // define the vectors psi and BW according to the current dimensions of lambda
    
    //Update number of Bernstein polynomials, Lambda
    List  Birth_death     = sample_Lambda(1, L, omegas, psi.row(iter-1).subvec(0,Lambda(iter-1)-1).t() , BW.row(iter-1).subvec(0,Lambda(iter-1)-1).t(), Pgram, P.row(iter-1).t(), V.row(iter-1).t() ,  Z.row(iter-1).t(), epsilon_BW,
                                          Lambda(iter-1), Lambda_prior, Lambda_init, Lambda_max, MaxFreq,  partition.row(iter-1).subvec(0,Lambda(iter-1)).t(), q);
    
    
    
    Lambda(iter)=as<int>( Birth_death[0]);
    
    partition.row(iter).subvec(0,Lambda(iter))   =as<vec>(Birth_death[1]).t();
    
    
    // for each value then make a for
    
    for(int l=0; l<Lambda(iter); l++){
      // Update psi
    psi.row(iter).subvec(0,Lambda(iter)-1) = sample_psi(l,  L,  as<vec>(Birth_death[2]), as<vec>(Birth_death[3]) , omegas, Pgram, P.row(iter-1).t(), Z.row(iter-1).t(),
            Lambda(iter)     ,  MaxFreq,partition.row(iter).subvec(0,Lambda(iter)).t() ).t();
      //Update Bandwidths
      BW.row(iter).subvec(0,Lambda(iter)-1) = sample_L(l,  L,   psi.row(iter).subvec(0,Lambda(iter)-1).t() , as<vec>(Birth_death[3]) , omegas, Pgram, P.row(iter-1).t(), Z.row(iter-1).t(),
             Lambda(iter), epsilon_BW,  MaxFreq, q).t();
      
          }
    
    
    
    

    
    
    
    // after sampling the number of compoentns the weights are obtained again using  the same atoms and bdp weigths
    
    //the partition is proposed selecting one subinterval at random
    // then a new point is drawn randomly on the interval selected
    //  a new location location is drawn in he subpartition where there is no location paramters  psi_jump
    // a new bandwidth is proposed uniformily
    // in the death process parameters are substitued by one value drawn randomly according to the transition probabilities
    //in location the mid point of the two parameters is select and from it a new location is proposed in
    // the selected  subinterval
    
    //Update weights V and P
    NewVP = sample_V(1, L, omegas, psi.row(iter).subvec(0,Lambda(iter)-1).t(), BW.row(iter).subvec(0,Lambda(iter)-1).t(), Pgram, V.row(iter-1).t(), P.row(iter-1).t(),
                     Z.row(iter-1).t(), Lambda(iter), epsilon, alpha_L(iter-1), MaxFreq, false);
    V.row(iter) = as<vec>(NewVP("V_new")).t();
    P.row(iter) = as<vec>(NewVP("P_new")).t();
    
    //Update Z
    Z.row(iter) = sample_Z(1, L, omegas, psi.row(iter).subvec(0,Lambda(iter)-1).t(),BW.row(iter).subvec(0,Lambda(iter)-1).t(),  Pgram, P.row(iter).t(), Z.row(iter-1).t(),
          Lambda(iter), epsilon, MaxFreq).t();
    
    //Dirichlet process Update alpha_L via Gibbs step
    
    // modofication gibss build from gamma prior west 1995 
    eta= R::rbeta( alpha_L(iter-1)+1, Tsize);
    
    
    if(alpha_const==false){
      //   alpha_L(iter)=my_rgamma(a_alpha_L+L-1, b_alpha_L - log(P(iter,L-1))    )[0] ;
      //this update is based on IOshwaran and James 2002 , 
      // based on the truncated DP approximation 
      // a\pha | V /sim GA( L +a -1, b- log pm)  pm is the last weight of the truncated atoms 
      
      // update below is based on slice sapling and the gibbs sampler from west and escobar  (1995 ), method from prof shahbaba
    //  alpha_L(iter)=pickAlpha(alpha_L(iter-1), Lambda(iter) ,  a_alpha_L  , b_alpha_L, Tsize  );
    alpha_L(iter)=alphamixgampost( eta, Lambda(iter), a_alpha_L, b_alpha_L, Tsize);
      
    }
    if(alpha_const==true){
      alpha_L(iter)=alpha_L(iter-1);
    }
    
   // if((iter-1)/Sup!=iter/Sup){cout << iter << " posterior samples complete" << endl;}
    
    // retrieve the weights
    w.row(iter).subvec(0,Lambda(iter)-1)= GetPhi(Z.row(iter).t(),P.row(iter).t(),L,Lambda(iter) ).t();
    
    
    loglik(iter)=Whittle( omegas,  psi.row(iter).subvec(0,Lambda(iter)-1).t() , BW.row(iter).subvec(0,Lambda(iter)-1).t() ,  Pgram,  P.row(iter).t(),  Z.row(iter).t(), Lambda(iter), L, MaxFreq);  
    
    
  }
  
  // change o return only the curves as a matrix 
  
  mat Curves=GetCurvesBDP(Nsamp, omegas,psi,BW,   Lambda, w);
  

  // return only effective samples
  return List::create(
    _["Lambda"] = Lambda,
    _["weights"] = w,
    _["psi"] = psi,
    _["BW"] = BW,
    _["CurveEsts"] = Curves
    ,_["omegas"] = omegas,
    _["Pgram"] = Pgram ,  
   _["alpha"] = alpha_L, 
   _["Whittle"] = loglik, 
   _["Tsize"] = Tsize,
   _["TruncL"] = L,
   _["epsilonvec"] = epsilon,
   _["epsilonBW"] = epsilon_BW,
   _["SamplingRate"] = SamplingRate,
   _["MaxFreq"] = MaxFreq,
   _["Lambda_max"] = Lambda_max,
   _["alpha_const"] = alpha_const,
   _["a_alpha_L"] = a_alpha_L,
   _["b_alpha_L"] = b_alpha_L,
   _["q"] = q
  );
  
  
}












