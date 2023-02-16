#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <progress.hpp>
#include <progress_bar.hpp>




// [[Rcpp::export]]
mat thetaUpdate(double rho, int corstr_, 
                vec z, vec id, vec id_unique, mat B, vec varY_sqrt, 
                int nthread = 1) {
  
  // vec id_unique = unique(id);
  int n = id_unique.size();
  int ntheta = B.n_cols;
  mat XtXY = zeros<mat>(ntheta, ntheta+1);
  //Rcout << "The sample size n : " << n << "\n";
  
  
  // loop over n subjects
  // Rcout << "Calculated theta update for " << n << " patients:\n";
  
// #ifdef _OPENMP
//   if ( nthread > 0 )
//     omp_set_num_threads( nthread );
// #endif  
//   
//   Progress p(n, true);
//   
// #pragma omp parallel for schedule(dynamic)  
  
  for(int i = 0; i < n; i++){
    
    // if ( ! Progress::check_abort() ) {
    //   p.increment();  // update progress
      
      //if (i % 10000 == 0) Rcout << "Calculated gradian and jacobian for " << i << " out of " << n << " patients" << '\n';
      
      uvec idxi = find(id == id_unique[i]); // Find indices
      int ni = idxi.size();
      vec zi = z.elem(idxi);
      mat Bi = B.rows(idxi);
      vec varYi_sqrt = varY_sqrt.elem(idxi);
      mat Ai_sqrt = diagmat(varYi_sqrt);
      
      mat Ri = eye<mat>(ni, ni);
      
      if(corstr_ == 2) { // exchangeable
        
        Ri = ones<mat>(ni, ni) * rho + eye<mat>(ni, ni) * (1-rho);
        
      }else if(corstr_ == 3) { // ar1
        
        vec v1 = linspace<vec>(1, ni, ni); // +1 to get the end value
        mat t1 = abs(log(exp(v1) * exp(-v1.t())));
        Ri = exp(log(rho) * t1);
        
      }
      mat Vi = Ai_sqrt % Ri % Ai_sqrt;
      //if (eig_sym(Vi).min()<0){
      //  stop("Vi is singular.");
      //}
      mat Vi_inv = inv(Vi + 1e-6 * eye<mat>(ni, ni) ); 
      
      // jacobiani[1,1] ----
      mat XtXYi = Bi.t() * Vi_inv * join_horiz(Bi, zi);
      //Rcout << gradn << "\n";
      //Rcout << gradni << "\n";
      XtXY += XtXYi;
      
    // }
    
  }// end loop over n subjects
  
  return XtXY;
}


// [[Rcpp::export]]
double thetaCorr(vec err, vec id, vec id_unique,
                 int corstr_, int nthread = 1) {
  
  // vec id_unique = unique(id);
  int n = id_unique.size();
  double numerator = 0;
  double denominator1 = 0;
  double denominator2 = mean(square(err));
  
  //Rcout << "The sample size n : " << n << "\n";
  
  for(int i = 0; i < n; i++){
    uvec idxi = find(id == id_unique[i]); // Find indices
    int ni = idxi.size();
    vec erri = err.elem(idxi);
    
    if(corstr_ == 2) { // exchangeable
      
      for (int j = 1; j < ni; j++) {
        for (int k = 0; k < j; k++) {
          numerator += erri[j] * erri[k];
          denominator1 += 1;
        }
      }
      
    }else if(corstr_ == 3) { // ar1
      
      for (int j = 1; j < ni; j++) {
        numerator += erri[j] * erri[j-1];
        denominator1 += 1;
      }
      
    }
  }// end loop over n subjects
  
  double rho = numerator / denominator1  / denominator2;
  return rho;
}

// [[Rcpp::export]]
mat thetaUpdate_cen(double rho_sts, double rho_lts, int corstr_, 
                    vec z, vec id, vec id_unique, mat B_sts, mat B_lts, 
                    vec varY_sts_sqrt, vec varY_lts_sqrt, int nthread = 1) {
  
  // vec id_unique = unique(id);
  int n = id_unique.size();
  int ntheta_sts = B_sts.n_cols;
  int ntheta_lts = B_lts.n_cols;
  
  mat XtXY = zeros<mat>(ntheta_sts+ntheta_lts, ntheta_sts+ntheta_lts+1);
  //Rcout << "The sample size n : " << n << "\n";
  
  // loop over n subjects
  //Rcout << "Calculated theta update for " << n << " patients:\n";
  
// #ifdef _OPENMP
//   if ( nthread > 0 )
//     omp_set_num_threads( nthread );
// #endif  
//   
//   Progress p(n, true);
//   
// #pragma omp parallel for schedule(dynamic)  
  
  for(int i = 0; i < n; i++){
    
    // if ( ! Progress::check_abort() ) {
    //   p.increment();  // update progress
      
      //if (i % 10000 == 0) Rcout << "Calculated gradian and jacobian for " << i << " out of " << n << " patients" << '\n';
      
      uvec idxi = find(id == id_unique[i]); // Find indices
      int ni = idxi.size();
      vec zi = z.elem(idxi);
      mat Bi_sts = B_sts.rows(idxi);
      mat Bi_lts = B_lts.rows(idxi);
      vec varYi_sts_sqrt = varY_sts_sqrt.elem(idxi);
      vec varYi_lts_sqrt = varY_lts_sqrt.elem(idxi);
      mat Ai_sts_sqrt = diagmat(varYi_sts_sqrt);
      mat Ai_lts_sqrt = diagmat(varYi_lts_sqrt);
      
      
      mat Ri_sts = eye<mat>(ni, ni);
      mat Ri_lts = eye<mat>(ni, ni);
      
      if(corstr_ == 2) { // exchangeable
        
        Ri_sts = ones<mat>(ni, ni) * rho_sts + eye<mat>(ni, ni) * (1-rho_sts);
        Ri_lts = ones<mat>(ni, ni) * rho_lts + eye<mat>(ni, ni) * (1-rho_lts);
        
      }else { // ar1
        
        vec v1 = linspace<vec>(1, ni, ni); // +1 to get the end value
        mat t1 = abs(log(exp(v1) * exp(-v1.t())));
        Ri_sts = exp(log(rho_sts) * t1);
        Ri_lts = exp(log(rho_lts) * t1);
        
      }
      mat Vi_sts = Ai_sts_sqrt % Ri_sts % Ai_sts_sqrt;
      mat Vi_lts = Ai_lts_sqrt % Ri_lts % Ai_lts_sqrt;
      //if (eig_sym(Vi).min()<0){
      //  stop("Vi is singular.");
      //}
      mat Vi_sts_inv = inv(Vi_sts + 1e-6 * eye<mat>(ni, ni)); 
      mat Vi_lts_inv = inv(Vi_lts + 1e-6 * eye<mat>(ni, ni)) ; 
      
      
      // jacobiani[1,1] ----
      mat XtXYi2 = join_vert(Bi_sts.t() * Vi_sts_inv, Bi_lts.t() * Vi_lts_inv) ;
      mat XtXYi1 = join_horiz(join_horiz(Bi_sts,Bi_lts), zi);
      
      mat XtXYi = XtXYi2 * XtXYi1;
      // Rcout << XtXYi.n_cols << XtXYi.n_rows << "\n";
      // mat XtXYi = join_vert(Bi_sts.t() * Vi_sts_inv, Bi_lts.t() * Vi_lts_inv) * join_horiz(join_horiz(Bi_sts,Bi_lts), zi);
      //Rcout << gradn << "\n";
      //Rcout << gradni << "\n";
      XtXY += XtXYi;
      
    // }
    
  }// end loop over n subjects
  
  return XtXY;
}


// [[Rcpp::export]]
mat thetaVariance(double rho_sts, double rho_lts, int corstr_,
                  vec id, vec id_unique, mat B_sts, mat B_lts, mat gradn_beta, mat deriv_beta, 
                  vec varY_sts_sqrt, vec varY_lts_sqrt,//mat surv_gradn, 
                  vec err, int ntheta, int nparam,//int ntheta_sts, int ntheta_lts, 
                  int nthread = 1) {
  
  // vec id_unique = unique(id);
  int n = id_unique.size();
  mat XtXY = zeros<mat>(ntheta*2+nparam,ntheta+nparam);
  uvec ii(1);
  //Rcout << "The sample size n : " << n << "\n";
  
  
  // loop over n subjects
  // Rcout << "Calculated variance estimate for " << n << " patients:\n";
  
// #ifdef _OPENMP
//   if ( nthread > 0 )
//     omp_set_num_threads( nthread );
// #endif
//   
//   Progress p(n, true);
//   
// #pragma omp parallel for schedule(dynamic)
  
  for(int i = 0; i < n; i++){
    
    // if ( ! Progress::check_abort() ) {
    //   p.increment();  // update progress
      
      //if (i % 10000 == 0) Rcout << "Calculated gradian and jacobian for " << i << " out of " << n << " patients" << '\n';
      
      uvec idxi = find(id == id_unique[i]); // Find indices
      int ni = idxi.size();
      vec err_i = err.elem(idxi);
      ii(0)=i;
      
      mat Bi_sts = B_sts.rows(idxi);
      mat Bi_lts = B_lts.rows(idxi);
      mat deriv_betai = deriv_beta.rows(idxi);
      vec varYi_sts_sqrt = varY_sts_sqrt.elem(idxi);
      vec varYi_lts_sqrt = varY_lts_sqrt.elem(idxi);
      mat Ai_sts_sqrt = diagmat(varYi_sts_sqrt);
      mat Ai_lts_sqrt = diagmat(varYi_lts_sqrt);
      
      mat Ri_sts = eye<mat>(ni, ni);
      mat Ri_lts = eye<mat>(ni, ni);
      
      if(corstr_ == 2) { // exchangeable
        
        Ri_sts = ones<mat>(ni, ni) * rho_sts + eye<mat>(ni, ni) * (1-rho_sts);
        Ri_lts = ones<mat>(ni, ni) * rho_lts + eye<mat>(ni, ni) * (1-rho_lts);
        
      }else { // ar1
        
        vec v1 = linspace<vec>(1, ni, ni); // +1 to get the end value
        mat t1 = abs(log(exp(v1) * exp(-v1.t())));
        Ri_sts = exp(log(rho_sts) * t1);
        Ri_lts = exp(log(rho_lts) * t1);
        
      }
      //Rcout << Ri_sts << '\n';
      mat Vi_sts = Ai_sts_sqrt % Ri_sts % Ai_sts_sqrt;
      mat Vi_lts = Ai_lts_sqrt % Ri_lts % Ai_lts_sqrt;
      mat Vi_sts_inv = inv(Vi_sts+ 1e-6 * eye<mat>(ni, ni)); // 
      mat Vi_lts_inv = inv(Vi_lts+ 1e-6 * eye<mat>(ni, ni));
      
      // // jacobiani[1,1] ----
      mat breadi2 = join_vert(Bi_sts.t() * Vi_sts_inv, Bi_lts.t() * Vi_lts_inv);
      mat tmp = join_horiz(Bi_sts,Bi_lts);
      mat breadi1 = join_horiz(tmp,deriv_betai);//
      mat breadi = breadi2 * breadi1;
      mat meati_sqrt = breadi2 * err_i; 
      mat gradn_betai = gradn_beta.cols(ii);
      meati_sqrt = join_cols(meati_sqrt, gradn_betai);
      mat meati = meati_sqrt * meati_sqrt.t();
      // Rcout << "meati " << meati.n_cols << " " << meati.n_rows << "\n";
      mat XtXYi = join_cols(breadi, meati);
      // Rcout << "XtXYi " << XtXYi.n_cols << " " << XtXYi.n_rows << "\n";
      XtXY += XtXYi;
    // }
    
    
  }// end loop over n subjects
  
  return XtXY;
}
// }



