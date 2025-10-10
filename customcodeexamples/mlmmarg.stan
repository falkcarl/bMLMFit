// MLM using marginal likelihood
/*
 Copyright 2024-2025 Carl F. Falk
 
 This program is free software: you can redistribute it and/or
 modify it under the terms of the GNU General Public License as
 published by the Free Software Foundation, either version 3 of
 the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.
 <http://www.gnu.org/licenses/>
*/

functions{
  
}

data {
  int n; // total number of observations
  int nclust; // number of clusters
  int nfixef; // number of fixed effects
  int nranef; // number of random effects parameters
  int p; // dimensions of random effects matrix
  array[nclust] int clustsizes; // size of each cluster
  
  row_vector[n] y; // response
  matrix[n, nfixef] Xmat; // Fixed effects design matrix
  matrix[n, p*nclust] Zmat; // Random effects design matrix, TODO: sparse matrix?
}

parameters {

  // possible TODO: cholesky decomp of random effects matrix
  real <lower=0> sig; // error sd
  matrix[nfixef,1] betas; // fixed effects
  corr_matrix[p] gsmat; // random effects correlation matrix
  vector<lower=0>[p] gssig; // scaling of random effects
}

transformed parameters {

 // actual random effects covariance matrix
 cov_matrix[p] Tau;
 Tau = quad_form_diag(gsmat, gssig);

}

model {
  
  // TODO: make some matrix computations more efficient?
  matrix[n, 1] XB;
  int idx = 1;
  int idxzmat = 1;
  
  // priors - try to mimic brms
  sig ~ student_t(3, 0, 2.5);
  betas[,1] ~ uniform(-10000, 10000); // is that flat enough?
  
  // random effects
  gsmat ~ lkj_corr(1); // LKJ prior for correlation matrix
  gssig ~ student_t(3, 0, 2.5); // half-t
  
  XB = Xmat*betas;

  for(j in 1:nclust){
    y[idx:(idx+clustsizes[j]-1)] ~ multi_normal(XB[idx:(idx+clustsizes[j]-1),1], quad_form_sym(Tau, Zmat[idx:(idx+clustsizes[j])-1,idxzmat:(idxzmat+p-1)]') + diag_matrix(rep_vector(square(sig),clustsizes[j])));  
    idx += clustsizes[j];
    idxzmat += p;
  }

}

generated quantities{
  
  // TODO: is there a way to avoid repeating some of these computations
  int idx = 1;
  int idxzmat = 1;
  matrix[n, 1] yhat;
  yhat = Xmat*betas;
  
  vector[nclust] log_lik;
  for (j in 1:nclust) {
    log_lik[j] = multi_normal_lpdf(y[idx:(idx+clustsizes[j]-1)] | yhat[idx:(idx+clustsizes[j]-1),1], quad_form_sym(Tau, Zmat[idx:(idx+clustsizes[j])-1,idxzmat:(idxzmat+p-1)]') + diag_matrix(rep_vector(square(sig),clustsizes[j])));
    idx += clustsizes[j];
    idxzmat += p;
  }

}
