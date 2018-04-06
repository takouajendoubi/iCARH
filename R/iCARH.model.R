#' @title Runs the integrative CAR Horseshoe model
#'
#' @description Infers treatment effects, association with heterogenous omic variables, pathway perturbation
#' among other parameters (e.g. time dependance). Regression coefficients (beta parameter) are initialised
#' using a univariate regression ignoring time and metabolite dependence.
#'
#' @param X the metabolomics time-course data
#' @param Y the additional omic time-course data
#' @param drug treatment effect (NA values not allowed in drug)
#' @param pathways pathway adjacency matrices
#' @param tau global sparsity parameter \eqn{\tau} as in Jendoubi, T., & Ebbels, T. (2018)
#' @param NA_value NA values are incompatible with stan.
#' NAs will be replaced by NA_value and will be inferred (only for X and Y data).
#' @param ... additional stan parameters
#'
#' @return stan object
#'
#' @examples data.sim = iCARH.simulate(4, 8, 10, 2, 2, path.probs=0.3, Zgroupeff=c(0,4),
#' beta.val=c(1,-1,0.5, -0.5))
#' XX = data.sim$XX
#' Y = data.sim$Y
#' Z = data.sim$Z
#' pathways = data.sim$pathways
#' fit = iCARH.model(XX, Y, Z, pathways, control = list(adapt_delta = 0.99, max_treedepth=10),
#' iter = 500, chains = 2)
#'
#' @export iCARH.model

iCARH.model = function(X, Y, drug, pathways, tau=1.2, NA_value=-99999, ...){
  x_i_mis = which(is.na(X), arr.ind = T);
  x_i_obs = which(!is.na(X), arr.ind = T);
  X[which(is.na(X))] =  NA_value;
  y_i_mis = which(is.na(Y), arr.ind = T);
  y_i_obs = which(!is.na(Y), arr.ind = T);
  Y[which(is.na(Y))] =  NA_value;
  adjmat = lapply(pathways, function(x) 1/(x+1)) #To check
  adjmat = lapply(adjmat, function(x) {diag(x)=0; 1/max(rowSums(x>0))*x})
  lambdas = lapply(adjmat,function(x) sort(eigen(x)$values)[c(1,nrow(x))]);
  regression.stan = "
  data {
  int<lower=0> J;
  int<lower=0> T;
  int<lower=0> N;
  int<lower=0> P;
  int<lower=0> K;
  int<lower=0> x_n_mis;
  int<lower=0> y_n_mis;
  matrix[N,J]  X[T];
  int<lower=1> x_i_mis[x_n_mis,3];
  int<lower=1> y_i_mis[y_n_mis,3];
  matrix[N,K]  Y[T];
  vector[N]  drug[T];
  matrix<lower=0>[J,J] adjmat[P];
  vector[2] lambdas[P];
  real<lower=1> nu;
  real  NA_value;
  }
  transformed data{
  matrix[J,J] I;
  I = diag_matrix(rep_vector(1,J));
  }
  parameters {
  real Xmis[x_n_mis];
  real Ymis[y_n_mis];
  real<lower=-1, upper=1> alpha[J];
  vector[K] gam_std[J];
  real beta[J];
  real inter[J]; // fixed intercept
  real<lower=0> sigma;
  real<lower=0> sigmae[J];
  real<lower=0> sigma_y;
  vector[N] err[J];
  vector<lower=0, upper=1>[2] phi_std[P];
  vector<lower=0>[J] nu1_g;
  vector<lower=0>[J] nu2_g;
  vector<lower=0>[K] nu1_l[J];
  vector<lower=0>[K] nu2_l[J];
  }
  transformed parameters{
  matrix[N,J] XX[T];
  matrix[N,J] Xm[T];
  matrix[N,K] YY[T];
  matrix[N,J] mu[T];
  vector[N] intern[J]; // random intercept
  vector[K] gam[J];
  vector[2] phi[P];
  cov_matrix[J] Sigma[2];
  vector<lower=0>[K] theta[J];
  vector<lower=0>[J] sigma_gam;
  for (i in 1:T){
  for(j in 1:J) {for(n in 1:N) if((X[i,n,j])!= NA_value) XX[i,n,j]=X[i,n,j];}
  }
  for(i in 1:T){
  for(k in 1:K) {for(n in 1:N) if((Y[i,n,k])!= NA_value) YY[i,n,k]=Y[i,n,k];}
  }
  for(l in 1:x_n_mis){
  XX[x_i_mis[l,1], x_i_mis[l,2], x_i_mis[l,3]] = Xmis[l];
  }
  for(l in 1:y_n_mis){
  YY[y_i_mis[l,1], y_i_mis[l,2], y_i_mis[l,3]] = Ymis[l];
  }
  Sigma[1] = rep_matrix(0,J,J);
  Sigma[2] = rep_matrix(0,J,J);
  for(p in 1:P){
  phi[p,1] = 1/(lambdas[p,1])+0.005 + (1/(lambdas[p,2]) - 1/(lambdas[p,1])-0.005) * phi_std[p,1];
  phi[p,2] = 1/(lambdas[p,1])+0.005 + (1/(lambdas[p,2]) - 1/(lambdas[p,1])-0.005) * phi_std[p,2];
  Sigma[1] = Sigma[1] + phi[p,1]*adjmat[p];
  Sigma[2] = Sigma[2] + phi[p,2]*adjmat[p];
  }
  Sigma[1] = (I-Sigma[1]/P)/sigma;
  Sigma[2] = (I-Sigma[2]/P)/sigma;
  sigma_gam = nu1_g .* sqrt(nu2_g);
  for(j in 1:J){
  intern[j] = inter[j] + sqrt(sigmae[j])*err[j];
  theta[j] = nu1_l[j] .* sqrt(nu2_l[j]);
  gam[j] = gam_std[j].*theta[j] * sigma_gam[j];
  }
  for (i in 1:(T)){ for(j in 1:J){ Xm[i,1:N,j] = intern[j] + YY[i]*(gam[j]) +beta[j]*drug[i];}}
  mu[1] = Xm[1];
  for (i in 2:(T)){ for(j in 1:J){ mu[i,1:N,j] = Xm[i,1:N,j] + alpha[j]*(XX[i-1,1:N,j]-Xm[i-1,1:N,j]);}}
  }
  model {
  //priors on pathway significance
  for(p in 1:P) phi_std[p,1] ~ beta(0.5, 0.5);
  for(p in 1:P) phi_std[p,2] ~ beta(0.5, 0.5);
  //variance
  sigma ~ inv_gamma(N*T/4,N*T/4-1);
  nu1_g ~ normal(0.0, 1.0);
  nu2_g ~ inv_gamma(0.5, 0.5);
  sigma_y ~ inv_gamma(1,1);
  sigmae ~ inv_gamma(1,0.1);
  for(j in 1:J){
  err[j] ~ normal(0,1);
  nu1_l[j] ~ normal(0.0, 1.0);
  nu2_l[j] ~ inv_gamma(0.5*nu, 0.5*nu);
  gam_std[j] ~ normal(0, 1);
  }
  // Inferring missing values in YY
  for(k in 1:K) YY[1,,k] ~ multi_normal(rep_vector(0,N),diag_matrix(rep_vector(1,N)));
  for(i in 2:T){
  for(k in 1:K) target += multi_normal_lpdf(YY[i,,k] | YY[i-1,,k], diag_matrix(rep_vector(sqrt(sigma_y),N)));
  }
  for (i in 1:(T)){
  for(n in 1:N){
  target += multi_normal_prec_lpdf(XX[i,n] | mu[i,n], Sigma[n<(N/2)?1:2]);
  }
  }
  }
  generated quantities {
  real log_lik[T,N];
  matrix[J,N] mupred[T];
  for (i in 1:T){
  for(n in 1:N) log_lik[i,n] = multi_normal_prec_lpdf(XX[i,n] | mu[i,n], Sigma[n<(N/2)?1:2]);
  for(n in 1:N) mupred[i,,n] = multi_normal_rng((mu[i,n])', Sigma[n<(N/2)?1:2]);
  }
  }"
  initf = function(){
    #if(perturb)
    X[which(X== NA_value)] = NA
    Y[which(Y== NA_value)] = NA
    X = X + array(rnorm(prod(dim(X)),0,1), dim=dim(X))
    coeff = array(dim=c(3,dim(X)[3], dim(Y)[3]))
    for(k in 1:dim(Y)[3]){
    for(i in 1:dim(X)[3]){
      dd = lm(X~.-1,data=data.frame(X=as.vector(X[2:nrow(X),,i]), as.vector(X[1:(nrow(X)-1),,i]),
                                    as.vector(Y[1:(nrow(X)-1),,k]), as.vector(drug[1:(nrow(X)-1)])))
      coeff[,i,k] = dd$coefficients
    }
    }
    return(list(alpha=rowMeans(coeff[1,,]),gam=coeff[2,,],beta=rowMeans(coeff[3,,])))
  }
dat = list(N=ncol(X), T=nrow(X), J=dim(X)[3], P=length(pathways), K=dim(Y)[3], X=X, Y=Y,
           drug=drug, lambdas=lambdas, adjmat=adjmat, x_i_mis=x_i_mis, y_i_mis=y_i_mis,
           x_n_mis=nrow(x_i_mis), y_n_mis=nrow(y_i_mis), nu=tau,  NA_value= NA_value)
fit = stan(model_name="iCARH",model_code = regression.stan, data = dat, init=initf, ...)
return(fit)
  }
