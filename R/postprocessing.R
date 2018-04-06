#' @title Postprocess and plot model parameters
#'
#'@description Group of functions to postprocess and plot model parameters of interest, compute WAIC
#' (Watanabe-Akaike Information Criterion) and MADs (Mean Absolute Deviation) for posterior predictive checks
#' and check normality assumptions.
#'
#' @describeIn iCARH.plotBeta Plot boxplots of posterior densities of \eqn{\beta} coefficients.
#'
#' @param X the metabolomics time-course data
#' @param fit stan object
#' @param path.names pathway names
#'
#' @return the \code{iCARH.plot[*]} functions return a ggplot graph object. \code{iCARH.checkNormality} returns the normalised data.
#' \code{iCARH.waic} and \code{iCARH.mad} return corresponding waic (scalar) and mad (vector of \eqn{J*(J+1)/2}) values
#'
#' @examples data.sim = iCARH.simulate(4, 20, 25, 3, 10, path.probs=0.3, Zgroupeff=c(0,4),
#' beta.val=c(1,-1,0.5, -0.5))
#' XX = data.sim$XX
#' Y = data.sim$Y
#' Z = data.sim$Z
#' pathways = data.sim$pathways
#' fit = iCARH.model(XX, Y, Z, pathways, control = list(adapt_delta = 0.99, max_treedepth=10),
#' iter = 1500, chains = 2)
#' gplot = iCARH.plotBeta(fit)
#'
#' @export iCARH.plotBeta
#'
#' @importFrom ggplot2 aes

iCARH.plotBeta = function(fit){
  X=value=NULL
  gam1= extract(fit, inc_warmup = FALSE, pars="gam")$gam
  attr(gam1,"dimnames")[[2]] = paste("X",1:dim(gam1)[2])
  attr(gam1,"dimnames")[[3]] = paste("Y",1:dim(gam1)[3])
  gamdf1 = melt(gam1, varnames = c("mcmc","X","Y"))
  gamdf1$X = as.factor(gamdf1$X)
  gg=ggplot(data=gamdf1, aes(y = value, x = X)) +
    facet_wrap(~Y,nrow = 2, scales="free_y") +
    geom_boxplot(outlier.size = NA, fill="grey") + geom_hline(yintercept = 0) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size=10, face="bold"),
          axis.text.y = element_text(size=15, face="bold"),
          axis.title=element_text(size=20,face="bold"),
          strip.text.x = element_text(size = 20, face="bold"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    ylab(expression(bold(beta)))
  return(gg)
}

#' @describeIn iCARH.plotBeta Plot boxplots of posterior densities of treatment effect coefficients.
#' @export iCARH.plotTreatmentEffect

iCARH.plotTreatmentEffect = function(fit){
  X=value=NULL
  gam1= extract(fit, inc_warmup = FALSE, pars="beta")$beta
  attr(gam1,"dimnames")[[2]] = paste("X",1:dim(gam1)[2])
  gamdf1 = melt(gam1, varnames = c("mcmc","X"))
  gamdf1$X = as.factor(gamdf1$X)
  gg=ggplot(data=gamdf1, aes(y = value, x = X)) +
    geom_boxplot(outlier.size = NA, fill="grey") + geom_hline(yintercept = 0) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size=10, face="bold"),
          axis.text.y = element_text(size=15, face="bold"),
          axis.title=element_text(size=20,face="bold"),
          strip.text.x = element_text(size = 20, face="bold"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    ylab(expression(bold(beta)))
  return(gg)
}

#' @describeIn iCARH.plotBeta Plot posterior densities of pathway perturbation paramaters
#' @export iCARH.plotPathwayPerturbation


iCARH.plotPathwayPerturbation = function(fit, path.names){
  loc=dens=treatment=pathway=NULL
  phi= extract(fit, inc_warmup = FALSE, pars="phi")$phi
  P = dim(phi)[2]
  phi.diff = apply(phi[,,1],2,function(x) density(x)$x)
  colnames(phi.diff) = c(1:P)
  phi.diff = melt(phi.diff, varnames = c("ind","pathway"), value.name = "loc")
  phi.diff$dens = melt(apply(phi[,,1],2,function(x) density(x)$y))$value
  phi.diff$treatment = "Controls"
  phi.diff.u = apply(phi[,,2],2,function(x) density(x)$x)
  colnames(phi.diff.u) = c(1:P)
  phi.diff.u = melt(phi.diff.u, varnames = c("ind","pathway"), value.name = "loc")
  phi.diff.u$dens = melt(apply(phi[,,2],2,function(x) -density(x)$y))$value
  phi.diff.u$treatment = "Cases"
  for(i in 1:nrow(phi.diff)){phi.diff$dens[i] <- phi.diff$dens[i] +phi.diff$pathway[i]}
  for(i in 1:nrow(phi.diff.u)){phi.diff.u$dens[i] <- phi.diff.u$dens[i] +phi.diff.u$pathway[i]}
  pathway.names = path.names
  names(pathway.names) = 1:P
  gg=ggplot(data=rbind(phi.diff, phi.diff.u), aes(y = loc, x = dens, fill=treatment,
                                                  group=interaction(as.factor(pathway),treatment))) +
    geom_polygon(colour="black") +
    facet_grid(~pathway, switch = "both", labeller = as_labeller(pathway.names)) +
    scale_fill_manual(values=c(Controls="white",Cases="grey"))+
    theme(axis.text.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y = element_text(size=15),
          axis.title.y = element_text(size=18),
          axis.title=element_text(size=14,face="bold"),
          strip.background = element_rect(fill="white"),
          strip.text.x = element_text(size = 12, hjust = 1, face="bold"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    xlab("") + ylab(expression(bold(phi^e)))
  return(gg)
}

#' @describeIn iCARH.plotBeta Check normality assumptions. Returns normalised data
#' and performs quantile-quantile plot
#' @export iCARH.checkNormality


iCARH.checkNormality = function(fit){
  mu = colMeans(extract(fit, inc_warmup=F, pars="mu")$mu)
  Sigma = colMeans(extract(fit, inc_warmup=F, pars="Sigma")$Sigma)
  XX = colMeans(extract(fit, inc_warmup=F, pars="XX")$XX)
  psi1 = chol(Sigma[1,,])
  psi2 = chol(Sigma[2,,])
  err=array(dim=dim(XX))
  N = dim(XX)[2]
  for(i in 1:nrow(XX)) for(n in 1:(N/2)) err[i,n,] = psi1%*%(XX[i,n,]-mu[i,n,])
  for(i in 1:nrow(XX)) for(n in (N/2+1):N) err[i,n,] = psi2%*%(XX[i,n,]-mu[i,n,])
  qqnorm(err[,1:(N/2),], cex=2.5, main="Controls")
  abline(0,1)
  qqnorm(err[,(N/2+1):N,], cex=2.5, main = "Cases")
  abline(0,1)
  return(err)
}

#' @describeIn iCARH.plotBeta Compute Watanabe-Akaike Information Criterion (WAIC)
#' @export iCARH.waic


iCARH.waic = function(fit){
  x = extract(fit, inc_warmup=F, pars="log_lik")$log_lik
  pwaic = sum(apply(x, c(2,3), var))
  lpd = sum(log(apply(exp(x), c(2,3), mean)))
  return(-2*(lpd-pwaic))
}

#' @describeIn iCARH.plotBeta Compute MADs (Mean Absolute Deviation) between true covariance matrix
#' and inferred covariance matrix for posterior predictive checks
#' @export iCARH.mad

iCARH.mad = function(fit, X){
  estimate = extract(fit, inc_warmup=F, pars="mupred")$mupred
  Tp = dim(estimate)[2]
  J = dim(estimate)[3]
  mad = array(0, dim = c(Tp, J*(J+1)/2))
  for(t in 1:Tp){
    true.cov = cov(X[t,,], use="complete.obs")[lower.tri(cov(X[t,,]),diag=T)]
    for(i in 1:nrow(estimate)){
      cov.pred = cov(t(estimate[i,t,,]))[lower.tri(cov(t(estimate[i,t,,])), diag=T)]
      mad[t,] = mad[t,] + abs( cov.pred- true.cov)/nrow(estimate)
    }
  }
  return(mad)
}
