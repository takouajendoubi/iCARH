#' Simulation inspired from the real provided data set


Tp=4 # timepoints
N=10 # number of samples
J=14 # number of metabolites
K=2  # number of bacteria species
P=8

set.seed(124573)

## For real data
## Build pathway matrices using GetPathwaysMat
# keggid = list("Unk1", "C03299","Unk2","Unk3",
#          c("C08363", "C00712"), ...  # allowing multiple ids per metabolite
#          )
# pathways = iCARH.getPathwaysMat(keggid, "rno")

## To simulate data use iCARH.simulate

# Manually choose pathways
path.names = c("path:map00564","path:map00590","path:map00061","path:map00591",
               "path:map00592","path:map00600","path:map01040","path:map00563")
# Elements in kegg id list may contain multiple kegg ids per metabolite
# If KEGG id unknown : "Unk[number]"

## Simulation
data.sim = iCARH.simulate(Tp, N, J, P, K, path.names=path.names, Zgroupeff=c(0,4),
                          beta.val=c(1,-1,0.5, -0.5))

XX = data.sim$XX
Y = data.sim$Y
Z = data.sim$Z
pathways = data.sim$pathways

## Check inaccuracies
pathways.bin = lapply(pathways, function(x) { y=1/(x+1); diag(y)=0; y})
adjmat = rowSums(abind::abind(pathways.bin, along = 3), dims=2)
cor.thresh = 0.7
# check percentage of metabolites correlated but not in same pathway
for(i in 1:Tp) print(sum(abs(cor(X[i,,])[which(adjmat==0)])>cor.thresh)/length(which(adjmat==0)))

## Run iCARH model

# takes about 12 mins
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
fit.sim = iCARH.model(XX, Y, Z, pathways, control = list(adapt_delta = 0.99, max_treedepth=10),
                      iter = 1100, chains = 2, pars=c("YY","Xmis","Ymis"), include=F)

## Processing results

# Bacteria effects
gplot = iCARH.plotBeta(fit.sim)
gplot

# treatments effects
gplot = iCARH.plotTreatmentEffect(fit.sim)
gplot

# Pathway analysis
gplot = iCARH.plotPathwayPerturbation(fit.sim, path.names)
gplot

# Normality assumptions
par(mfrow=c(1,2))
iCARH.checkNormality(fit.sim)

# WAIC
waic = iCARH.waic(fit.sim)
waic

# Posterior predictive checks
# MAD : mean absolute deviation between covariance matrices
mad = iCARH.mad(fit.sim, XX)
quantile(mad)

