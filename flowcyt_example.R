
setwd("~/Dropbox/Duke/Thesis/locally_tied_stick_breaking/examples/eqapol/")
library(MPG)

#####################################################
### same subject different stimulation conditions ###
### to see changes with stimulation               ###
#####################################################

data.1 = read.csv("~/Dropbox/Duke/Thesis/Flowcytometry/Equapol/EQAPOL_101_EP4_K69038GP_Unstim_02064.csv",  header = TRUE)
data.2 = read.csv("~/Dropbox/Duke/Thesis/Flowcytometry/Equapol/EQAPOL_101_EP4_K69038GP_CEF_02072.csv",  header = TRUE)
data.3 = read.csv("~/Dropbox/Duke/Thesis/Flowcytometry/Equapol/EQAPOL_101_EP4_K69038GP_CMVpp65_02069.csv",  header = TRUE)

# function to remove "extreme" data points 
remove.extremes <- function(data)
{
  for(i in 1:ncol(data))
  {
    idx.max = which( data[,i] == max(data[,i],na.rm = TRUE) )
    if(length(idx.max)>1)
    {
      data[idx.max, ] = NA
    }
    idx.min = which( data[,i] == min(data[,i],na.rm = TRUE) )
    if(length(idx.min)>1)
    {
      data[idx.min, ] = NA
    }  
  }    
  return(na.omit(data))
}

data.1 = remove.extremes(data.1)
data.2 = remove.extremes(data.2)
data.3 = remove.extremes(data.3)

data = rbind(data.1, data.2, data.3)

# subsample 5000 observations from each sample
nobs = 5000
set.seed(1)
dim = seq(5,10)
Y_tot = scale( as.matrix( rbind(  data.1[sample(nrow(data.1),nobs),], 
                                  data.2[sample(nrow(data.2),nobs),],
                                  data.3[sample(nrow(data.3),nobs),] )  ) )[,dim]

C =  c(rep(0,nobs),rep(1,nobs), rep(2,nobs) )
idx.order = sample(nrow(Y_tot),nrow(Y_tot))
Y_tot = Y_tot[idx.order,]
C = C[idx.order]



mcmc = list(nburn = 50000, nsave = 500, nskip = 100, ndisplay = 100)



p = ncol(Y_tot)

# Here K is equal to the number of shared weight clusters as well as the number of weight variation clusters
# The model has in total 2K clusters

# epsilons_range defines the extremes of the uniform prior on epsilon_1 (first row) and epsilon_2 (second row)

# tau_varphi are the dirichlet hyperparameters controlling the probability:
# (no local perturbation, small perturbation, large perturbation)
# if you fix one of those to zero, than that component is removed.
# Example:
# tau_varphi = (a,b,0) the model has no large perturbation and the prior 
# on (no local perturbation, small perturbation) is a Beta(a,b)

prior = list( K = 75,
              epsilons_range =  matrix( c(10^-10, 1, 1, 2), ncol = 2),   
              m_2 = rep(0,p), S_2 = 1000*diag(1,p),
              nu_1 = p+2, 
              nu_2 = p+2, Psi_2 = diag(1,p),        
              tau_k0 = c(4,4),
              tau_alpha = c(1,1),
              tau_rho = c(1,1), point_masses_rho = c(0.0, 0.0),
              tau_varphi = c(1,1,0), point_masses_varphi = c(0.0, 0.0, 0.0)
)


ans = mpg(Y_tot, C, prior, mcmc)


## posterior of rho, varphi and epsilon
par(mfrow=c(1,3))
x = seq(0,1,length=1000)
hist(ans$chain$rho, prob=TRUE, xlim = c(0,1), 
     main = "", xlab = expression(rho), cex.lab = 2, ylab = "")
lines(x, dbeta(x, ans$prior$tau_rho[1],ans$prior$tau_rho[2]), col="red")
hist(ans$chain$varphi[,1], prob=TRUE, xlim = c(0,1), 
     main = "", xlab = expression(varphi), cex.lab = 2, ylab = "")
lines(x, dbeta(x, ans$prior$tau_varphi[1],ans$prior$tau_varphi[2]), col="red")
hist(ans$chain$epsilon[,1], prob=TRUE, xlim = c(0, max (ans$prior$epsilons_range[1,])), 
     main = "", xlab = expression(epsilon), cex.lab = 2, ylab = "")
lines(x, dunif(x, ans$prior$epsilons_range[1,1], ans$prior$epsilons_range[1,2]), col="red")


par(mfrow=c(1,1))
hist(colMeans(ans$chain$Z>=ans$prior$K), xlim = c(0,1), nclass = 20, 
     main="Model with local perturbations", xlab = expression(paste( "P(", Z[i][j]>K[0], "|data)")), 
     col = c(rep("grey",19),rep("darkred",1)) )

## different ways to plot the differences across the samples

dim = c(4,5)
plotDiff2(ans, dim=dim, main = c("Unstim", "CEF", "CMVpp65"))


dim = c(4,5)
plotDiff(ans, dim=dim, type="misalignment")

plotDiff(ans, dim=dim, type="weight")


## predict 

pred = predict.MPG(object = ans, remove_misalignment = TRUE)

par(mfrow=c(1,3))
obs = 1000
hist(pred$class_probabilities[obs,1,])
hist(pred$class_probabilities[obs,2,])
hist(pred$class_probabilities[obs,3,])

par(mfrow=c(1,1))
plotPred(pred, q = .5)


## remove small perturbation from the data

cal = calibrate(ans)


par(mfrow=c(1,2))
dim = c(1,2)
plot(ans$data$Y[,dim], col = ans$data$C+1, main = "original", pch = ans$data$C+1 )
plot(cal$Ycal[,dim], col=ans$data$C+1, main = "calibrated", pch = ans$data$C+1)



