library(MASS)
library(MPG)

### generate data ###

K = 4
p = 4
n = 1000
G = 3

mu = array(NA, dim=c(K,G,p))
Sigma = array(NA, dim=c(K,p,p) )

probs.1 = c(0.08, 0.9, 0.01, 0.01 ) 
probs.2 = c(0.06, 0.9, 0.03, 0.01 ) 
probs.3 = c(0.04, 0.9, 0.05, 0.01 ) 

mu[1,,] = matrix( rep(c(1,8), G*p/2), ncol = p, byrow=TRUE )
mu[2,,] = matrix( rep(c(8,8), G*p/2), ncol = p, byrow=TRUE )
mu[3,,] = matrix( rep(c(1,1), G*p/2), ncol = p, byrow=TRUE )
mu[4,,] = matrix( rep(c(8,1), G*p/2), ncol = p, byrow=TRUE )

mu[1,2,2] = mu[1,1,2] - 1
mu[1,3,2] = mu[1,1,2] - 2

mu[4,2,1] = mu[4,1,1] + 1
mu[4,3,1] = mu[4,1,1] + 2

Sigma[1,,] = diag(rep(1,p))
Sigma[2,,] = diag(rep(2,p))
Sigma[3,,] = diag(rep(.2,p))
Sigma[4,,] = diag(rep(.1,p))

mixture.components.1 = table( factor(sample(K, n, replace=TRUE, prob=probs.1),levels=1:K) )
mixture.components.2 = table( factor(sample(K, n, replace=TRUE, prob=probs.2),levels=1:K) )
mixture.components.3 = table( factor(sample(K, n, replace=TRUE, prob=probs.3),levels=1:K) )

mix.com.sum.1 = c(0, cumsum(mixture.components.1))
mix.com.sum.2 = c(0, cumsum(mixture.components.2))
mix.com.sum.3 = c(0, cumsum(mixture.components.3))

data.1 = matrix(NA, nrow=n, ncol=p)
data.2 = matrix(NA, nrow=n, ncol=p)
data.3 = matrix(NA, nrow=n, ncol=p)

class.1 = rep(NA, n)
class.2 = rep(NA, n)
class.3 = rep(NA, n)

for(k in 1:K)
{
  data.1[(mix.com.sum.1[k]+1):mix.com.sum.1[k+1],] = 
    mvrnorm(mixture.components.1[k], mu[k,1,], Sigma[k,,])
  
  class.1[(mix.com.sum.1[k]+1):mix.com.sum.1[k+1]] = k 
  
  data.2[(mix.com.sum.2[k]+1):mix.com.sum.2[k+1],] = 
    mvrnorm(mixture.components.2[k], mu[k,2,], Sigma[k,,])  
  
  class.2[(mix.com.sum.2[k]+1):mix.com.sum.2[k+1]] = k
  
  data.3[(mix.com.sum.3[k]+1):mix.com.sum.3[k+1],] = 
    mvrnorm(mixture.components.3[k], mu[k,3,], Sigma[k,,])  
  
  class.3[(mix.com.sum.3[k]+1):mix.com.sum.3[k+1]] = k
  
}


idx = sample(3*n,3*n)
Y = rbind(data.1, data.2, data.3)[idx,]
C = c(rep(0,n), rep(1,n), rep(2, n))[idx]
colnames(Y) = rep("",p)

### plot data ###
## projection on the first two dimensions 

par(mfrow=c(1,3))
par(oma=c(.1,.1,.1,.1))
par(mar=c(3.1, 3.1, 3.1, 2.1))
plot(Y[C==0,], xlim = range(Y[,1]), ylim = range(Y[,2]), main="sample 1")
plot(Y[C==1,], xlim = range(Y[,1]), ylim = range(Y[,2]), main="sample 2")
plot(Y[C==2,], xlim = range(Y[,1]), ylim = range(Y[,2]), main="sample 3")

### define prior ###

prior = list( K = 10 , 
              epsilons_range = matrix( c(10^-10, 1, 1, 2), ncol = 2, byrow=TRUE),
              tau_varphi = c(0.5,0.5,0))

### run model ###

ans = mpg(Y,C, prior = prior)

### plots ###

## color the data-points if the posterior weight difference is larger than "conf"
## "dim" is a vector indicating the pair of dimensions on which projecting the data
plotDiff(ans, type = "weight", conf = 0.9, dim = c(1,2))
## some of the mass moves from the top left cluster to the bottom left one.

# color the data-points if the posterior local perturbation is larger than "conf"
plotDiff(ans, type = "misalignment", con = 0.9, dim = c(1,2))
## the mean of the top left cluster is moving down 
## the mean of the bottom right cluster is moving to the right

### DO NOT LOOK AT SHIFT ###
## color the data-points if the posterior local shift is larger than "conf"
## I am putting zero prior mass on the local shift, so also the posterior is zero.
plotDiff(ans, type = "shift", con = 0.9, dim = c(1,2))
### -------------------- ###

