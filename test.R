library(MPG)


n = 500
p = 2
J = 2

Y = matrix(  rnorm(n*p, mean=0) , ncol= p, byrow=TRUE )
C = sample(0:(J-1), n, replace = TRUE)
Y[C==1] = Y[C==1] + 0.5
colnames(Y) = c("X1","X2")

mcmc = list(nburn = 5000, nsave = 1000, nskip = 1, ndisplay = 1000)

prior = list( K = 10,
              epsilons_range =  matrix( c(10^-10, 10^-1, 10^-1, 2), ncol = 2),   
              m_2 = rep(0,p),
              nu_2 = p+2, 
              nu_1 = p+2,
              Psi_2 = diag(1,p),
              S_2 = 1000*diag(1,p),
              tau_k0 = rep(4,2),
              tau_alpha = rep(1,2),
              tau_rho = c(9,1),
              point_masses_rho = c(0.0, 0.0),
              tau_varphi = c(4.5,4.5,1),
              point_masses_varphi = c(0.0, 0.0, 0.0)
            )

ans = mpg(Y, C, prior, mcmc)




plotDiff(ans, type = "weight", conf = 0.9)
plotDiff(ans, type = "misalignment", con = 0.9)
plotDiff(ans, type = "shift", con = 0.9)


pred = predict(ans)
plotPred(pred)




par(mfrow=c(1,1))
hist(ans$chain$rho, nclass = 20)

par(mfrow=c(3,1))
plot(ans$chain$alpha[,1], type = "l")
plot(ans$chain$alpha[,2], type = "l")
plot(ans$chain$alpha[,3], type = "l")


par(mfrow=c(1,3))
plot(ans$chain$varphi[,1], type = "l")
plot(ans$chain$varphi[,2], type = "l")
plot(ans$chain$varphi[,3], type = "l")


par(mfrow=c(1,2))
plot(ans$chain$epsilon[,1], type = "l")
plot(ans$chain$epsilon[,2], type = "l")

par(mfrow=c(1,2))
hist(ans$chain$epsilon[,1])
hist(ans$chain$epsilon[,2])


par(mfrow=c(2,1))
plot(colMeans(ans$chain$Z<ans$prior$K))
plot(rowMeans(ans$chain$Z<ans$prior$K))

S = matrix(NA, nrow=nrow(ans$chain$Z),  ncol=ncol(ans$chain$Z))
for( it in 1:nrow(S))
  S[it,] = ans$chain$T[it, ans$chain$Z[it,]+1]

par(mfrow=c(2,3))
plot(colMeans(S==0), ylim=c(0,1))
plot(colMeans(S==1), ylim=c(0,1))
plot(colMeans(S==2), ylim=c(0,1))
plot(rowMeans(S==0), ylim=c(0,1))
plot(rowMeans(S==1), ylim=c(0,1))
plot(rowMeans(S==2), ylim=c(0,1))

par(mfrow=c(2,1))
plot(colMeans(ans$chain$T==0))
plot(rowMeans(ans$chain$T==0))

  