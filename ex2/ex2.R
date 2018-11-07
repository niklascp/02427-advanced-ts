run_sim <- function(N, a, sigma_v_sq, sigma_e_sq) {
    x = numeric(length = N)
    v = rnorm(N, sd = sqrt(sigma_v_sq))
    e = rnorm(N, sd = sqrt(sigma_e_sq))
    x[1] <- v[1] 
    for (t in 2:N) {
        x[t] = a * x[t - 1] + v[t]
    }
    y <- x + e
    list(x = x, y = y)
}

N = 1000
a = .4
sigma_v_sq = 1
sigma_e_sq = 1
K = 20

# Run K simulation each with N time steps.
ress <- list()
for (s in 1:K) {
    ress[[s]] <- run_sim(N, a, sigma_v_sq, sigma_e_sq)    
}

library('ggplot2')
library('RColorBrewer')

getPalette = colorRampPalette(brewer.pal(9, "Set3"))

gg <- ggplot()
for (s in 1:K) {
    gg <- gg + geom_line(data = data.frame(ress[[s]]), aes(x = seq(N), y = y), color = getPalette(20)[s])
}
gg <- gg + xlab('t')
ggsave('plots/4_a.pdf', plot = gg, width = 16, height = 9)

y <- ress[[1]]$y

run_ekf <- function(y, aInit, aVarInit, sigma.v) {
    ## Initialize
    ## Init the state vector estimate
    zt <- c(0,aInit)
    ## Init the variance matrices
    Rv <- matrix(c(sigma.v^2,0,0,0), ncol=2)
    ## sigma.e : Standard deviation of the measurement noise in the filter
    Re <- 1 

    ## Init the P matrix, that is the estimate of the state variance
    Pt <- matrix(c(Re,0,0,aVarInit), nrow=2, ncol=2)
    ## The state is [X a] so the differentiated observation function is
    Ht <- t(c(1,0))
    ## Init a vector for keeping the parameter a variance estimates
    aVar <- rep(NA,length(y))
    ## and keeping the states
    Z <- matrix(NA, nrow=length(y), ncol=2)
    Z[,1] <- zt

    ## The Kalman filtering
    for(t in 1:(length(y)-1))
    {
    ## Derivatives (Jacobians)
    Ft <- matrix(c(zt[2],0,zt[1],1), ncol=2)  # F_t-1
    # Ht does not change 
    
    ## Prediction step
    zt = c(zt[2]*zt[1],zt[2]) #z_t|t-1 f(z_t-1|t-1)
    Pt = Ft %*% Pt %*% t(Ft) + Rv #P_t|t-1
    
    ## Update step
    res = y[t] - zt[1] # the residual at time t
    St =  Ht %*% Pt %*% t(Ht) + Re # innovation covariance
    Kt = Pt %*% t(Ht) %*% St^-1 # Kalman gain
    zt = zt + Kt * res # z_t|t
    Pt = (diag(2) - Kt%*%Ht)%*%Pt #P_t|t
    
    ## Keep the state estimate
    Z[t+1,] <- zt
    ## Keep the P[2,2], which is the variance of the estimate of a
    aVar[t+1] <- Pt[2,2]
    
    }
    Z
}

run_ekf_all <- function(sims, aInit, aVarInit, sigma.v) {
    ekf_ress <- list()
    for (i in 1:length(sims)) {
        ekf_ress[[i]] <- run_ekf(sims[[i]]$y, aInit, aVarInit, sigma.v)
    }
    ekf_ress
}

plot_convergance <- function(sims, aInit, aVarInit, sigma.v, filename) {
    Zs <- run_ekf_all(sims, aInit, aVarInit,  sigma.v)
    gg <- ggplot()
    est = c()
    for (s in 1:length(Zs)) {
        dfZ <- data.frame(Zs[[s]][2:N, 2])
        colnames(dfZ) <- c('a')
        dfZ$t <- seq(N-1)
        gg <- gg + geom_line(data = dfZ, aes(x = t, y = a))
        est <- c(est, Zs[[s]][N,2])
    }
    gg <- gg + xlab('t')
    ggsave(filename, plot = gg, width = 5, height = 3)

    print(mean(est))
    print(quantile(est, probs = c(0.025, 0.975)))
}

plot_convergance(ress, aInit = .5, aVarInit = 1, sigma.v = 10, filename = 'plots/convergance01a.pdf')
plot_convergance(ress, aInit = -.5, aVarInit = 1, sigma.v = 10, filename = 'plots/convergance01b.pdf')
plot_convergance(ress, aInit = .5, aVarInit = 1, sigma.v = 1, filename = 'plots/convergance02a.pdf')
plot_convergance(ress, aInit = -.5, aVarInit = 1, sigma.v = 1, filename = 'plots/convergance02b.pdf')
plot_convergance(ress, aInit = .5, aVarInit = 10, sigma.v = 10, filename = 'plots/convergance03a.pdf')
plot_convergance(ress, aInit = -.5, aVarInit = 10, sigma.v = 10, filename = 'plots/convergance03b.pdf')
plot_convergance(ress, aInit = .5, aVarInit = 10, sigma.v = 1, filename = 'plots/convergance04a.pdf')
plot_convergance(ress, aInit = -.5, aVarInit = 10, sigma.v = 1, filename = 'plots/convergance04b.pdf')

gg
