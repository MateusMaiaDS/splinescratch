rm(list=ls())
# Building up the scratch functions
pos <- function(x, x_i){
     
     # Getting the positive part only
     return(ifelse((x-x_i)>0,(x-x_i),0))
}

# Function to generate B
bspline <- function(x,x_obs, knots_ = nrow(x_obs)){
     
     B <- matrix(NA, nrow = nrow(x),ncol = nrow(x_obs)+2)
     
     error <- .Machine$double.eps 
     
     B[,1] <- 1     
     B[,2] <- x
     x_n <- max(x_obs)
     x_n_1 <- max(x_obs[-which.max(x_obs)])
     
     for(i in 1:nrow(x_obs)){
          B[,i+2] <- (pos(x = x,x_i = x_obs[i,1])^3-pos(x = x,x_i = x_n)^3)/(x_n-x_obs[i,1]+error) - (pos(x = x,x_i = x_n_1)^3-pos(x = x,x_i = x_n)^3)/(x_n-x_n_1)
     }
     
     return(B)
}

# Creating some simulated data
n <- 100
x <- matrix(runif(n = n,min = -5,max = 5))
y <- sin(x) + rnorm(n = n,sd = 0.1)

plot(x,y)

# Running function arguments
bspline_ <- bspline(x = x,x_obs = x,knots_ = x)

# Sample beta coefficients
beta_sample <- function(bspline_, y,tau, tau_b){
     
     mvn_mean <- solve(crossprod(bspline_)+(tau_b/tau)*diag(nrow = ncol(bspline_)),crossprod(bspline_,y))
     mvn_var <- solve(crossprod(bspline_)+(tau_b/tau)*diag(nrow = ncol(bspline_)))
     
     beta_post <- mvnfast::rmvn(n = 1,mu = mvn_mean,sigma = mvn_var)
     
     return(beta_post)
}

tau_sample <- function(y,y_hat, a_tau, d_tau){
     
     rgamma(n = 1,shape = a_tau + 0.5*length(y),rate = d_tau + 0.5*crossprod((y-y_hat)))
}

# Initalising the MCMC post sample
mcmc_sampler <- function(x,y,n_post,n_burn, tau = 1){
     
     
     # Getting the bspline object
     B <- bspline(x = x,x_obs = x,knots_ = x)
     
     # Gettin the n_mcmc samples
     n_mcmc <- n_post + n_burn
     y_hat_post <- matrix(NA, nrow = n_post,ncol = length(y))
     beta_post <- matrix(NA, nrow = n_post, ncol = ncol(B))
     tau_post <- numeric(n_post)
     curr <- 0
     for(i in 1:n_mcmc){
          
          beta <- beta_sample(bspline_ = B,y = y,tau = tau,tau_b = 1)
          y_hat <- tcrossprod(B,beta)
          tau <- tau_sample(y = y,y_hat = y_hat,a_tau = 0.001,d_tau = 0.001)
          
          if(i > n_burn){
               curr <- curr + 1
               y_hat_post[curr,] <- unlist(y_hat)
               tau_post[curr] <- tau
               beta_post[curr,] <- beta
          }
     }
     
     return( list(y_hat_post = y_hat_post,
                  tau_post = tau_post,
                  beta_post = beta_post))
     .
}


sample_obj <- mcmc_sampler(x = x,y = y,n_post = 2000,n_burn = 500)

plot(x,sin(x), col = "orange")
points(x,y)
points(x,sample_obj$y_hat_post %>% colMeans(), col = "blue",pch = 20)

plot(sample_obj$tau_post,type = "l")

# Getting for a x_new
x_new <- seq(min(x)-1,max(x)+1,length.out = 1000) %>% as.matrix()
new_b <- bspline(x = x_new,x_obs = x,knots_ = x)
y_hat_new <- tcrossprod(x = new_b,t(matrix(apply(sample_obj$beta_post,2,mean))))

plot(x_new,y_hat_new)
points(x,y)
points(x_new,sin(x_new), col = "orange")
