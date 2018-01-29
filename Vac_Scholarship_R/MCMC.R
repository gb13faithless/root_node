library("mvtnorm")
library("coda")
library("tidyverse")



# Target parameters for univariate normal distributions
rho <- 0.9
mu1 <- 1; s1 <- 1
mu2 <- 2; s2 <- 1

# Parameters for bivariate normal distribution
mu <- c(mu1,mu2) # Mean 
sigma <- matrix(c(s1^2, s1*s2*rho, s1*s2*rho, s2^2), 2) # Covariance matrix

n <- 10000

  vec1 <- vector("numeric", n)
  vec2 <- vector("numeric", n)
  x1 <- 10 # theta1 initial
  x2 <- 10 # theta2 initial
  vec1[1] <- x1    
  vec2[1] <- x2
  delta <- 0.3
  
  for (i in 2:n) {
    can1 <- vec1[i-1]+runif(1, -delta, delta)  # candidate point for theta1
    can2 <- vec2[i-1]+runif(1, -delta, delta)  # candidate point for theta2
    aprob <- min(1, dmvnorm(c(can1,can2),mu,sigma)/dmvnorm(c(vec1[i-1],vec2[i-1]), mu,sigma)) # uniform is symmetric
    u2 <- runif(1) #makes uniform 1 dist
    if(u1 < aprob){vec1[i] <- can1} else{vec1[i] <- vec1[i-1]} # applies the algorithmic tests
    if(u2 < aprob){vec2[i] <- can2} else{vec2[i] <- vec2[i-1]}
  }
  vec1
  vec2

  
  
v <- cbind(vec1,vec2) # makes a vector of theta 1, 2
i <-  data.frame(Iterations = 1:n)

v <- cbind(v,i)
v <- as.data.frame(v)

ggplot(v, aes(i,vec1)) + geom_line()
ggplot(v, aes(i,vec2)) + geom_line()

install.packages("elllipse")
library("ellipse")
#scatterplot
v2 <- v[1000:10000,]
ggplot(v2, aes(vec1,vec2)) + geom_point(size=0.1)+stat_ellipse
  
tail(v)

vec1 <- as.mcmc(vec1)

traceplot(vec1)





