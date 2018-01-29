
n_flips <- 20
n_heads <- 2
n_tails <- n_flips-n_heads
p_true <- 0.1

prior_alpha <- 5
prior_beta <- 5

prop_head <- rbeta(1, prior_alpha, prior_beta)
head <- rbinom(1, 20, prop_head)

n = 10000

sim=data.frame(prop = rep(NA, n), heads=rep(NA, n), prediction=rep(2,n), deviation=rep(NA, n),  accept=rep(NA, n), distribution=rep(NA,n) )
temp <- data.frame(prop = rep(NA, n),head=rep(NA,n))

for (i in 1:n){
  temp[i,1] <- rbeta(1, prior_alpha,prior_beta )
  temp[i,2] <- rbinom(1, 20, temp[i,1])
  sim[i,1] <- temp[i,1]
  sim[i,2] <- temp[i,2]
  sim[i,4] <- abs(sim[i,3]-sim[i,2])
  if(sim[i,4]/n_flips>0.2){sim[i,5] <- "no"} else{sim[i,5] <- "yes"}
  if(sim[i,4]/n_flips>0.2){sim[i,6] <- "prior"} else{sim[i,6] <- "posterior"}
}

#e specified to be 0.2

n_post <- sum(sim$accept == "yes")

head(sim)

#create data frame for 
dist <- split(sim,sim$accept)
posterior <- dist$yes
sim[,6] <- "prior"
dist <- rbind(posterior, sim)


p1<- ggplot(sim, aes(prop, heads,colour=accept ))+ geom_point(size=0.1)
p1

p2 <- ggplot(dist, aes(x=prop,fill=distribution)) + geom_density(position = "identity", alpha = 0.3)
p2

#stat_function

