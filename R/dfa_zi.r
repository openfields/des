# DFA ZI
library(jagsUI)

#dfa = read.csv('/home/will/Documents/Stream/Abundance/qryDFA_counthabitat.csv', header=TRUE)
dfa = read.table('C:/Users/wfields/Documents/Stream/Abundance/qryDFA_counthabitat.csv', header=TRUE, sep=',')
y <- dfa[,14:16]
y[is.na(y)]<-0
nsites <- dim(dfa)[1]
ncap<-apply(y,1,sum)
ymax<-ncap

# set nsites, ncaps, ymax
nsites<-nrow(y)
ncap<-apply(y,1,sum)
ymax<-apply(y,1,sum)

# create data list
tokunaga(dfa[,12],dfa[,10])->tok

# scale catchment area and betweenness
betsc <- scale(dfa$bet)
betsc[is.na(betsc)]<-0
fasc <- scale(dfa$facc)
fasc[is.na(fasc)] <- 0


data.dfa.zi <- list(y=y,nsites=nsites,ncap=ncap,cov1=scale(dfa[,4]),cov2=scale(dfa[,5]),cov3=scale(dfa[,6]),cov4=scale(dfa[,11]),cov5=scale(tok$t1),
             cov6=scale(dfa[,12]),cov7=betsc, cov8=fasc)
# initial values
inits <- function(){
  list (p0=runif(1),beta0=runif(1,-1,1),N=ymax+1,z=rep(1,59))
}

# parameters to monitor
parameters <- c("N","p0","beta0","beta1","beta2","beta3","beta4","beta5","beta6","beta7","beta8","omega","b4","b5","fit","fit.new")

# mcmc settings
nthin<-150
nc<-3
nb<-75000
ni<-500000

system.time(dfa.zi2<-jags(data.dfa.zi, inits, parameters, "zi3.txt", n.chains=nc, n.iter=ni, n.burnin=nb, parallel=TRUE))
print(dfa.zi2)


dfj = read.csv('/home/will/Documents/Stream/Abundance/qryDFJ_stats.csv', header=TRUE)
dfj <- read.table('C:/Users/wfields/Documents/Stream/Abundance/qryDFJ_stats.csv', header=TRUE, sep=',')
head(dfj)
y <- dfj[,17:19]
y[is.na(y)]<-0
nsites <- dim(dfj)[1]
ncap<-apply(y,1,sum)
ymax<-ncap

# set nsites, ncaps, ymax
nsites<-nrow(y)
ncap<-apply(y,1,sum)
ymax<-apply(y,1,sum)

# create data list
tokunaga(dfa[,12],dfa[,10])->tok
data.dfj <- list(y=y,nsites=nsites,ncap=ncap,cov1=scale(dfa[,4]),cov2=scale(dfa[,5]),cov3=scale(dfa[,6]),cov4=scale(dfa[,11]),
                 cov5=scale(tok$t1),cov6=scale(dfa[,12]),cov7=betsc,cov8=fasc)

# initial values
inits <- function(){
  list (p0=runif(1),beta0=runif(1,-1,1),N=ymax+1,z=rep(1,59))
}

# parameters to monitor
parameters <- c("N","p0","beta0","beta1","beta2","beta3","beta4","beta5","beta6","beta7","beta8","omega","b4","b5","fit","fit.new")

# mcmc settings
nthin<-200
nc<-3
nb<-75000
ni<-500000
# fit model
system.time(dfj.zi2<-jags(data.dfj, inits, parameters, "zi3b.txt", n.chains=nc, n.iter=ni, n.burnin=nb, parallel=TRUE))
print(dfj.zi2)

