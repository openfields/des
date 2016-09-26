source("F:/Stream/Occupancy/tokunaga.r")

cat("
model {
    beta0 ~ dnorm(0,.1)
    beta1 ~ dnorm(0,.1) #pc1
    beta2 ~ dnorm(0,.1) #pc2
    beta3 ~ dnorm(0,.1) #pc3
    beta4 ~ dnorm(0,.1) # total stream
    beta5 ~ dnorm(0,.1) # T1
    beta6 ~ dnorm(0,.1) # length ratio
        
    b4 ~ dnorm(0,.1)
    b5 ~ dnorm(0,.1)
    p0~dnorm(0,.1)
    
    for(i in 1:nsites){
    p[i]<- p0 + b4*cov1[i] + b5*cov2[i]# could have covariates here
    mu[i,1]<- p[i]
    mu[i,2]<- p[i]*(1-p[i])
    mu[i,3]<- p[i]*(1-p[i])*(1-p[i])
    pi0[i]<- 1 - mu[i,1]-mu[i,2]-mu[i,3]
    pcap[i]<-1-pi0[i]
    
    for(j in 1:3){
    muc[i,j]<-mu[i,j]/pcap[i]
    }
    
    # 1. model part 1: the conditional multinomial
    y[i,1:3] ~ dmulti(muc[i,1:3],ncap[i])
    
    # 2. model for the observed count of uniques
    ncap[i] ~ dbin(pcap[i],N[i])
    
    # 3. abundance model
    N[i] ~ dpois(lambda[i])
    log(lambda[i])<- beta0 + beta1*cov1[i] + beta2*cov2[i] + beta3*cov3[i] + beta4*cov4[i] + beta5*cov5[i]+ beta6*cov6[i]

    }

    for(i in 1:nsites){
      ncap.fit[i] ~ dbin(pcap[i],N[i])
      y.fit[i,1:3] ~ dmulti(muc[i,1:3],ncap[i]) # note this is conditional on ncap[i]
      for(t in 1:3){
        e1[i,t]<- muc[i,t]*ncap[i]
        resid1[i,t]<- pow(pow(y[i,t],0.5)-pow(e1[i,t],0.5),2)
        resid1.fit[i,t]<- pow(pow(y.fit[i,t],0.5) - pow(e1[i,t],0.5),2)
      }
      e2[i]<- N[i]*lambda[i]
      resid2[i]<- pow( pow(ncap[i],0.5) - pow(e2[i],0.5),2)
      resid2.fit[i]<- pow( pow(ncap.fit[i],0.5) - pow(e2[i],0.5),2)
    }
    ft1.data<- sum(resid1[,])
    ft1.post<- sum(resid1.fit[,])
    ft2.data<- sum(resid2[])
    ft2.post<- sum(resid2.fit[])

    }
    ",fill=TRUE,file="rem3p.txt")

# get observed counts
# y<- as.matrix(getY(ovenFrame))
dfa = read.csv('F:/Stream/Abundance/qryDFA_counthabitat.csv', header=TRUE)
dfa <- read.table('C:/Users/wfields/Documents/Stream/Abundance/qryDFA_counthabitat.csv', header=TRUE, sep=',')
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
data.dfa <- list(y=y,nsites=nsites,ncap=ncap,cov1=scale(dfa[,4]),cov2=scale(dfa[,5]),cov3=scale(dfa[,6]),cov4=scale(dfa[,11]),cov5=scale(tok$t1),cov6=scale(dfa[,12]))

# initial values
inits <- function(){
  list (p0=runif(1),beta0=runif(1,-1,1),N=ymax+2 )
}

# parameters to monitor
parameters <- c("N","p0","beta0","beta1","beta2","beta3","beta4","beta5","beta6","b4","b5","ft1.data","ft1.post","ft2.data","ft2.post")

# mcmc settings
nthin<-1
nc<-3
nb<-2000
ni<-10000

# call jags
dfa.fittest<-jags(data.dfa, inits, parameters, "rem3p.txt", n.chains=nc, n.iter=ni, n.burnin=nb, parallel=TRUE)
print(dfa.fittest)







#-------------------------------------------------------------------
#-------------------------------------------------------------------

# read data for df j
dfj = read.csv('F:/Stream/Abundance/qryDFJ_stats.csv', header=TRUE)
head(dfj)
y <- dfj[,17:19]
y[is.na(y)]<-0
nsites <- dim(dja)[1]
ncap<-apply(y,1,sum)
ymax<-ncap

# set nsites, ncaps, ymax
nsites<-nrow(y)
ncap<-apply(y,1,sum)
ymax<-apply(y,1,sum)

# create data list
tokunaga(dfa[,12],dfa[,10])->tok
data <- list(y=y,nsites=nsites,ncap=ncap,cov1=scale(dfa[,4]),cov2=scale(dfa[,5]),cov3=scale(dfa[,6]),cov4=scale(dfa[,11]),cov5=scale(tok$t1),cov6=scale(dfa[,12]))

# initial values
inits <- function(){
  list (p0=runif(1),beta0=runif(1,-1,1),N=ymax+2 )
}

# parameters to monitor
parameters <- c("N","p0","beta0","beta1","beta2","beta3","beta4","beta5","beta6","b4","b5","ft1.data","ft1.post","ft2.data","ft2.post")

# mcmc settings
nthin<-1
nc<-3
nb<-10000
ni<-60000
# fit model
dfj.fit1<-jags(data, inits, parameters, "rem3p.txt", n.chains=nc, n.iter=ni, n.burnin=nb)
print(dfj.fit1)



#---------------------------------------------------
#---------------------------------------------------
# read data for dm a
dma = read.csv('F:/Stream/Abundance/qryDMA_stats.csv', header=TRUE)
head(dma)
y<-dma[,14:16]
y[is.na(y)]<-0
nsites <- dim(dma)[1]
ncap<-apply(y,1,sum)
ymax<-ncap
data <- list(y=y,nsites=nsites,ncap=ncap,cov1=scale(dfa[,4]),cov2=scale(dfa[,5]),cov3=scale(dfa[,6]),cov4=scale(dfa[,11]),cov5=scale(tok$t1),cov6=scale(dfa[,12]))

# initial values
inits <- function(){
  list (p0=runif(1),beta0=runif(1,-1,1),N=ymax+2)
}

# parameters to monitor
parameters <- c("N","p0","beta0","beta1","beta2","beta3","beta4","beta5","beta6","b4","b5","ft1.data","ft1.post","ft2.data","ft2.post")

# mcmc settings
nthin<-1
nc<-3
nb<-2000
ni<-20000
# fit model
dma.fit1<-jags(data, inits, parameters, "rem3p.txt", n.chains=nc, n.iter=ni, n.burnin=nb)
print(dma.fit1)




#---------------------------------------------------
# read data for dm sa
dmsa = read.csv('F:/Stream/Abundance/qryDMSA_stats.csv', header=TRUE)
head(dmsa)
y<-dmsa[,12:14]
y[is.na(y)]<-0
nsites <- dim(dfa)[1]
ncap<-apply(y,1,sum)
ymax<-ncap

# set nsites, ncaps, ymax
nsites<-nrow(y)
ncap<-apply(y,1,sum)
ymax<-apply(y,1,sum)

# create data list
data <- list(y=y,nsites=nsites,ncap=ncap,cov1=scale(dfa[,4]),cov2=scale(dfa[,5]),cov3=scale(dfa[,6]),cov4=scale(dfa[,11]),cov5=scale(tok$t1),cov6=scale(dfa[,12]))

# initial values
inits <- function(){
  list (p0=runif(1),beta0=runif(1,-1,1),N=ymax+2 )
}

# parameters to monitor
parameters <- c("N","p0","beta0","beta1","beta2","beta3","beta4","beta5","beta6","b4","b5","ft1.data","ft1.post","ft2.data","ft2.post")

# mcmc settings
nthin<-1
nc<-3
nb<-2000
ni<-10000

dmsa.fit1<-jags(data, inits, parameters, "rem3p.txt", n.chains=nc, n.iter=ni, n.burnin=nb)
print(dmsa.fit1)

#---------------------------------------------------


# read data for dm j
dmj = read.csv('F:/Stream/Abundance/qryDMJ_stats.csv', header=TRUE)
head(dmj)
y<-dmj[,12:14]
y[is.na(y)]<-0
nsites <- dim(dmj)[1]
ncap<-apply(y,1,sum)
ymax<-ncap

# set nsites, ncaps, ymax
nsites<-nrow(y)
ncap<-apply(y,1,sum)
ymax<-apply(y,1,sum)

# create data list
data <- list(y=y,nsites=nsites,ncap=ncap,cov1=scale(dfa[,4]),cov2=scale(dfa[,5]),cov3=scale(dfa[,6]),cov4=scale(dfa[,11]),cov5=scale(tok$t1),cov6=scale(dfa[,12]))

# initial values
inits <- function(){
  list (p0=runif(1),beta0=runif(1,-1,1),N=ymax+2 )
}

# parameters to monitor
parameters <- c("N","p0","beta0","beta1","beta2","beta3","beta4","beta5","beta6","b4","b5","ft1.data","ft1.post","ft2.data","ft2.post")

# mcmc settings
nthin<-1
nc<-3
nb<-5000
ni<-20000

# fit model
dmj.fit1<-jags(data, inits, parameters, "rem3p.txt", n.chains=nc, n.iter=ni, n.burnin=nb)
print(dmj.fit1)




#------------------------------------------------
# zero inflated model  
#------------------------------------------------

cat("
model {
    beta0 ~ dnorm(0,.1)
    beta1 ~ dnorm(0,.1) #pc1
    beta2 ~ dnorm(0,.1) #pc2
    beta3 ~ dnorm(0,.1) #pc3
    beta4 ~ dnorm(0,.1) # total stream
    beta5 ~ dnorm(0,.1) # T1
    beta6 ~ dnorm(0,.1) # length ratio
    
    b4 ~ dnorm(0,.1)
    b5 ~ dnorm(0,.1)
    p0~dnorm(0,.1)

    omega ~ dunif(0,1)
    
    for(i in 1:nsites){
    p[i]<- p0 + b4*cov1[i] + b5*cov2[i]# could have covariates here
    mu[i,1]<- p[i]
    mu[i,2]<- p[i]*(1-p[i])
    mu[i,3]<- p[i]*(1-p[i])*(1-p[i])
    pi0[i]<- 1 - mu[i,1]-mu[i,2]-mu[i,3]
    pcap[i]<-1-pi0[i]
    
    for(j in 1:3){
    muc[i,j]<-mu[i,j]/pcap[i]
    }
    
    # 1. model part 1: the conditional multinomial
    y[i,1:3] ~ dmulti(muc[i,1:3],ncap[i])
    
    # 2. model for the observed count of uniques
    ncap[i] ~ dbin(pcap[i],N[i])
    
    # 3. abundance model
    z[i] ~ dbern(omega)
    N[i] ~ dpois(lam.eff[i])
    lam.eff[i] <- z[i] * lambda[i]
    log(lambda[i])<- beta0 + beta1*cov1[i] + beta2*cov2[i] + beta3*cov3[i] + beta4*cov4[i] + beta5*cov5[i]+ beta6*cov6[i]
    
    }
    
    for(i in 1:nsites){
    ncap.fit[i] ~ dbin(pcap[i],N[i])
    y.fit[i,1:3] ~ dmulti(muc[i,1:3],ncap[i]) # note this is conditional on ncap[i]
    for(t in 1:3){
    e1[i,t]<- muc[i,t]*ncap[i]
    resid1[i,t]<- pow(pow(y[i,t],0.5)-pow(e1[i,t],0.5),2)
    resid1.fit[i,t]<- pow(pow(y.fit[i,t],0.5) - pow(e1[i,t],0.5),2)
    }
    e2[i]<- N[i]*lambda[i]
    resid2[i]<- pow( pow(ncap[i],0.5) - pow(e2[i],0.5),2)
    resid2.fit[i]<- pow( pow(ncap.fit[i],0.5) - pow(e2[i],0.5),2)
    }
    ft1.data<- sum(resid1[,])
    ft1.post<- sum(resid1.fit[,])
    ft2.data<- sum(resid2[])
    ft2.post<- sum(resid2.fit[])
    
    }
    ",fill=TRUE,file="rem3zi.txt")

dfa = read.csv('F:/Stream/Abundance/qryDFA_counthabitat.csv', header=TRUE)
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
data <- list(y=y,nsites=nsites,ncap=ncap,cov1=scale(dfa[,4]),cov2=scale(dfa[,5]),cov3=scale(dfa[,6]),cov4=scale(dfa[,11]),cov5=scale(tok$t1),cov6=scale(dfa[,12]))

# initial values
inits <- function(){
  list (p0=runif(1),beta0=runif(1,-1,1),beta1=runif(1,-1,1),beta2=runif(1,-1,1),beta3=runif(1,-1,1),beta4=runif(1,-1,1),
        beta5=runif(1,-1,1),beta6=runif(1,-1,1),N=ymax,omega=.54)
}

# parameters to monitor
parameters <- c("N","p0","beta0","beta1","beta2","beta3","beta4","beta5","beta6","omega","b4","b5","ft1.data","ft1.post","ft2.data","ft2.post")

# mcmc settings
nthin<-5
nc<-3
nb<-15000
ni<-80000

# call jags
dfa.zi1<-jags(data, inits, parameters, "rem3zi.txt", n.chains=nc, n.iter=ni, n.burnin=nb)
print(dfa.zi1)

#-----------------------------------------------------
#-----------------------------------------------------
#JUVENILE FUSCUS

# read data for df j
dfj = read.csv('F:/Stream/Abundance/qryDFJ_stats.csv', header=TRUE)
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
data <- list(y=y,nsites=nsites,ncap=ncap,cov1=scale(dfa[,4]),cov2=scale(dfa[,5]),cov3=scale(dfa[,6]),cov4=scale(dfa[,11]),cov5=scale(tok$t1),cov6=scale(dfa[,12]))

# initial values
inits <- function(){
  list (p0=runif(1),beta0=runif(1,-1,1),N=ymax,omega=.69)
}

# parameters to monitor
parameters <- c("N","p0","omega","beta0","beta1","beta2","beta3","beta4","beta5","beta6","b4","b5","ft1.data","ft1.post","ft2.data","ft2.post")

# mcmc settings
nthin<-5
nc<-3
nb<-20000
ni<-60000
# fit model
dfj.zi1<-jags(data, inits, parameters, "rem3zi.txt", n.chains=nc, n.iter=ni, n.burnin=nb)
print(dfj.fit1)


