library(doParallel)

source("simfun.R")

uicc=0.02 # ICC of the outcome model
taut=(pi^2*uicc/3)/(1-uicc) # corresponding variance parameter

vicc=0.05 # ICC of the membership model
gat=(pi^2*vicc/3)/(1-vicc) # corresponding variance parameter

nc = 20 # number of clusters in each arm
msize = 25 # average cluster size
sdsize = 3 # standard deviation

nsim=5 # number of simulations

#### simulate cluster level random intercepts

su1=matrix(rnorm(nc*nsim,mean=0,sd=sqrt(taut)),nc,nsim)
su0=matrix(rnorm(nc*nsim,mean=0,sd=sqrt(taut)),nc,nsim)

sv1=matrix(rnorm(nc*nsim,mean=0,sd=sqrt(gat)),nc,nsim)
sv0=matrix(rnorm(nc*nsim,mean=0,sd=sqrt(gat)),nc,nsim)

rn = simn(msize, sdsize, nc, nsim) # simulate cluster size

# simulate predictors

cc=detectCores()
cl=makeCluster(cc)
registerDoParallel(cl)

sx1=foreach(kk=1:nsim,.combine='rbind') %:%
foreach(ii=1:nc,.combine='rbind') %dopar%{
	source("simfun.R")
			simx(rn[ii,kk],ii,kk)

			}

sx0=foreach(kk=1:nsim,.combine='rbind') %:%
foreach(ii=1:nc,.combine='rbind') %dopar%{
	source("simfun.R")
            simx(rn[ii,kk],ii,kk)	
			}

stopCluster(cl)
            
#########

xnames=c("x0","x1","x2")
xm1=data.matrix(sx1[,xnames])
xm0=data.matrix(sx0[,xnames])

## regression coefficients

tbsn=c(0.8,-0.4,0.4)
tbss =c(1.2,-0.3,0.3) 
tbss0=c(1,-0.5,0.5)

tass=c(1.6,0.2,0.1)
tasn=c(-0.1,-0.1,-0.2)

mss=xm1%*%tbss
msn=xm1%*%tbsn
mss0=xm0%*%tbss0
xss=xm1%*%tass
xsn=xm1%*%tasn
xss0=xm0%*%tass
xsn0=xm0%*%tasn

### simulate the outcome variable

cc=detectCores()
cl=makeCluster(cc)
registerDoParallel(cl)

sda1=foreach(kk=1:nsim,.combine='rbind') %:%
foreach(ii=1:nc,.combine='rbind') %dopar%{

indk=(sx1$simid==kk)

sind=which(indk&(sx1$cid==ii))

mssi=mss[sind,]
msni=msn[sind,]
xssi=xss[sind,]
xsni=xsn[sind,]

xi=sx1[sind,xnames]
ui=su1[ii,kk]
vi=sv1[ii,kk]
simby1(xi,ui,vi,mssi,msni,xssi,xsni,ii,kk)
}

sda0=foreach(kk=1:nsim,.combine='rbind') %:%
foreach(ii=1:nc,.combine='rbind') %dopar%{

indk=(sx0$simid==kk)
sind0=which(indk&(sx0$cid==ii))

mssi0=mss0[sind0,]
xssi0=xss0[sind0,]
xsni0=xsn0[sind0,]

xi0=sx0[sind0,xnames]
ui0=su0[ii,kk]
vi0=sv0[ii,kk]
simby0(xi0,ui0,vi0,mssi0,xssi0,xsni0,ii,kk)
}

simo = list(sda1=sda1,sda0=sda0)
save(simo, file = "simdata.RData")

### simulation summary

cm0=colMeans(sda0[,c("yss","ss","sn","nn")])
cm1=colMeans(sda1[,c("yss","ysn","ss","sn","nn")])

sace=rep(0,nsim)

for(kk in 1:nsim){
	da1=subset(sda1,simid==kk)
	da0=subset(sda0,simid==kk)
sace[kk]=mean(da1$y[which(da1$ss==1)])-mean(da0$y[which(da0$ss==1)]) 
}

tsace=mean(sace) # true SACE value

ssim=list(cm0=cm0,cm1=cm1,tsace=tsace)

save(ssim,file="ssim.RData")

