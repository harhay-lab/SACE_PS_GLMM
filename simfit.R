library(doParallel)
source("bme.R")
source("mefit.R")
source("simfun.R")
nbp=10 # number of bootstrap samples

########

load("simdata.RData")

sda0=simo$sda0
sda1=simo$sda1

nc=max(sda0$cid)
nsim=max(sda0$simid)

xnames=c("x0","x1","x2")

cc=detectCores()
cl=makeCluster(cc)
registerDoParallel(cl)

esace=foreach(kk=1:nsim,.combine='c',.multicombine=T,
.maxcombine=nsim,.errorhandling ="remove") %dopar%{

da0=subset(sda0,simid==kk)
da1=subset(sda1,simid==kk)
io=ini(da1,da0,xnames)
bfit(da1,da0,io)}

stopCluster(cl) 

#######

cl=makeCluster(cc)
registerDoParallel(cl)

bsace=foreach(kk=1:nsim) %:%
foreach(icount(nbp),.combine='c',.multicombine=T,.maxcombine=nbp,
.errorhandling ="pass") %dopar% {

da0=subset(sda0,simid==kk)
da1=subset(sda1,simid==kk)

s1=sample(nc,nc,replace=T)
s0=sample(nc,nc,replace=T)

bda1=bda0=NULL

for(ii in 1:nc){
da1j=subset(da1,cid==s1[ii])
da1j$cid=ii
bda1=rbind(bda1,da1j)}

for(ii in 1:nc){
da0k=subset(da0,cid==s0[ii])
da0k$cid=ii
bda0=rbind(bda0,da0k)}

bio=ini(bda1,bda0,xnames)
bfit(bda1,bda0,bio)
}

stopCluster(cl) 

######## bias, MSE and coverage

load("ssim.RData")

trsace=ssim$tsace # true SACE value

bias=abs(mean(esace)-trsace)
mse=mean((esace-trsace)^2)

cover=rep(0,nsim)
for(kk in 1:nsim){
bbi=bsace[[kk]]
cover[kk]=bcover(bbi,trsace)
}
mean(cover) # proportion of coverage