trpmod=function(sindi,xss,xsn,gat,ni,nms){

svi=rnorm(nms,0,sqrt(gat))

pss=psn=pnn=sdess=sdesn=matrix(0,ni,nms)
fz1vi=rep(0,nms)

indi=rep(0,ni)
indi[sindi]=1 # survivor
nindi=1-indi

for(k in 1:nms){

ess=exp(xss+svi[k])
esn=exp(xsn+svi[k])
dsn=(1+ess+esn)

pnn[,k]=pnnk=1/dsn
fz1vi[k]=prod(((1-pnnk)^indi)*(pnnk^nindi))
pss[,k]=pssk=ess*pnnk
psn[,k]=psnk=esn*pnnk	

sdess[,k]=pssk*(1-pssk)
sdesn[,k]=psnk*(1-psnk)
}

mpss=rowMeans(pss)
mpsn=rowMeans(psn)
mpnn=rowMeans(pnn)

fz1=mean(fz1vi)
if(fz1==0){ev1i2=0
msdess1=rowMeans(sdess)
msdesn1=rowMeans(sdesn)
} else{ev1i2=mean(svi^2*fz1vi)/fz1
msdess1=rowMeans(sdess*fz1vi)/fz1
msdesn1=rowMeans(sdesn*fz1vi)/fz1}

return(list(mpss=mpss,mpsn=mpsn,mpnn=mpnn,
ev1i2=ev1i2,msdesn1=msdesn1,msdess1=msdess1))	
}

conpmod=function(sindi,xss,xsn,gat,ni,nms){

svi=rnorm(nms,0,sqrt(gat))

pss=psn=pnn=sdess=sdesn=matrix(0,ni,nms)
fz0vi=rep(0,nms)

indi=rep(0,ni)
indi[sindi]=1
nindi=1-indi

for(k in 1:nms){

ess=exp(xss+svi[k])
esn=exp(xsn+svi[k])
dsn=(1+ess+esn)

pnn[,k]=pnnk=1/dsn
pss[,k]=pssk=ess*pnnk
psn[,k]=psnk=esn*pnnk	

fz0vi[k]=prod(((1-pssk)^nindi)*(pssk^indi))

sdess[,k]=pssk*(1-pssk)
sdesn[,k]=psnk*(1-psnk)
}

mpss=rowMeans(pss)
mpsn=rowMeans(psn)
mpnn=rowMeans(pnn)

fz0=mean(fz0vi)
if(fz0==0){
ev0i2=0
msdess0=rowMeans(sdess)
msdesn0=rowMeans(sdesn)
} else{
ev0i2=mean(svi^2*fz0vi)/fz0
msdess0=rowMeans(sdess*fz0vi)/fz0
msdesn0=rowMeans(sdesn*fz0vi)/fz0}

meta0=mpsn/(mpsn+mpnn)

return(list(mpss=mpss,mpsn=mpsn,mpnn=mpnn,meta0=meta0,
ev0i2=ev0i2,msdesn0=msdesn0,msdess0=msdess0))	
}

alp=function(xi,txi,mpss,msdss,mpsn,msdsn,mi,p){

mpssm=matrix(mpss,nrow=p,ncol=mi,byrow=T)
msdssm=matrix(msdss,nrow=p,ncol=mi,byrow=T)

dass1=-rowSums(txi*mpssm)
dass2=-(txi*msdssm)%*%xi

mpsnm=matrix(mpsn,nrow=p,ncol=mi,byrow=T)
msdsnm=matrix(msdsn,nrow=p,ncol=mi,byrow=T)
	
dasn1=-rowSums(txi*mpsnm)
dasn2=-(txi*msdsnm)%*%xi

return(list(dass1=dass1,dass2=dass2,
dasn1=dasn1,dasn2=dasn2))
}

ini=function(da1,da0,xnames){
 # cid: cluster name
n1=length(unique(da1$cid))
m1v=as.vector(table(da1$cid))
m1=sum(m1v)
m1v1=as.vector(table(da1$cid,da1$yind)[,2])

xm1=da1[,c("cid",xnames),drop=F]
p=ncol(xm1)-1
y1=da1[,"y"]
sind1=tapply(da1$yind,da1$cid,function(x){which(x==1)})	    

####

n0=length(unique(da0$cid))
m0v=as.vector(table(da0$cid))
m0=sum(m0v)
m0v1=as.vector(table(da0$cid,da0$yind)[,2])

xm0=da0[,c("cid",xnames),drop=F]
y0=da0[,"y"]
sind0=tapply(da0$yind,da0$cid,function(x){which(x==1)})	 
	    
alphass=t(t((1:p)/p))/20
alphasn=-t(t((1:p)/p))/10

s1=lm(y1~as.matrix(xm1[,-1])-1)
s0=lm(y0~as.matrix(xm0[,-1])-1)

betass1=s1$coef
betass0=s0$coef
betasn=(betass1+betass0)/2

sat=(summary(s1)$sigma^2+summary(s0)$sigma^2)/2
taut=sat/5
gat=sat/10

return(list(alphass=alphass,alphasn=alphasn,betass1=betass1,
betasn=betasn,betass0=betass0,p=p,sat=sat,gat=gat,taut=taut,
xm1=xm1,sind1=sind1,y1=y1,n1=n1,m1v=m1v,m1v1=m1v1,m1=m1,m11=sum(m1v1),
xm0=xm0,sind0=sind0,y0=y0,n0=n0,m0v=m0v,m0v1=m0v1,m0=m0,m01=sum(m0v1)))
}

dff=function(to,co,betass1,betasn,betass0,alphass,alphasn){
	
dif.betass1=betass1-to$nbetass1
dif.betasn=betasn-to$nbetasn
dif.betass0=betass0-co$nbetass0
dif.beta=c(dif.betass1,dif.betasn,dif.betass0)

betass1=to$nbetass1
betasn=to$nbetasn
betass0=co$nbetass0

#####

nalphass=alphass-solve(to$dass2+co$dass2,to$dass1+co$dass1+to$aess+co$aess)

nalphasn=alphasn-solve(to$dasn2+co$dasn2,to$dasn1+co$dasn1+to$aesn+co$aesn)

dif.alphass=alphass-nalphass
dif.alphasn=alphasn-nalphasn
dif.alpha=c(dif.alphass,dif.alphasn)

alphass=nalphass
alphasn=nalphasn

return(list(alphass=alphass,alphasn=alphasn,betass1=betass1,
betasn=betasn,betass0=betass0,difab=c(dif.alpha,dif.beta)))
}

sace=function(da,y.name,tr.name,xnames,cluster.name,dist){

da$y=da[,y.name]
da$yind=1
da$yind[is.na(da$y)]=0

da0=da[which(da[,tr.name]==0),]
da1=da[which(da[,tr.name]==1),]

da0$cid=da0[,cluster.name]
ugp0=sort(unique(da0[,cluster.name]))
lugp0=length(ugp0)
for(i in 1:lugp0){
da0$cid[which(da0[,cluster.name]==ugp0[i])]=i
}

da1$cid=da1[,cluster.name]
ugp1=sort(unique(da1[,cluster.name]))
lugp1=length(ugp1)
for(i in 1:lugp1){
da1$cid[which(da1[,cluster.name]==ugp1[i])]=i
}

io=ini(da1,da0,xnames)

if(dist=="binomial") fobj=bfit(da1,da0,io)
if(dist=="poisson") fobj=pfit(da1,da0,io)
if(dist=="normal") fobj=nfit(da1,da0,io)

return(fobj)
}

