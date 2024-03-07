btreu=function(y1i,m1i1,mss1,msn1,taut,mpss1,mpsn1,nms){
	
d1=d2=matrix(0,m1i1,nms)
fyij=rep(0,nms)
sui=rnorm(nms,0,sqrt(taut))
pbss=pbss2=pbsn=pbsn2=matrix(0,m1i1,nms)

for(k in 1:nms){
	uk=sui[k]
pyss=exp(mss1+uk)
pbss[,k]=pyss/(1+pyss)
pbss2[,k]=pbss[,k]/(1+pyss)

pysn=exp(msn1+uk)
pbsn[,k]=pysn/(1+pysn)
pbsn2[,k]=pbsn[,k]/(1+pysn)

	d1[,k]=dbinom(y1i,1,prob=pbss[,k])*mpss1
	d2[,k]=dbinom(y1i,1,prob=pbsn[,k])*mpsn1
	fij=d1[,k]+d2[,k]
	fyij[k]=prod(fij)
	}

fyi=mean(fyij)
sym=matrix(fyij,nrow=m1i1,ncol=nms,byrow=T)

if(fyi==0) {eui2=eui=0
            epuiss=rowMeans(pbss)
            epuiss2=rowMeans(pbss2)
            epuisn=rowMeans(pbsn)
            epuisn2=rowMeans(pbsn2)  } else {
eui=mean(sui*fyij)/fyi
eui2=mean(sui^2*fyij)/fyi
epuiss=rowMeans(pbss*sym)/fyi
epuiss2=rowMeans(pbss2*sym)/fyi
epuisn=rowMeans(pbsn*sym)/fyi
epuisn2=rowMeans(pbsn2*sym)/fyi
}

rd1=rowMeans(d1)
etass=rd1/(rd1+rowMeans(d2))
etass[is.na(etass)]=1/2
                              
return(list(eui=eui,eui2=eui2,epuiss=epuiss,epuiss2=epuiss2,
epuisn=epuisn,epuisn2=epuisn2,etass=etass))
}

btr=function(xm1,sind1,y1,n1,m1v,m1v1,
alphass,alphasn,betass1,betasn,taut,gat,p,nms,xnames){
	
bssp=bsnp=matrix(0,p,p)
bssy=bsny=matrix(0,p,1)
dass1=dasn1=matrix(0,p,1)
dass2=dasn2=matrix(0,p,p)
aess=aesn=matrix(0,p,1)

nbyi1=byi1=ta1=pssi1=psni1=ga1=0

ru=rnorm(nms,0,sqrt(taut))

for(i in 1:n1){
	
cind=which(xm1$cid==i)
sind1i=sind1[[i]]
	
x1i=as.matrix(xm1[cind,xnames,drop=F]) # m_{1,i} x p
tx1i=t(x1i) # p x m_{1,i}

m1i1=m1v1[i]
m1i=m1v[i]

xss=x1i%*%alphass
xsn=x1i%*%alphasn

######### E step

pmod1=trpmod(sind1i,xss,xsn,gat,m1i,nms)
mpss=pmod1$mpss
mpsn=pmod1$mpsn

msdss=pmod1$msdess1
msdsn=pmod1$msdesn1
#pnn=pmod1$pnn

ga1=ga1+pmod1$ev1i2
mss=x1i%*%betass1

emss=exp(outer(mss,ru,"+"))
eemss=rowMeans((emss/(1+emss)))

nbyi1=nbyi1+sum(mpss*eemss)
pssi1=pssi1+sum(mpss)
psni1=psni1+sum(mpsn)

if(m1i1>0){

x1i1=x1i[sind1i,,drop=F] # m_{1,i,1} x p
tx1i1=t(x1i1)  # p x m_{1,i,1}

mss1=mss[sind1i]
msn1=x1i1%*%betasn
	
mpss1=mpss[sind1i]
mpsn1=mpsn[sind1i]

y1i=y1[cind][sind1i]

#######

euo=btreu(y1i,m1i1,mss1,msn1,taut,mpss1,mpsn1,nms)

########### survivors 

#### beta
etass=euo$etass
etasn=1-etass
metai=matrix(etass,nrow=p,ncol=m1i1,byrow=T)

metaiss=metai*matrix(euo$epuiss2,nrow=p,ncol=m1i1,byrow=T)
etaix1i1ss=tx1i1*metaiss

bssp=bssp+etaix1i1ss%*%x1i1
bssy=bssy+tx1i1%*%(etass*(y1i-euo$epuiss))

metaisn=matrix(etasn*euo$epuisn2,nrow=p,ncol=m1i1,byrow=T)
etai1x1i1sn=tx1i1*metaisn

bsnp=bsnp+etai1x1i1sn%*%x1i1
bsny=bsny+tx1i1%*%(etasn*(y1i-euo$epuisn))

##### alpha

etaix1i1=tx1i1*metai
etai1x1i1=tx1i1-etaix1i1
aess=aess+rowSums(etaix1i1)
aesn=aesn+rowSums(etai1x1i1)

#####

ta1=ta1+euo$eui2
esu=exp(mss+euo$eui)
byi1=byi1+sum(mpss*(esu/(1+esu)))

} else {bssp=bssp
        bssy=bssy
        bsnp=bsnp
        bsny=bsny
        aess=aess
        aesn=aesn
        byi1=byi1
        ta1=ta1}
        
############### all

ao=alp(x1i,tx1i,mpss,msdss,mpsn,msdsn,m1i,p)
dass1=dass1+ao$dass1
dasn1=dasn1+ao$dasn1
dass2=dass2+ao$dass2
dasn2=dasn2+ao$dasn2
}

nbetass1=betass1+solve(bssp,bssy)
nbetasn=betasn+solve(bsnp,bsny)

return(list(nbetass1=nbetass1,nbetasn=nbetasn,dass1=dass1,
dass2=dass2,dasn1=dasn1,dasn2=dasn2,ta1=ta1,
aess=aess,aesn=aesn,s1=byi1/pssi1,ns1=nbyi1/pssi1,
pssi1=pssi1,psni1=psni1,ga1=ga1))
}

bconeu=function(y0i,m0i1,mss01,taut,nms){

sui=rnorm(nms,mean=0,sd=sqrt(taut))
fyij0=rep(0,nms)

pbss0=pbss02=matrix(0,m0i1,nms)

for(k in 1:nms){
	uk=sui[k]
pyss0=exp(mss01+uk)
pbss0[,k]=pyss0/(1+pyss0)
pbss02[,k]=pbss0[,k]/(1+pyss0)
	d0=dbinom(y0i,1,prob=pbss0[,k])
	fyij0[k]=prod(d0)
	}

fyi0=mean(fyij0)
sym=matrix(fyij0,nrow=m0i1,ncol=nms,byrow=T)

if(fyi0==0) {eui02=eui0=0
            epuiss0=rowMeans(pbss0)
            epuiss02=rowMeans(pbss02)
 } else {
eui0=mean(sui^2*fyij0)/fyi0
eui02=mean(sui^2*fyij0)/fyi0
epuiss0=rowMeans(pbss0*sym)/fyi0
epuiss02=rowMeans(pbss02*sym)/fyi0
}

return(list(eui02=eui02,eui0=eui0,epuiss0=epuiss0,epuiss02=epuiss02))
}

bcon=function(xm0,sind0,y0,n0,m0v,m0v1,
alphass,alphasn,betass0,taut,gat,p,nms,xnames){
	
bssp0=matrix(0,p,p)
bssy0=matrix(0,p,1)
dasn1=dass1=matrix(0,p,1)
dasn2=dass2=matrix(0,p,p)
aess=aesn=matrix(0,p,1)

ta0=nbyi0=ga0=byi0=pssi0=psni0=0

ru=rnorm(nms,0,sqrt(taut))

for(i in 1:n0){

cind=which(xm0$cid==i)	
sind0i=sind0[[i]]

x0i=as.matrix(xm0[cind,xnames,drop=F]) # m_{0,i} x p
tx0i=t(x0i) # p x m_{0,i}

m0i=m0v[i]
m0i1=m0v1[i]
m0i0=m0i-m0i1

xss=x0i%*%alphass
xsn=x0i%*%alphasn

pmod0=conpmod(sind0i,xss,xsn,gat,m0i,nms)
mpss0=pmod0$mpss
mpsn0=pmod0$mpsn

msdss0=pmod0$msdess
msdsn0=pmod0$msdesn
#pnn0=pmod0$pnn

ga0=ga0+pmod0$ev0i2
mss0=x0i%*%betass0

emss=exp(outer(mss0,ru,"+"))
eemss=rowMeans((emss/(1+emss)))

nbyi0=nbyi0+sum(mpss0*eemss)
pssi0=pssi0+sum(mpss0)
psni0=psni0+sum(mpsn0)

########### non-survivors

if(m0i0>0) {
	
##### alpha

if(m0i1>0){
	
x0i0=x0i[-sind0i,,drop=F] # m_{0,i,0} x p
etai=pmod0$meta0[-sind0i]} else { x0i0=x0i
                                  etai=pmod0$meta0}

tx0i0=t(x0i0)  # p x m_{0,i,0}
metai=matrix(etai,nrow=p,ncol=m0i0,byrow=T)
p0=rowSums(tx0i0*metai)
aesn=aesn+p0} else {aesn=aesn}

########### survivors 

if(m0i1>0){
	
x0i1=x0i[sind0i,,drop=F] # m_{0,i,1} x p
tx0i1=t(x0i1)  # p x m_{0,i,1}

mss01=mss0[sind0i]

y0i=y0[cind][sind0i]

eu0=bconeu(y0i,m0i1,mss01,taut,nms)

etaix0i1ss0=tx0i1*matrix(eu0$epuiss02,nrow=p,ncol=m0i1,byrow=T)

#### beta

bssp0=bssp0+etaix0i1ss0%*%x0i1
bssy0=bssy0+tx0i1%*%(y0i-eu0$epuiss0)

##### alpha

aess=aess+rowSums(tx0i1)

ta0=ta0+eu0$eui02
esu=exp(mss0+eu0$eui)
byi0=byi0+sum(mpss0*(esu/(1+esu)))

} else {bssp0=bssp0
        bssy0=bssy0
        aess=aess
        byi0=byi0
        ta0=ta0
}

############### all

ao=alp(x0i,tx0i,mpss0,msdss0,mpsn0,msdsn0,m0i,p)
dass1=dass1+ao$dass1
dasn1=dasn1+ao$dasn1
dass2=dass2+ao$dass2
dasn2=dasn2+ao$dasn2

}

nbetass0=betass0+solve(bssp0,bssy0)

return(list(nbetass0=nbetass0,dass1=dass1,dass2=dass2,
dasn1=dasn1,dasn2=dasn2,aess=aess,aesn=aesn,ta0=ta0,
s0=byi0/pssi0,ns0=nbyi0/pssi0,pssi0=pssi0,psni0=psni0,
ga0=ga0))
}

#######################

bfit=function(da1,da0,io){

nms=100

alphass=io$alphass
alphasn=io$alphasn
betass1=io$betass1
betasn=io$betasn
betass0=io$betass0

sat=io$sat
gat=io$gat
taut=io$taut

dif=1
nit=35
thr=1e-3
nuit=0

while(any(dif>thr)&(nuit<nit)){

to=btr(io$xm1,io$sind1,io$y1,io$n1,io$m1v,io$m1v1,
alphass,alphasn,betass1,betasn,taut,gat,io$p,nms,xnames)

co=bcon(io$xm0,io$sind0,io$y0,io$n0,io$m0v,io$m0v1,
alphass,alphasn,betass0,taut,gat,io$p,nms,xnames)

dfo=dff(to,co,betass1,betasn,betass0,alphass,alphasn)

alphass=dfo$alphass
alphasn=dfo$alphasn
betass1=dfo$betass1
betasn=dfo$betasn
betass0=dfo$betass0
difab=dfo$difab

ntaut=(to$ta1+co$ta0)/(io$n1+io$n0)
ngat=(to$ga1+co$ga0)/(io$n1+io$n0)

dif.v=c(taut-ntaut,ngat-gat)

taut=ntaut
gat=ngat

dif=abs(c(dif.v,difab))
nuit=nuit+1
}

nsace=to$ns1-co$ns0

pss=(to$pssi1+co$pssi0)/(io$m1+io$m0)
psn=(to$psni1+co$psni0)/(io$m1+io$m0)
pnn=1-pss-psn

return(nsace)
}


