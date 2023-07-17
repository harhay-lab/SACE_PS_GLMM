ntreta=function(y1i,mss1,msn1,taut,sat,mpss1,mpsn1){

jointsd=sqrt(sat+taut)

nume=dnorm(y1i,mean=mss1,sd=jointsd)*mpss1
dsn=dnorm(y1i,mean=msn1,sd=jointsd)*mpsn1
eta1i=nume/(nume+dsn)
eta1i[is.na(eta1i)]=1/2

return(eta1i)
}

ntreu=function(y1i,mss1,msn1,taut,sat,mpss1,mpsn1,nms){

sdtau=sqrt(taut)
sdsa=sqrt(sat) 

sui=rnorm(nms,mean=0,sd=sdtau)
fyij=rep(0,nms)

for(k in 1:nms){
	uk=sui[k]
	d1=dnorm(y1i,mean=mss1+uk,sd=sdsa)*mpss1
	d2=dnorm(y1i,mean=msn1+uk,sd=sdsa)*mpsn1
	fij=d1+d2
	fyij[k]=prod(fij)}

fyi=mean(fyij)

if(fyi==0){
eui=eui2=vui=0
} else{
eui=mean(sui*fyij)/fyi
eui2=mean(sui^2*fyij)/fyi
vui=eui2-(eui)^2}

return(list(vui=vui,eui2=eui2,eui=eui))
}

ntr=function(xm1,sind1,y1,n1,m1v,m1v1,
alphass,alphasn,betass1,betasn,taut,sat,gat,p,nms,xnames){
	
bssp=bsnp=matrix(0,p,p)
bssy=bsny=matrix(0,p,1)
dass1=dasn1=matrix(0,p,1)
dass2=dasn2=matrix(0,p,p)
aess=aesn=matrix(0,p,1)

den1=0
ta1=ga1=0
pssi1=nfyi1=fyi1=psni1=0

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
#mpnn=pmod1$mpnn

mpss1=mpss[sind1i]
mpsn1=mpsn[sind1i]

msdss=pmod1$msdess1
msdsn=pmod1$msdesn1

ga1=ga1+pmod1$ev1i2
mss=x1i%*%betass1

nfyi1=nfyi1+sum(mpss*mss)
pssi1=pssi1+sum(mpss)
psni1=psni1+sum(mpsn)

if(m1i1>0){

x1i1=x1i[sind1i,,drop=F] # m_{1,i,1} x p
tx1i1=t(x1i1)  # p x m_{1,i,1}

mss1=mss[sind1i]
msn1=x1i1%*%betasn

y1i=y1[cind][sind1i]

etai=ntreta(y1i,mss1,msn1,taut,sat,mpss1,mpsn1)

########### survivors 

euo=ntreu(y1i,mss1,msn1,taut,sat,mpss1,mpsn1,nms)
eu1i=euo$eui
cyi=y1i-eu1i

#### beta

metai=matrix(etai,nrow=p,ncol=m1i1,byrow=T)
etaix1i1=tx1i1*metai
etai1x1i1=tx1i1-etaix1i1

bssp=bssp+etaix1i1%*%x1i1
bssy=bssy+etaix1i1%*%cyi

bsnp=bsnp+etai1x1i1%*%x1i1
bsny=bsny+etai1x1i1%*%cyi

##### alpha

aess=aess+rowSums(etaix1i1)
aesn=aesn+rowSums(etai1x1i1)

######

ta1=ta1+euo$eui2

ty1ss=etai*(y1i-mss1-eu1i)^2
ty1sn=(1-etai)*(y1i-msn1-eu1i)^2
den1=den1+sum(ty1ss+ty1sn+euo$vui)
fyi1=fyi1+sum(mpss*(mss+eu1i))

} else{bssp=bssp
       bssy=bssy
       bsnp=bsnp
       bsny=bsny
       aess=aess
       aesn=aesn
       fyi1=fyi1
       ta1=ta1}

############### all

ao=alp(x1i,tx1i,mpss,msdss,mpsn,msdsn,m1i,p)
dass1=dass1+ao$dass1
dasn1=dasn1+ao$dasn1
dass2=dass2+ao$dass2
dasn2=dasn2+ao$dasn2
}

nbetass1=lm.fit(bssp,bssy)$coef
nbetasn=lm.fit(bsnp,bsny)$coef

return(list(nbetass1=nbetass1,nbetasn=nbetasn,dass1=dass1,
dass2=dass2,dasn1=dasn1,dasn2=dasn2,aess=aess,aesn=aesn,
s1=fyi1/pssi1,ns1=nfyi1/pssi1,pssi1=pssi1,psni1=psni1,
ta1=ta1,ga1=ga1,den1=den1))
}

nconeu=function(y0i,m0i1,mss01,taut,sat){

coef=taut/(m0i1*taut+sat)
eui=coef*sum(y0i-mss01)
vui=coef*sat
eui2=vui+eui^2
return(list(vui=vui,eui2=eui2,eui=eui))
}

ncon=function(xm0,sind0,y0,n0,m0v,m0v1,
alphass,alphasn,betass0,taut,sat,gat,p,nms,xnames){
	
bssp0=matrix(0,p,p)
bssy0=matrix(0,p,1)
dasn1=dass1=matrix(0,p,1)
dasn2=dass2=matrix(0,p,p)
aess=aesn=matrix(0,p,1)

den0=0
ta0=ga0=0
pssi0=nfyi0=fyi0=psni0=0

for(i in 1:n0){

cind=which(xm0$cid==i)	
sind0i=sind0[[i]]

x0i=as.matrix(xm0[cind,xnames,drop=F]) # m_{0,i} x p
tx0i=t(x0i) # p x m_{0,i}
m0i1=m0v1[i]
m0i=m0v[i]
m0i0=m0i-m0i1

xss=x0i%*%alphass
xsn=x0i%*%alphasn

pmod0=conpmod(sind0i,xss,xsn,gat,m0i,nms)
mpss0=pmod0$mpss
mpsn0=pmod0$mpsn
#mpnn0=pmod0$mpnn
	
msdss0=pmod0$msdess0
msdsn0=pmod0$msdesn0

ga0=ga0+pmod0$ev0i2
mss0=x0i%*%betass0

nfyi0=nfyi0+sum(mpss0*mss0)
pssi0=pssi0+sum(mpss0)
psni0=psni0+sum(mpsn0)

########### non-survivors

if(m0i0>0){

##### alpha

if(m0i1>0){
x0i0=x0i[-sind0i,,drop=F] # m_{0,i,0} x p
etai=pmod0$meta0[-sind0i]} else {
x0i0=x0i
etai=pmod0$meta0
}

tx0i0=t(x0i0)  # p x m_{0,i,0} 
metai=matrix(etai,nrow=p,ncol=m0i0,byrow=T)
p0=rowSums(tx0i0*metai)
aesn=aesn+p0 } else {aesn=aesn}

########### survivors 

if(m0i1>0){

x0i1=x0i[sind0i,,drop=F] # m_{0,i,1} x p
tx0i1=t(x0i1)  # p x m_{0,i,1}

mss01=mss0[sind0i]

y0i=y0[cind][sind0i]

eu0=nconeu(y0i,m0i1,mss01,taut,sat)
eu0i=eu0$eui
cyi=y0i-eu0i

#### beta

bssp0=bssp0+tx0i1%*%x0i1
bssy0=bssy0+tx0i1%*%cyi

##### alpha

aess=aess+rowSums(tx0i1)

ta0=ta0+eu0$eui2
ty0ss=(y0i-mss01-eu0i)^2
den0=den0+sum(ty0ss+eu0$vui)
fyi0=fyi0+sum(mpss0*(mss0+eu0i))

} else{bssp0=bssp0
       bssy0=bssy0
       aess=aess
       fyi0=fyi0
       ta0=ta0
}

############### all

ao=alp(x0i,tx0i,mpss0,msdss0,mpsn0,msdsn0,m0i,p)
dass1=dass1+ao$dass1
dasn1=dasn1+ao$dasn1
dass2=dass2+ao$dass2
dasn2=dasn2+ao$dasn2
}

nbetass0=lm.fit(bssp0,bssy0)$coef

return(list(nbetass0=nbetass0,dass1=dass1,dass2=dass2,
dasn1=dasn1,dasn2=dasn2,aess=aess,aesn=aesn,
s0=fyi0/pssi0,ns0=nfyi0/pssi0,pssi0=pssi0,psni0=psni0,
den0=den0,ta0=ta0,ga0=ga0))
}


#######################

nfit=function(da1,da0,io){

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

to=ntr(io$xm1,io$sind1,io$y1,io$n1,io$m1v,io$m1v1,
alphass,alphasn,betass1,betasn,taut,sat,gat,io$p,nms,xnames)

co=ncon(io$xm0,io$sind0,io$y0,io$n0,io$m0v,io$m0v1,
alphass,alphasn,betass0,taut,sat,gat,io$p,nms,xnames)

dfo=dff(to,co,betass1,betasn,betass0,alphass,alphasn)

alphass=dfo$alphass
alphasn=dfo$alphasn
betass1=dfo$betass1
betasn=dfo$betasn
betass0=dfo$betass0
difab=dfo$difab

nsat=(to$den1+co$den0)/(io$m11+io$m01)
ntaut=(to$ta1+co$ta0)/(io$n1+io$n0)
ngat=(to$ga1+co$ga0)/(io$n1+io$n0)

dif.v=c(sat-nsat,taut-ntaut,ngat-gat)

sat=nsat
taut=ntaut
gat=ngat

dif=abs(c(dif.v,difab))
nuit=nuit+1
}

sace=to$s1-co$s0
nsace=to$ns1-co$ns0

pss=(to$pssi1+co$pssi0)/(io$m1+io$m0)
psn=(to$psni1+co$psni0)/(io$m1+io$m0)
pnn=1-pss-psn

return(obj=c(nsace,sace,pss,psn,pnn,taut,
gat,sat,betass1,betasn,betass0,alphass,alphasn))
}
