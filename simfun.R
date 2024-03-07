library(MASS)

# simulate the cluster size
# msize: mean
# sdsize: standard deviation
# nc: number of clusters
# nsim: number of simulations

simn=function(msize,sdsize,nc,nsim){
	rn=matrix(round(rnorm(nc*nsim,msize,sdsize)),nc,nsim)
	return(rn)
}

# simulate the predictors
# ni: cluster size
# cid: cluster id
# simid: simulation index

simx=function(ni,cid,simid){
	x1=rbinom(ni,1,1/2)
	x2=rnorm(ni,0,1)
	xi=data.frame(x0=1,x1=x1,x2=x2,id=1:ni,cid=cid,simid=simid)	
	return(xi)
	}

#  principal stratum membership model
# vi: random intercept
# ni: cluster size
# xssi: xi%*%alpha_ss
# xsni: xi%*%alpha_sn

mlogit=function(xssi,xsni,vi,ni){

        ess=exp(xssi+vi)
        esn=exp(xsni+vi)
        pnni=1/(1+ess+esn)
        pssi=ess*pnni
        psni=esn*pnni
        sta=matrix(0,ni,3)
        for(j in 1:ni){
        sta[j,]=rmultinom(n=1,size=1,prob=c(pssi[j],psni[j],pnni[j]))
		}
		
	return(list(sta=sta,pssi=pssi,psni=psni,pnni=pnni))
		}

# simulate outcome in the treatment group
# xi: predictors 
# ui: random intercept in the outcome model
# vi: random intercept in the membership model
# mssi: xi%*%beta_ss
# msni: xi%*%beta_sn
# xssi: xi%*%alpha_ss
# xsni: xi%*%alpha_sn
# cid: cluster id
# simid: simulation index
simby1=function(xi,ui,vi,mssi,msni,xssi,xsni,cid,simid){

        ni=nrow(xi)
        pyess=exp(mssi+ui)
	yiss=rbinom(ni,1,prob=pyess/(1+pyess))

        pyesn=exp(msni+ui)
	yisn=rbinom(ni,1,prob=pyesn/(1+pyesn))

        mo=mlogit(xssi,xsni,vi,ni)
        sta=mo$sta

	yi=rowSums(cbind(yiss,yisn)*sta[,1:2])
	yind=rowSums(sta[,1:2])

	dai=data.frame(y=yi,yss=yiss,ysn=yisn,pss=mo$pssi,psn=mo$psni,
	pnn=mo$pnni,ss=sta[,1],sn=sta[,2],nn=sta[,3],id=1:ni,cid=cid,
	simid=simid,x0=1,x1=xi[,"x1"],x2=xi[,"x2"],yind=yind,u=ui,v=vi)		

return(dai)
}

# simulate outcome in the control group

simby0=function(xi0,ui0,vi0,mssi0,xssi0,xsni0,cid,simid){

        ni=nrow(xi0)
        pyb0=exp(mssi0+ui0)
	yi=rbinom(ni,1,prob=pyb0/(1+pyb0))			

        mo=mlogit(xssi0,xsni0,vi0,ni)
        sta=mo$sta
		
	dai=data.frame(y=yi*sta[,1],yss=yi,pss=mo$pssi,psn=mo$psni,
	pnn=mo$pnni,ss=sta[,1],sn=sta[,2],nn=sta[,3],id=1:ni,cid=cid,
	simid=simid,x0=1,x1=xi0[,"x1"],x2=xi0[,"x2"],yind=sta[,1],u=ui0,v=vi0)		

return(dai)
}

##### calculate the frequency of coverage
# tv: true value of the SACE
# bobj: SACE estimates from bootstrap samples

bcover=function(bobj,tv){

pb=c(0.025,0.975)

if(is.list(bobj)) { lb=length(bobj)
                    nbo=NULL
                    for(j in 1:lb){
                    bj=unlist(bobj[j],use.names=F)
                    indj=is.numeric(bj)
                    if(indj) nbo=c(nbo,bj)
                    }
                    bobj=nbo
                       }
ind=which(is.na(bobj))
if(length(ind)>0) bobj=bobj[-ind]

bq=quantile(bobj,probs=pb)

cds=((bq[1]<=tv)&(bq[2]>=tv))
return(cds)
}
