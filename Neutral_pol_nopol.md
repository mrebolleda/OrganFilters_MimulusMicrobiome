##Paper:  Floral organs act as environmental filters and interact with pollinators to structure the yellow monkeyflower *Mimulus guttatus* floral microbiome.
###Code by: María Rebolleda Gómez
####Last updated: April 2018

Set directory, load packages and export data. 
```{r}
library(VennDiagram)
library(ggplot2)


setwd(directory)
OTU=read.table("16s/OTU_rarefied_organs.txt",sep="\t",check.names = FALSE)
metadata=read.table("16s/Metadata_subset.txt",sep="\t")

OTU=OTU[rownames(OTU)!="306"&rownames(OTU)!="311",]
metadata=metadata[metadata$SampleID!="306"&metadata$SampleID!="311",]

metadata$SampleID=as.character(metadata$SampleID)
dim(OTU)

tapply(metadata$Loc,metadata$Organ, length) #sample sizes 

#Color palettes:
cols=c("#9c755f","#f28e2b","#76b7b2")
c_neu=c("#f27a6a","#465b6c","#bab0ac")
```


```{r}
#Remove transplant and sample for each organ
subset=metadata$SampleID[metadata$Pollinator=="WP"&metadata$Organ=="Petals"]
OTU_WP_P=OTU[match(as.character(subset),rownames(OTU)),]
OTU_WP_P=as.matrix(OTU_WP_P)

subset=metadata$SampleID[metadata$Pollinator=="NP"&metadata$Organ=="Petals"]
OTU_NP_P=OTU[match(as.character(subset),rownames(OTU)),]
OTU_NP_P=as.matrix(OTU_NP_P)

subset=metadata$SampleID[metadata$Pollinator=="WP"&metadata$Organ=="Anthers"]
OTU_WP_A=OTU[match(as.character(subset),rownames(OTU)),]
OTU_WP_A=as.matrix(OTU_WP_A)

subset=metadata$SampleID[metadata$Pollinator=="NP"&metadata$Organ=="Anthers"]
OTU_NP_A=OTU[match(as.character(subset),rownames(OTU)),]
OTU_NP_A=as.matrix(OTU_NP_A)


subset=metadata$SampleID[metadata$Pollinator=="WP"&metadata$Organ=="Style"]
OTU_WP_S=OTU[match(as.character(subset),rownames(OTU)),]
OTU_WP_S=as.matrix(OTU_WP_S)


subset=metadata$SampleID[metadata$Pollinator=="NP"&metadata$Organ=="Style"]
OTU_NP_S=OTU[match(as.character(subset),rownames(OTU)),]
OTU_NP_S=as.matrix(OTU_NP_S)
```


Neutral fitting after Burns et al. 2016 aplication of Sloan et al. 2006. 
```{r}
#Adam Burns - 2/10/2015
#adburns@stanford.edu
#spp: A community table for communities of interest with local communities/samples as rows and taxa as columns. All samples must be rarefied to the same depth.
#pool: A community table for defining source community (optional; Default=NULL).
#taxon: A table listing the taxonomic calls for each otu, with OTU ids as row names and taxonomic classifications as columns.
#If stats=TRUE the function will return fitting statistics.
#If stats=FALSE the function will return a table of observed and predicted values for each otu.

sncm.fit <- function(spp, pool=NULL, stats=TRUE, taxon=NULL){
	require(minpack.lm)
	require(Hmisc)
	require(stats4)
	
	options(warn=-1)

	#Calculate the number of individuals per community
	N <- mean(apply(spp, 1, sum))
	
	#Calculate the average relative abundance of each taxa across communities
	if(is.null(pool)){
		p.m <- apply(spp, 2, mean)
		p.m <- p.m[p.m != 0]
		p <- p.m/N
	} else {
		p.m <- apply(pool, 2, mean)
		p.m <- p.m[p.m != 0]
		p <- p.m/N
	}

	#Calculate the occurrence frequency of each taxa across communities
	spp.bi <- 1*(spp>0)
	freq <- apply(spp.bi, 2, mean)
	freq <- freq[freq != 0]

	#Combine
	C <- merge(p, freq, by=0)
	C <- C[order(C[,2]),]
	C <- as.data.frame(C)
	C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),] #Removes rows with any zero (absent in either source pool or local communities)
	p <- C.0[,2]
	freq <- C.0[,3]
	names(p) <- C.0[,1]
	names(freq) <- C.0[,1]

	#Calculate the limit of detection
	d = 1/N

	##Fit model parameter m (or Nm) using Non-linear least squares (NLS)
	m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE), start=list(m=0.1))
	m.ci <- confint(m.fit, 'm', level=0.95)
	
	##Fit neutral model parameter m (or Nm) using Maximum likelihood estimation (MLE)
	sncm.LL <- function(m, sigma){
		R = freq - pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE)
		R = dnorm(R, 0, sigma)
		-sum(log(R))
	}
	m.mle <- mle(sncm.LL, start=list(m=0.1, sigma=0.1), nobs=length(p))
	
	##Calculate Akaike's Information Criterion (AIC)
	aic.fit <- AIC(m.mle, k=2)
	bic.fit <- BIC(m.mle)

	N=as.integer(N) #For binomial N needs to be integer
	
	##Calculate goodness-of-fit (R-squared and Root Mean Squared Error)
	freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1-p), lower.tail=FALSE)
	Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
	RMSE <- sqrt(sum((freq-freq.pred)^2)/(length(freq)-1))
	
	pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=F)

	
	##Calculate AIC for binomial model
	bino.LL <- function(mu, sigma){
		R = freq - pbinom(d, N, p, lower.tail=FALSE)
		R = dnorm(R, mu, sigma)
		-sum(log(R))
	}
	bino.mle <- mle(bino.LL, start=list(mu=0.1, sigma=0.1), nobs=length(p))
	
	aic.bino <- AIC(bino.mle, k=2)
	bic.bino <- BIC(bino.mle)
	
	##Goodness of fit for binomial model
	bino.pred <- pbinom(d, N, p, lower.tail=FALSE)
	Rsqr.bino <- 1 - (sum((freq - bino.pred)^2))/(sum((freq - mean(freq))^2))
	RMSE.bino <- sqrt(sum((freq - bino.pred)^2)/(length(freq) - 1))

	bino.pred.ci <- binconf(bino.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
	
	##Calculate AIC for Poisson model
	pois.LL <- function(mu, sigma){
		R = freq - ppois(d, N*p, lower.tail=FALSE)
		R = dnorm(R, mu, sigma)
		-sum(log(R))
	}
	pois.mle <- mle(pois.LL, start=list(mu=0, sigma=0.1), nobs=length(p))
	
	aic.pois <- AIC(pois.mle, k=2)
	bic.pois <- BIC(pois.mle)
	
	##Goodness of fit for Poisson model
	pois.pred <- ppois(d, N*p, lower.tail=FALSE)
	Rsqr.pois <- 1 - (sum((freq - pois.pred)^2))/(sum((freq - mean(freq))^2))
	RMSE.pois <- sqrt(sum((freq - pois.pred)^2)/(length(freq) - 1))

	pois.pred.ci <- binconf(pois.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)

	##Results
	if(stats==TRUE){
		fitstats <- data.frame(m=numeric(), m.ci=numeric(), m.mle=numeric(), maxLL=numeric(), binoLL=numeric(),
		                       poisLL=numeric(), Rsqr=numeric(), Rsqr.bino=numeric(), Rsqr.pois=numeric(), RMSE=numeric(), RMSE.bino=numeric(), RMSE.pois=numeric(), AIC=numeric(), BIC=numeric(), AIC.bino=numeric(), BIC.bino=numeric(), AIC.pois=numeric(), BIC.pois=numeric(), N=numeric(), Samples=numeric(), Richness=numeric(), Detect=numeric())
		fitstats[1,] <- c(coef(m.fit), coef(m.fit)-m.ci[1], m.mle@coef['m'], m.mle@details$value, bino.mle@details$value, pois.mle@details$value, Rsqr, Rsqr.bino, Rsqr.pois, RMSE, RMSE.bino, RMSE.pois, aic.fit, bic.fit, aic.bino, bic.bino, aic.pois, bic.pois, N, nrow(spp), length(p), d)
		return(fitstats)
	} else {
		A <- cbind(p, freq, freq.pred, pred.ci[,2:3], bino.pred, bino.pred.ci[,2:3])
		A <- as.data.frame(A)
		colnames(A) <- c('p', 'freq', 'freq.pred', 'pred.lwr', 'pred.upr', 'bino.pred', 'bino.lwr', 'bino.upr')
		if(is.null(taxon)){
			B <- A[order(A[,1]),]
		} else {
			B <- merge(A, taxon, by=0, all=TRUE)
			row.names(B) <- B[,1]
			B <- B[,-1]
			B <- B[order(B[,1]),]
		}
		return(B)
	}
}


list.all=list(OTU_WP_P,OTU_NP_P,OTU_WP_A,OTU_NP_A,OTU_WP_S,OTU_NP_S)
BL=lapply(list.all,  sncm.fit,pool = OTU,stats = F)


B=data.frame(rbind(BL[[1]],BL[[2]],BL[[3]],BL[[4]],BL[[5]],BL[[6]]),
             Organ=c(rep("Petals",nrow(BL[[1]])+nrow(BL[[2]])),
             rep("Anthers",nrow(BL[[3]])+nrow(BL[[4]])),
             rep("Style",nrow(BL[[5]])+nrow(BL[[6]]))),
             Pollinator=c(rep("WP",nrow(BL[[1]])),rep("NP",nrow(BL[[2]])),
                          rep("WP",nrow(BL[[3]])),rep("NP",nrow(BL[[4]])),
                          rep("WP",nrow(BL[[5]])),rep("NP",nrow(BL[[6]]))),
             OTU=c(rownames(BL[[1]]),rownames(BL[[2]]),
             rownames(BL[[3]]),rownames(BL[[4]]),
             rownames(BL[[5]]),rownames(BL[[6]])))
head(B)

B$color="Neutral"
B$color[B$freq>B$pred.upr]="Above"
B$color[B$freq<B$pred.lwr]="Below"



ggplot(data=B, aes(x=log(p),y=freq, colour=color))+
  facet_wrap(~Pollinator*Organ)+
  geom_point()+
  geom_line(aes(y=freq.pred),color="black")+
  xlab("log(Mean relative abundance)")+
  ylab("Ocurrence frequency")+
  geom_line(aes(y=pred.lwr), linetype=2,color="black")+
  geom_line(aes(y=pred.upr), linetype=2,color="black")+
  scale_color_manual(values=c_neu)+
  theme_bw()


write.csv(B,"Neutral_OTU_all_ByPollinator.csv")

StatsL=lapply(list.all,  sncm.fit,pool = OTU,stats = T)
Stats=data.frame(rbind(StatsL[[1]],StatsL[[2]],StatsL[[3]],StatsL[[4]],StatsL[[5]],StatsL[[6]]),Organ=rep(c("Petals","Anthers","Style"), each=2),Pollinator=rep(c("WP","NP"),3))

write.csv(Stats,"Neutral_stats_all.csv")

```

Compare differences in proportions for treatment. 
```{r}
Part_bypol=paste(B$color,B$Pollinator)
tpol=table(Part_bypol)
countspol=matrix(tpol,nrow=2)
chisq.test(countspol)
```

Get proportion of over/under per organ and treatment
```{r}
Part_org=paste(B$Organ,B$Pollinator,B$color)
torg=table(Part_org)
countsorg=matrix(torg,nrow=3)
percent=countsorg/rep(colSums(countsorg),each=3)

stamens=rbind(countsorg[,1],countsorg[,2])
```

Get the proportion of OTUs that is shared across at least two organs in the different partitions and compare against a null model assuming that each partition is just a random sample of OTUs associated with that organ of the size of the partition. 
```{r}
####Observed values
#Create two data sets
Bpol=B[B$Pollinator=="WP",]
Bnpol=B[B$Pollinator=="NP",]

#1. Separate data according to each partition
splt_pol=split(Bpol,Bpol$color)
splt_npol=split(Bnpol,Bnpol$color)
Organs=c("Petals","Anthers","Style")

#2. Get the list of OTUs for each organ in each partition
sub_org=function(x){
  OTU_list=vector(mode="list",length = 3L)
  for (i in 1:3){
    OTU_list[[i]]=x$OTU[x$Organ==Organs[i]]
  }
  OTU_list
}

part_org_wp=lapply(splt_pol,sub_org)
part_org_np=lapply(splt_npol,sub_org)

#3. Get the partition
getparts=function(x){
  part_above=get.venn.partitions(x$Above)
  part_bel=get.venn.partitions(x$Below)
  part_neu=get.venn.partitions(x$Neutral)
  obserM=matrix(c(sum(part_above$..count..[c(4,6,7)]),
                    sum(part_above$..count..[c(1,2,3,5)]),
                    sum(part_bel$..count..[c(4,6,7)]),
                    sum(part_bel$..count..[c(1,2,3,5)]),
                    sum(part_neu$..count..[c(4,6,7)]),
                    sum(part_neu$..count..[c(1,2,3,5)])),ncol=3)
  obs=obserM[2,]/colSums(obserM)
  observed_data=data.frame(partition=c("Above","Below","Neutral"),obs_values=obs)
  observed_data
}

obs_pol=getparts(part_org_wp)
obs_npol=getparts(part_org_np)


####Null model
#1. Function to sample each partition 
sample_neutr=function (x,organ,n){
  rdm=vector("list",length(n))
  for (i in 1:length(n)){
    rdm[[i]]=sample(x$OTU[x$Organ==organ],size = n[i])
  }
  
  
  rdm
}

#2. Get sample sizes for sample
expect=function(x){
  N=table(paste(x$Organ,x$color))
  exp_null=matrix(ncol=3,nrow=10000)
  
  for (j in 1:10000){
    null=matrix(nrow=3,ncol=2)
    pet=sample_neutr(x,"Petals",N[4:6])
    stam=sample_neutr(x,"Anthers",N[1:3])
    styles=sample_neutr(x,"Style",N[7:9])
    for (i in 1:3){
      tmpl=list(pet[[i]],stam[[i]],styles[[i]])
      vpart=get.venn.partitions(tmpl)
      null[i,1]=sum(vpart$..count..[c(4,6,7)])
      null[i,2]=sum(vpart$..count..[c(1,2,3,5)])
      }
  exp_null[j,]=null[,2]/rowSums(null)
  }
  
expec_data=data.frame(simul=as.vector(exp_null),partition=rep(c("Above","Below","Neutral"),each=10000))
expec_data
}

exp_pol=expect(Bpol)
exp_npol=expect(Bnpol)

ggplot(exp_pol, aes(x=simul, fill=partition))+
  facet_grid(.~ partition)+
  geom_density(aes(y=..scaled..))+ 
  scale_fill_manual(values=c_neu)+
  geom_vline(data=obs_pol, aes(xintercept=obs_values),linetype="dashed")+
  ylab("Frequency (null distribution)")+
  xlab("Proportion of OTUs shared across floral organs")+
  theme_bw()

ggplot(exp_npol, aes(x=simul, fill=partition))+
  facet_grid(.~ partition)+
  geom_density(aes(y=..scaled..))+ 
  scale_fill_manual(values=c_neu)+
  geom_vline(data=obs_npol, aes(xintercept=obs_values),linetype="dashed")+
  ylab("Frequency (null distribution)")+
  xlab("Proportion of OTUs shared across floral organs")+
  theme_bw()


all=data.frame(rbind(exp_pol,exp_npol),Pollinator=rep(c("wp","np"),each=30000))
all_obs=data.frame(rbind(obs_pol,obs_npol),Pollinator=rep(c("wp","np"),each=3))

ggplot(all, aes(x=simul, fill=partition))+
  facet_grid(Pollinator~ partition)+
  geom_density(aes(y=..scaled..))+ 
  scale_fill_manual(values=c_neu)+
  geom_vline(data=all_obs, aes(xintercept=obs_values),linetype="dashed")+
  ylab("Frequency (null distribution)")+
  xlab("Proportion of OTUs shared across floral organs")+
  theme_bw()



```


```{r}
M_exp_pol=tapply(exp_pol$simul,exp_pol$partition,mean)
M_exp_npol=tapply(exp_npol$simul,exp_npol$partition,mean)

dist=data.frame(dist=c(obs_pol$obs_values-M_exp_pol,obs_npol$obs_values-M_exp_npol), partition=rep(c("Above","Below","Neutral"),2), Pollinator=rep(c("wp","np"),each=3))

ggplot(dist, aes(x=Pollinator, y=dist, fill=partition))+
  geom_point(shape=21,size=2)+
  scale_fill_manual(values=c_neu)+
  ylab("Distance from expectation")+
  xlab("Pollinator treatment")+
  theme_bw()

```

