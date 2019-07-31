##Paper:  Floral organs act as environmental filters and interact with pollinators to structure the yellow monkeyflower *Mimulus guttatus* floral microbiome.
###Code by: María Rebolleda Gómez
####Last updated: April 2018



Set directory, load packages and export data. 
```{r, eval=FALSE}
library("vegan")
library("ggplot2")
library("tidyr")
library("phyloseq")
library("betapart")
library("multcomp")


treefile = "16s/rep_set.tre"
mymap <- import_qiime_sample_data("../Raw/MIGU_metadata.txt")
class(mymap)
dim(mymap) #199 12
biom_file = import_biom("16s/otu_table_filtered_no_euks.biom", treefile)
biom_file

#Save OTU table and taxonomic table as their own tables
OTU_all=otu_table(biom_file)
OTU_all=t(as.matrix(OTU_all@.Data))

#Remove samples with less than 1100 reads (See notebook Oct?)
keep.row=rowSums(OTU_all)>1200
OTU_min=OTU_all[keep.row,]
dim(OTU_min)
```

Once excluding chloroplast sequences, data numbers are low and can exagerate differences between communities. To get a better estimate of rarefied community I re-sampled 100 times. 
```{r, eval=FALSE}
#Remove tomato seed|samples for a different project and leave only washes. 
subset=mymap$SampleID[mymap$Organ!="Seed"&mymap$Type=="Wash"]
OTU_sub=OTU_min[match(as.character(subset),rownames(OTU_min)),]
remove=which(is.na(rownames(OTU_sub)))
OTU_sub=OTU_sub[-remove,]


#Create empty vector for simulation of rare
Rarefied_data=array(data=NA,dim=c(dim(OTU_sub),1000))

for (i in 1:1000){
  Rarefied_data[,,i]=rrarefy(OTU_sub,1200)
}

Rarefied_OTU=apply(Rarefied_data,c(1,2),mean)
dim(Rarefied_OTU)==dim(OTU_sub)
rownames(Rarefied_OTU)=rownames(OTU_sub)
colnames(Rarefied_OTU)=colnames(OTU_sub)

```


Test differentiation between soil and organs
```{r, eval=FALSE}
#Remove OTUs not present  
not_0=colSums(Rarefied_OTU)>0
ncol(Rarefied_OTU)-sum(not_0)
OTU_wash_clean=Rarefied_OTU[,not_0]
dim(OTU_wash_clean)

#Get metadadata in same order for correct labeling 
m=match(rownames(OTU_wash_clean),as.character(mymap$SampleID))
metadata=(mymap[m,])
nrow(metadata)
head(metadata)


#write.table(OTU_wash_clean,file="16s/OTU_rarefied_allwash.txt",sep="\t")
```


Make PCoA with all organs
```{r}
#Remove soil
subset=metadata$SampleID[metadata$Organ!="Soil"]
OTU_nosoil=OTU_wash_clean[match(as.character(subset),rownames(OTU_wash_clean)),]

#Remove OTUs not present  
not_0=colSums(OTU_nosoil)>0
ncol(OTU_nosoil)-sum(not_0)
OTU_nosoil_clean=OTU_nosoil[,not_0]
dim(OTU_nosoil_clean)

#Get metadadata in same order for correct labeling 
m=match(rownames(OTU_nosoil_clean),as.character(metadata$SampleID))
metadata_ns=(metadata[m,])
head(metadata_ns)


dist=function(OTU,tree){
  #Get Sorensen and Bray-Curtis distances
  distances=list(dist_soren=vegdist(OTU,method="bray",binary=T),
               dist_bray=vegdist(OTU,method="bray",binary=F))
  
  #Get unifrac distances:
  #1. Load tree and removed branches out of subset
  tree=read_tree(treefile)
  drop=tree$tip.label[!(tree$tip.label%in%colnames(OTU))]
  tree.new=drop.tip(tree,drop)
  if (dim(OTU)[2]!=length(tree.new$tip.label)){
    warning("Check OTU table and treefile")
  }
  #2. Load tree, OTU and taxa as phyloseq data
  taxa_s=taxa[match(colnames(OTU),rownames(taxa)),]
  phylo_data=phyloseq(sample_data(metadata),
                      otu_table(t(OTU),taxa_are_rows=T),
                      phy_tree(tree.new),tax_table(taxa_s))
  #3. Calculate distances and append to distances list
  distances$unifrac=phyloseq::distance(phylo_data,method="unifrac")
  distances$wunifrac=phyloseq::distance(phylo_data,method="wunifrac")
  distances
}


distances=dist(OTU_nosoil_clean,treefile)


#PCoA
PCoA_all=function(x){
  pcoa.x=cmdscale(x,eig=T)
  #add percentage of variation explained by each of the first two axis
  pcoa.x$v=c(round(pcoa.x$eig[1]/sum(pcoa.x$eig)*100),
              round(pcoa.x$eig[2]/sum(pcoa.x$eig)*100))
  pcoa.x
}

ord=lapply(distances,PCoA_all)
str(ord)

#Plot ordinations
shapes=c(19,1)
par(oma=c(1, 1, 1, 8), xpd=TRUE)
par(mar=c(4,4,4,4))
par(mfrow=c(2,2))

col_all=c("#e59423","#F77F71","#91B24D","#e3072a","gray20","#10a8b0")

for (i in 1:4){
    plot(ord[[i]]$points,col=col_all[metadata$Organ],pch=19,
       xlab=paste("PCoA1(",ord[[i]]$v[1],"%)",sep=""),
       ylab=paste("PCoA2 (",ord[[i]]$v[2],"%)",sep=""), main=names(ord)[i])
}

#Add legend to the 4-way plot
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend(0.58,0.85,legend=levels(metadata$Organ),pch=19,col=col_all,bty ="n",cex=1.2)


```



The *Distance_pool* function calculates the mean distance (Sorensen) and 95% CI between samples within a location (Seep)
and between these samples and the potential source pool of microbes. 
To run the formula requires a list of abundance matrices for each location with samples as rows and OTUs as columns, the Organ or community type within a location, the puttative source (vector of OTU abundances) and the name of the source pool. 

```{r, eval=FALSE}
pres=(OTU_wash_clean>0)*1

#Function to sample each organ of interest and pollinator or not
sampleOTU=function(df,Org,otu,is_p=F,pol){
  if (is_p==T){
  sub=df$SampleID[df$Organ==Org&df$Pollinator==pol]} else {
  sub=df$SampleID[df$Organ==Org]}
  OTU=otu[rownames(otu)%in%sub,]
  OTU
}

Distance_Pool=function(List,Org,is_p=F,pol,Pool,nameP){
  #subset OTU by organ 
  org_OTU=lapply(List,sampleOTU,Org=Org,otu=pres,is_p=is_p,pol=pol)
  s_org=mapply(rbind,Pool,org_OTU, SIMPLIFY = F)
  soren=lapply(s_org,beta.pair,index.family = "sorensen")
  rows=3
  seepdata=data.frame(seep=rep("NA",length(soren)*rows),
                      type=rep("NA",length(soren)*rows),
                      div_in=rep("NA",length(soren)*rows),
                      dist=rep(0,length(soren)*rows),
                      organ=rep(Org,length(soren)*rows), stringsAsFactors = F)
  for (i in 1:length(soren)){
    #Get seep name
    seep=names(soren)[i]
    n=nrow(as.matrix(soren[[seep]][[1]]))
    
    #Get distances (from soil/and samples)
    sim=c(mean(as.matrix(soren[[seep]][[1]])[-1,1]))
    sne=c(mean(as.matrix(soren[[seep]][[2]])[-1,1]))
    sor=c(mean(as.matrix(soren[[seep]][[3]])[-1,1]))
    
    seepdata$dist[((i-1)*3+1):(((i-1)*3)+3)]=c(sim,sne,sor)
    seepdata$seep[((i-1)*3+1):(((i-1)*3)+3)]=seep
    seepdata$type[((i-1)*3+1):(((i-1)*3)+3)]=rep(nameP,3)
    seepdata$div_in[((i-1)*3+1):(((i-1)*3)+3)]=c("turnover","nestedness","total")
  }
  seepdata
}

#1. Separate in a list the data.frame by seep or location.
seep_list=split(metadata,metadata$Seep)

#2. Sample soils,community and leaves. 
soils=lapply(seep_list,sampleOTU,Org="Soil",otu=pres)
leaves=lapply(seep_list,sampleOTU,Org="Leaves",otu=OTU_wash_clean)
leaves=lapply(leaves,colMeans)
leaves=lapply(leaves,function(x) (x>0)*1)


community=lapply(seep_list,sampleOTU,Org="Community",otu=OTU_wash_clean)
community=lapply(community,colMeans)
community=lapply(community,function(x) (x>0)*1)


sample_neg_OTU=function(df,org,otu){
  sub=df$SampleID[df$Organ!=org]
  OTU=otu[rownames(otu)%in%sub,]
  OTU
}

no_ant=lapply(seep_list,sample_neg_OTU,org="Anthers",otu=OTU_wash_clean)
no_ant=lapply(no_ant,colMeans)
no_ant=lapply(no_ant,function(x) (x>0)*1)

no_pet=lapply(seep_list,sample_neg_OTU,org="Petals",otu=OTU_wash_clean)
no_pet=lapply(no_pet,colMeans)
no_pet=lapply(no_pet,function(x) (x>0)*1)

no_st=lapply(seep_list,sample_neg_OTU,org="Styles",otu=OTU_wash_clean)
no_st=lapply(no_st,colMeans)
no_st=lapply(no_st,function(x) (x>0)*1)



#3. Apply function for each organ and pool for all the CONTROL treatments use-(List,Org,is_p=F,pol,Pool,nameP)
antherswp_soil=Distance_Pool(seep_list,"Anthers",is_p=T,pol="WP",soils,"soil")
petalswp_soil=Distance_Pool(seep_list,"Petals",is_p=T,pol="WP",soils,"soil")
stylewp_soil=Distance_Pool(seep_list,"Style",is_p=T,pol="WP",soils,"soil")


antherswp_leaves=Distance_Pool(seep_list,"Anthers",is_p=T,pol="WP",leaves,"leaves")
petalswp_leaves=Distance_Pool(seep_list,"Petals",is_p=T,pol="WP",leaves,"leaves")
stylewp_leaves=Distance_Pool(seep_list,"Style",is_p=T,pol="WP",leaves,"leaves")


community=community[c(1,3,5)] #For community remove empty seeps RH1 and RHB
seep_list_com=seep_list[c(1,3,5)]
antherswp_com=Distance_Pool(seep_list_com,"Anthers",is_p=T,pol="WP",community,"neighbors")
petalswp_com=Distance_Pool(seep_list_com,"Petals",is_p=T,pol="WP",community,"neighbors")
stylewp_com=Distance_Pool(seep_list_com,"Style",is_p=T,pol="WP",community,"neighbors")

antherswp_flower=Distance_Pool(seep_list,"Anthers",is_p=T,pol="WP",no_ant,"flower")
petalswp_flower=Distance_Pool(seep_list,"Petals",is_p=T,pol="WP",no_pet,"flower")
stylewp_flower=Distance_Pool(seep_list,"Style",is_p=T,pol="WP",no_st,"flower")

#3. Apply function for each organ and pool for all of the POLLINATOR EXCLUSION/ use-(List,Org,is_p=F,pol,Pool,nameP)

anthersnp_soil=Distance_Pool(seep_list,"Anthers",is_p=T,pol="NP",soils,"soil")
petalsnp_soil=Distance_Pool(seep_list,"Petals",is_p=T,pol="NP",soils,"soil")
stylenp_soil=Distance_Pool(seep_list,"Style",is_p=T,pol="NP",soils,"soil")


anthersnp_leaves=Distance_Pool(seep_list,"Anthers",is_p=T,pol="NP",leaves,"leaves")
petalsnp_leaves=Distance_Pool(seep_list,"Petals",is_p=T,pol="NP",leaves,"leaves")
stylenp_leaves=Distance_Pool(seep_list,"Style",is_p=T,pol="NP",leaves,"leaves")


anthersnp_com=Distance_Pool(seep_list_com,"Anthers",is_p=T,pol="NP",community,"neighbors")
petalsnp_com=Distance_Pool(seep_list_com,"Petals",is_p=T,pol="NP",community,"neighbors")
stylenp_com=Distance_Pool(seep_list_com,"Style",is_p=T,pol="NP",community,"neighbors")

anthersnp_flower=Distance_Pool(seep_list,"Anthers",is_p=T,pol="NP",no_ant,"flower")
petalsnp_flower=Distance_Pool(seep_list,"Petals",is_p=T,pol="NP",no_pet,"flower")
stylenp_flower=Distance_Pool(seep_list,"Style",is_p=T,pol="NP",no_st,"flower")

#4. Merge data, summarize and plot
all_data=data.frame(rbind(antherswp_soil,antherswp_leaves,antherswp_com,antherswp_flower,
               petalswp_soil,petalswp_leaves,petalswp_com,petalswp_flower,
               stylewp_soil,stylewp_leaves,stylewp_com,stylewp_flower,
               anthersnp_soil,anthersnp_leaves,anthersnp_com,anthersnp_flower,
               petalsnp_soil,petalsnp_leaves,petalsnp_com,petalsnp_flower,
               stylenp_soil,stylenp_leaves,stylenp_com,stylenp_flower),
               pollinator=rep(c("WP","NP"),each=((15*3)+9)*3))


no_tot=all_data[all_data$div_in!="total",]

sumdf<-ddply(all_data,.(type,organ, div_in,pollinator),summarise,
            Mean=mean(dist),
            stdev=2*sd(dist))
sumdf$div_in=factor(sumdf$div_in,levels=c("total","turnover","nestedness"))

ggplot(sumdf, aes(y=Mean,x=type,colour=type,shape=pollinator))+
  facet_wrap(~div_in*organ)+
  ylab("Beta diversity (Sorensen)")+
  xlab("Potential sources")+
  scale_color_manual(values=c("#FFBB51","#90b86a","#f95252","gray20"))+
  geom_errorbar(position=position_dodge(0.6), color="black",width=0.2,
                aes(ymax=Mean+stdev,ymin=Mean-stdev))+
  geom_point(position=position_dodge(0.6), size=2)+
  scale_shape_manual(values=c(19,1))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45,hjust=1))
  
  
#Data overlaps perfectly for no polllinator and pollinator treatment. So we analysed them together. 
a=aov(dist~div_in*type*organ,data=no_tot)
summary(a)

mu=mean(all_data$dist[all_data$div_in=="total"])
se=sd(all_data$dist[all_data$div_in=="total"])/sqrt(length(all_data$dist[all_data$div_in=="total"]))

g1=interaction(no_tot$div_in,no_tot$type)
m2=lm(dist~g1*organ,data=no_tot)
#PLANNED CONTRASTS
cont=matrix(c(0,0,0,0,0,0,
              1,1,1,0,0,-1,
              0,0,0,0,0,0,
              -1,0,0,0,-1,0,
              0,0,0,0,0,0,
              0,-1,0,-1,0,0,
              0,0,0,0,0,0),ncol=8)

G = glht(m2,linfct = mcp(g1 = cont))

G$linfct

summary(G,test=adjusted("fdr"))

```

To evaluate likelihood of sources I used source tracker (Kights et al. 2011). 
```{r}
#Convert OTU table to integers
otus_int=round(OTU_wash_clean*1000)

source=which(metadata$Organ=="Community"|metadata$Organ=="Leaves"|
                  metadata$Organ=="Soil")
sink=which(metadata$Organ=="Anthers"|metadata$Organ=="Petals"|
                  metadata$Organ=="Style")
envs=metadata$Organ

#

# load SourceTracker package
source('sourcetracker-1.0.1/src/SourceTracker.r') #Can be found in https://github.com/danknights/sourcetracker/releases

#tune the alpha values using cross-validation (this is slow!)
tune.results <- tune.st(OTU_sub[source,], envs[source])
alpha1 <- tune.results$best.alpha1
alpha2 <- tune.results$best.alpha2

# train SourceTracker object on sources
st <- sourcetracker(OTU_sub[source,], envs[source])

# Estimate source proportions in test data
results <- predict(st,OTU_sub[sink,], alpha1=alpha1, alpha2=alpha2)

# Estimate leave-one-out source proportions in training data 
results.train <- predict(st, alpha1=alpha1, alpha2=alpha2)

#Plot results
Proportions=as.data.frame(results$proportions)
Proportions$Organ=metadata$Organ[match(rownames(Proportions),metadata$SampleID)]
Proportions$Pollinator=metadata$Pollinator[match(rownames(Proportions),metadata$SampleID)]
Proportions$Seep=metadata$Seep[match(rownames(Proportions),metadata$SampleID)]
Proportions=Proportions[Proportions$Pollinator!="Transp",]
Proportions$Pollinator=droplevels(Proportions$Pollinator)

#For barplot assign an index to have all plots distributed in the same scale and 
# in order of community as source
Proportions=Proportions[order(Proportions$Community),]

#Add an organ 
Organ=c("Anthers","Petals","Style")
tab=table(paste(Proportions$Pollinator,Proportions$Seep))
Seep=c("BNS","RH1","RHA","RHB","TP9")
  
for (i in 1:5){
  org=Seep[i]
  Proportions$index[Proportions$Pollinator=="NP"&Proportions$Seep==org]=1:tab[i]
  Proportions$index[Proportions$Pollinator=="WP"&Proportions$Seep==org]=1:tab[i+5]
}

Prop_long=gather(Proportions,Source,Proportion,Community:Unknown,factor_key = T)


ggplot(Prop_long, aes(x=index,y=Proportion,fill=Source))+
  geom_bar(stat="identity")+
  facet_grid(Pollinator~Seep)+
  scale_fill_manual(values=c("#F77F71","#91B24D","gray20","gray75"))+
  theme_bw()

sumdf_prop<-ddply(Prop_long,.(Source,Organ,Pollinator),summarise,
            Mean=mean(Proportion))

ggplot(sumdf_prop, aes(x=Organ,y=Mean,fill=Source))+
  geom_bar(stat="identity")+
  facet_wrap(~Pollinator, nrow = 2)+
  scale_fill_manual(values=c("#F77F71","#91B24D","gray20","gray75"))+
  theme_bw()

```

