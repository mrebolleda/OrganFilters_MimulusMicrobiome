##Paper: Floral organs act as environmental filters and interact with pollinators to structure the yellow monkeyflower *Mimulus guttatus* floral microbiome.
###Code by: María Rebolleda Gómez
####Last updated: July 2019

Set directory, load packages and export data (data will be available in Dryad upon publication)
```{r}
library("vegan")
library("phyloseq")
library("ggplot2")
library("tidyr")
library(plyr)
library(VennDiagram)
library(picante)
library(ape)

directory<- #write here in quotation marks the path to files
setwd(directory)

treefile = "16s/rep_set.tre"
taxa=as.matrix(read.table("16s/Taxa.txt",sep="\t"))
colnames(taxa)=c("Kingdom","Phyllum","Class","Order","Family","Genus","Species")

#Open processed reads from OTUtable_Processing
OTU_tabl=read.table(file="16s/OTU_rarefied_organs.txt",sep="\t",check.names = F)

#Open sample information for all samples (mymap) and only organ washes (metadata)
mymap <- import_qiime_sample_data("MIGU_metadata.txt")
metadata=read.table(file="16s/Metadata_subset.txt",sep="\t")
metadata$SampleID=as.character(metadata$SampleID)
rownames(metadata)=metadata$SampleID


#Open flower data
fl1=read.csv("FlowerAbundance_170510.csv")
fl2=read.csv("FlowerAbundance_170516_Block2.csv")

#Open pollinator data
pol=read.csv("../RAW/EnvironmentalCovariates/PollinatorData.csv")

#Color palettes
cols=c("#e59423","#e3072a","#10a8b0")

#I removed samples that could have been contaminated by leaves- results are the same. 
OTU_tabl=OTU_tabl[rownames(OTU_tabl)!="306"&rownames(OTU_tabl)!="311",]
metadata=metadata[metadata$SampleID!="306"&metadata$SampleID!="311",]
```

Test for variables shaping differentiation across communities with different beta diversity indeces. This step needs to be performed for most of the analyses. 
```{r}
#This fuction calculates distances using Sorensen, Bray-Curtis, Unifrac and wunifrac with OTU table, tree, taxa and metadata.
dist=function(OTU,tree,meta,taxa){
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
  phylo_data=phyloseq(sample_data(meta),
                      otu_table(t(OTU),taxa_are_rows=T),
                      phy_tree(tree.new),tax_table(taxa_s))
  #3. Calculate distances and append to distances list
  distances$unifrac=phyloseq::distance(phylo_data,method="unifrac")
  distances$wunifrac=phyloseq::distance(phylo_data,method="wunifrac")
  distances
}


distances=dist(OTU_tabl,treefile,metadata,taxa)


#Formula to obtain PCoA with percentage variation from an OTU table
PCoA_all=function(x){
  pcoa.x=cmdscale(x,eig=T)
  #add percentage of variation explained by each of the first two axis
  pcoa.x$v=c(round(pcoa.x$eig[1]/sum(pcoa.x$eig)*100),
              round(pcoa.x$eig[2]/sum(pcoa.x$eig)*100))
  pcoa.x
}

#Do PCoA for each distance table obtained
ord=lapply(distances,PCoA_all)
str(ord)

#Check that order of points is the same for the different measures
head(rownames(ord$dist_soren$points))
head(rownames(ord$unifrac$points))

#Plot ordinations
shapes=c(19,1)
par(oma=c(1, 1, 1, 8), xpd=TRUE)
par(mar=c(4,4,4,4))
par(mfrow=c(2,2))

for (i in 1:4){
    plot(ord[[i]]$points,pch=shapes[metadata$Pollinator],col=cols[metadata$Organ],
       xlab=paste("PCoA1(",ord[[i]]$v[1],"%)",sep=""),
       ylab=paste("PCoA2 (",ord[[i]]$v[2],"%)",sep=""), main=names(ord)[i])
}

#Add legend to the 4-way plot
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend(0.58,0.85,legend=c("Anthers","Petals","Style"),pch=19,col=cols,bty ="n",cex=1.2)
legend(0.58,0.55, legend=levels(metadata$Pollinator),pch=c(1,19),bty ="n",cex=1.2)

dev.off()
par(mar=c(4,4,1,8))
plot(ord[[3]]$points,pch=shapes[metadata$Pollinator],col=cols[metadata$Organ],
       xlab=paste("PCoA1(",ord[[3]]$v[1],"%)",sep=""),
       ylab=paste("PCoA2 (",ord[[3]]$v[2],"%)",sep=""))
par(xpd=NA)
legend(0.3,0.32,legend=c("Anthers","Petals","Style"),pch=15,
       col=cols,bty ="n",cex=1.2)
legend(0.3,0.2, legend=levels(metadata$Pollinator),pch=c(1,19),bty ="n",cex=1.2)

dev.off()
par(mar=c(4,4,1,8))
seep_col=c("#f69053","#def2b4","#91cba9","#2b83ba","#d7191c")
plot(ord[[3]]$points,pch=shapes[metadata$Pollinator],col=seep_col[metadata$Seep],
       xlab=paste("PCoA1(",ord[[3]]$v[1],"%)",sep=""),
       ylab=paste("PCoA2 (",ord[[3]]$v[2],"%)",sep=""))
par(xpd=NA)
legend(0.3,0.32,legend=levels(metadata$Seep),pch=15,
       col=seep_col,bty ="n",cex=1.2)
legend(0.3,0.1, legend=levels(metadata$Pollinator),pch=c(1,19),bty ="n",cex=1.2)

#############################################################
########## Test differences between organs and seep##########
#############################################################
#Permanovas
permanovas=vector("list",4L)
for (i in 1:4){
  permanovas[[i]]=adonis(distances[[i]]~metadata$Organ*metadata$Pollinator+metadata$Organ*metadata$Seep+metadata$Pollinator*metadata$Seep)
}

permanovas_nopol=vector("list",4L)
permanovas_nopol=vector("list",4L)
for (i in 1:4){
  permanovas_nopol[[i]]=adonis(distances[[i]]~metadata$Organ*metadata$Seep)
}


#Beta-disperser (test assumptions for permanova)
betadisp=vector("list",4L)
for (i in 1:4){
  betadisp[[i]]=betadisper(distances[[i]], metadata$Organ, bias.adjust = FALSE, sqrt.dist = FALSE, add = FALSE)
}

#Plot
shapes=c(0,16,17,15,1)
par(oma=c(1, 1, 1, 8), xpd=TRUE)
par(mar=c(4,4,4,4))
par(mfrow=c(2,2))

for (i in 1:4){
  boxplot(betadisp[[i]], col=cols, main=names(ord)[i])
}

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend(0.6,0.8,legend=c("Anthers","Petals","Style"),pch=15,col=cols,bty ="n",cex=1.2)

lapply(betadisp,permutest)


metadata$bray_betadisp=betadisp[[3]]$distances

```


Obtain distances and model fit for pollinator no pollinator
```{r}

pol.meta=metadata[metadata$Pollinator=="WP",]
sub.pol=pol.meta$SampleID
OTU.pol=OTU_tabl[match(as.character(sub.pol),rownames(OTU_tabl)),]

dist.pol=dist(OTU.pol,treefile)

perm.pol=vector("list",4L)
for (i in 1:4){
  perm.pol[[i]]=adonis(dist.pol[[i]]~pol.meta$Organ*pol.meta$Seep)
}

npol.meta=metadata[metadata$Pollinator=="NP",]
sub.npol=npol.meta$SampleID
OTU.npol=OTU_tabl[match(as.character(sub.npol),rownames(OTU_tabl)),]

dist.npol=dist(OTU.npol,treefile)

perm.npol=vector("list",4L)
for (i in 1:4){
  perm.npol[[i]]=adonis(dist.npol[[i]]~npol.meta$Organ*npol.meta$Seep)
}

Rpol=unlist(lapply(perm.pol, function(x) {x$aov.tab$R2[c(1:3)]*100}))
Rnpol=unlist(lapply(perm.npol, function(x) {x$aov.tab$R2[c(1:3)]*100}))

r_comp=data.frame(R=c(Rpol,Rnpol),Pol=rep(c("WP","NP"),each=12),
                  Factor=rep(c("Organ","Seep","Organ*Seep"),8),
                  index=rep(rep(c("Sorensen","Bray-Curtis",
                                  "Unifrac","Weighted Unifrac"),each=3),2))

ggplot(r_comp[r_comp$Factor=="Organ",],aes(x=index,y=R,colour=Pol,shape=Pol))+
  geom_point(size=3)+
  xlab("Beta diversity index")+
  ylab("% variation explained by 'Organ'")+
  scale_colour_manual(values=c("gray50","gray20"))+
  scale_shape_manual(values=c(1,19))+
  theme_bw()

```


Plot relative abundance of most abundant taxa (order level) by organ and pollinator treatment. 
```{r}
#Open taxa table as data.frame
taxa_s_df=as.data.frame(taxa)

#Calculate relative abundances for each OTU
OTU_relab=OTU_tabl/rowSums(OTU_tabl)
#sanity check
rowSums(OTU_relab)

# Merge taxa and OTU table 
OTU_tax=merge(taxa,t(OTU_relab),by=0)


#Get most abundant genus for each organ
organ=c("Anthers", "Style","Petals")
L=vector("list", length(sub)) 

By_Genus=function(x){
  tapply(x,OTU_tax$Genus,sum)
}

for (i in 1:length(organ)){
  org=organ[i]
  s=c("Genus",metadata$SampleID[metadata$Organ==org])
  otu_sub=OTU_tax[,match(s,colnames(OTU_tax))]
  L[[i]]=as.data.frame(lapply(otu_sub[,2:ncol(otu_sub)],By_Genus))
  }

#Order by most abundant 
mostab=function(x){
  x[order(rowSums(x),decreasing=T),]
}

otu_order_byorgan=lapply(L, mostab)

mean_abundance=function(x){
  mean=rowMeans(x,na.rm=T)
  UP95=apply(x,1,quantile,probs=0.975,na.rm=T)
  LOW95=apply(x,1,quantile,probs=0.025,na.rm=T)
  data.frame(mean,UP95,LOW95)
}

otu_list_sum=lapply(otu_order_byorgan,mean_abundance)


L2=vector("list", length(sub)) 
for (i in 1:length(organ)){
  names_ab=rownames(otu_list_sum[[i]][1:30,])
    L2[[i]]=data.frame(genus=substr(names_ab,4,30),
    rbind(otu_list_sum[[1]][match(names_ab,rownames(otu_list_sum[[1]])),],
    otu_list_sum[[2]][match(names_ab,rownames(otu_list_sum[[2]])),],
    otu_list_sum[[3]][match(names_ab,rownames(otu_list_sum[[3]])),]),
    org=rep(organ,each=30))
}
  
ggplot(L2[[1]],aes(x=genus,y=mean,fill=org))+
  geom_point(shape=21)+
  scale_fill_manual(values=cols)+
  coord_flip()+
  theme_bw()
  
ggplot(L2[[2]],aes(x=genus,y=mean,fill=org))+
  geom_point(shape=21)+
  scale_fill_manual(values=cols)+
  coord_flip()+
  theme_bw()

ggplot(L2[[3]],aes(x=genus,y=mean,fill=org))+
  geom_point(shape=21)+
  scale_fill_manual(values=cols)+
  coord_flip()+
  theme_bw()
  

#Sum relative abundances of OTUs within each order to get an order level relative abundance. 
By_Order=function(x){
  tapply(x,OTU_tax$Order,sum)
}

avg_relab_ord=as.data.frame(lapply(OTU_tax[,9:ncol(OTU_tax)],By_Order))

#Obtain order names and assign to new column
avg_relab_ord$Order=substr(rownames(avg_relab_ord),4,stop=20)

#Re-fromat data frame for substet and plotting
av_rab=gather(avg_relab_ord,SampleID,avg_rab,X294:X540)
av_rab$SampleID=substr(av_rab$SampleID,2,4)

#Subset only orders of more than 25% across samples
x=tapply(av_rab$avg_rab,av_rab$Order,sum)
av_rab=av_rab[av_rab$Order%in%names(x)[x>=0.25],]

phyl_data=merge(metadata,av_rab)

ggplot(data=phyl_data, aes(x=reorder(Order,avg_rab),y=avg_rab, fill=Organ))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(shape=21, alpha=0.6)+
  facet_grid(Organ~Pollinator)+
  stat_smooth(method = "lm", se=F)+
  ylab("Relative abundance")+
  xlab(NULL)+
  scale_fill_manual(values=cols)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45,hjust=1))


```


Beta diversity differences across pollinator and organ treatment, and across open styles
```{r}
dev.off()
beta.pol=betadisper(distances[[2]], paste(metadata$Organ,metadata$Pollinator), bias.adjust = FALSE, sqrt.dist = FALSE, add = FALSE)

# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
#           Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)   
# Groups     5 0.23652 0.047303 3.8721    999  0.003 **
# Residuals 86 1.05062 0.012217                        
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

col.new=c("gray30","#FF7573")
boxplot(beta.pol, col=col.new)
permutest(beta.pol)
TukeyHSD(beta.pol)

#  Tukey multiple comparisons of means
#     95% family-wise confidence level
# 
# Fit: aov(formula = distances ~ group, data = df)
# 
# $group
#                               diff          lwr          upr     p adj
# Anthers WP-Anthers NP  0.073403908  0.004815252  0.141992563 0.0287529
# Petals NP-Anthers NP   0.021854744 -0.046733911  0.090443399 0.9379943
# Petals WP-Anthers NP  -0.006281573 -0.077201195  0.064638049 0.9998387
# Style NP-Anthers NP    0.038209525 -0.032710097  0.109129147 0.6199702
# Style WP-Anthers NP    0.048595835 -0.019009655  0.116201324 0.2995271
# Petals NP-Anthers WP  -0.051549164 -0.119022482  0.015924154 0.2363343
# Petals WP-Anthers WP  -0.079685481 -0.149527000 -0.009843961 0.0158649
# Style NP-Anthers WP   -0.035194383 -0.105035902  0.034647137 0.6847507
# Style WP-Anthers WP   -0.024808073 -0.091281731  0.041665585 0.8849481
# Petals WP-Petals NP   -0.028136317 -0.097977837  0.041705203 0.8478077
# Style NP-Petals NP     0.016354781 -0.053486739  0.086196300 0.9834373
# Style WP-Petals NP     0.026741091 -0.039732567  0.093214749 0.8485859
# Style NP-Petals WP     0.044491098 -0.027640913  0.116623109 0.4725170
# Style WP-Petals WP     0.054877408 -0.013998832  0.123753648 0.1964757
# Style WP-Style NP      0.010386310 -0.058489930  0.079262550 0.9978568




```

Geography-Bray Curtis: Chose Bray-Curtis because under neutral model of dispersal we would not expect particularly strong signal of phylogeny and we do not expect 
```{r}
#Create pairwise indexes of sample and location
id.grid=expand.grid(metadata$SampleID,metadata$SampleID)
id.new=paste(id.grid$Var1,id.grid$Var2, sep="_")
loc.grid=expand.grid(paste(metadata$Seep,metadata$Loc, sep="."),paste(metadata$Seep,metadata$Loc,sep="."))
loc.new=paste(loc.grid$Var1,loc.grid$Var2, sep="_")

############## Subset by pollinator 
Organs=c("Anthers","Petals","Style")
  
bray.all=function(x){
  sub.pol=metadata[metadata$Pollinator==x,]
  bray.pol <- vector("list", 3)
  for (i in 1:3){
  #Subset by organ
  subset=sub.pol$SampleID[sub.pol$Organ==Organs[i]]
  OTU=OTU_tabl[match(as.character(subset),rownames(OTU_tabl)),]
  #Remove OTUs not present  
  not_0=colSums(OTU)>0
  OTU_clean=OTU[,not_0]
  #Calculate distance matrix and save as matrix for easier re-formating
  matrices = as.matrix(vegdist(OTU_clean))
    #Create data.frames
    bray.pol[[i]]=data.frame(Dist=as.vector(matrices),
                                      Id.row=rep(rownames(matrices),ncol(matrices)),
                                      Id.col=rep(colnames(matrices),each=nrow(matrices)),
                                      id.new=paste(rep(rownames(matrices),ncol(matrices)),
                                                   rep(colnames(matrices),each=nrow(matrices)),sep="_"),
                                      Organ=rep(Organs[i],ncol(matrices)*nrow(matrices)), stringsAsFactors = F)
  }
bray.data=data.frame(rbind(bray.pol[[1]],bray.pol[[2]],bray.pol[[3]]))
bray.data$Pollinator=rep(x,nrow(bray.data))
bray.data
}

bray.pol=bray.all("WP")
bray.npol=bray.all("NP")
bray.pol.orgs=data.frame(rbind(bray.pol,bray.npol))

bray.pol.orgs$loc.id=loc.new[match(bray.pol.orgs$id.new,id.new)]

#OPEN GEOGRAPHICAL, MERGE DATA AND PLOT!
geo=read.csv("Geographic/DistanceMatLong.csv")

geo$loc.id=paste(paste(geo$Seep1,geo$Loc1, sep="."),paste(geo$Seep2,geo$Loc2, sep="."), sep="_")
dist.pol.data=merge(bray.pol.orgs,geo)

bray.data.no0=dist.pol.data[dist.pol.data$Dist!=0,]


#Plot geographic distance by organ
ggplot(data=bray.data.no0, aes(x=Distance,y=Dist, colour=Pollinator))+
  geom_point(alpha=0.6,size=0.5)+
  facet_wrap(~Organ)+
  stat_smooth(method = "lm", se=F)+
  ylab("Bray-curtis")+
  xlab("Distance (m)")+
  scale_color_manual(values=c("gray30","#FF7573"))+
  theme_bw()

```

Mantel tests
```{r}
#Change format of geographic data to matrix for each sample
Geo_batch=function(df,x){
  sub.pol=df[df$Pollinator==x,]
  geo=sub.pol[,c(3,4,12)]
  temp.list=split(geo,sub.pol$Organ)
  geo_wide=lapply(temp.list,spread,key=Id.col,Distance)
  geo_wide_sm=lapply( geo_wide, function(x) { x["Id.row"] <- NULL; x })
  geo_mat=lapply(geo_wide_sm,as.matrix,dimnames=list(geo_wide_sm[1],geo_wide_sm[1]))
  dist=sub.pol[,c(2,3,4)]
  temp.list2=split(dist,sub.pol$Organ)
  dist_wide=lapply(temp.list2,spread,key=Id.col,Dist)
  dist_wide_sm=lapply(dist_wide, function(x) { x["Id.row"] <- NULL; x })
  dist_mat=lapply(dist_wide_sm,as.matrix,dimnames=list(geo_wide_sm[1],geo_wide_sm[1]))
  mapply(mantel,geo_mat,dist_mat,SIMPLIFY = F)
}


Geo_batch(dist.pol.data,"WP")
Geo_batch(dist.pol.data,"NP")


```


Repeat distance analysis and mantel test for leaves as reference
```{r}
#Open processed reads from OTUtable_Processing
OTU_rel_leaves=read.table(file="16s/OTU_leaves.txt",sep="\t",check.names = F)


#Get metadadata in same order for correct labeling 
ml=match(rownames(OTU_rel_leaves),as.character(mymap$SampleID))
meta_leaf=(mymap[ml,])
nrow(meta_leaf)
head(meta_leaf)

id.grid_l=expand.grid(meta_leaf$SampleID,meta_leaf$SampleID)
id.new_l=paste(id.grid_l$Var1,id.grid_l$Var2, sep="_")
loc.grid_l=expand.grid(paste(meta_leaf$Seep,meta_leaf$Loc, sep="."),paste(meta_leaf$Seep,meta_leaf$Loc,sep="."))
loc.new_l=paste(loc.grid_l$Var1,loc.grid_l$Var2, sep="_")  

leaf_dist=as.matrix(vegdist(OTU_rel_leaves))
bray.leaf=data.frame(Dist=as.vector(leaf_dist),
                     Id.row=rep(rownames(leaf_dist),ncol(leaf_dist)),
                     Id.col=rep(colnames(leaf_dist),each=nrow(leaf_dist)),
                     id.new=paste(rep(rownames(leaf_dist),ncol(leaf_dist)),
                                  rep(colnames(leaf_dist),each=nrow(leaf_dist)),sep="_"),
                     Organ=rep("Leaf",ncol(leaf_dist)*nrow(leaf_dist)),
                     Pollinator=rep("NP",ncol(leaf_dist)*nrow(leaf_dist)), stringsAsFactors = F)


bray.leaf$loc.id=loc.new_l[match(bray.leaf$id.new,id.new_l)]
dist.leaf.data=merge(bray.leaf,geo)

leaf_bray.no0=dist.leaf.data[dist.leaf.data$Dist!=0,]

new_geo=rbind(bray.data.no0,leaf_bray.no0)
new_geo$Organ <- factor(new_geo$Organ, levels = c("Anthers", "Petals", "Style","Leaf"))

#Plot geographic distance by organ
ggplot(data=new_geo, aes(x=Distance,y=Dist, colour=Pollinator))+
  facet_wrap(~Organ,ncol=4)+
  geom_point(alpha=0.6,size=0.5)+
  stat_smooth(method = "lm", se=F)+
  scale_color_manual(values=c("gray30","#FF7573"))+
  ylab("Bray-curtis")+
  xlab("Distance (m)")+
  theme_bw()


geoL=dist.leaf.data[,c(3,4,12)]
geo_wide=spread(geoL,key=Id.col,Distance)
geo_wide["Id.row"]=NULL
geo_mat=as.matrix(geo_wide,dimnames=list(geo_wide[1],geo_wide[1]))
dist=dist.leaf.data[,c(2,3,4)]
dist_wide=spread(dist,key=Id.col,Dist)
dist_wide["Id.row"]=NULL
dist_mat=as.matrix(dist_wide,dimnames=list(geo_wide[1],geo_wide[1]))
mantel(geo_mat,dist_mat)


```


Distance between pollinator and no-pollinator samples
```{r}
pairs=read.csv("../Metadata/PollinatorPairs_v2_OTU.csv")
pairs$Pollinator=as.character(pairs$Pollinator)
pairs$No.Pollinator=as.character(pairs$NoPollinator)

#Make distance into matrix object
M=as.matrix(distances)

matrices=function(x,pairs){
  #Make distance into matrix object and select only the distances 
  #between pollinator/no pollinator pairs 
  M=as.matrix(x)
  pol.pairs=M[rownames(M)%in%pairs$No.Pollinator,rownames(M)%in%pairs$Pollinator]
  pairs_new=pairs
  pairs_new$beta=diag(pol.pairs)
  pairs_new
}

pairs_list=lapply(distances, matrices, pairs=pairs)

#####Make data frame with all indeces
#Create burner row
pairs_df=data.frame(pairs_list[[1]][1,],index="remove")

for (i in 1:4){
  index=names(pairs_list)[i]
  pairs_tmp=data.frame(pairs_list[[i]],index=rep(index,nrow(pairs)))
  pairs_df=rbind(pairs_df,pairs_tmp)
}

#Remove burner
pairs_df=pairs_df[-1,]

ggplot(pairs_df, aes(x=Organ, y=beta, fill=Organ))+
  facet_wrap(~index)+
  geom_boxplot()+
  geom_jitter(shape=21)+
  xlab("Flower structure")+
  ylab("Pair-wise beta diversity")+
  scale_fill_manual(values=cols)+
  theme_bw()

ax=lm(beta~Organ*index,data=pairs_df)
anova(ax)

```

To process pollinator data:
```{r}
pol_wide <- spread(pol, Mimulus, Visits)
pol_wide$In[is.na(pol_wide$In)]=0
pol_wide$No[is.na(pol_wide$No)]=0
pol_wide$Out[is.na(pol_wide$Out)]=0
pol_wide$tot=pol_wide$In+pol_wide$No+pol_wide$Out
pol_wide$mimul=pol_wide$In+pol_wide$Out

pol_wide$time_ID=paste(pol_wide$Block,pol_wide$Time, sep="_")
pol_wide=na.omit(pol_wide)


pol.sums=ddply(pol_wide,c("Location","Seep","time_ID"),summarise,tot=sum(tot),
              In=sum(In),Out=sum(Out),mimul=sum(mimul))

#To make figure S6 
totalOut=tapply(pol_wide$Out, pol_wide$FuncionalGroup, sum)
totalIn=tapply(pol_wide$In, pol_wide$FuncionalGroup, sum)

table_cont=cbind(totalIn,totalOut)
chisq.test(table_cont)
fisher.test(table_cont,simulate.p.value=TRUE,B=10000)

Percentage=c((totalOut/sum(totalOut)*100), (totalIn/sum(totalIn)*100))

data_prop=data.frame(Group=names(totalIn),Percentage,
                     type=rep(c("Out","In"),each=14))

data_prop$type <- factor(data_prop$type, levels =c("Out","In"))
data_prop$Group <- factor (data_prop$Group, 
                           levels=c("XS","LS","SB","MB","LB","LL",
                                    "BF","SF","LF","FL","LE",
                                    "BE","VE","OT"))
col_pol=c("#d73027","#a50026","#f46d43","#fdae61","#fee090","#ffffbf",
          "#e0f3f8","#abd9e9","#4575b4","#313695","#F5F1DE",
          "gray20","gray50","gray80")

ggplot(data_prop,aes(y=Percentage,x=type,fill=Group))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values=col_pol)+
  theme_bw()



#Multiply times 4 to get hour rate. 
pol.new=ddply(pol.sums,c("Location","Seep"),summarise,mean_visits=mean(tot)*4,
              mean_in=mean(In)*4,mean_out=mean(Out)*4,
              mean_mimul=mean(mimul)*4)

pol.new$Loc=as.factor(pol.new$Location)
pairs_pol=merge(pairs_df,pol.new)

#Make long format for plot
pairs_pol_long=gather(pairs_pol,pol_type,visits,c(mean_visits,mean_in,mean_mimul))

pairs_pol_long$pol_type=factor(pairs_pol_long$pol_type, levels=c("mean_visits","mean_mimul","mean_in"))

ggplot(pairs_pol_long, aes(x=visits, y=beta, color=Organ, shape=Organ))+
  facet_grid(index~pol_type, scales="free")+
  geom_point()+
  stat_smooth(method="lm", se=F)+
  ylab("Difference")+
  xlab("Mean number of visits")+
  scale_color_manual(values=cols)+
  theme_bw()


outside=pairs_pol_long[pairs_pol_long$pol_type=="mean_mimul",]
inside=pairs_pol_long[pairs_pol_long$pol_type=="mean_in",]
inside=droplevels(inside)
outside=droplevels(outside)

####Function to obtain correlation (visitation-pairwise beta) coeficients and p-values for each diversity index and across organs####

correlations=function (df) {
  #Create empty matrices for for loops
  corr_coef=matrix(nrow=4,ncol=3)
  pvalues=matrix(nrow=4,ncol=3)
  
  #For each div. index and organ obtain p-value and corr-coef
  for (i in 1:4){
    ind=levels(df$index)[i]
    for (j in 1:3){
      org=levels(df$Organ)[j]
      sub=df[df$index==ind&df$Organ==org,]
      corr=cor.test(sub$beta,sub$visits)
      corr_coef[i,j]=corr$estimate
      pvalues[i,j]=corr$p.value
    }
  }
  
  #Adjust p-values for multiple comparisons
  p_adj=p.adjust(pvalues,method = "BH")
  
  #Organize data output
  p_data=data.frame(Index=rep(levels(inside$index),3),
                  Organ=rep(levels(inside$Organ),each=4),p_adj)
  corr_data=data.frame(Index=rep(levels(inside$index),3),
                  Organ=rep(levels(inside$Organ),each=4),
                  cor=as.vector(corr_coef))
  corr_data$Index = factor(corr_data$Index,
                    levels = c("dist_soren", "unifrac", "dist_bray","wunifrac"))
  corr_data$Organ= factor(corr_data$Organ,
                        levels=c("Anthers","Style","Petals"))
  p_data$Index = factor(p_data$Index,
                        levels = c("dist_soren", "unifrac", "dist_bray","wunifrac"))
  p_data$Organ= factor(p_data$Organ, levels=c("Anthers","Style","Petals"))
  list(p_data,corr_data)
}

#Apply function to different pollinator groups
cor_inside=correlations(inside)
cor_migu=correlations(outside)


#Plot data
#1. p-value inside
ggplot(data=cor_inside[[1]],aes(x=Organ,y=Index,fill=p_adj))+
  geom_tile(color="gray20",size=0.2)+
  scale_fill_gradient2(low ="#003c30", mid ="#f6e8c3",
  high ="#543005", midpoint = 0.2)+
  #scale_fill_distiller(type="div",palette ="BrBG",direction=1)+
  xlab("Organ")+
  ylab("Diversity index")+
  theme_bw()

#2. p-value mimulus
ggplot(data=cor_migu[[1]],aes(x=Organ,y=Index,fill=p_adj))+
  geom_tile(color="gray20",size=0.2)+
  scale_fill_gradient2(low ="#003c30", mid ="#f6e8c3",
  high ="#543005", midpoint = 0.2)+
  #scale_fill_distiller(type="div",palette ="BrBG",direction=1)+
  xlab("Organ")+
  ylab("Diversity index")+
  theme_bw()


#Plot correlations side by side. 
#Create new data for ploting
cor_all=data.frame(rbind(cor_migu[[2]],cor_inside[[2]]),
                      type=rep(c("all","inside"),each=12))

ggplot(data=cor_all,aes(x=Organ,y=Index,fill=cor))+
  facet_wrap(~type)+
  geom_tile(color="gray20",size=0.2)+
  #scale_fill_gradient2(low ="#003c30", mid ="#f6e8c3",
  #high ="#543005", midpoint = 0.3)+
  scale_fill_distiller(type="div",palette ="BrBG",direction=1)+
  xlab("Organ")+
  ylab("Diversity index")+
  theme_bw()

```


Make Venn-Diagrams with core microbiome and calculate phylogenetic distance
```{r}
# 1. Unique samples
Organs=c("Anthers","Petals","Style")

uniqueOTU=function(x,p){
  OTU_list=vector(mode="list",length = 3L)
  sub.pol=metadata[metadata$Pollinator==x,]
  for (i in 1:3){
    subset=sub.pol$SampleID[sub.pol$Organ==Organs[i]]
    OTU=OTU_tabl[match(as.character(subset),rownames(OTU_tabl)),]
    OTU_pa=(OTU>0)*1
    OTU_list[[i]]=colnames(OTU_pa)[colSums(OTU_pa)/nrow(OTU_pa)>p]
  }
  OTU_list
}

prop=seq(0,1,by=0.1)

Shared=function(p,pol){
  M=matrix(ncol=2,nrow=length(p), dimnames=list(p,c("N","NShared")))
  for (i in 1:length(p)){
    uni=uniqueOTU(pol,p[i])
    if (length(uni[[1]])==0&length(uni[[2]])==0&length(uni[[3]])==0){
      M[i,1]=0
      M[i,2]=0
    } else {
      over=get.venn.partitions(uni)
      M[i,1]=sum(over$..count..)
      M[i,2]=sum(over$..count..[c(4,6,7)])/M[i,1]
    }
  }
  M
}

pol_uniqueOTU=Shared(prop,"WP")
npol_uniqueOTU=Shared(prop,"NP")

prop_shared=data.frame(rbind(pol_uniqueOTU,npol_uniqueOTU))
prop_shared$p=rep(prop,2)
prop_shared$Pollinator=rep(c("WP","NP"),each=length(prop))



ggplot(prop_shared,aes(x=p,y=N,colour=Pollinator))+
  geom_point()+
  geom_line()+
  xlab("Minimal proportion of samples with OTU present")+
  ylab("Total number of OTUs")+
  scale_color_manual(values=c("gray80","gray20"))+
  theme_bw()

ggplot(prop_shared,aes(x=p,y=log(N),colour=Pollinator))+
  geom_point()+
  stat_smooth(method="lm")+
  xlab("Minimal proportion of samples with OTU present")+
  ylab("log(Total number of OTUs)")+
  scale_color_manual(values=c("gray80","gray20"))+
  theme_bw()

ggplot(prop_shared,aes(x=p,y=NShared,colour=Pollinator))+
  geom_point()+
  geom_line()+
  ylab("OTUs present in only one floral organ")+
  xlab("Minimal proportion of samples with OTU present")+
  scale_color_manual(values=c("gray80","gray20"))+
  theme_bw()
   
#Calculate core microbiome at 20% shared across samples of a determinate organ and pollinator treatment. 
pol_uniqueOTU=uniqueOTU("WP",0.20)
npol_uniqueOTU=uniqueOTU("NP",0.20)

#Assign genus to unique OTU
OTU_core=OTU_tax[match(unlist(pol_uniqueOTU),OTU_tax$Row.names),]
abund=data.frame(genus=OTU_core$Genus,mean_ab=rowMeans(OTU_core[,9:ncol(OTU_core)]))
tapply(abund$mean_ab, abund$genus, sum)

pol_part=get.venn.partitions(pol_uniqueOTU)
npol_part=get.venn.partitions(npol_uniqueOTU)

count.mat=matrix(c(sum(pol_part$..count..[c(4,6,7)]),
                   sum(npol_part$..count..[c(4,6,7)]),
                   sum(pol_part$..count..[c(1:3,5)]),
                   sum(npol_part$..count..[c(1:3,5)])),ncol=2)

fisher.test(count.mat)


#
tree=read_tree(treefile)
drop=tree$tip.label[!(tree$tip.label%in%colnames(OTU_tabl))]
tree.new=drop.tip(tree,drop)



#Calculate expected phylogenetic distance and variance
E_pd=expected.pd(tree.new)
V_pd=variance.pd(tree.new,upper.bound = F)
#Variance takes a long time to compute so it is better to save for future analyses
#write.csv(E_pd_all, "expected_pd_all.csv")
#write.csv(V_pda_all, "variance_expectedpd_all.csv")

z_pd=function (x){
  #1. Make new OTU table sampling
  mat=matrix(rep(1,length(x)),nrow=1,
           dimnames=list("samp",x))
  #2. Calculate observed pd
  PD=pd(mat,tree.new,include.root = F)
  #3. Calculate standarized deviation
  ex=E_pd$expected.pd[PD$SR]
  sd=sqrt(V_pd$variance.pd[PD$SR])
  Z=(PD$PD-ex)/sd
  Z
}

z_pol=unlist(lapply(pol_uniqueOTU,z_pd))
z_npol=unlist(lapply(npol_uniqueOTU,z_pd))
part_data=data.frame(Organ=rep(c("Stamens","Petals","Style"),2),
                     Pollinator=rep(c("WP","NP"),each=3),
                     z=c(z_pol,z_npol))

ggplot(part_data,aes(x=Organ,y=z,colour=Organ,shape=Pollinator))+
  geom_point(size=3)+
  geom_hline(aes(yintercept=-1.96),linetype="dashed")+
  geom_hline(aes(yintercept=1.96),linetype="dashed")+
  geom_hline(aes(yintercept=0))+
  ylab("Deviation from expectation (standarized effect sizes)")+
  scale_color_manual(values=cols)+
  scale_shape_manual(values=c(1,19))+
  theme_bw()


#Repeat the same analyses but for the whole data
#We can use the same function to subsaple per organ, just setting p to 0
pol_OTU=uniqueOTU("WP",0)
npol_OTU=uniqueOTU("NP",0)

z_pol_all=unlist(lapply(pol_OTU,z_pd))
z_npol_all=unlist(lapply(npol_OTU,z_pd))
part_data_all=data.frame(Organ=rep(c("Anthers","Petals","Style"),2),
                     Pollinator=rep(c("WP","NP"),each=3),
                     z=c(z_pol_all,z_npol_all))

ggplot(part_data_all,aes(x=Organ,y=z,colour=Organ,shape=Pollinator))+
  geom_point(size=3)+
  geom_hline(aes(yintercept=-1.96),linetype="dashed")+
  geom_hline(aes(yintercept=1.96),linetype="dashed")+
  geom_hline(aes(yintercept=0))+
  ylab("Deviation from expectation (standarized effect sizes)")+
  scale_color_manual(values=cols)+
  scale_shape_manual(values=c(1,19))+
  theme_bw()

```


```{r}
#Soil samples
subsoil=mymap$SampleID[mymap$Organ=="Soil"]
OTU_soil=Rarefied_OTU[match(as.character(subsoil),rownames(Rarefied_OTU)),]

#Remove OTUs not present  
not_0s=colSums(OTU_soil)>0
ncol(OTU_soil)-sum(not_0s)
OTU_soil_clean=OTU_soil[,not_0s]
dim(OTU_soil_clean)

write.table(OTU_soil_clean,file="16s/OTU_soil.txt",sep="\t")

#Leaves samples
subleave=mymap$SampleID[mymap$Organ=="Leaves"]
OTU_leaves=Rarefied_OTU[match(as.character(subleave),rownames(Rarefied_OTU)),]

#Remove OTUs not present  
not_0l=colSums(OTU_leaves)>0
ncol(OTU_leaves)-sum(not_0l)
OTU_leaves_clean=OTU_leaves[,not_0l]
dim(OTU_leaves_clean)

write.table(OTU_leaves_clean,file="16s/OTU_leaves.txt",sep="\t")


#Community samples
subcom=mymap$SampleID[mymap$Organ=="Community"]
OTU_com=Rarefied_OTU[match(as.character(subcom),rownames(Rarefied_OTU)),]
remove=which(is.na(rownames(OTU_com)))
OTU_com=OTU_com[-remove,]

#Remove OTUs not present  
not_0c=colSums(OTU_com)>0
ncol(OTU_com)-sum(not_0c)
OTU_com_clean=OTU_com[,not_0c]
dim(OTU_com_clean)

write.table(OTU_com_clean,file="16s/OTU_community.txt",sep="\t")


```

Procrsutes analysis to test impact of flower community on microbial community
```{r}
#Average flower community and format matrix
L_fl=list(as.matrix(fl1[,-1]),as.matrix(fl2[,-1]))
M_fl=do.call(cbind,L_fl)
M_fl=array(M_fl, dim=c(dim(L_fl[[1]]), length(L_fl)))
fl_mean=apply(M_fl,c(1,2),mean)
dimnames(fl_mean)=list(fl1$X,colnames(as.matrix(fl1[,-1])))
fl_mean=t(fl_mean)

#Get matrix OTU for each organ/treatment
Organs=c("Anthers","Petals","Style")

#Make list of OTU tables by organ following format of fl data
organ_list=function(x,meta){
  OTU_list=vector(mode="list",length = 3L)
  for (i in 1:3){
    subset=meta$SampleID[meta$Organ==Organs[i]]
    OTU_list[[i]]=x[match(as.character(subset),rownames(x)),]
    loc=paste(meta$Seep,meta$Loc,sep="_")
    rownames(OTU_list[[i]])=loc[match(rownames(OTU_list[[i]]),meta$SampleID)]
  }
  OTU_list
}

otus.organ.pol=organ_list(OTU.pol,pol.meta)
otus.organ.npol=organ_list(OTU.npol,npol.meta)

#Get same number of rows for matrices to compare and run procrustes. 
proc=function(x,t){
  Xfl=fl_mean[match(rownames(x),rownames(fl_mean)),]
  mds_fl=metaMDSdist(Xfl,distance="bray")
  mds_x=metaMDSdist(x,distance="bray")
  protest(mds_fl,mds_x)
}

proc_pol=lapply(otus.organ.pol, proc)
lapply(proc_pol, summary)
lapply(proc_pol, plot, kind=1)
lapply(proc_pol, plot, kind=2)


proc_npol=lapply(otus.organ.npol, proc)
lapply(proc_npol, summary)
lapply(proc_npol, plot, kind=1)
lapply(proc_npol, plot, kind=2)
```

